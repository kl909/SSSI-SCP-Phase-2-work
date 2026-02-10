library(galah)
library(rgbif)              #for downloading datasets from gbif
library(countrycode)        #for getting country names based on countryCode important
library(rnaturalearth)      #for downloading maps
library(sf)                 #for manipulating downloaded maps
library(tidyverse)          #for tidy analysis
library(CoordinateCleaner)  #for quality checking of occurrence data
library(dplyr)
library(stringr)
library(here)

# Pantheon species list
pantheon_data <- read_csv(here("pantheon_w2_species_list.csv"))


# Pull species from Pantheon dataframe
species_list <- pantheon_data %>%
  pull(Species) %>%
  unique()

get_species_occ <- function(species_name, country = "GB", year = "1970,*") {
  
  # 1. Get the taxon key
  key <- name_backbone(species_name)$usageKey
  if (is.null(key)) {
    warning(paste("No taxonKey found for", species_name))
    return(NULL)
  }
  
  # 2. CHECK THE COUNT (whether under 100k or not)
  total_count <- occ_count(taxonKey = key, country = country, year = year)
  message(paste(species_name, "has", total_count, "records."))
  
  # 3. DECIDE METHOD
  if (total_count < 100000) {
    # --- METHOD A: Your existing loop for small/medium datasets ---
    message("Using occ_search method...")
    start <- 0
    limit <- 50000 
    all_data <- list()
    
    repeat {
      res <- occ_search(taxonKey = key, country = country, year = year,
                        limit = limit, start = start)
      if (length(res$data) == 0 || nrow(res$data) == 0) break
      all_data[[length(all_data) + 1]] <- res$data
      if (nrow(res$data) < limit) break 
      start <- start + limit
    }
    return(do.call(rbind, all_data))
    
  } else {
    # --- METHOD B: The "Download" method for >100k (Colleague's logic) ---
    message("Switching to Download method (this will take a few minutes)...")
    
    # This triggers the request on GBIF servers
    dl_request <- occ_download(
      pred("taxonKey", key),
      pred("country", country),
      pred("year", year),
      user = Sys.getenv("GBIF_USER"), 
      pwd = Sys.getenv("GBIF_PWD"), 
      email = Sys.getenv("GBIF_EMAIL")
    )
    
    # Wait for the file to be prepared (Important!)
    occ_download_wait(dl_request)
    
    # Download and Import
    res <- occ_download_get(dl_request) %>% 
      occ_download_import()
    
    return(res)
  }
}

##### end of function

# loop over all species
all_occurrences <- lapply(species_list, get_species_occ) %>%
  bind_rows(.id = "species_index")

# add species names 
all_occurrences$species_name <- species_list[all_occurrences$species_index]

# save to csv
#write.csv(all_occurrences, "test_occurances.csv", row.names = FALSE)

######## start data cleaning
# Removes:
# - missing coordinates
# - duplicates
# - zero / invalid / centroid / institution / GBIF HQ points
# Removes records with >1 km spatial uncertainty
data <- all_occurrences %>%
  filter(
    !is.na(decimalLatitude),
    !is.na(decimalLongitude)
  )

# run coordinate cleaner
clean <- data %>%
  filter(!str_detect(license, regex("by-nc|non-commercial", ignore_case = TRUE))) %>% # filters non-commercial licenses
  cc_dupl(lon = "decimalLongitude", lat = "decimalLatitude", species = "scientificName", additions = c("eventDate")) %>%
  cc_zero(lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_equ(lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_val(lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_sea(lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_cap(buffer = 2000, lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_cen(buffer = 2000, lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_gbif(buffer = 2000, lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_inst(buffer = 2000, lon = "decimalLongitude", lat = "decimalLatitude")

# message to say what was done
message(
  nrow(data) - nrow(clean), " records deleted by coordinate cleaning; ",
  nrow(clean), " records remaining."
)

# filter out coordinate uncertainty
clean_2 <- clean %>%
  filter(
    is.na(coordinateUncertaintyInMeters) |
      coordinateUncertaintyInMeters < 1000
  )

message(
  nrow(clean) - nrow(clean_2), " records deleted by uncertainty filter; ",
  nrow(clean_2), " records remaining."
)  

# make sure species identities are retained
clean_2 <- clean_2 %>%
  mutate(species_name = as.character(species_name))

### remove entries without full date but record how many are removed
# 1. Record the starting number of rows
rows_before <- nrow(clean_2)

# 2. Apply the Date Filter
# We keep records that:
# a) Are NOT NA
# b) Match the pattern YYYY-MM-DD (Start with 4 digits, hyphen, 2 digits, hyphen, 2 digits)
clean_3 <- clean_2 %>%
  filter(!is.na(eventDate)) %>% 
  # Step A: Must start with YYYY-MM-DD (Removes "2004" or "2004-06")
  filter(str_detect(eventDate, "^\\d{4}-\\d{2}-\\d{2}")) %>%
  
  # Step B: Must NOT contain a slash (Removes ranges like "2004-06-15/2014-...")
  filter(!str_detect(eventDate, "/"))

# 3. Calculate and Print the removed count
rows_removed <- rows_before - nrow(clean_3)

message(paste(rows_removed, "records removed due to missing/incomplete dates; ", 
              nrow(clean_3), "records remaining."))


###### convert coordinates to grid references
#1. Ensure we only have valid coordinates
gbif_to_convert <- clean_3 %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude))

# 2. Convert to Spatial Object (WGS84 -> OSGB36)
gbif_sf <- st_as_sf(gbif_to_convert, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  st_transform(crs = 27700) # Convert to British National Grid

# 3. Define the Grid Reference Lookup (Same matrix as before)
grid_letters <- matrix(
  c("SV", "SW", "SX", "SY", "SZ", "TV", "TW",
    "SQ", "SR", "SS", "ST", "SU", "TQ", "TR",
    "SL", "SM", "SN", "SO", "SP", "TL", "TM",
    "SF", "SG", "SH", "SJ", "SK", "TF", "TG",
    "SA", "SB", "SC", "SD", "SE", "TA", "TB",
    "NV", "NW", "NX", "NY", "NZ", "OV", "OW",
    "NQ", "NR", "NS", "NT", "NU", "OQ", "OR",
    "NL", "NM", "NN", "NO", "NP", "OL", "OM",
    "NF", "NG", "NH", "NJ", "NK", "OF", "OG",
    "NA", "NB", "NC", "ND", "NE", "OA", "OB",
    "HV", "HW", "HX", "HY", "HZ", "JV", "JW"),
  ncol = 7, byrow = TRUE
)


# 4. Generate the 1km Grid Ref
gbif_gridref <- gbif_sf %>%
  mutate(
    # Get raw Easting/Northing (in meters)
    east = st_coordinates(.)[, 1],
    north = st_coordinates(.)[, 2],
    
    # Calculate indices (100km tiles)
    idx_e = floor(east / 100000) + 1,
    idx_n = floor(north / 100000),
    
    # Calculate Matrix Row (South to North)
    mat_row = idx_n + 1,
    
    # Lookup the Letters (e.g., "TL")
    grid_let = mapply(function(r, c) tryCatch(grid_letters[r, c], error=function(e) NA), mat_row, idx_e),
    
    # Calculate the 1km digits (e.g., "12" and "34")
    e_1km = sprintf("%02d", floor((east %% 100000) / 1000)),
    n_1km = sprintf("%02d", floor((north %% 100000) / 1000)),
    
    # Combine WITHOUT spaces (e.g., "TL1234")
    gridReference = paste0(grid_let, e_1km, n_1km) 
  ) %>%
  # Remove helper columns to keep it clean
  select(-east, -north, -idx_e, -idx_n, -mat_row, -grid_let, -e_1km, -n_1km)

write_csv(gbif_gridref, here("gbif_occurrences_cleaned.csv"))

##### remove unnecessary columns
final_data <- st_drop_geometry(gbif_gridref) %>%
  select(species, gridReference, eventDate)

write_csv(final_data, here("gbif_final_data.csv"))

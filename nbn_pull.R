library(galah)
library(rgbif)
library(dplyr)
library(CoordinateCleaner)
#library(lubridate)
library(sf)
library(stringr)

setwd("/home/kl909/Documents/NE_postdoc/my_scripts")

########### 1 species test ################

# galah_config(atlas = "United Kingdom",
#              email = "kl909@york.ac.uk",
#              expected_licences_nbn = c("OGL", "CC0", "CC-BY", "CC-BY-NC"),
#              download_reason_id = 17)
# 
# 
# species_name <- "Haliplus immaculatus"
# 
# data <- galah_call() |>
#   galah_identify(species_name) |>
#   galah_filter(year >= 1970) |>
#   # Manually select standard Darwin Core columns to avoid the error
#   galah_select(scientificName, 
#                decimalLatitude, 
#                decimalLongitude,
#                coordinateUncertaintyInMeters,
#                eventDate, 
#                year) |> 
#   atlas_occurrences()
# 
# head(data)
# 
# 
# # clean data
# #Flagging records with problematic occurrence information using functions of the coordinatecleaner package.
# test_clean <- data %>%
#   filter(!is.na(decimalLatitude),
#          !is.na(decimalLongitude)) %>%
#   cc_dupl(species = "scientificName") %>%
#   cc_zero() %>%
#   cc_equ() %>%
#   cc_val() %>%
#   cc_sea() %>%
#   cc_cap(buffer = 2000) %>%
#   cc_cen(buffer = 2000) %>%
#   cc_gbif(buffer = 2000) %>%
#   cc_inst(buffer = 2000)
# print(paste0(nrow(data)-nrow(test_clean), " records deleted; ",
#              nrow(test_clean), " records remaining."))
# 
# # Removing records with coordinate uncertainty and precision issues
# test_clean_2 <- test_clean %>%
#   filter(is.na(coordinateUncertaintyInMeters) |
#            coordinateUncertaintyInMeters < 1000) # removing anything over 1km 
# 
# print(paste0(nrow(data)-nrow(test_clean_2), " records deleted; ",
#              nrow(test_clean_2), " records remaining." ))


################ 5 species test ###########
galah_config(atlas = "United Kingdom",
             email = "kl909@york.ac.uk",
             download_reason_id = 17)

# Species list
species_list <- c(
  "Haliplus immaculatus",
  "Notonecta maculata",
  "Plea leachii", #Plea minutissima (NBN doesn't have records under minutissima)
  "Ischnura pumilio",
  "Somatochlora metallica"
)

# 2. DEFINE FUNCTION ---------------------------------------------------------
# Function to get all NBN occurrences since 1970 for one species
get_nbn_occ <- function(species_name) {
  
  message(paste("Downloading records for:", species_name))
  
  # Pause for 5 seconds to be polite to the server (prevents some 403s)
  Sys.sleep(5)
  
  tryCatch({
    df <- galah_call() |>
      galah_identify(species_name) |>
      galah_filter(year >= 1970,
                   license != "CC-BY-NC"  # <--- THIS excludes Non-Commercial data
                   ) |>
      galah_select(scientificName, 
                   decimalLatitude, 
                   decimalLongitude, 
                   coordinateUncertaintyInMeters,
                   gridReference,
                   eventDate, 
                   year,
                   license) |> 
      atlas_occurrences()
    
    # Check if empty
    if (is.null(df) || nrow(df) == 0) {
      warning(paste("No records found for", species_name))
      return(NULL)
    }
    
    # Force column types to prevent bind_rows() errors
    df <- df |>
      mutate(
        decimalLatitude = as.numeric(decimalLatitude),
        decimalLongitude = as.numeric(decimalLongitude),
        year = as.integer(year),
        coordinateUncertaintyInMeters = as.numeric(coordinateUncertaintyInMeters)
      )
    
    return(df)
    
  }, error = function(e) {
    warning(paste("Failed to download data for", species_name, ":", e$message))
    return(NULL)
  })
}


# 3. EXECUTE DOWNLOAD --------------------------------------------------------
# Loop over all species and combine into one dataframe
# We use bind_rows straight away; 'scientificName' is already in the data
all_occurrences <- lapply(species_list, get_nbn_occ) %>%
  bind_rows()

# Quick check of what we got
print(table(all_occurrences$scientificName))


# 4. DATA CLEANING -----------------------------------------------------------
# Filter 1: Basic Coordinates
data_initial <- all_occurrences %>%
  filter(
    !is.na(decimalLatitude),
    !is.na(decimalLongitude)
  )

# Run CoordinateCleaner
# Note: We specify species="scientificName" to match galah's output
clean <- data_initial %>%
  cc_dupl(lon = "decimalLongitude", lat = "decimalLatitude", species = "scientificName", additions = c("eventDate")) %>%
  cc_zero(lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_equ(lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_val(lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_sea(lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_cap(buffer = 2000, lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_cen(buffer = 2000, lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_gbif(buffer = 2000, lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_inst(buffer = 2000, lon = "decimalLongitude", lat = "decimalLatitude")

# Message 1
message(
  nrow(data_initial) - nrow(clean), " records deleted by coordinate cleaning; ",
  nrow(clean), " records remaining."
)

######## Filter 2: Coordinate Uncertainty (< 1km)
clean_2 <- clean %>%
  filter(
    is.na(coordinateUncertaintyInMeters) |
      coordinateUncertaintyInMeters < 1000
  )

# Message 2
message(
  nrow(clean) - nrow(clean_2), " records deleted by uncertainty filter; ",
  nrow(clean_2), " records remaining."
)


########### Filter 3: remove entries without full date but record how many are removed
# A. Record the starting number of rows
rows_before <- nrow(clean_2)

# B. Apply the Date Filter
# We keep records that:
# a) Are NOT NA
# b) Match the pattern YYYY-MM-DD (Start with 4 digits, hyphen, 2 digits, hyphen, 2 digits)
clean_3 <- clean_2 %>%
  filter(!is.na(eventDate)) %>% 
  # Step A: Must start with YYYY-MM-DD (Removes "2004" or "2004-06")
  filter(str_detect(eventDate, "^\\d{4}-\\d{2}-\\d{2}")) %>%
  
  # Step B: Must NOT contain a slash (Removes ranges like "2004-06-15/2014-...")
  filter(!str_detect(eventDate, "/"))

# C. Calculate and Print the removed count
rows_removed <- rows_before - nrow(clean_3)

message(paste(rows_removed, "records removed due to missing/incomplete dates; ", 
              nrow(clean_3), "records remaining."))

##### convert all ultra high reso gridreference codes to 1km
nbn_gridref <- clean_3 %>%
  mutate(
    # Remove spaces first (e.g. "TL 1234 5678" -> "TL12345678")
    gr_clean = str_replace_all(gridReference, " ", ""),
    letters = str_sub(gr_clean, 1, 2),
    digits  = str_sub(gr_clean, 3),
    n_digits = nchar(digits),
    # Calculate how many digits represent the Easting (half the total digits)
    half_len = n_digits / 2,
    
    # Extract the 1km component (The 1st and 2nd digit of Easting, 1st and 2nd of Northing)
    easting_1km = str_sub(digits, 1, 2), 
    northing_1km = str_sub(digits, half_len + 1, half_len + 2),
    
    # overwrite the existing gridReference column
    gridReference = paste0(letters, easting_1km, northing_1km)
  ) %>%
  # Filter out vague records (e.g. 10km squares which only had 2 digits total)
  filter(n_digits >= 4) %>%
  
  # 5. Remove the temporary columns so the dataset stays clean
  select(-gr_clean, -letters, -digits, -n_digits, -half_len, -easting_1km, -northing_1km)


head(nbn_gridref$gridReference)

# 5. SAVE
write.csv(nbn_gridref, "nbn_occurrences_cleaned.csv", row.names = FALSE)

##### remove unnecessary columns
final_data <- nbn_gridref %>%
  select(species = scientificName, gridReference, eventDate)

write.csv(final_data, "nbn_final_data.csv", row.names = FALSE)


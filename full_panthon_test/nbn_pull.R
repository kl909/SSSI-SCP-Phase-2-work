library(galah)
library(rgbif)
library(dplyr)
library(CoordinateCleaner)
#library(lubridate)
library(sf)
library(stringr)
library(here)
library(readr)

###### 1. get species list
# load Pantheon species csv
pantheon_data <- read_csv(here("pantheon_w2_species_list.csv"))

# Pull species from Pantheon dataframe
species_list <- pantheon_data %>%
  pull(Species) %>%
  gsub("\\s*\\([^\\)]+\\)", "", .) %>%  # Removes "(Ranatra)" etc.
  unique() %>%
  na.omit


##### 2. create directory 
dir.create(here("data", "nbn_temp"), recursive = TRUE, showWarnings = FALSE)

##### 3. load configurations
galah_config(atlas = "United Kingdom", 
             email = "kl909@york.ac.uk")

# set reason
galah_config(download_reason_id = 17)

# 4. THE LOOP
for (sp in species_list) {
  
  # Define the checkpoint file path
  file_path <- here("data", "nbn_temp", paste0(gsub(" ", "_", sp), ".csv"))
  
  # STEP A: Check if we already have it (The "Colleague" Strategy)
  if (file.exists(file_path)) {
    message(paste("skipping (already have data for):", sp))
    next
  }
  
  # STEP B: Try to download with retries
  success <- FALSE
  attempt <- 1
  max_retries <- 3
  
  while (!success && attempt <= max_retries) {
    message(paste("Checking NBN for:", sp, "| Attempt:", attempt))
    
    dat <- tryCatch({
      galah_call() %>%
        galah_identify(sp) %>%
        galah_filter(year >= 1970, license != "CC-BY-NC") %>%
        galah_select(scientificName, decimalLatitude, decimalLongitude, 
                     coordinateUncertaintyInMeters, gridReference, eventDate, year) %>%
        atlas_occurrences()
    }, error = function(e) {
      if (grepl("403", e$message)) return("403_ERROR")
      return(NULL)
    })
    
    # NEW LOGIC: If it's a 403, it's likely sensitive. Skip it after 1 try.
    if (identical(dat, "403_ERROR")) {
      message(paste(" 403 Forbidden for", sp, "- likely sensitive/protected. Skipping species."))
      
      # Save an empty file so the loop skips it next time you run the script
      write_csv(data.frame(noted_error = "403_Forbidden_Sensitive"), file_path)
      
      success <- TRUE # We set success to TRUE just to break the WHILE loop
      break           # Break the while loop immediately
      
    } else if (is.null(dat) || nrow(dat) == 0) {
      write_csv(data.frame(), file_path)
      message(paste("No records found for:", sp))
      success <- TRUE 
    } else {
      write_csv(dat, file_path)
      message(paste("Successfully saved:", sp))
      success <- TRUE
    }
  }
  
  # STEP C: Mandatory pause between different species
  if (success) {
    message("Waiting 10 seconds before next species...")
    Sys.sleep(10)
  }
}

##### 5. combine all files
all_occurances <- list.files(here("data", "nbn_temp"), full.names = TRUE) %>%
  lapply(read_csv, show_col_types = FALSE) %>%
  bind_rows()

# 4. DATA CLEANING -----------------------------------------------------------
# Filter 1: Basic Coordinates
data_initial <- all_occurances %>%
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



# 5. SAVE
write.csv(nbn_gridref, "nbn_occurrences_cleaned.csv", row.names = FALSE)

##### remove unnecessary columns
final_data <- nbn_gridref %>%
  select(species = scientificName, gridReference, eventDate)

write.csv(final_data, "nbn_final_data.csv", row.names = FALSE)


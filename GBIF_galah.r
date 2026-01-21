library(galah)
library(rgbif)              #for downloading datasets from gbif
library(countrycode)        #for getting country names based on countryCode important
library(rnaturalearth)      #for downloading maps
library(sf)                 #for manipulating downloaded maps
library(tidyverse)          #for tidy analysis
library(CoordinateCleaner)  #for quality checking of occurrence data
library(dplyr)

setwd("/home/kl909/Documents/NE_postdoc/my_scripts")
# galah wasn't working due to a version gap so using rgbif to pull data instead
# species:
# Haliplus immaculatus - least concern - + 500 obs
# Notonecta maculata - least concern - + 500 obs
# Plea minutissima - least concern - + 500  obs
# Ischnura pumilio - near threatened - + 500 obs
# Somatochlora metallica - vulnerable - + 500 obs

##################### TEST DATA (1 SPECIES)#########
# USING RGBIF
test_results <- occ_search(
  scientificName = "Haliplus immaculatus",
  country = "GB",
  year = "1970,*")

test_data <- test_results$test_data

# Flagging records with problematic occurrence information using functions of the coordinatecleaner package.
test_clean <- test_data %>%
  filter(!is.na(decimalLatitude),
         !is.na(decimalLongitude)) %>%
  cc_dupl() %>%
  cc_zero() %>%
  cc_equ() %>%
  cc_val() %>%
  cc_sea() %>%
  cc_cap(buffer = 2000) %>%
  cc_cen(buffer = 2000) %>%
  cc_gbif(buffer = 2000) %>%
  cc_inst(buffer = 2000)
print(paste0(nrow(test_data)-nrow(test_clean), " records deleted; ",
             nrow(test_clean), " records remaining."))

# Removing records with coordinate uncertainty and precision issues
test_clean_2 <- test_clean %>%
  filter(is.na(coordinateUncertaintyInMeters) |
           coordinateUncertaintyInMeters < 1000) # removing anything over 1km 

print(paste0(nrow(test_data)-nrow(test_clean_2), " records deleted; ",
             nrow(test_clean_2), " records remaining." ))



################### 5 SPECIES PIPELINE ############

# species list
species_list <- c(
  "Haliplus immaculatus",
  "Notonecta maculata",
  "Plea minutissima",
  "Ischnura pumilio",
  "Somatochlora metallica"
)

# create function to get all GB ocurrences since 1970 for one species
get_species_occ <- function(species_name, country = "GB", year = "1970,*") {
  # get the taxon key from GBIF
  key <- name_backbone(species_name)$usageKey
  
  if (is.null(key)) {
    warning(paste("No taxonKey found for", species_name))
    return(NULL)
  }
  
  # get the data for one
  start <- 0
  limit <- 50000  # maximum per occ_search
  all_data <- list()
  
  repeat {
    res <- occ_search(
      taxonKey = key,
      country = country,
      year = year,
      limit = limit,
      start = start
    )
    
    if (length(res$data) == 0) break
    
    all_data[[length(all_data) + 1]] <- res$data
    
    if (nrow(res$data) < limit) break  # last page
    
    start <- start + limit
  }
  
  # Combine all pages
  do.call(rbind, all_data)
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
# emoves:
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
  cc_dupl(lon = "decimalLongitude", lat = "decimalLatitude") %>%
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

write.csv(clean_2, "gbif_occurrences_cleaned.csv", row.names = FALSE)


  
  
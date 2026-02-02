library(readr)
library(dplyr)
library(here)


# read in gbif and nbn filtered dataframes
gbif_data <- read_csv(here("gbif_final_data.csv"))
nbn_data  <- read_csv(here("nbn_final_data.csv"))

# merge datasets
all_records <- bind_rows(gbif_data, nbn_data)


# remove exact duplicates
final_unique_data <- all_records %>%
  distinct(species, gridReference, eventDate)

# save to csv
write.csv(final_unique_data, "merged_data.csv", row.names = FALSE)

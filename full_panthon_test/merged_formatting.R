library(dplyr)
library(readr)
library(here)
library(tidyr)

# read in the data
data <- read_csv(here("merged_data.csv"))

# generate new dataframe for presences
formatted_presence <- data %>%
  # format date
  mutate(date = as.Date(eventDate)) %>%
  # create the unique visit ID (grid + underscore + date)
  mutate(visit = paste(gridReference, eventDate, sep = "_")) %>%
  # calculate "umRecords (total sightings of this species in the whole of data)
  add_count(visit, name = "numRecords") %>%
  
  # Add "presence" - binary
  mutate(presence = 1) %>%
  
  # select columns
  select(
    visit,
    gridRef = gridReference,
    date,
    species,
    presence,
    numRecords
  ) %>%
  distinct() # stops duplicates

#### step 2: absences of data
final_dataset <- formatted_presence %>% #This expands the data so every Survey (Grid+Date) has a row for every species
  # for every unique grid+date, make sure every species is listed
  # if species isn't there, give 0
  complete(nesting(visit, gridRef, date, numRecords), species, fill = list(presence = 0)) %>%
  
  # get rid of NAs
  group_by(visit) %>%
  mutate(numRecords = max(numRecords, na.rm = TRUE)) %>%
  ungroup()

saveRDS(final_dataset, "final_data.rds", compress = FALSE)
write.csv(final_dataset, "final_data_2.csv", row.names = FALSE)



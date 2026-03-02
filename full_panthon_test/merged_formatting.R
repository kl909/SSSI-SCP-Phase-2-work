library(dplyr)
library(readr)
library(here)
library(tidyr)
library(rnrfa)
library(lubridate)

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

#### sorting out gridrefs and xy coordinates
# remove Northern Ireland entries (only one letter at beginning of GridRef)
final_dataset <- final_dataset %>%
  rename(monad = gridRef) %>%
  filter(grepl("^[A-Z]{2}", monad))

# convert grid reference to xy coordinates
gb_coords <- osg_parse(final_dataset$monad)

final_dataset$x <- gb_coords$easting
final_dataset$y <- gb_coords$northing

# shift to centroid
final_dataset <- final_dataset %>%
  mutate(
    x = x + 500,
    y = y + 500
  )

# final formatting to match Charles'
# get year column from the date column
final_dataset <- final_dataset %>%
  mutate(year = year(date))

# get visitLength from numRecords
final_dataset <- final_dataset %>%
  mutate(visitLength = case_when(
    numRecords == 1           ~ "single",
    numRecords >= 2 & numRecords <= 3 ~ "short",
    numRecords >= 4           ~ "long"
  ))

# reorder columns
final_dataset <- final_dataset %>%
  select(monad, date, year, visit, species, presence, numRecords, visitLength, x, y)

saveRDS(final_dataset, "data/species_data/final_data.rds", compress = FALSE)
write.csv(final_dataset, "final_data_2.csv", row.names = FALSE)



# This script is the final stage to create the species list required by INLA. It has been set up
# to match the dataset used in INLA by Charles Cunningham and needs to be saved as a spatVector

library(dplyr)
library(readr)
library(here)
library(tidyr)
library(rnrfa)
library(lubridate)
library(terra)

# read in the data
data <- read_csv(here("data/species_csv_files/merged_data.csv"))

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

# save as a csv for easy checking
write.csv(final_dataset, here("data/species_csv_files/final_data_2.csv"), row.names = FALSE)

# convert dataframe to a SpatVector for INLA
final_vector <- vect(final_dataset, 
                     geom = c("x", "y"), 
                     crs = "EPSG:27700",
                     keepgeom = TRUE)

# save as a SpatVector
saveRDS(wrap(final_vector), here("data/species_data/final_vector.rds"))





# These are the scripts involved in the set up and running of INLA_SDM.R for wetland species in the UK.


# ------------------------------------ Required files before running scripts (check you have them):

# pantheon_w2_species_list.csv <- data/species_csv_files/ <- this data was pulled from Pantheon and contains wetland species from the filter W2. Currently we only have invertebrate data, this list will need to # be updated eventually to contain wetland vertebrates too. 

# WatercourseLink.shp <- ../../data/os_rivers/data/ <- this is from OS open data https://osdatahub.os.uk/data/downloads/open/OpenRivers and is the river data for the entire UK

# 1km_grid.shp <- data/spatial_data/river_data <- 1km grid I made in QGIS over the entire UK used for creating river rasters


# ------------------------------------ These scripts are for prepping the species data required to run 10_INLA_SDM.R

# 01_GBIF_pull.R: pulls pantheon_w2_species_list.csv species from GBIF and creates gbif_final_data.csv

# 02_NBN_pull.R: pulls pantheon_w2_species_list.csv species from NBN and creates nbn_final_data.csv

# 03_critchlow_pull.R: DON'T NEED TO RUN IF OUTPUT FILE (critchlow_species.csv) EXISTS. This script creates a list of species already analysed in phase 1 from a google drive of files. This output is then used # in the next script to delete species which have already been done.

# 04_gbif_nbn_combined.R: combines the data pulled from GBIF and NBN into a single dataframe merged_data.csv

# 05_merged_formatting.R: the final stage to create the species list required by INLA. It has been set up to match the dataset used in INLA by 
# Charles Cunningham and needs to be saved as a spatVector "data/species_data/final_vector.rds"

# ------------------ These scripts are for creating the 3 river order rasters required to run 10_INLA_SDM.R. Order_10 contains lake data and was 
# created in QGIS.

# 06_step_by_step_river_ordering.R: This extracts river orders 1 - 3 from WatercourseLink.shp and saves them as separate shp files. Order 1 rivers are small streams which aren't fed by any other stream. Order 2 rivers are larger and are created by at least two order 1 rivers feeding into them. Order 3 rivers are the largest order and are created by at least two order 2 rivers feeding into it. Order10 represents lakes.

# 07_order1_rivers.R: This generates the 1km resolution raster of all order 1 rivers with a 250 m buffer around each

# 08_order2_rivers.R: This generates the 1km resolution raster of all order 2 rivers with a 250 m buffer around each

# 09_order3_rivers.R: This generates the 1km resolution raster of all order 3 rivers with a 250 m buffer around each

# -------------------

# 10_INLA_SDM.R: This is the INLA species distribution model which has been adapted from Charles'"11a_INLABRU_ST_SDM_cluster_main.R". Charles' script is a spatio-temporal model, whereas 10_INLA_SDM.R is only spatial so all temporal arguments have been removed. 




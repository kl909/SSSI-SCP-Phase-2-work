library(INLA) 
library(inlabru)
library(sf)
library(terra)
library(stars)
library(tidyverse)
library(BRCmap)
library(ggplot2)
library(grid)
library(gridExtra)
library(here)

# SET PARAMETERS ---------------------------------------

# Organise raw data taxa groups
dataWET <- c( "Invertebrates", "Vertebrates", "Plants" ) # uncertain what this is doing for my data


# List range of years used for species records
range = c(1970,2026)


# Set batch number/species (as.numeric(args[1])) and taxa group (arg[2]), specified array in job script
# This is for viking & when if I were to have a file for every single species separately
# e.g. RScript MyScript.R 104 birds
args <- commandArgs(trailingOnly = TRUE)
batchN <- as.numeric(args[1]) # e.g. will grab "104"
taxaGroup <- args[2] # e.g. will grab "birds"

# Set number of threads (cores)
numThreads = 4

# Estimated range of spatial effect in km (determines mesh)
# mesh = triangles, counts how many triangle corners for distance
estimated_range <- 50 # nearby areas expected to be similar community. Further away = less similar (related to effort)

# Set bitmap type for ggplot2::ggsave to work on cluster
options(bitmapType='cairo') # for plotting

# SET UP INLA -----------------------------------------

# Set number of threads
inla.setOption("num.threads" = numThreads)
# ...and print out
paste(numThreads, "threads used") %>%
  print

# DATA FILES ------------------------------------------

### SPATIAL DATA
# SpatRasters
#-- connW.tif   : woodland connectivity in Amperes - DON'T NEED
#-- coverBF.tif : broadleaf woodland proportion cover of each 1x1km cell - DON'T NEED
#-- coverCF.tif : coniferous woodland proportion cover of each 1x1km cell - DON'T NEED
#-- GDD5.tif    : growing degree days - ???
#-- RAIN.tif    : annual precipitation - NEED
#-- soilM.tif   : soil moisture - NEED
#-- tasCV.tif   : temperature seasonality - NEED
#-- UK_R.tif    : rasterised UK boundary - NEED
#-- WMIN.tif    : winter cold - ???

# Scaling parameters
load("data/spatial_data/scalingParams.RData")

# SpatRasters
for (i in list.files("data/spatial_data/rasters/",
                     pattern =  "\\.tif$")) {
  
  assign(gsub(".tif", "", i),
         rast(paste0("data/spatial_data/rasters/",
                     i))) 
}

# SpatVectors
for (i in list.files("data/spatial_data/vectors/",
                     pattern =  "\\.shp$")) {
  
  assign(gsub(".shp", "", i),
         vect(paste0("data/spatial_data/vectors/",
                     i)))
}

## [KL] my raster projections need to match Charles' unknown CRS
crs(order1_length) <- crs(GDD5_grp)
crs(order2_length) <- crs(GDD5_grp)
crs(order3_length) <- crs(GDD5_grp)
crs(order10_area) <- crs(GDD5_grp)

order1_length <- terra::resample(order1_length, GDD5_grp, method = "bilinear")
order2_length <- terra::resample(order2_length, GDD5_grp, method = "bilinear")
order3_length <- terra::resample(order3_length, GDD5_grp, method = "bilinear")
order10_area <- terra::resample(order10_area, GDD5_grp, method = "bilinear")

### Download BNG WKT string
# N.B. Individual filename needed for each task to prevent
# different tasks tring to read/write at same time 

# Specify temporary file to download 'bng' CRS wkt to
tempFile <- paste0(taxaGroup, batchN, "bng.prj")

# Download file
download.file(url = "https://epsg.io/27700.wkt2?download=1",
              destfile = tempFile)

# Assign wkt string
bng <- sf::st_crs(tempFile)$wkt

# Remove temporary file
unlink(tempFile)

# PROCESS COVARIATES -----------------------------------

# Load example species if using
visitDataSpatial <- readRDS("data/species_data/final_vector.rds") # MY FILTERED SPECIES DATA
# remember Charles loads in one species at a time, whereas this file is for all species

# CREATE WEEK COVARIATE

### Add week of year column ( will be included in model as covariate)
# N.B. Week of the year as decimal number (01--53) as defined in ISO 8601.
# If the week (starting on Monday) containing 1 January has four or more 
# days in the new year, then it is considered week 1. Otherwise, 
# it is the last week of the previous year, and the next week is week 1.
visitDataSpatial$week <- visitDataSpatial$date %>%
  strftime(., format = "%V") %>%
  as.numeric(.)

# CREATE EFFORT COVARIATE

# Convert factor levels to dummy variables
visitDataSpatial <- visitDataSpatial %>%
  model.matrix(object = ~visitLength) %>%
  as.data.frame() %>%
  dplyr::select(-1) %>%
  cbind(visitDataSpatial, .)

# EXTRACT COVARIATES FOR EFFECTS PLOT -------------------------------
# Need to extract spatial covariates over species records for plots

### Create data frame of covariate values at visit locations

# Create a base data frame with iYear, presence and week (non-spatial) to build up from
covarValues <- dplyr::select(as.data.frame(visitDataSpatial), presence, week)


## [KL] extract xy from visitDataSpatial so there is just a 2 column dataframe for the following loop
#coords_matrix <- as.matrix(visitDataSpatial[, c("x", "y")]) 

# Loop through spatial (random) variables
for (i in c( "GDD5_grp", "WMIN_grp", "tasCV_grp", "RAIN_grp", "soilM_grp",
             "order1_length", "order2_length", "order3_length", "order10_area")) {
  # LENGTH OF ORDER 1, LENGTH OF ORDER 2, LENGTH OF ORDER 3, AREA OF STANDING WATER
  
  # Get covariate i spatRaster (each layer is for time period iYear)
  cov_R <- get(i)
  
  # Extract values for all iYear layers and add to beginning of data frame using species records
  #covarValues <- terra::extract(cov_R, visitDataSpatial[, c("x", "y")]) %>% # [KL] I removed ,ID = FALSE from end
    #cbind(covarValues) # Add existing data frame onto end
  
  # Drop now redundant columns from first columns
  #covarValues <- covarValues[, -1] # [KL] changed from covarValues <- covarValues[, -c(1,2)]
  
  pts <- terra::vect(visitDataSpatial[, c("x", "y")], 
                     geom = c("x", "y"), 
                     crs = bng)
  
  # 2. Extract (Using the spatial object 'pts' is much safer than a matrix)
  # This returns a data frame with [ID, values]
  temp_extract <- terra::extract(cov_R, pts)
  
  # 3. Create a clean 1-column data frame with the correct name
  # This stops the "NA.1", "NA.2" naming issue
  new_col <- data.frame(temp_extract[, 2])
  colnames(new_col) <- i
  
  # 4. Bind it to your master table
  covarValues <- cbind(new_col, covarValues)
  
  
}

### Process climate covariate values for rug plot

# Remove cover and climate variables and convert to long format
climCovarValues <-  dplyr::select(covarValues,
                                  -c(coverBF, coverCF, connW)) %>%
  gather(randomEff, value, 3:NCOL(.))

# Count number of each quantile value for each variable for presence/absence separately
climCovarValues <- climCovarValues %>%
  group_by(randomEff, presence, value) %>%
  mutate(count = n()) %>%
  ungroup %>%
  distinct

# TIDY MEMORY BEFORE MODEL RUN --------------------------------------

# Remove data frames no longer needed 
rm(visitXY, rawDataUK, taxaData, visitData)

# Garbage clean
gc()

# CREATE MESH -------------------------------------------------------

# Max edge is as a rule of thumb (range/3 to range/10)
maxEdge <- estimated_range/5

# Find record locations to build mesh from
recordCoords <- crds(visitDataSpatial) %>% 
  unique(.)

# Create mesh (convert boundary to sp object as leads to best convergence)
mesh <- inla.mesh.2d(boundary = st_as_sf(smoothUK) %>% as("Spatial"),
                     loc = recordCoords,
                     max.edge = c(1,5) * maxEdge,
                     offset = c(1,2) * maxEdge, 
                     cutoff = maxEdge/2,
                     min.angle = 26,
                     crs = gsub( "units=m", "units=km", st_crs(bng)$proj4string ))

# FIT SPATIO-TEMPORAL MODEL ---------------------------------

# Create indices -- NO TEMPORAL DIMENSION FOR US - MAKE SURE ONLY ONE YEAR
iYear <- visitDataSpatial$iYear
nYear <- length(unique(iYear))

# Define spatial SPDE priors
mySpace <- inla.spde2.pcmatern(
  mesh,
  prior.range = c(1 * maxEdge, 0.5),
  prior.sigma = c(1, 0.5))

# Priors for fixed effects
fixedHyper <- list( mean = 0,
                    prec = 1 ) # Precision for all fixed effects except intercept

# Priors for random effects
randomHyper <- list(theta = list(prior="pc.prec",
                                 param=c(0.5, 0.01)))

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

# change resolution of river rasters using climate rasters
ext(order1_length)
ext(order2_length)
ext(order3_length)
ext(order10_area)
ext(GDD5_grp)

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

# Load species data
master_data <- readRDS("data/species_data/final_vector.rds") # MY FILTERED SPECIES DATA
# remember Charles loads in one species at a time, whereas this file is for all species

# [KL] identify which species (mimicking Charles' 'batchN')
all_species_list <- unique(master_data$species)

# [KL] loop over covariate script for each species
for (target_species in all_species_list) {
  
  temp_df <- as.data.frame(master_data)
  
  visitDataSpatial <- temp_df %>%
    filter(species == target_species) %>% 
    group_by(visit) %>% 
    # If a visit has both a 1 and a 0 for some reason, keep the 1
    slice(which.max(presence)) %>% 
    ungroup()
   
  # generate week covariate 
  visitDataSpatial$week <- visitDataSpatial$date %>%
    strftime(., format = "%V") %>%
    as.numeric(.)
  
  # convert back to spatvector
  visitDataSpatial_Vect <- vect(visitDataSpatial, 
                                geom = c("x", "y"), 
                                crs = crs(master_data))
  
  covarValues <- visitDataSpatial
  
  clim_vars <- c("GDD5_grp", "WMIN_grp", "tasCV_grp", "RAIN_grp", "soilM_grp",
                 "order1_length", "order2_length", "order3_length", "order10_area"
  )
  
  for (i in clim_vars) {
    
    cov_R <- get(i)
    
    # Extract values from the SpatRaster (cov_R) at the species locations
    temp_ext <- terra::extract(cov_R, visitDataSpatial_Vect)
    
    # Ensure we take the column that matches the data, not the ID column
    new_data <- temp_ext[, 2, drop = FALSE]
    colnames(new_data) <- i
    
    covarValues <- cbind(new_data, covarValues)
  }
  
  # tidy up columns
  covarValues <- covarValues %>%
    dplyr::select(presence,
                  week,
                  all_of(clim_vars) 
                  )
  
  # Process climate covariate values for rug plot
  climCovarValues <- covarValues %>%
    dplyr::select(any_of(names(.)), -any_of(-c(order1_length, order2_length, order3_length, order10_area))) %>%
    gather(randomEff, value, 3:NCOL(.))
  
  # count number of each quantile value
  climCovarValues <- climCovarValues %>%
    group_by(randomEff, presence, value) %>%
    mutate(count = n()) %>%
    ungroup() %>%
    distinct()
  
  # save results
  tidy_name <- gsub(" ", "_", target_species)
  
  # save wide data
  saveRDS(covarValues, paste0("data/output/covars_", tidy_name, ".rds"))
  
  # save the long summary data
  saveRDS(climCovarValues, paste0("data/output/clim_summary_", tidy_name, ".rds"))
}

-----------------------------------------------------------------------

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

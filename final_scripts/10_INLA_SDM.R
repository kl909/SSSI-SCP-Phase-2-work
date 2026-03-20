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
#dataWET <- c( "Invertebrates", "Vertebrates", "Plants" ) # uncertain what this is doing for my data


# List range of years used for species records
#range = c(1970,2026)


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
# SpatRasters Charles made
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

# scale and match projection to RAIN raster
order1_length <- project(order1_length, RAIN)
order1_length <- scale(order1_length)
order2_length <- project(order2_length, RAIN)
order2_length <- scale(order2_length)
order3_length <- project(order3_length, RAIN)
order3_length <- scale(order3_length)
order10_area <- project(order10_area, RAIN)
order10_area <- scale(order10_area)

#-----------------------------------
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


# SET UP FOLDER FOR PLOTS -----------------------------
dir.create("data/output/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("data/output/models", recursive = TRUE, showWarnings = FALSE)

# PROCESS COVARIATES -----------------------------------

# Load species data
master_data <- readRDS("data/species_data/final_vector.rds") # MY FILTERED SPECIES DATA
# remember Charles loads in one species at a time, whereas this file is for all species

# [KL] identify which species (mimicking Charles' 'batchN')
all_species_list <- unique(master_data$species)

# --- TEST SETTINGS ---- 
test_mode <- FALSE  # Set to FALSE when you are ready for the full species run
if (test_mode) {
  species_to_run <- all_species_list[1:3] # Just run the first 3 species
} else {
  species_to_run <- all_species_list      # Run all species
}

# [KL] loop over covariate script for each species
for (target_species in species_to_run) {
  
  # skip species if already done
  tidy_name <- gsub(" ", "_", target_species)
  covar_file <- paste0("data/output/covars/covars_", tidy_name, ".rds")
  summary_file <- paste0("data/output/clim_sums/clim_summary_", tidy_name, ".rds")
  
  # If BOTH files already exist, skip to the next species
  if (file.exists(covar_file) & file.exists(summary_file)) {
    message("Skipping covariate extraction for: ", target_species, " (Already exists)")
    next
  }
  
  # Rest of loop
  
  temp_df <- as.data.frame(master_data)
  
  visitDataSpatial_df <- temp_df %>%
    filter(species == target_species) %>% 
    group_by(visit) %>% 
    # If a visit has both a 1 and a 0 for some reason, keep the 1
    slice(which.max(presence)) %>% 
    ungroup()
   
  # generate week covariate 
  visitDataSpatial_df$week <- visitDataSpatial_df$date %>%
    strftime(., format = "%V") %>%
    as.numeric(.)
  
  # convert back to spatvector
  visitDataSpatial <- vect(visitDataSpatial_df, 
                                geom = c("x", "y"), 
                                crs = crs(master_data))
  
  # CREATE EFFORT COVARIATE
  new_cols <- visitDataSpatial %>%
    as.data.frame() %>%                          # Convert to DF for math
    model.matrix(object = ~visitLength) %>%      # Create dummies
    as.data.frame() %>%
    dplyr::select(-1)
  
  # Safely put the new columns back into the spatial object
  visitDataSpatial$visitLengthShort  <- new_cols[,1]
  visitDataSpatial$visitLengthSingle <- new_cols[,2]
  
  # EXTRACT COVARIATES FOR EFFECTS PLOT -------------------------------
  
  covarValues <- visitDataSpatial_df
  
  clim_vars <- c("GDD5_grp", "WMIN_grp", "tasCV_grp", "RAIN_grp", "soilM_grp",
                 "order1_length", "order2_length", "order3_length", "order10_area"
  )
  
  for (i in clim_vars) {
    
    cov_R <- get(i)
    
    # Extract values from the SpatRaster (cov_R) at the species locations
    temp_ext <- terra::extract(cov_R, visitDataSpatial)
    
    # Ensure we take the column that matches the data, not the ID column
    new_data <- temp_ext[, 2, drop = FALSE]
    colnames(new_data) <- i
    
    covarValues <- cbind(new_data, covarValues)
  }
  
  # Remove rows where GDD5_grp is NA (outside the raster)
  #covarValues <- covarValues %>%
   # filter(!is.na(GDD5_grp))
  
  # Remove rows where GDD5_grp is NA (outside the raster)
  # This uses GDD5 as the indicator for 'offshore' or 'missing' data
  keep_idx <- which(!is.na(covarValues$GDD5_grp))
  
  # Apply to the data frame
  covarValues <- covarValues[keep_idx, ]
  
  # Apply to the SpatVector (This stops the inlabru warning!)
  visitDataSpatial <- visitDataSpatial[keep_idx, ]
  
  # Convert NAs in the river ordr columns to 0
  covarValues <- covarValues %>%
    mutate(across(starts_with("order"), ~replace_na(., 0)))
  
  ##### tidy up columns
  covarValues <- covarValues %>%
    dplyr::select(presence,
                  week,
                  all_of(clim_vars) 
                  )
  
  # Process climate covariate values for rug plot
  clim_groups <- c("GDD5_grp", "WMIN_grp", "tasCV_grp", "RAIN_grp", "soilM_grp")
  
  climCovarValues <- covarValues %>%
    dplyr::select(presence, week, all_of(clim_groups)) %>%
    gather(randomEff, value, all_of(clim_groups))
  
  # count number of each quantile value
  climCovarValues <- climCovarValues %>%
    group_by(randomEff, presence, value) %>%
    mutate(count = n()) %>%
    ungroup() %>%
    distinct()
  
  # save results
  tidy_name <- gsub(" ", "_", target_species)
  
  # save wide data
  saveRDS(covarValues, paste0("data/output/covars/covars_", tidy_name, ".rds"))
  
  # save the long summary data
  saveRDS(climCovarValues, paste0("data/output/clim_sums/clim_summary_", tidy_name, ".rds"))
  
  # save the spatial data
  saveRDS(visitDataSpatial, paste0("data/output/spatial_objs/sp_", tidy_name, ".rds"))
}


# TIDY MEMORY BEFORE MODEL RUN --------------------------------------

# Remove data frames no longer needed 
#rm(new_data, visitDataSpatial_df)

# Garbage clean
gc()

##### check data and load in covarValues and climCoverValues for loop
for (target_species in all_species_list[1:2]) { # remove [1:2] to run script on all species, or just delete [1:2] if using the test settings from
  # line 125. This loop will be very memory intensive so don't run on many species without Viking
  
  tidy_name <- gsub(" ", "_", target_species)
  model_file <- paste0("data/output/models/model_", tidy_name, ".rds")
  
  # CHECKPOINT: Skip if this species has already been done in a previous run
  if (file.exists(model_file)) {
    message(">>> Skipping ", target_species, " (Model already exists)")
    next 
  }
  
  message(">>> Processing Model for: ", target_species)
  
  # LOAD THE OUTPUTS from your first loop
  # We use try() here in case a specific RDS file is missing or corrupt
  covarValues <- try(readRDS(paste0("data/output/covars/covars_", tidy_name, ".rds")))
  climCovarValues <- try(readRDS(paste0("data/output/clim_sums/clim_summary_", tidy_name, ".rds")))
  visitDataSpatial <- try(readRDS(paste0("data/output/spatial_objs/sp_", tidy_name, ".rds")))
  
  if (inherits(covarValues, "try-error")) {
    message("!!! Could not find data for ", target_species, ". Skipping.")
    next
  }
  
  try({
    
    

  # CREATE MESH -------------------------------------------------------
  
  # [KL] create simplifieed boundary of UK:
  uk_boundary_sf <- smoothUK %>%
    aggregate(fact = 10) %>% 
    as.polygons() %>%
    st_as_sf() %>%
    st_union()
  
  # [KL] convert to INLA format
  uk_boundary_sp <- as(uk_boundary_sf, "Spatial")
  
  # [KL] define km based BNG projection once
  crs_km <- st_crs(27700)$proj4string %>% 
    gsub("units=m", "units=km", .)
  
  # Max edge is as a rule of thumb (range/3 to range/10)
  maxEdge <- estimated_range/5
  
  # Find record locations to build mesh from
  recordCoords <- crds(visitDataSpatial) %>% 
    unique(.)
  
  # [KL] convert recordCoords to km
  recordCoords_km <- recordCoords / 1000
  
  # Create mesh (convert boundary to sp object as leads to best convergence)
  mesh <- inla.mesh.2d(boundary = uk_boundary_sp,
                       loc = recordCoords_km,
                       max.edge = c(1,5) * maxEdge,
                       offset = c(1,2) * maxEdge, 
                       cutoff = maxEdge/2,
                       min.angle = 26,
                       crs = crs_km)
  
  
  # FIT SPATIO-TEMPORAL MODEL ---------------------------------
  
  # Create indices -- NO TEMPORAL DIMENSION FOR US - MAKE SURE ONLY ONE YEAR
  
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
  
  # Set components
  inlabruCmp  <-  presence ~ 0 + Intercept(1) +
    
    GDD5(main = GDD5_grp,
         #main_layer = iYear,
         model = "rw2",
         scale.model = TRUE,
         hyper = randomHyper) +
    WMIN(main = WMIN_grp,
         #main_layer = iYear,
         model = "rw2",
         scale.model = TRUE,
         hyper = randomHyper) +
    tasCV(main = tasCV_grp,
          #main_layer = iYear,
          model = "rw2",
          scale.model = TRUE,
          hyper = randomHyper) +
    RAIN(main = RAIN_grp,
         #main_layer = iYear,
         model = "rw2",
         scale.model = TRUE,
         hyper = randomHyper) +
    soilM(main = soilM_grp,
          #main_layer = iYear,
          model = "rw2",
          scale.model = TRUE,
          hyper = randomHyper) +
    order1(main = order1_length,
            #main_layer = iYear,
            model = "linear") +
    order2(main = order2_length,
            #main_layer = iYear,
            model = "linear") +
    order3(main = order3_length,
                 #main_layer = iYear,
                 model = "linear") +
    order10(main = order10_area,
              #main_layer = iYear,
              model = "linear") +
    
    visitLengthSingle(main = visitLengthSingle, model = "linear") +
    visitLengthShort(main = visitLengthShort, model = "linear") +
    
    week(main = week,
         model = "rw2",
         cyclic = TRUE,
         hyper = randomHyper) +
    spaceTime(main = geometry, # this only space (no Time)
              #group = iYear,
              #ngroup = nYear,
              model = mySpace)
              #control.group = list(model = "ar1"), # don't need as only using space not time
                                   #hyper = ar1Hyper))
  
  # Fit model
  model <- bru(components = inlabruCmp,
               family = "binomial",
               control.family = list(link = "cloglog"),
               data = st_as_sf(visitDataSpatial),
               options=list(control.fixed = fixedHyper,
                            control.inla= list(int.strategy='eb'),
                            control.compute = list(waic = TRUE, dic = FALSE, cpo = TRUE),
                            verbose = TRUE))
  
  # Assign model summary object and output
  modelSummary <- summary(model)
  
  # save model
  saveRDS(model, paste0("data/output/models/model_", tidy_name, ".rds"))
  
  # PREDICT -----------------------------------------
  # Create grid prediction pixels
  ppxl <- mask(UK_R, smoothUK) %>%
    crop(.,smoothUK ) %>%
    as.points %>%
    st_as_sf
  
  
  
  # Predict using spatio model
  # ( N.B. excluding VISIT_LENGTH means including the reference factor level- long- which is what we want!)
  modelPred <- predict(model, 
                       ppxl, 
                       ~ data.frame(
                         lambda =  1 - exp( -exp( spaceTime + # cloglog back transform
                                                    soilM +
                                                    WMIN +
                                                    tasCV +
                                                    GDD5 +
                                                    RAIN +
                                                    order1 +
                                                    order2 +
                                                    order3 +
                                                    order10 +
                                                    # Max value for week to predict over (removed later)
                                                    max(model$summary.random$week$mean) + 
                                                    Intercept ))),
                       exclude = c("week")) 
  
  # MODEL EVALUATION --------------------------
  
  # SET PARAMETERS
  
  # Labels
  randomEffLabels <- c('GDD5' = "Growing degree days", 
                       'RAIN' = "Annual precipiation",
                       'soilM' = "Soil moisture", 
                       'tasCV' = "Temperature seasonality",
                       'week' = "Week of year" , 
                       'WMIN' = "Winter minimum temperature"
                       )
  
  linearEffLabels <- c('order1' = "order 1 river length",
                       'order2' = "order 2 river length",
                       'order3' = "order 3 river length",
                       'order10' = "order 10 lake area",
                       'visitLengthSingle' = "Single record visit",
                       'visitLengthShort' = "Short visit (2-3 records)")
  
  
  # Template raster for converting from sf to terra raster objects
  template_R <- st_as_stars(UK_R)
  template_R[[1]][1:ncell(template_R)] <- NA
  
  # CALCULATE LOGCPO
  
  logCPO_vect = log(model$cpo$cpo[model$cpo$cpo != 0])
  logCPO_vect = logCPO_vect[is.finite(logCPO_vect)]
  logCPO = round(-sum(logCPO_vect, na.rm = T), digits = 2)
  
  # NON-SPATIAL RANDOM EFFECTS PLOT
  ### Create effects data frame
  
  # Extract random effects from model, and exclude spatial
  randomEff_df <- model$summary.random
  randomEff_df["spaceTime"] <- NULL
  
  # Add name of random effect to each dataframe in list
  randomEff_df <- imap(randomEff_df, ~mutate(.x, randomEff = .y))
  
  # Unlist, then rename and select quantile columns
  randomEff_df <- do.call(rbind, randomEff_df)%>%
    rename("q0.025" = "0.025quant",
           "q0.5" = "0.5quant",
           "q0.975" = "0.975quant")
  
  ### Back scale non-spatial random covariate values
  
  # Apply 'unscaling' function to every row of effects data frame
  randomEff_df$unscaledID <- sapply(1:NROW(randomEff_df), function(x) {
    
    # If week, just use ID as not scaled
    if (randomEff_df$randomEff[x] == "week") {
      
      return(randomEff_df$ID[x])
      
    } else { # Otherwise, 'unscale'!
      
      # Extract covariate mean and sd for scaling function
      randomEffMean <-
        scalingParams[scalingParams$variable == randomEff_df$randomEff[x],
                      "variableMean"]
      randomEffSD <-
        scalingParams[scalingParams$variable == randomEff_df$randomEff[x],
                      "variableSD"]
      
      # Unscale using xSCALED = (x - xbar)/sd --> x = (xSCALED * sd) + xbar principle
      unscaledID <- (randomEff_df$ID[x] * randomEffSD) + randomEffMean
      
      return(unscaledID) # Return unscaled value
      
    }})
  
  # Apply 'unscaling' function to every row of record locations
  climCovarValues$unscaledValue <- sapply(1:NROW(climCovarValues), function(x) {
    
    # If week, just use value as not scaled
    if (climCovarValues$randomEff[x] == "week") {
      
      return(climCovarValues$value[x])
      
    } else { # Otherwise, 'unscale'!
      
      # Extract covariate mean and sd for scaling function
      randomEffMean <- scalingParams[scalingParams$variable == climCovarValues$randomEff[x],
                                     "variableMean"]
      randomEffSD <- scalingParams[scalingParams$variable == climCovarValues$randomEff[x],
                                   "variableSD"]
      
      # Unscale using xSCALED = (x - xbar)/sd --> x = (xSCALED * sd) + xbar principle
      unscaledValue <- (climCovarValues$value[x] * randomEffSD) + randomEffMean
      
      return(unscaledValue)
      
    }})
  
  ### Plot
  
  randomEffPlot <- ggplot(randomEff_df) +
    
    # Random effect size
    geom_line(aes(x = unscaledID, y = q0.5)) +
    geom_line(aes(x = unscaledID, y = q0.025), lty = 2, alpha = .5) +
    geom_line(aes(x = unscaledID, y = q0.975), lty = 2, alpha = .5) +
    
    facet_wrap(~ randomEff, scale = 'free_x', labeller = as_labeller(randomEffLabels)) +
    ggtitle("Non-linear random effects") + theme_minimal()
  
  # display plot
  randomEffPlot
  
  # save plot
  ggsave(filename = paste0("data/output/plots/random_eff_", tidy_name, ".png"),
         plot = randomEffPlot, width = 10, height = 8, bg = "white")
  
  # FIXED EFFECTS PLOT
  # Loop through covariates and extract estimates
  fixed_effects_table <- modelSummary$inla$fixed
  
  # 2. Verify it's not NULL anymore (this should finally work!)
  #print(head(fixed_effects_table))
  
  # 3. Run the loop
  effectSizeAll <- data.frame()
  
  for (i in names(linearEffLabels)) {
    if (i %in% rownames(fixed_effects_table)) {
      effectSize <- as.data.frame(fixed_effects_table[i, , drop = FALSE])
      effectSize$Covariate <- i
      effectSizeAll <- rbind(effectSizeAll, effectSize)
    }
  }
  
  #rownames(effectSizeAll) <- NULL
  
  # Plot fixed effects
  fixedEffPlot <- ggplot(effectSizeAll, 
                         aes(y = `0.5quant`, x = Covariate,
                             ymin = `0.025quant`, ymax=`0.975quant`, 
                             col = Covariate, fill = Covariate)) + 
    #specify position here
    geom_linerange(linewidth=4, colour = "lightblue") +
    ggtitle("Linear effects") +
    geom_hline(yintercept=0, lty=2) +
    geom_point(size=2, shape=21, colour="white", fill = "black", stroke = 0.1) +
    scale_x_discrete(name="",
                     limits = rev(names(linearEffLabels)),
                     labels = as_labeller(linearEffLabels)) +
    scale_y_continuous(name="Effect size") +
    coord_flip() +
    theme_minimal() + 
    guides(colour = "none") +
    theme(axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          legend.text = element_text(size = 16),
          plot.title = element_text(hjust = 0.5, vjust = -0.5))
  
  # display plot
  fixedEffPlot
  
  ggsave(filename = paste0("data/output/plots/fixed_eff_", tidy_name, ".png"),
         plot = fixedEffPlot, width = 8, height = 6, bg = "white")
  
  ### [KL] adding in prediction map plot - this was just for Colin 
  
  # prob_raster <- st_rasterize(modelPred["mean"], template_R)
  # prob_raster_terra <- rast(prob_raster)
  # 
  # # Plot
  # ggplot() +
  #   geom_stars(data = prob_raster) +
  #   scale_fill_viridis_c(option = "magma", name = "Prob. of Presence") +
  #   coord_equal() +
  #   theme_minimal() +
  #   labs(title = paste("Predicted Distribution:", target_species))
  
  
  #####################################################
  
  # SPDE PARAMETER POSTERIORS
  
  # Extract range and variance of space-time SPDE, and plot
  range.plot <- plot( spde.posterior(model, "spaceTime", what = "range")) +
    ggtitle("SPDE range") +
    theme(plot.title = element_text(hjust = 0.5))
  var.plot <- plot(spde.posterior(model, "spaceTime", what = "log.variance")) +
    ggtitle("SPDE log variance") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ########################
  # MATERN CORRELATION AND COVARIANCE
  
  corplot <- plot(spde.posterior(model, "spaceTime", what = "matern.correlation")) +
    ggtitle("Matern correlation") +
    theme(plot.title = element_text(hjust = 0.5))
  covplot <- plot(spde.posterior(model, "spaceTime", what = "matern.covariance")) +
    ggtitle("Matern covariance") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste0("data/output/plots/cor_plot_", tidy_name, ".png"),
         plot = corplot, width = 8, height = 6, bg = "white")
  ggsave(filename = paste0("data/output/plots/cov_plot_", tidy_name, ".png"),
         plot = covplot, width = 8, height = 6, bg = "white")
  
  # MEDIAN AND SD PREDICTION PLOTS
  ### Plot median
  
  # Convert to medians to spatRast for saving/plotting
  median_R <- st_rasterize(modelPred[, "median"],
                             template = template_R,
                             options = c("a_nodata = NA")) %>%
    rast
  
  # Add names, i.e. iYear
  names(median_R) <- "Spatial_Prediction" # chage the name of this to "median"
  
  # Convert to data frame for plotting (single layer now)
  median_df <- as.data.frame(median_R, xy = TRUE) 
  colnames(median_df)[3] <- "median" # Ensure column is named 'median'
  
  # Plot posterior median
  predMedian <- ggplot(data = median_df) +
    ggtitle("Median posterior occupancy") +
    coord_fixed() +
    geom_tile(aes(x=x, y=y, fill = median, colour = median)) +
    scale_fill_distiller(palette = "BuGn",
                         direction = 1,
                         limits = c(0,1),
                         guide = guide_colourbar(title = "Occupancy\nprobability")) +
    scale_colour_distiller(palette = "BuGn",
                           direction = 1,
                           limits = c(0,1),
                           guide = "none") +
    # Removed facet_wrap and geom_text(timeLabels) as there is only one map
    theme_void() + 
    theme(plot.title = element_text(hjust = 0.5, vjust = -1)) +
    geom_sf(data = st_as_sf(smoothUK), fill = NA, colour = "black")
  
  ggsave(filename = paste0("data/output/plots/pred_median_", tidy_name, ".png"),
         plot = predMedian, width = 8, height = 6, bg = "white")
  
  ### Plot posterior sd
  # Convert SD to spatRast (single layer)
  sd_R <- st_rasterize(modelPred[, "sd"],
                       template = template_R,
                       options = c("a_nodata = NA")) %>%
    rast()
  
  # Convert to data frame for plotting
  names(sd_R) <- "sd"
  sd_df <- as.data.frame(sd_R, xy = TRUE)
  
  
  predSD <- ggplot(data = sd_df) +
    ggtitle("Posterior standard deviation") +
    geom_tile(aes(x=x, y=y, fill = sd, colour = sd)) +
    scale_fill_distiller(palette = "BuGn",
                         direction = 1) +
    # Removed facet_wrap and geom_text
    theme_void() + 
    #theme(plot.title = element_text(hjust = 0.5, vjust = -1)) +
    coord_fixed() +
    geom_sf(data = st_as_sf(smoothUK), fill = NA, colour = "black")
  
  ggsave(filename = paste0("data/output/plots/pred_SD_", tidy_name, ".png"),
         plot = predSD, width = 8, height = 6, bg = "white")
  
  # Predict spatial field only (no time groups)
  spaceTimePred <- predict(model, 
                           ppxl, # Predicts for just the grid once
                           ~ data.frame(effectSize = spaceTime), # Looks for 'field'
                           include = c("spaceTime"))
  
  # Convert to a single-layer raster
  spaceTime_R <- st_rasterize(spaceTimePred[, "median"],
                              template = template_R,
                              options = c("a_nodata = NA")) %>%
    rast()
  
  
  # Convert to data frame for plotting
  spaceTime_df <- as.data.frame(spaceTime_R, xy = TRUE) # Convert to data frame
  colnames(spaceTime_df)[3] <- "median"
  
  # Plot - THIS IS ONLY SPATIAL
  spaceTimePlot <- ggplot(data = spaceTime_df) +
    ggtitle("Spatio-temporal field")  +
    geom_tile(aes(x=x, y=y, fill = median, colour = median)) +
    scale_fill_distiller(palette = "RdYlBu",
                         direction = 1,
                         limits = c(-1,1) * max(abs(spaceTime_df$median))) +
    theme_void() + 
    theme(plot.title = element_text(hjust = 0.5, vjust = 1)) +
    coord_fixed() +
    geom_sf(data = st_as_sf(smoothUK), fill = NA, colour = "black", inherit.aes = FALSE)
  
  ggsave(filename = paste0("data/output/plots/space_time_", tidy_name, ".png"),
         plot = spaceTimePlot, width = 8, height = 6, bg = "white")


  })
}

###################### [KL] works up to here
###################################################################################
# # below is all for treescape connectivity - change for river connectivity?
# 
# # COVER-CONNECTIVITY INTERACTION PLOTS
# # COVER-CONNECTIVITY INTERACTION PLOTS
# 
# ### Set up prediction data frames
# 
# # How many prediction steps?
# nSamp <- 100
# 
# # Extract max scaled value
# maxConnectivity <- global(connW, fun = "max", na.rm = TRUE) %>% 
#   max
# 
# # Create unscaled data frame of cover and connectivity values to predict over
# # Separate broadleaf and coniferous data frames
# BF_pred_df <-  expand.grid(BF_pred = seq(0, 1, by = 1/nSamp),
#                            conn_pred = seq(0, maxConnectivity, by = maxConnectivity/nSamp))
# CF_pred_df <- expand.grid(CF_pred = seq(0, 1, by = 1/nSamp),
#                           conn_pred = seq(0, maxConnectivity, by = maxConnectivity/nSamp))
# 
# ### Scale covariates
# 
# # Scale the prediction steps for broadleaf and coniferous woodland separately, and connectivity
# # N.B. Have to name columns the same as the original datasets!
# BF_pred_df$coverBF_scaled <- ( BF_pred_df$BF_pred - 
#                                  scalingParams[scalingParams$variable == "coverBF", "variableMean"] ) /
#   scalingParams[scalingParams$variable == "coverBF", "variableSD"]
# 
# CF_pred_df$coverCF_scaled <- ( CF_pred_df$CF_pred - 
#                                  scalingParams[scalingParams$variable == "coverCF", "variableMean"] ) /
#   scalingParams[scalingParams$variable == "coverCF", "variableSD"]
# 
# BF_pred_df$connW_scaled <- CF_pred_df$connW_scaled <- # N.B. Connectivity is the same for both cover types
#   (BF_pred_df$conn_pred - scalingParams[scalingParams$variable == "connW", "variableMean"]) /
#   scalingParams[scalingParams$variable == "connW", "variableSD"]
# 
# # Calculate scaled interaction terms for prediction
# BF_pred_df$coverBF_connW <- BF_pred_df$coverBF_scaled * BF_pred_df$connW_scaled
# CF_pred_df$coverCF_connW <- CF_pred_df$coverCF_scaled * CF_pred_df$connW_scaled
# 
# ### Create 'unscaled' vectors of covariate values where species is present for plot
# 
# # Subset covarValues to presence records, and then unscale for 
# # broadleaf and coniferous woodland, and connectivity
# presentCoverBF <- subset(covarValues, presence == 1)$coverBF
# presentCoverBF <-
#   ((presentCoverBF * scalingParams[scalingParams$variable == "coverBF", "variableSD"]) +
#      scalingParams[scalingParams$variable == "coverBF", "variableMean"])
# 
# presentCoverCF <- subset(covarValues, presence == 1)$coverCF
# presentCoverCF <-
#   ((presentCoverCF * scalingParams[scalingParams$variable == "coverCF", "variableSD"]) +
#      scalingParams[scalingParams$variable == "coverCF", "variableMean"])
# 
# presentConnW <- subset(covarValues, presence == 1)$connW
# presentConnW <-
#   ((presentConnW * scalingParams[scalingParams$variable == "connW", "variableSD"]) +
#      scalingParams[scalingParams$variable == "connW", "variableMean"])
# 
# # Join unscaled covariate values together, and take unique values to speed up later steps
# presentFixedEff <- data.frame(presentCoverBF, presentCoverCF, presentConnW) %>%
#   distinct
# 
# # Remove obsolete objects
# rm(presentCoverBF, presentCoverCF, presentConnW)
# 
# ### Predict
# 
# # Predict broadleaf cover and connectivity interaction at link scale
# BFconnINTpred <- predict(model,
#                          BF_pred_df,
#                          formula = ~ coverBF_eval(coverBF_scaled) +
#                            connectivity_eval(connW_scaled) +
#                            BFconnINT_eval(coverBF_connW),
#                          exclude = c("spaceTime", "week",
#                                      "soilM",  "WMIN", "tasCV", "GDD5", "RAIN",
#                                      "coverCF", "CFconnINT"))
# 
# # Predict coniferous cover and connectivity interaction at link scale
# CFconnINTpred <- predict(model,
#                          CF_pred_df,
#                          formula = ~ coverCF_eval(coverCF_scaled) +
#                            connectivity_eval(connW_scaled) +
#                            CFconnINT_eval(coverCF_connW),
#                          exclude = c("spaceTime", "week",
#                                      "soilM",  "WMIN", "tasCV", "GDD5", "RAIN",
#                                      "coverBF", "BFconnINT"))
# 
# # Rename prediction columns needed for plot(median)
# BFconnINTpred <- rename(BFconnINTpred, BFmedian = median)
# CFconnINTpred <- rename(CFconnINTpred, CFmedian = median)
# 
# # Join into a single dataframe
# allPred <- cbind(BFconnINTpred[, c("BF_pred", "conn_pred", "BFmedian" )], 
#                  CFconnINTpred[, c("CF_pred", "CFmedian" )])
# 
# ### Filter predictions by cover and connectivity values which species is present in
# 
# # Find prediction values (discrete) which species fixed effect values overlap with
# # Using both broadleaf and connectivity columns of unscaled covariate values where the species is present,
# # find if there are any of these values that is within the prediction bin (i.e. x +- max/nSamp)
# # for the respective covariate
# allPred$BFpresent <- mapply( x = allPred$BF_pred, y = allPred$conn_pred,
#                              FUN = function(x, y)
#                                
#                                any((x - 1 / nSamp) < presentFixedEff$presentCoverBF & 
#                                      (x + 1 / nSamp) > presentFixedEff$presentCoverBF  &
#                                      (y - maxConnectivity / nSamp) < presentFixedEff$presentConnW  &
#                                      (y + maxConnectivity / nSamp) > presentFixedEff$presentConnW ))
# 
# allPred$CFpresent <- mapply( x = allPred$CF_pred, y = allPred$conn_pred,
#                              FUN = function(x, y)
#                                
#                                any((x - 1 / nSamp) < presentFixedEff$presentCoverCF &
#                                      (x + 1 / nSamp) > presentFixedEff$presentCoverCF  &
#                                      (y - maxConnectivity / nSamp) < presentFixedEff$presentConnW  &
#                                      (y + maxConnectivity / nSamp) > presentFixedEff$presentConnW ))
# 
# # Create separate column where prediction is replaced with 'NA' if no presence records from the cell
# allPred <-  allPred %>%
#   mutate(BFmedianPres = if_else(BFpresent == TRUE, BFmedian, NA)) %>%
#   mutate(CFmedianPres = if_else(CFpresent == TRUE, CFmedian, NA))
# 
# # Convert to long data format
# allPred <- allPred %>% 
#   gather(coverType, cover_pred, "BF_pred" , "CF_pred") %>%
#   mutate(median = if_else(coverType == "BF_pred", BFmedian , CFmedian )) %>%
#   mutate(medianPres = if_else(coverType == "BF_pred", BFmedianPres , CFmedianPres ))
# 
# ### Plot
# 
# # Plot broadleaf and coniferous cover-connectivity interaction
# # (contours based on all predictions, fill only uses prediction space where species is present)
# intPlot <- ggplot(allPred, aes(x = cover_pred, y = conn_pred, z = median)) +
#   ggtitle("Woodland cover - connectivity interaction") +
#   facet_wrap( ~coverType, labeller = as_labeller(intLabels)) +
#   geom_tile(aes(fill = medianPres, colour = medianPres)) +
#   stat_contour(bins = 50, colour = "black") +
#   scale_fill_distiller(na.value = NA,
#                        palette = "RdYlBu",
#                        direction = 1,
#                        guide = guide_colourbar(title = "Posterior\nmedian"),
#                        limits = c(-1,1) * max(abs(na.omit(allPred$medianPres)))) +
#   scale_colour_distiller(na.value = NA,
#                          palette = "RdYlBu",
#                          direction = 1,
#                          guide = "none",
#                          limits = c(-1,1) * max(abs(na.omit(allPred$medianPres)))) +
#   scale_x_continuous(name="Proportion woodland cover") +
#   scale_y_continuous(name="Connectivity (A)") +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5, vjust = -1),
#         panel.grid.major = element_line(colour = "darkgrey"),
#         strip.text.x = element_text(size = 12))
# 
# # JOINT EVALUATION PLOT
# evalPlot <- arrangeGrob(predMedian, predSD, fixedEffPlot, intPlot, randomEffPlot, spaceTimePlot,
#                         nrow = 3, ncol = 2, 
#                         layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 6)),
#                         top = grid::textGrob(paste0(iSpecies, ", logCPO = ", logCPO), gp = grid::gpar(fontsize=20)))
# 
# # Compose matern plot
# spdeAndMaternPlot <- arrangeGrob(range.plot, covplot,
#                                  var.plot, corplot, ncol = 2)
# 
# # SAVE -----------------------------------------
# 
# # Define species directory
# iSpeciesDir <- paste0("Data/SDM_output/Output/", 
#                       taxaGroup, "/", 
#                       iSpeciesTidy)
# 
# # Create species directory
# if (!dir.exists(iSpeciesDir)) {
#   
#   dir.create(iSpeciesDir,
#              recursive = TRUE)
# }
# 
# ### Save objects
# 
# # Model
# save(model,
#      file = paste0(iSpeciesDir,
#                    "/modelFit.RData"))
# # Model summary
# save(modelSummary,
#      file = paste0(iSpeciesDir,
#                    "/modelSummary.RData"))
# 
# # Model prediction object
# save(modelPred,
#      file = paste0(iSpeciesDir,
#                    "/modelPred.RData"))
# 
# # Posterior median relative occurrence probability prediction
# writeRaster(median_R,
#             file = paste0(iSpeciesDir,
#                           "/medianPred.tif"),
#             overwrite= TRUE)
# 
# # Evaluation plot
# ggsave(paste0(iSpeciesDir,
#               "/effectsPlot_range_", estimated_range,
#               "_logCPO_", logCPO, ".png"),
#        evalPlot,
#        width = 6000, height = 6000, 
#        units = "px", dpi = 400,
#        limitsize = FALSE)
# 
# # SPDE parameter posterior (range and variance) and
# # matern correlation and covariance plot
# ggsave(paste0(iSpeciesDir,
#               "/MatCorCovPlot_range_", estimated_range,
#               "_logCPO_", logCPO, ".png"),
#        spdeAndMaternPlot,
#        width = 6000, height = 3000,
#        units = "px", dpi = 400,
#        limitsize = FALSE)

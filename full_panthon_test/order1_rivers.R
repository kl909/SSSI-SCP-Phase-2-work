library(sf)
library(terra)
library(dplyr)

#load data
rivers <- st_read('data/spatial_data/river_data/order1/order_1.gpkg')
grid <- st_read('data/spatial_data/river_data/1km_grid.shp')

# assign grid same projection as river
st_crs(grid) <- 27700

# create grid_id column
grid <- grid %>% mutate(grid_id = row_number())

# generate 250 m buffer of rivers (dissolved)
river_buffer <- st_buffer(rivers, dist = 250) %>% 
  st_union() # st_union performs the 'dissolve'

# save buffer
#st_write(river_buffer, "data/spatial_data/river_data/order2/river_buffer_250m.shp", delete_dsn = TRUE)

# intersect and calculate length of rivers
# 'st_intersection' splits the river lines at the grid boundaries
intersections <- st_intersection(rivers, grid) %>%
  mutate(river_length = st_length(.)) # Calculate length of segments

# dissolve by grid_id (summing lengths per pixel)
river_stats <- intersections %>%
  st_drop_geometry() %>% # Drop geometry for faster aggregation
  group_by(grid_id) %>%
  summarize(total_length = sum(river_length))

# save layer
#grid_with_lengths <- grid %>%
#  inner_join(river_stats, by = "grid_id")

#st_write(grid_with_lengths, "data/spatial_data/river_data/order2/grid_river_lengths_dissolved.shp", delete_dsn = TRUE)

# Filter grid by buffer and join data (keeping only pixels within 250 m of rivers)
grid_filtered <- grid[st_intersects(grid, river_buffer, sparse = FALSE), ] %>%
  left_join(river_stats, by = "grid_id")

# handle NULLS (Coalesce to 0)
grid_filtered <- grid_filtered %>%
  mutate(final_length = ifelse(is.na(total_length), 0, total_length))

# rasterize
# create template raster based on the grid extent and 1km reso
template_raster <- rast(grid, res = 1000)

# burn the final length of each cell into the raster
river_raster <- rasterize(vect(grid_filtered), 
                          template_raster, 
                          field = "final_length", 
                          background = NA)

# save output
writeRaster(river_raster, "data/spatial_data/river_data/order1/order1.tif", overwrite = TRUE)


library(sf)
library(dplyr)
library(tidygraph)
library(sfnetworks)
library(here)

# 1. Import OS river data
rivers <- st_read("../../data/os_rivers/data/WatercourseLink.shp")

# 2. Convert to a network object
net <- as_sfnetwork(rivers, directed = TRUE, from = startNode, to = endNode)

# 3. Get node information
node_status <- net %>%
  activate("nodes") %>%
  mutate(
    in_degree = centrality_degree(mode = "in"),
    out_degree = centrality_degree(mode = "out")
  ) %>%
  as_tibble()

# 4. Apply Order 1 and Order 2
net_ordered <- net %>%
  activate("edges") %>%
  mutate(
    # Find the status of the 'from' node for every edge
    in_at_start = node_status$in_degree[from],
    
    # Logic: 
    # - If in_degree is 0, it's a headwater (Order 1)
    # - If in_degree is 1, it's a continuation of the previous stream (Stay Order 1)
    # - If in_degree >= 2, it's a junction (Order 2)
    order = ifelse(in_at_start >= 2, 2, 1)
  )

#------------------------------------------------
rivers_final2 <- net_ordered %>%
  activate("edges") %>%
  st_as_sf()

# 8. Clean up columns and geometry
rivers_final2 <- rivers_final2 %>%
  select(-from, -to) %>%
  st_zm() # Make 2D

# 9. Save as shapefile
st_write(
  rivers_final2, 
  "../../data/os_rivers/data/OS_Watercourses_Ordered_2.shp", 
  delete_layer = TRUE
)

library(sf)
library(dplyr)
library(tidygraph)
library(sfnetworks)
library(here)
library(purrr)

# 1. Import OS river data
rivers <- st_read("data/spatial_data/river_data/os_rivers/data/WatercourseLink.shp")

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

# step for continuing order 2 rivers
for (i in 1:50) {
  # Get the 'to' node of every edge that is already Order 2
  # this is vectorized and very fast
  order2_nodes <- net_ordered %>%
    activate("edges") %>%
    as_tibble() %>%
    filter(order == 2) %>%
    pull(to) %>%
    unique()
  
  # Update the edges: If an edge starts (from) at one of these nodes, 
  # it becomes Order 2.
  net_ordered <- net_ordered %>%
    activate("edges") %>%
    mutate(
      order = ifelse(from %in% order2_nodes, 2, order)
    )
}


####### order 3
# Identify nodes where two or more Order 2 rivers meet
order3_birth_nodes <- net_ordered %>%
  activate("edges") %>%
  as_tibble() %>%
  group_by(to) %>%
  summarise(count_order2 = sum(order == 2)) %>%
  filter(count_order2 >= 2) %>%
  pull(to)

# Assign Order 3 to edges starting at these junctions
net_ordered <- net_ordered %>%
  activate("edges") %>%
  mutate(
    order = ifelse(from %in% order3_birth_nodes, 3, order)
  )

# make sure order 3 rivers stay order 3
# Push the Order 3 status downstream through simple nodes
for (i in 1:50) {
  # Find nodes that have an Order 3 river flowing INTO them
  order3_nodes <- net_ordered %>%
    activate("edges") %>%
    as_tibble() %>%
    filter(order == 3) %>%
    pull(to) %>%
    unique()
  
  # Update the edges: if it starts at an Order 3 node, it stays Order 3
  net_ordered <- net_ordered %>%
    activate("edges") %>%
    mutate(
      order = ifelse(from %in% order3_nodes, 3, order)
    )
}


#------------------------------------------------
rivers_final2 <- net_ordered %>%
  activate("edges") %>%
  st_as_sf()

# 8. Clean up columns and geometry
rivers_final2 <- rivers_final2 %>%
  select(-from, -to) %>%
  st_zm() # Make 2D

# 9. Save as whole shapefile
st_write(
  rivers_final2, 
  "data/spatial_data/river_data/OS_Watercourses_Ordered_2.shp", 
  delete_layer = TRUE
)

# save order 1:
order1_rivers <- rivers_final2 %>%
  filter(order == 1)

st_write(
  order1_rivers, 
  "data/spatial_data/river_data/order1/order_1.gpkg", 
  delete_dsn = TRUE
)

# save order 2:
order2_rivers <- rivers_final2 %>%
  filter(order == 2)

st_write(
  order2_rivers, 
  "data/spatial_data/river_data/order2/order_2.gpkg", 
  delete_dsn = TRUE
)

# save order 3:
order3_rivers <- rivers_final2 %>%
  filter(order == 3)

st_write(
  order3_rivers, 
  "data/spatial_data/river_data/order3/order_3.gpkg", 
  delete_dsn = TRUE
)
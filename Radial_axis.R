rm(list=ls())
gc()
.rs.restartR()
set.seed(67)

library(Seurat)
library(RANN)
library(tidyr)
library(magrittr)
library(ggplot2)

# Load object
obj <- readRDS(spatial_atlas_path)

unique(obj$cell_annotations)
obj$cell_annotations[which(obj$cell_annotations %in% c('Pan-AC','Grey_AC','White_AC','Pial_AC'))] <- 'AC'

# Extract coordinates and cell names
coordinates <- obj@meta.data[, c("compactx", "compacty")]
cell_names <- obj@meta.data$Cell_ID
cell_types <- obj$cell_annotations

# Identify endothelial cells
ec_indices <- which(cell_types %in% c("Arteriole",'Capillary','Venule'))
ec_coordinates <- coordinates[ec_indices, , drop = FALSE]

# Set maximum search radius
max_radius <- 50

unique(cell_types)
search_CT <- c('Arteriole','Capillary','Venule','Matrix_PC','Transport_PC','aaSMC','Parenchymal_FB','Parenchymal_PVM','MG',
               'L2/3','L4','L5','L6',
               'VIP','PAX6/SNCG','LAMP5','SST_CHODL','SST','PVALB',
               'OL','OPC','AC')

# Perform kNN search up to max_radius
nn <- nn2(coordinates, ec_coordinates, searchtype = "radius", radius = max_radius)

get_neighbors_in_bin <- function(i, nn, cell_names, cell_types, min_dist, max_dist) {
  distances <- nn$nn.dists[i, ]
  valid_indices <- which(distances >= min_dist & distances < max_dist)
  neighbor_ids <- nn$nn.idx[i, valid_indices]
  
  valid_neighbors <- cell_names[neighbor_ids]
  valid_types <- cell_types[neighbor_ids]
  return(sum(valid_types == CT, na.rm = TRUE))
}

# Store results in bins
neighbor_bins <- expand.grid(Bin = 1:max_radius, CT = search_CT)
neighbor_bins$Count <- 0

for (CT in search_CT) {
  for (bin in 1:max_radius) {
    min_dist <- bin - 1
    max_dist <- bin
    neighbor_bins$Count[neighbor_bins$Bin == bin & neighbor_bins$CT == CT] <- sum(
      sapply(1:length(ec_indices), get_neighbors_in_bin, nn, cell_names, cell_types, min_dist, max_dist)
    )
  }
}

# Identify the peak bin for each CT
peak_info <- neighbor_bins %>%
  group_by(CT) %>%
  summarise(max_bin = Bin[which.max(Count)], max_count = max(Count)) %>%
  arrange(max_bin)

# Convert CT to factor with levels ordered by ascending peak bin
neighbor_bins <- left_join(neighbor_bins, peak_info, by = "CT") %>%
  mutate(CT = factor(CT, levels = peak_info$CT))

# Use ggplot for plotting hereafter. 


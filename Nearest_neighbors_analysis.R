rm(list=ls())
gc()
.rs.restartR()
set.seed(67)

library(Banksy)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(gridExtra)
library(pals)
library(tibble)
library(dplyr)
library(RANN)
library(tidyr)
library(purrr)

#Start KNN
obj <- readRDS(spatial_atlas_path)
coordinates <- obj@meta.data[, c("compactx", "compacty")]
radius <- 10 # mural cells
radius <- 20 # perivascular cells
cell_names <- obj@meta.data$Cell_ID
nn <- nn2(coordinates, coordinates, searchtype = "radius", radius = radius)
get_neighbors_n <- function(i, nn) {
  nn$nn.idx[i, nn$nn.idx[i,] != 0]
}
get_neighbor_names <- function(i, nn, cell_names) {
  neighbor_indices <- nn$nn.idx[i, nn$nn.idx[i,] != 0]
  return(cell_names[neighbor_indices])
}
neighbors_name_list <- lapply(1:nrow(coordinates), get_neighbor_names, nn, cell_names)
names(neighbors_name_list) <- cell_names

neighbors_n_list <- lapply(1:nrow(coordinates), get_neighbors_n, nn)
obj$neighbors <- sapply(neighbors_n_list, length)

summary(obj$neighbors)
FeaturePlot(obj,features = 'neighbors')

cell_types <- obj$cell_annotations
perform_hypergeometric_test <- function(neighbors, all_cells, cell_type_A, cell_type_B) {
  n_neighbors <- length(neighbors)
  n_type_B_neighbors <- sum(cell_types[neighbors] == cell_type_B)
  n_type_B_total <- as.numeric(table(cell_types == cell_type_B)['TRUE'])
  n_total <- length(all_cells)
  
  p_value <- phyper(n_type_B_neighbors - 1, n_type_B_total, n_total - n_type_B_total, n_neighbors, lower.tail = FALSE)
  odds_ratio <- (n_type_B_neighbors / (n_neighbors - n_type_B_neighbors)) / 
    (n_type_B_total / (n_total - n_type_B_total))
  
  return(c(p_value = p_value, odds_ratio = odds_ratio))
}
pairwise_enrichment <- function(neighbors_list, cell_types) {
  unique_cell_types <- unique(cell_types)
  n_types <- length(unique_cell_types)
  p_value_matrix <- matrix(NA, nrow = n_types, ncol = n_types, 
                           dimnames = list(unique_cell_types, unique_cell_types))
  odds_ratio_matrix <- matrix(NA, nrow = n_types, ncol = n_types, 
                              dimnames = list(unique_cell_types, unique_cell_types))
  
  for (i in 1:n_types) {
    type_A_cells <- names(cell_types)[cell_types == unique_cell_types[i]]
    type_A_neighbors <- unlist(neighbors_list[type_A_cells])
    
    for (j in 1:n_types) {
      if (i != j) {
        test_result <- perform_hypergeometric_test(type_A_neighbors, names(cell_types), 
                                                   unique_cell_types[i], unique_cell_types[j])
        p_value_matrix[i, j] <- test_result["p_value"]
        odds_ratio_matrix[i, j] <- test_result["odds_ratio"]
      }
    }
  }
  
  return(list(p_values = p_value_matrix, odds_ratios = odds_ratio_matrix))
}
enrichment_results <- pairwise_enrichment(neighbors_name_list, cell_types)

# Convert matrices to data frames
p_value_df <- as.data.frame(enrichment_results$p_values)
odds_ratio_df <- as.data.frame(enrichment_results$odds_ratios)

# Add cell type column to each data frame
p_value_df$cell_type_A <- rownames(p_value_df)
odds_ratio_df$cell_type_A <- rownames(odds_ratio_df)

p_value_long <- tidyr::pivot_longer(p_value_df, 
                                    cols = -cell_type_A, 
                                    names_to = "cell_type_B", 
                                    values_to = "p_value")

odds_ratio_long <- tidyr::pivot_longer(odds_ratio_df, 
                                       cols = -cell_type_A, 
                                       names_to = "cell_type_B", 
                                       values_to = "odds_ratio")

# Remove rows where cell_type_A == cell_type_B
p_value_long <- p_value_long[p_value_long$cell_type_A != p_value_long$cell_type_B, ]
odds_ratio_long <- odds_ratio_long[odds_ratio_long$cell_type_A != odds_ratio_long$cell_type_B, ]

# Apply FDR correction
p_value_long$p_value_fdr <- p.adjust(p_value_long$p_value, method = "BH")

# Apply the cap to FDR-corrected p-values
min_p_value <- 1e-300
cap_p_value <- function(p) {
  return(max(p, min_p_value))
}
p_value_long$p_value_fdr_capped <- sapply(p_value_long$p_value_fdr, cap_p_value)

# Merge with odds ratio data
combined_df <- merge(p_value_long, odds_ratio_long, by = c("cell_type_A", "cell_type_B"))



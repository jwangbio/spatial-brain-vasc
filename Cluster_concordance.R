rm(list=ls())
gc()
.rs.restartR()
set.seed(67)

library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(magrittr)

spatial_ref <- readRDS(spatial_atlas_path)
cell_ref <- readRDS(cell_atlas_path)
spatial <- spatial_ref
cell <- cell_ref

spatial$cell_annotations[which(spatial$cell_annotations %in% c('Grey_AC','Pan-AC','White_AC','Pial_AC'))] <- 'AC'
spatial$cell_annotations[which(spatial$cell_annotations %in% c('Pan-Ex','VIP','','L6','L5/6','L4/5','L2/3','L5','PVALB','SST','LAMP5','L4','SST_CHODL','PAX6/SNCG'))] <- 'Neu'
spatial$cell_annotations[which(spatial$cell_annotations %in% c('Pial_PVM','Parenchymal_PVM'))] <- 'PVM'
spatial$cell_annotations[which(spatial$cell_annotations %in% c('Artery','Arteriole','Capillary','Venule','Vein'))] <- 'EC'
spatial$cell_annotations[which(spatial$cell_annotations %in% c('Matrix_PC','Transport_PC'))] <- 'PC'
spatial$cell_annotations[which(spatial$cell_annotations %in% c('aSMC','aaSMC','FBMC'))] <- 'SMC'
spatial$cell_annotations[which(spatial$cell_annotations %in% c('MG','PVM','Mono'))] <- 'Myeloid'
spatial$cell_annotations[which(spatial$cell_annotations %in% c('Dural_Border_FB','Arachnoid_Barrier_FB','Inner_Arachnoid_FB','Pial_FB','Parenchymal_FB'))] <- 'FB'

cell$reconciliation[which(cell$reconciliation %in% c('Artery','Arteriole','Capillary','Venule','Vein'))] <- 'EC'
cell$reconciliation[which(cell$reconciliation %in% c('Matrix_PC','Transport_PC'))] <- 'PC'
cell$reconciliation[which(cell$reconciliation %in% c('aSMC','aaSMC','FBMC'))] <- 'SMC'
cell$reconciliation[which(cell$reconciliation %in% c('MG','PVM','Mono'))] <- 'Myeloid'
cell$reconciliation[which(cell$reconciliation %in% c('Dural_Border_FB','Arachnoid_Barrier_FB','Inner_Arachnoid_FB','Pial_FB','Parenchymal_FB'))] <- 'FB'

spatial_aggregate <- AggregateExpression(spatial,return.seurat = T,group.by = 'cell_annotations')
cell_aggregate <- AggregateExpression(cell,return.seurat = T,group.by = 'reconciliation')

CT_intersect <- intersect(unique(cell$reconciliation),unique(spatial$cell_annotations))
CT_intersect <- gsub('_','-',CT_intersect)
CT_interest1 <- CT_intersect[1]
CT_interest2 <- CT_intersect[2]

print(CT_interest1 %in% unique(cell_aggregate$reconciliation) & CT_interest1 %in% unique(spatial_aggregate$cell_annotations))
print(CT_interest2 %in% unique(cell_aggregate$reconciliation) & CT_interest2 %in% unique(spatial_aggregate$cell_annotations))

spatial_subset <- spatial_aggregate
cell_subset <- cell_aggregate

cell_matrix <- cell_subset@assays$RNA$data
cell_matrix <- cell_matrix[which(rownames(cell_matrix) %in% genes),which(colnames(cell_matrix) %in% CT_interest1)]
cell_matrix <- cell_matrix[genes]
cell_matrix <- na.omit(cell_matrix)
spatial_matrix <- spatial_subset@assays$SCT$counts
spatial_matrix <- spatial_matrix[which(rownames(spatial_matrix) %in% genes),colnames(spatial_matrix) %in% CT_interest2]
spatial_matrix <- spatial_matrix[genes]
spatial_matrix <- na.omit(spatial_matrix)
correlation <- cor(cell_matrix, spatial_matrix, method = "spearman")
print(correlation)

# Choose which one to fine-tune.
Idents(spatial) <- spatial$cell_annotations
markers <- FindAllMarkers(spatial,logfc.threshold = 0.1,min.diff.pct = 0.1,only.pos = TRUE,assay = 'SCT',test.use = 'wilcox')
Idents(cell) <- cell$reconciliation
markers2 <- FindAllMarkers(cell,logfc.threshold = 0.1,min.diff.pct = 0.1,only.pos = TRUE,assay = 'RNA',test.use = 'wilcox')
markers$merged <- paste0(markers$gene,'_',markers$cluster)
markers2$merged <- paste0(markers2$gene,'_',markers2$cluster)

genes <- intersect(markers$merged,markers2$merged)
split_genes <- strsplit(genes, "_")
split_genes_df <- do.call(rbind, split_genes)
genes <- split_genes_df[,1]

# Get intersected cell types
CT_intersect <- intersect(unique(cell$reconciliation), unique(spatial$cell_annotations))
CT_intersect <- gsub('_', '-', CT_intersect)

# Ensure that all intersected cell types are present in both datasets
CT_intersect <- CT_intersect[
  CT_intersect %in% unique(cell_aggregate$reconciliation) &
    CT_intersect %in% unique(spatial_aggregate$cell_annotations)
]

CT_intersect <- c('EC','PC','FB','SMC','Myeloid','Neu','AC','OPC','OL')

# Initialize an empty correlation matrix
correlation_matrix <- matrix(
  NA, 
  nrow = length(CT_intersect), 
  ncol = length(CT_intersect), 
  dimnames = list(CT_intersect, CT_intersect)
)

# Iterate over all combinations of CT_intersect
for (CT_interest1 in CT_intersect) {
  for (CT_interest2 in CT_intersect) {
    cell_matrix <- cell_subset@assays$RNA$data
    cell_matrix <- cell_matrix[which(rownames(cell_matrix) %in% genes),which(colnames(cell_matrix) %in% CT_interest1)]
    cell_matrix <- cell_matrix[genes]
    cell_matrix <- na.omit(cell_matrix)
    spatial_matrix <- spatial_subset@assays$SCT$counts
    spatial_matrix <- spatial_matrix[which(rownames(spatial_matrix) %in% genes),colnames(spatial_matrix) %in% CT_interest2]
    spatial_matrix <- spatial_matrix[genes]
    spatial_matrix <- na.omit(spatial_matrix)
    
    correlation <- cor(as.numeric(cell_matrix), as.numeric(spatial_matrix), method = "pearson")
    correlation_matrix[CT_interest1,CT_interest2] <- correlation
  }
}

formatted_correlation_matrix <- correlation_matrix
formatted_correlation_matrix_text <- matrix(
  "", 
  nrow = nrow(correlation_matrix), 
  ncol = ncol(correlation_matrix),
  dimnames = dimnames(correlation_matrix)
)

for (i in 1:nrow(correlation_matrix)) {
  max_col <- which.max(correlation_matrix[i, ])
  formatted_correlation_matrix_text[i, max_col] <- formatC(
    correlation_matrix[i, max_col], 
    format = "f", digits = 2
  )
}
correlation_long <- as.data.frame(as.table(correlation_matrix)) %>%
  rename(Cell = Var1, Spatial = Var2, Correlation = Freq)
formatted_text_long <- as.data.frame(as.table(formatted_correlation_matrix_text)) %>%
  rename(Cell = Var1, Spatial = Var2, Label = Freq)
heatmap_data <- correlation_long %>%
  left_join(formatted_text_long, by = c("Cell", "Spatial"))

# Create the heatmap with ggplot2
ggplot(heatmap_data, aes(x = Spatial, y = Cell, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Correlation"
  ) +
  geom_text(aes(label = Label), size = 5) +
  scale_y_discrete(limits = rev(unique(heatmap_data$Cell))) + 
  labs(
    title = "Correlation Heatmap",
    x = "Spatial",
    y = "Cell"
  ) +
  theme_classic(base_size = 25) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20),
    plot.title = element_text(hjust = 0.5, size = 30),
    legend.position = "right"
  )
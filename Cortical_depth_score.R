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
library(tidyr)
library(harmony)
library(monocle3)
library(viridis)
library(slingshot)
library(rstatix)
library(RANN)
library(magrittr)
library(R.cache)

obj <- readRDS(spatial_atlas_path)
keepers <- unique(obj$cell_annotations[which(! obj$cell_annotations %in% c('Admixture','Pan-AC','Pan-Ex','Non-specific'))])
obj <- subset(obj,subset = cell_annotations %in% keepers)
BANKSY = RunBanksy(obj, lambda = 0.95, assay = 'SCT', slot = 'data',
                   dimx = 'x', dimy = 'y', features = 'all',
                   group = 'Section_ID', split.scale = TRUE, k_geom = 50)
BANKSY = RunPCA(BANKSY, assay = 'BANKSY', features = rownames(BANKSY), npcs = 30)
BANKSY = RunUMAP(BANKSY, dims = 1:30)
BANKSY = FindNeighbors(BANKSY, dims = 1:30)
BANKSY = FindClusters(BANKSY, resolution = 0.25,algorithm = 2)

DimPlot(BANKSY, pt.size = 0.75, label = TRUE, label.size = 5,split.by = 'Slide_ID')
FeaturePlot(BANKSY, features = c("AQP4", "TREM2", "MOG", "SNAP25")) #Parenchyma
FeaturePlot(BANKSY, features = c('CUX2','TSHZ2','RORB','GRIK3')) # Cortical Lamination
FeaturePlot(BANKSY, features = c("MYH11",'CLDN5','PDGFRB','FBLN1')) #Vascular and Perivascular
FeaturePlot(BANKSY, features = c('LTBP4','ALPL','MFSD2A','IL1R1')) # AV Zonations
FeaturePlot(BANKSY, features = c('P2RY12','SPP1','NAV3','TREM2')) # PVM

#Manually group domains
BANKSY <- subset(BANKSY,subset =  BANKSY_snn_res.0.25 %in% names(which(table(BANKSY$BANKSY_snn_res.0.25)>100))) #keep clusters with >100 cells
Idents(BANKSY) <- BANKSY$BANKSY_snn_res.0.25
DimPlot(BANKSY, pt.size = 0.75, label = TRUE, label.size = 5)
as.data.frame(table(BANKSY$BANKSY_snn_res.0.25,BANKSY$cell_annotations)) %>% 
  group_by(Var1) %>%
  mutate(scale_1 = Freq / max(Freq)) %>% 
  ungroup() %>% 
  group_by(Var2) %>%
  mutate(scale_2 = Freq / max(Freq)) %>% 
  ungroup() %>% 
  filter(Var1 %in% names(which(table(BANKSY$BANKSY_snn_res.0.25)>100))) %>% 
  ggplot(aes(x = Var2,y = Var1,fill=scale_2)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
manual <- c( #domains
  "L2/3", #0
  "L2/3",
  "L5/6",
  "L4/5",
  "Sulci",
  "L4/5",
  "WM", #6
  "WM",
  "WM",
  "L4/5",
  "Sulci", #10
  "WM",
  "Sulci",
  "Non-specific",
  "Subarachnoid",
  "Subarachnoid",
  "Non-specific",
  "Non-specific", #17
  "Non-specific",
  "Non-specific",
  "Non-specific", #20
  "Non-specific",
  "Non-specific",
  "Non-specific",
  "Non-specific", #24
  "Non-specific",
  "Non-specific"
)
annot_tbl <- cbind(0:(length(manual)-1),as.character(manual))
colnames(annot_tbl) <- c("seurat_clusters","domains")
annot_tbl <- as_tibble(annot_tbl)
BANKSY@meta.data <- left_join(BANKSY@meta.data,annot_tbl,by = "seurat_clusters")
rownames(BANKSY@meta.data) <- BANKSY$Cell_ID
Idents(BANKSY) <- BANKSY$domains
DimPlot(BANKSY, pt.size = 0.75, label = TRUE, label.size = 5)
table(BANKSY$domains)
FeatureScatter(BANKSY, 'compactx', 'compacty', pt.size = 1,group.by = 'domains') + coord_equal()

#Pseudotime
obj <- BANKSY
obj <- subset(obj,subset = domains %in% c('L2/3','L4/5','L5/6','Sulci'))
cds <- SeuratWrappers::as.cell_data_set(obj)
cds <- cluster_cells(cds,partition_qval = 0.05,random_seed = 67)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
cds <- learn_graph(cds,learn_graph_control = list(minimal_branch_len = 5))
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE,cell_size = 5,group_label_size = 5)
max1 <- which.max(unlist(FetchData(obj, "CD22")))
max5 <- which.max(unlist(FetchData(obj, "COL1A1")))
max_all <- colnames(obj)[c(max1,max5)]
cds <- order_cells(cds, root_cells = max_all)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE,show_trajectory_graph = F,rasterize = F,cell_size = 1) + 
  scale_color_viridis()  
traj.plot <- plot_cells(cds, color_cells_by = "pseudotime")
point.data <- ggplot_build(traj.plot)[["plot"]][["data"]] %>% 
  as_tibble()
ggplot(point.data,aes(x = data_dim_1,data_dim_2,color = cell_color)) +
  geom_point() +
  scale_color_viridis() +
  theme_classic()

# Step 4: Create the plot
color_key_full <- readRDS(color_key_path)
color_key_full
celltypes <- c('aSMC','Pial_AC','L2/3','L4/5','L5/6','L6')
point.data <- subset(point.data,subset = cell_annotations %in% celltypes)
ggplot(point.data, aes(x = cell_color, fill = cell_annotations)) +
  geom_histogram(position = "identity", binwidth = 0.25, alpha = 0.6) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 14, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 27, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 35, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 42, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 49, linetype = "dashed", color = "black") +
  labs(x = "Pseudotime (Cortical Depth)",
       y = "Frequency of Cell Types",
       fill = "Neighbor Cell Type") +
  facet_wrap(~factor(cell_annotations,levels = celltypes),
             ncol = 6,
             scales = "free_x"
  ) +
  theme_classic() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_reverse() +
  coord_flip() +
  scale_fill_manual(values = color_key_full)



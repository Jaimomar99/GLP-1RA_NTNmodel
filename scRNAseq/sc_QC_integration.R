suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(cowplot)))
suppressWarnings(suppressMessages(library(patchwork)))
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(cluster)))
suppressWarnings(suppressMessages(library(harmony)))

manifest <- read.csv("/Data/scRNA/kidney/mouse_NTN/manifest_snRNAseq.csv")

list_NTN <- lapply(1:nrow(manifest), function(i) {
  row <- manifest[i,]
  seurat_obj <- CreateSeuratObject(counts = Read10X_h5(row$path_h5),
                                 project = "mouse_NTN_snRNAseq")
  seurat_obj$alias_id <- as.factor(row$alias_id)
  seurat_obj$Treatment <- as.factor(row$Treatment)
  
  return(seurat_obj)
})

names(list_NTN) <- manifest$alias_id

####QC
merged_seurat <- merge(x = list_NTN[[1]],
                       y = list_NTN[2:length(list_NTN)],
                       add.cell.ids = names(list_NTN),
                       project = "mouse_NTN_snRNAseq")

merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)


mito_genes <- rownames(merged_seurat)[grep("^mt-", rownames(merged_seurat))]
total_counts_per_cell <- colSums(merged_seurat@assays$RNA@counts)
merged_seurat$percent_mito <- colSums(merged_seurat@assays$RNA@counts[mito_genes, ])/total_counts_per_cell

# Spot level
cells <- WhichCells(merged_seurat, expression = nCount_RNA > 500 &
                      nFeature_RNA > 200 &
                      percent_mito < 30)
# Gene level
genes <- rownames(merged_seurat)[Matrix::rowSums(merged_seurat) > 10]

# Subset the Seurat object based on the filtering criteria
filtered_seurat <- subset(merged_seurat, features = genes, cells = cells)

###Clustering
res_range <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
normalized_seurat <- SCTransform(filtered_seurat,
                                 vars.to.regress = "percent_mito",
                                 variable.features.n = 3000,
                                 vst.flavor = "vst2",
                                 return.only.var.genes = FALSE) %>% 
                    RunPCA(.) %>% 
                    RunHarmony(., c("alias_id", "Treatment"), #default parameters
                               plot_convergence=TRUE,
                               lambda = rep(1, length(c("alias_id", "Treatment"))),
                               max.iter = 20) %>% #lambda is importance to each 
                    RunUMAP(., reduction = "harmony", assay = "SCT", dims = 1:50) %>% 
                    FindNeighbors(., dims = 1:50) %>% 
                    FindClusters(., resolution  = res_range)

#silohuette score
distance_matrix_harmony <- dist(Embeddings(normalized_seurat[["harmony"]]))
  
silhouette_scores_list <- list()
  
for (i in seq_along(res_range)) {
    res <- paste0("SCT_snn_res.", res_range[i])
    print(res)
    clusters <- normalized_seurat@meta.data[res] %>% pull()
    silhouette_scores_list[[i]] <- silhouette(as.numeric(clusters), dist = distance_matrix_harmony) %>%
      as_tibble() %>%
      mutate(res=res_range[i])
  }
  #plot
sil_plot <- do.call(rbind, silhouette_scores_list) %>% 
          mutate(cluster = factor(cluster)) %>% 
          mutate(res = factor(res)) %>% 
          group_by(cluster,res) %>% 
          summarise("mean_silhoute" = mean(sil_width)) %>% 
          ggplot(aes(x=res, y= mean_silhoute)) +
          geom_violin(draw_quantiles=c(0.5), scale="width") + 
          geom_jitter(
            shape = 21, 
            alpha = .8, 
            size = 2,
            aes(fill = cluster)) +
          labs(x = "Cluster Resolution", y = "Cell Silhouette Score", title= "Silhoutte score clustering") +
          theme_classic() +
          theme(legend.position = "none")

#maximum
max <- lapply(silhouette_scores_list, function(df) {
  tryCatch(
      {
      df[, 3] %>% pull() %>% mean()
      },
      error = function(e) {
      NA
      }
  )
  }) %>% unlist() %>% which.max()
best_res <- paste0("SCT_snn_res.", res_range[max])
print(best_res)

Idents(normalized_seurat) <- best_res
normalized_seurat$clustering_harmony <- Idents(normalized_seurat)

saveRDS(normalized_seurat, "Data/scRNA/kidney/mouse_NTN/processed/normalized_v2.rds")

print("Script done, object saved as normalized_v2.rds")
#input dataset information by developers
plot_dir <- paste0("./projects/kidney_mouse_ntn/reports/")
project_dir <- "./Data/Spatial/visium/kidney/mouse_NTN/"
processed_dir <- "./Data/Spatial/visium/kidney/mouse_NTN/processed/ultimate_processed/"
pipeline_dir <- "./software/project_repos/visium_pipeline"
manifest_path <- paste0(project_dir, "/manifest.csv")
group_variable <- "diagnosis"

suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(patchwork)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(gridExtra)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(scales)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(cowplot)))
suppressMessages(source(paste0(pipeline_dir, "/R/utils.R")))

###READING DATA
sample_metadata <- read.csv(manifest_path)

# by sample
sample_list <- sample_metadata$alias_id
object_list <- vector('list', length(sample_list))
for (i in 1:length(sample_list)) {
    sample_name <- paste0("sample", sample_list[[i]])
    curr_obj <- readRDS(file = file.path(processed_dir, "by_sample/rds/", paste0(sample_list[[i]], ".rds")))
    
    if (length(curr_obj@assays$Spatial@scale.data) == 0){
        # normalize and scale data so we can plot and color code
        # curr_obj <- NormalizeData(curr_obj, assay = "Spatial")
        # curr_obj <- ScaleData(curr_obj, assay = "Spatial", verbose = F)
        message("Norm & scaling sample")
    }

    object_list[[i]] <- curr_obj
}

#integrated
integrated_obj <- readRDS(file.path(processed_dir, "integrated/rds/integrated.rds"))

map_meta_integrated <- function(integrated_obj, object_list, columns) {
    """
    This function maps the metadata from the integrated object to the individual samples.
    """ 

  object_list <- lapply(1:length(object_list), FUN = function(x) {
    subset_integrated <- integrated_obj@meta.data[, columns] %>% 
      filter(alias_id == x) %>% 
      mutate(seq = rownames(.)) %>% 
      mutate(seq = gsub("_\\d+$", "", seq))
    
    rownames(subset_integrated) <- subset_integrated$seq
    object_list[[x]]@meta.data <- merge(object_list[[x]]@meta.data, subset_integrated, by = 0, all.x = TRUE, suffixes = "") %>% 
      column_to_rownames("Row.names") %>% 
      select(!c("seq", "alias_idNA"))
      
    # Convert specified columns to factors
    for (column in columns) {
      object_list[[x]]@meta.data[[column]] <- as.factor(object_list[[x]]@meta.data[[column]])
    }
    return(object_list[[x]])
  })
  return(object_list)
}


#fix columns to factor
object_list <- lapply(1:length(object_list), FUN = function(x) {
    object_list[[1]]$Compositional_Cluster <- NULL
    return(object_list[[x]])
  })

integrated_obj$Compositional_Cluster <-  as.factor(integrated_obj$Compositional_Cluster)

object_list <- map_meta_integrated(integrated_obj, object_list, c("alias_id", "Compositional_Cluster"))

#calculate other clusters to annotate regions
integrated_obj <- integrated_obj %>% FindClusters(., resolution = c(0.05,0.1, 0.15))

# Rename the columns in the Seurat object
cluster_cols <- grep("^SCT_snn_res", colnames(integrated_obj@meta.data), value = TRUE)
new_colnames <- sub("^SCT_snn_res", "clustering_harmony", cluster_cols)

current_colnames <- colnames(integrated_obj@meta.data)
cols_to_change <- current_colnames %in% cluster_cols
colnames(integrated_obj@meta.data)[cols_to_change] <- new_colnames

object_list <- map_meta_integrated(integrated_obj, object_list, c("alias_id", 'clustering_harmony.0.05','clustering_harmony.0.1','clustering_harmony.0.15', "clustering_harmony.0.2"))

#looking at plots, we can see that the 0.05 and 0.1 resolutions are the best approach to annotate the regions
integrated_obj$annotation <- integrated_obj$clustering_harmony.0.1 %>% as.character() 
ids <- integrated_obj@meta.data %>% filter(clustering_harmony.0.05 =="3") %>% rownames()
integrated_obj@meta.data[ids,"annotation"] <- "5"
object_list <- map_meta_integrated(integrated_obj, object_list, c("alias_id", 'annotation'))

Idents(integrated_obj) <- "annotation"
integrated_obj <- RenameIdents(object = integrated_obj, 
                                 "0" = "Cortex",
                                 "1" = "OSOM",
                                 "2" = "ISOM",
                                 "3" = "Cortex",
                                 "4" = "Inner_Medulla",
                                 "5" = "Renal Pelvis")

integrated_obj$Region <- Idents(integrated_obj)
object_list <- map_meta_integrated(integrated_obj, object_list, c("alias_id", 'Region'))

saveRDS(integrated_obj,"./mdu/Data/Spatial/visium/kidney/mouse_NTN/processed/ultimate_processed/integrated/rds/integrated_annotated.rds")

lapply(1:length(sample_list), FUN = function(x) {
    saveRDS(object_list[[x]], paste0("./mdu/Data/Spatial/visium/kidney/mouse_NTN/processed/ultimate_processed/by_sample/rds/", x, "_annotated.rds"))
    })


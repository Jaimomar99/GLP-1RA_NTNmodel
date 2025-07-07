require(Libra)
require(Seurat)

DE_function <- function(seurat_object, celltype_column, cells, library_id_column, disease_status_column, cond1, cond2, de_model = "wilcox",min.pct= 0, min.logFC=0, pval_filter=0, assay="RNA"){
    DefaultAssay(seurat_object) <- assay
    Idents(seurat_object) <- celltype_column
    seurat_object <- subset(seurat_object, idents = cells)

    Idents(seurat_object) <- disease_status_column
    seurat_object <- subset(seurat_object, idents = c(cond1,cond2))
    seurat_object[[disease_status_column]] <- factor(pull(seurat_object[[disease_status_column]]), levels = c(cond1,cond2))
    message("cond1 vs cond2: ", levels(pull(seurat_object[[disease_status_column]])))

  if(assay== "SCT"){seurat_object <- PrepSCTFindMarkers(seurat_object)}
  
  message("Running DEA with model: ", de_model, " for cell types: ", cells, " and conditions: ", cond1, " vs ", cond2)
  if(de_model == "wilcox"){
    markers <- FindMarkers(seurat_object, 
                           logfc.threshold = min.logFC,
                           test.use = "wilcox",
                           ident.1 =  cond1,
                           ident.2 = cond2,
                           assay= assay,
                           min.pct=min.pct,
                           slot="data") %>% 
      tibble::rownames_to_column(var = "gene") %>% 
      dplyr::rename("avg_logFC" = avg_log2FC) %>%
      mutate(cell_type = cells) %>%
      mutate(model=de_model) %>%
      mutate(cond1= cond1) %>%
      mutate(cond2= cond2) %>%
      dplyr::select(cell_type, gene, avg_logFC, p_val, p_val_adj, pct.1, pct.2, model, cond1, cond2)
  }

 
  if(!is.null(pval_filter)){
    markers <- markers %>% 
      dplyr::filter(p_val_adj < pval_filter) 
  }
  markers <- markers %>%
      arrange(desc(avg_logFC))
  
  message("DEA done")

  return (markers)
}
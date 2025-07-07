suppressWarnings(library(cluster))
suppressWarnings(library(gridExtra))
suppressWarnings(library(scellpam))
suppressWarnings(library(R.utils))
suppressWarnings(library(Matrix))
suppressWarnings(library(jsonlite))


.get_legend <- function(a.gplot){
  # Function to extracts the legend from a ggplot object
  #' @param a.gplot A ggplot object from which the legend is to be extracted
  #' @return A gtable object containing the legend of the input ggplot

  # Build the ggplot object and convert it to a gtable
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  
  # Find the index of the legend ("guide-box") in the list of grobs
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  
  # Extract the legend grob using the index
  legend <- tmp$grobs[[leg]]
  
  # Return the legend
  return(legend)
}

QC_ST <- function(seurat_obj, nCount = 500, nFeature=250, threshold_mt= 30, times_genes = 10){
  # Function to performs quality control on a Seurat object
  #' @param seurat_obj A Seurat object to be processed
  #' @param nCount Minimum number of counts for a cell to be considered valid
  #' @param nFeature Minimum number of features for a cell to be considered valid
  #' @param threshold_mt Maximum percentage of mitochondrial genes for a cell to be considered valid
  #' @param times_genes Minimum number of cells in which a gene must be present to be considered valid
  #' @return A Seurat object after quality control
    # Calculate the total number of spots before QC
    seurat_obj$n_spots_before_QC <- Cells(seurat_obj) %>% length() %>% as.factor()
    
    # Calculate the percentage of mitochondrial genes
    total_counts_per_cell <- colSums(seurat_obj@assays$Spatial@counts)
    mito_genes <- rownames(seurat_obj)[grep("^MT-", rownames(seurat_obj))]
    seurat_obj$percent_mito <- colSums(seurat_obj@assays$Spatial@counts[mito_genes, ])/total_counts_per_cell

    # Filtering
    # Spot level
    cells <- WhichCells(seurat_obj, expression = nCount_Spatial > nCount &
                          nFeature_Spatial > nFeature &
                          percent_mito < threshold_mt)
    # Gene level
    genes <- rownames(seurat_obj)[Matrix::rowSums(seurat_obj) > times_genes]
    
    # Subset the Seurat object based on the filtering criteria
    seurat_obj <- subset(seurat_obj, features = genes, cells = cells)
    
    # Calculate the total number of spots after QC
    seurat_obj$n_spots_after_QC <- Cells(seurat_obj) %>% length() %>% as.factor()

    # Calculate the percentage of spots lost during QC
    pct <- (as.numeric(levels(seurat_obj$n_spots_before_QC)) - as.numeric(levels(seurat_obj$n_spots_after_QC))) / as.numeric(levels(seurat_obj$n_spots_before_QC)) * 100 
    seurat_obj$pct_spots_lost <- round(pct,3) %>% as.factor()

  return(seurat_obj)
}


filter_spatial <- function(seurat_obj, min_spot_fraction=0.01, min_count_quantile=0.02, min_gene_count_quantile=0.02, max_mt_quantile=0.98, assay = "Spatial"){
    #' Apply spot and gene filters 
    #' @param seurat_obj Seurat object 
    #' @param min_spot_fraction Filter genes that are detected in less than this fraction of spots
    #' @param min_count_quantile Use quantile to calcualte minimum spot total UMI count threshold; ex. if 2%
    #' of spots have total count < 100, remove spots that fall below this 100 count threshold
    #' @param min_gene_count_quantile Use quantile to calcualte minimum number of genes detected per spot; filter
    #' out spots that have less genes that this threshold
    #' @returns seurat_obj with spot and gene filters applied

    n_spot <- Cells(seurat_obj) %>% length()
    min_spot <- round(n_spot * min_spot_fraction)

    keep.genes <- rowSums(GetAssayData(seurat_obj, assay = assay, slot = "counts") > 0) >= min_spot

    count_threshold = quantile(seurat_obj$nCount_Spatial, min_count_quantile)
    gene_count_threshold = quantile(seurat_obj$nFeature_Spatial, min_count_quantile)
    percent_mito_threshold = 0

    # check if there are mito genes (MT-, mt-, or Mt-)
    mito_genes <- rownames(seurat_obj)[grep("^MT-", rownames(seurat_obj), ignore.case = TRUE)]
    if (length(mito_genes) > 0) {
        print("Filtering based on percent mitochondrial transcripts")
        total_counts_per_cell <- colSums(seurat_obj@assays$Spatial@counts)
        seurat_obj$percent_mito <- colSums(seurat_obj@assays$Spatial@counts[mito_genes, ]) / total_counts_per_cell
        print("percent mito NA:")
        print(sum(is.na(seurat_obj$percent_mito)))
        percent_mito_threshold = quantile(seurat_obj$percent_mito, max_mt_quantile, na.rm = TRUE)

        # keep spots that have minimum # of total counts and minimum number of genes detected
        keep.spots <- (seurat_obj[[]][, 'nCount_Spatial'] >= count_threshold) & (seurat_obj[[]][, 'nFeature_Spatial'] >= gene_count_threshold) & (seurat_obj[[]][, 'percent_mito'] <= percent_mito_threshold)
    } else {
        seurat_obj$percent_mito <- NA
        keep.spots <- (seurat_obj[[]][, 'nCount_Spatial'] >= count_threshold) & (seurat_obj[[]][, 'nFeature_Spatial'] >= gene_count_threshold)
    }

    # apply filter
    seurat_obj <- seurat_obj[keep.genes, keep.spots]

    return(seurat_obj)
}


Normalization_clustering <- function(seurat_obj, resolution_range = c(0.2, 0.3, 0.4, 0.6, 0.8), reduction ="pca", sample=TRUE){
  # Function performs normalization and clustering on a Seurat object.
  #' @param seurat_obj A Seurat object to be normalized and clustered.
  #' @param resolution_range A numeric vector specifying the resolutions at which to find clusters.
  #' @param reduction The dimensionality reduction technique to be used. Default is "pca".
  #' @param sample A logical value indicating whether to perform sampling. Default is TRUE.
  #' @return A Seurat object that has been normalized and clustered.

  # If sampling is TRUE, perform SCTransform normalization and PCA
  if (sample){
    seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = TRUE, vst.flavor = "vst2", return.only.var.genes = FALSE)
    DefaultAssay(seurat_obj) <- "SCT"
    seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  }
  
  # Determine the number of principal components to use
  pcs <- .Elbow_plot(seurat_obj)
  
  # Perform UMAP dimensionality reduction, find neighbors, and find clusters at specified resolutions
  seurat_obj <- RunUMAP(seurat_obj, reduction = reduction, assay = "SCT", dims = 1:pcs)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:pcs)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution_range) 
  
  # Set the identity of the Seurat object to the clustering result
  # best_res <- silhouette_score(seurat_obj, resolution_range, pcs, embeddings = reduction)
  # Idents(seurat_obj) <- best_res
  Idents(seurat_obj) <- "SCT_snn_res.0.2"

  # If sampling is TRUE, assign the clustering result to the 'clustering_by_sample' column, else assign it to the 'clustering_harmony' column

    # Get the column names that start with "SCT_snn_res"
    cluster_cols <- grep("^SCT_snn_res", colnames(seurat_obj@meta.data), value = TRUE)

  if (sample){
    seurat_obj$clustering_by_sample <- Idents(seurat_obj)

    new_colnames <- sub("^SCT_snn_res", "clustering_by_sample", cluster_cols)
    # Rename the columns in the Seurat object
    current_colnames <- colnames(seurat_obj@meta.data)
    cols_to_change <- current_colnames %in% cluster_cols
    colnames(seurat_obj@meta.data)[cols_to_change] <- new_colnames

  } else {
    seurat_obj$clustering_harmony <- Idents(seurat_obj)

    new_colnames <- sub("^SCT_snn_res", "clustering_harmony", cluster_cols)
    # Rename the columns in the Seurat object
    current_colnames <- colnames(seurat_obj@meta.data)
    cols_to_change <- current_colnames %in% cluster_cols
    colnames(seurat_obj@meta.data)[cols_to_change] <- new_colnames
  }

  # Remove the 'seurat_clusters' column from the Seurat object
  seurat_obj$seurat_clusters <- NULL

  # Return the Seurat object
  return(seurat_obj)
}



silhouette_score <- function(object, resolution_range, pcs, embeddings, plot = FALSE){
  # Code: https://github.com/satijalab/Integration2019/blob/master/analysis_code/integration/integration_metrics.R#L36
  
  distance_matrix <- dist(Embeddings(object[[embeddings]])[, pcs])
  
  silhouette_scores_list <- list()
  
  for (i in seq_along(resolution_range)) {
    res <- paste0("SCT_snn_res.", resolution_range[1])
    print(res)
    clusters <- object@meta.data[res] %>% pull()
    silhouette_scores_list[[i]] <- CalculateSilhouette(as.numeric(clusters), fdist = distance_matrix, nthreads = 16) %>%
      as_tibble() %>%
      mutate(res=resolution_range[i])
#     silhouette_scores_list[[i]] <- silhouette(as.numeric(clusters), dist = distance_matrix) %>%
#       as_tibble() %>%
#       mutate(res=resolution_range[i])
  }
  #plot
  if (plot){
      print(do.call(rbind, silhouette_scores_list) %>% 
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
  )
  } 
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
    best_res <- paste0("SCT_snn_res.", resolution_range[max])
    return(best_res)
}

QC_plot <- function(harmonized_ST, subject_id_column, group_column, count_column, feature_column, percent.mt = NULL, spots=NULL){
  
  #get common features
  vars_to_fetch <- c(subject_id_column, group_column, count_column, feature_column)
  column_names <- c("subject_id", "group_column", "nCount", "nFeature")
  
  if (!is.null(percent.mt)) {
    vars_to_fetch <- c(vars_to_fetch, percent.mt)
    column_names <- c(column_names, "percent.mt")
  }
  
  if (!is.null(spots)) {
    vars_to_fetch <- c(vars_to_fetch, spots)
    column_names <- c(column_names, "spots")
  }
  
  meta <- FetchData(harmonized_ST, vars = vars_to_fetch)
  meta <- setNames(meta, column_names)
  
  n_subjects <- meta$subject_id %>% unique() %>% length()
  n_classes <- meta$group_column %>% unique() %>% length()
  
  #start plotting
  plots <- list()
  
  plots[["p1"]] <- ggplot(meta, aes(x = subject_id, color = group_column, y = nCount, fill=group_column)) + 
    geom_boxplot(alpha = 0.2) +
    ggtitle("Transcripts count") +
    labs(fill = group_column, color=group_column) + 
    theme(legend.position="bottom") # +
    # scale_color_NN() +
    # scale_fill_NN()
  
  plots[["p2"]] <- ggplot(meta, aes(x = subject_id, color = group_column, y = nFeature, fill=group_column)) + 
    geom_boxplot(alpha = 0.2) + 
    ggtitle("Number of Features ") 
  
  if (!is.null(percent.mt)) {
    plots[["p3"]] <- ggplot(meta, aes(x = subject_id, color = group_column, y = percent.mt, fill=group_column)) + 
      geom_boxplot(alpha = 0.2) + 
      ggtitle("Mitochondrial percent")
  }
  
  if (!is.null(spots)) {
    plots[["p4"]] <- ggplot(meta, aes(x=subject_id, color=group_column, y=spots, fill=group_column)) +
      geom_col(width=0.5) +
      ggtitle("Number spots")
  }
  
  mylegend <- .get_legend(plots[[1]])
  
  plots <- lapply(plots, function(plot) {
    if (n_subjects > 50) {
      plot <- plot + theme(axis.text.x = element_text(size = 4))
    }
    
    plot <- plot + 
      scale_x_discrete(labels = as.character(1:n_subjects))  +
      theme(
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none"
      ) # +
      # scale_color_NN() +
      # scale_fill_NN()
    return(plot)
  })
  
  final_plot <-  grid.arrange(arrangeGrob(grobs= plots,
                                          ncol=1),
                              mylegend, nrow=2,
                              heights = unit(c(14, 2), c("null", "cm")))
  return(final_plot)
}

.Elbow_plot <-  function(object, plot = FALSE){
  # Function generates an elbow plot for a given object and returns the number of principal components (PCs).
  #' @param object An object for which the elbow plot is to be generated. It should contain PCA results.
  #' @param plot A logical value indicating whether to generate the elbow plot. Default is FALSE.
  #' @return The number of principal components (PCs) to be used.

  # Calculate the percentage of standard deviation contributed by each PC
  pct <- object[["pca"]]@stdev / sum(object[["pca"]]@stdev) * 100
  
  # Calculate the cumulative sum of the percentages
  cumu <- cumsum(pct)
  
  # Find the first PC where the cumulative sum is greater than 90% and the percentage is less than 5%
  co1 <- which(cumu > 90 & pct < 5)[1]
  
  # Find the first PC where the percent change in variation between consecutive PCs is greater than 0.1%
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  
  # Choose the smaller of the two PCs
  pcs <- min(co1, co2)
  
  # If plot is TRUE, generate the elbow plot
  if (plot){
    print(data.frame(pct = pct[1:25], 
                     cumu = cumu[1:25], 
                     rank = 1:length(pct[1:25])) %>%
        ggplot(aes(cumu, pct, label = rank, color = rank > pcs)) + 
        geom_text(size=5) +
        ggtitle("Elbow plot") +
        theme(legend.position = "none",
              plot.title = element_text(size = 12)))
  }

  # Return the number of PCs
  return(pcs)
}

plot_raw_data <- function(list_ST, title_size=9){
  # Function generates spatial dimension plots for a list of Seurat objects.
  #' @param list_ST A list of Seurat objects for which the spatial dimension plots are to be generated.
  #' @param title_size The size of the title text in the plots. Default is 9.
  #' @return A grid of spatial dimension plots for the input Seurat objects.

  # Get the layout for the grid of plots
  layout <- .get_layout(list_ST)
  
  # Generate a list of spatial dimension plots for each Seurat object in the input list
  raw_plot <- lapply(1:length(list_ST), function(i) {
    SpatialDimPlot(list_ST[[i]],  crop=FALSE, pt.size.factor = 1, alpha=c(0,0)) +
      NoLegend() +
      ggtitle(names(list_ST[i])) +
      theme(plot.title = element_text(size = title_size))
  })
  
  # Set equal width for all columns and equal height for all rows
  col_widths <- rep(1, layout$n_cols)  
  row_heights <- rep(1, layout$n_rows)  
  
  # Arrange the plots in a grid
  p <- grid.arrange(grobs = raw_plot, ncol = layout$n_cols, nrow = layout$n_rows,
               widths = col_widths, heights = row_heights)
  
  # Return the grid of plots
  return(p)
}

Harmony_integration <- function(list_ST, batch_cond_vector = c("seq_sample_id"), apply_filter=TRUE){
  # Function performs Harmony integration on a list of Seurat objects.
  #' @param list_ST A list of Seurat objects to be integrated.
  #' @param batch_cond_vector A character vector specifying the batch conditions to be used for integration. Default is "seq_sample_id".
  #' @param apply_filter A logical value indicating whether to apply filtering. Default is TRUE.
  #' @return A Seurat object that has been integrated using Harmony.

  # Select features for integration
  integ_features <- SelectIntegrationFeatures(object.list = list_ST, nfeatures = 3000) 
  
  # Merge the Seurat objects in the input list
  merge_ST <- merge(x = list_ST[[1]],
                    y = list_ST[2:length(list_ST)],
                    merge.data = TRUE)

  # Set the default assay to "SCT"
  DefaultAssay(merge_ST) <- "SCT"

  # Set the variable features of the merged Seurat object
  VariableFeatures(merge_ST) <- integ_features

  # Perform PCA on the merged Seurat object
  merge_ST <- RunPCA(merge_ST, assay = "SCT")

  # Set the lambda parameter for Harmony integration
  lambda_list = rep(1, length(batch_cond_vector)) 

  # Perform Harmony integration on the merged Seurat object
  harmonized_ST <- RunHarmony(merge_ST, 
                              group.by.vars = batch_cond_vector, lambda = lambda_list, 
                              reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
  
  # Return the integrated Seurat object
  return(harmonized_ST)
}

harmony_checkpoint_plot <- function(harmonized_ST, batch_cond_vector = c("orig.ident")){
  # Function generates dimension plots for a Harmony-integrated Seurat object, grouped by specified batch conditions.
  #' @param harmonized_ST A Seurat object that has been integrated using Harmony.
  #' @param batch_cond_vector A character vector specifying the batch conditions to be used for grouping in the dimension plots. Default is "orig.ident".
  #' @return A grid of dimension plots for the input Seurat object, grouped by the specified batch conditions.

  # Initialize a list to store the dimension plots
  list_plots <- list()
  
  # Generate a dimension plot for each batch condition in the input vector
  for (variable in batch_cond_vector) {
    list_plots[[variable]] <- DimPlot(object = harmonized_ST, reduction = "harmony", pt.size = .1, group.by = variable) +
      ggtitle(variable) +
      theme_bw() +
      theme(legend.title = element_blank(),
            legend.text = element_text(size = 8),
            legend.position = "bottom", 
            axis.title = element_blank())
  }
  
  # Get the layout for the grid of plots
  layout <- .get_layout(list_plots)
  
  # Set equal width for all columns and equal height for all rows
  col_widths <- rep(1, layout$n_cols)  
  row_heights <- rep(1, layout$n_rows)  
  
  # Arrange the plots in a grid
  grid.arrange(grobs = list_plots, ncol = layout$n_cols, nrow = layout$n_rows,
               widths = col_widths, heights = row_heights)
}

.get_layout <- function(list_ST) {
  # Function calculates the layout for a grid of plots based on the number of Seurat objects in a list.
  #' @param list_ST A list of Seurat objects for which the grid layout is to be calculated.
  #' @return A list containing the number of columns and rows for the grid layout.

  # Calculate the number of Seurat objects in the list
  n <- length(list_ST)
  
  # Calculate the number of columns as the square root of the number of Seurat objects
  columns <- as.integer(sqrt(n))
  
  # Calculate the number of rows as the ceiling of the number of Seurat objects divided by the number of columns
  rows <- as.integer(ceiling(n / columns))
  
  # Return the layout as a list
  return(list(n_cols = columns, n_rows = rows))
}


DE_function <- function(seurat_object, cluster_column, assay, min.pct= 0.2, ntop = 20){
    # Function finds all markers (differentially expressed genes) for each cluster in a Seurat object.
    #' @param seurat_object A Seurat object that contains the data to be analyzed.
    #' @param cluster_column A character string that specifies the column in the Seurat object that contains the cluster identities.
    #' @param assay A character string that specifies the assay to be used.
    #' @param min.pct A numeric value that specifies the minimum percentage of cells where a gene must be detected in order to be considered. Default is 0.2.
    #' @param ntop A numeric value that specifies the number of top markers to return per cluster. Default is 20.
    #' @return A data frame containing the top markers for each cluster.

    # Set the cluster identities in the Seurat object
    Idents(seurat_object) <- cluster_column
  
    # Find all markers for each cluster in the Seurat object
    markers <- FindAllMarkers(seurat_object, 
                              logfc.threshold = 0.25,
                              test.use = "wilcox",
                              assay = assay,
                              min.pct = min.pct,
                              slot = "data") 

    # Filter and sort the markers
    markers <- markers %>%  
        mutate(cluster = as.integer(as.character(cluster)), na.rm=TRUE) %>% 
        select(cluster, gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2) %>%
        dplyr::filter(p_val_adj < 0.05) %>%
        dplyr::filter(abs(avg_log2FC) > 0.25) %>%
        arrange(cluster, desc(avg_log2FC)) %>% 
        mutate_if(is.numeric,round, 3) %>%
        group_by(cluster) %>% 
        top_n(n=ntop, wt=avg_log2FC)

    # Return the top markers for each cluster
    return (markers)
}


infer_CC <- function(integrated_obj, cell_matrix){
  # Function infers compositional clusters (CC) in a Seurat object using k-means clustering.
  #' @param integrated_obj A Seurat object to which the inferred clusters are to be added.
  #' @param cell_matrix A matrix of cell data to be used for k-means clustering.
  #' @return A Seurat object with the inferred clusters added as metadata.

  # Initialize a vector to store the within-cluster sum of squares (WCSS) for each number of clusters
  wcss <- NULL
  
  # Perform k-means clustering for 1 to 10 clusters and calculate the WCSS for each. 
  #Watch out k-means assign different CC each time, order would be different
  for (i in 1:15) {
    sample.fit = kmeans(cell_matrix, centers = i)
    wcss <- c(wcss, sample.fit$tot.withinss)
  }
  
  # Calculate the difference in WCSS between consecutive numbers of clusters
  diff_wcss <- diff(wcss)
  
  # Calculate the cumulative sum of the WCSS
  cumulative_wcss <- cumsum(wcss)
  
  # Determine the optimal number of clusters as the smallest number for which the difference in cumulative WCSS is less than the mean difference
  optimal_k <- which(diff(cumulative_wcss) < mean(diff(cumulative_wcss))) %>% min()
  
  # Perform k-means clustering with the optimal number of clusters
  comp_clusters  <- kmeans(cell_matrix, optimal_k)
  
  # Add a "CC" prefix to the cluster labels
  comp_clusters$cluster <- paste("CC", comp_clusters$cluster, sep = "")
  
  # Add the inferred clusters to the Seurat object as metadata
  integrated_obj <- AddMetaData(integrated_obj, comp_clusters$cluster, "Compositional_Cluster")
  
  # Return the Seurat object with the added clusters
  return(integrated_obj)
}


Composition_cluster_enrichedCelltypes <- function(Spatial.object, abund_df){
  # Function generates a boxplot of cell type proportions for a Seurat object.
  #' @param Spatial.object A Seurat object for which the boxplot is to be generated.
  #' @param abund_df A data frame containing the cell type proportions for each spot in the Seurat object.
  #' @return A ggplot object representing the boxplot of cell type proportions.

  # Suppress warnings from the matrixStats and reshape2 packages
  suppressWarnings(library(matrixStats))
  suppressWarnings(library(reshape2))
  
  # Convert the cell type proportions to a matrix
  COI.topic.probs <- abund_df %>% as.matrix()
  
  # Calculate the median cell type proportion for each cell type and sort in decreasing order
  COI.topic.probs.medians <- sort(colMedians(COI.topic.probs),decreasing = T)
  
  # Get the names of the two cell types with the highest median proportions
  COI.leading.topics <- names(COI.topic.probs.medians[1:2])
  
  # Melt the matrix of cell type proportions into a data frame
  COI.topic.probs.melt <- melt(COI.topic.probs)
  colnames(COI.topic.probs.melt) <- c("spot_id", "Celltype", "prob")
  
  # Create a data frame of cell type proportions
  Topic.prob.all.df <- COI.topic.probs.melt
  Topic.prob.all.df$Topic <- factor(Topic.prob.all.df$Celltype)
  
  # Generate a boxplot of cell type proportions
  p <- ggplot(Topic.prob.all.df, aes(x=Topic, y=prob, fill=Topic)) + 
    geom_boxplot(show.legend = FALSE) + 
    theme_bw() + 
    theme (
      axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1,size = 10, face = "bold", colour = "black"),
      axis.text.y = element_text( hjust = 1, size = 17, face = "bold", colour = "black"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")) + 
    ggtitle("Deconvolution proportions")
  
  # Return the boxplot
  return(p)
}


.BlendMap <- function(color.matrix, xlab = NULL, ylab = NULL) {
  # BlendMap => copied from Seurat for use in custom spatial blended plots
  color.heat <- matrix(
    data = 1:prod(dim(x = color.matrix)) - 1,
    nrow = nrow(x = color.matrix),
    ncol = ncol(x = color.matrix),
    dimnames = list(
      1:nrow(x = color.matrix),
      1:ncol(x = color.matrix)
    )
  )
  xbreaks <- seq.int(from = 0, to = nrow(x = color.matrix), by = 2)
  ybreaks <- seq.int(from = 0, to = ncol(x = color.matrix), by = 2)
  color.heat <- .Melt(x = color.heat)
  color.heat$rows <- as.numeric(x = as.character(x = color.heat$rows))
  color.heat$cols <- as.numeric(x = as.character(x = color.heat$cols))
  color.heat$vals <- factor(x = color.heat$vals)
  plot <- ggplot(
    data = color.heat,
    mapping = aes_string(x = 'rows', y = 'cols', fill = 'vals')
  ) +
    geom_raster(show.legend = FALSE) +
    theme(plot.margin = unit(x = rep.int(x = 0, times = 4), units = 'cm')) +
    scale_x_continuous(breaks = xbreaks, expand = c(0, 0), labels = xbreaks) +
    scale_y_continuous(breaks = ybreaks, expand = c(0, 0), labels = ybreaks) +
    scale_fill_manual(values = as.vector(x = color.matrix)) +
    theme_cowplot()

    if (!is.null(xlab)) {
        plot <- plot + xlab(xlab)
    }
    if (!is.null(ylab)) {
        plot <- plot + ylab(ylab)
    }
  return(plot)
}

.Melt <- function(x) {
  # .Melt => copied from Seurat for use in custom spatial blended plots
  if (!is.data.frame(x = x)) {
    x <- as.data.frame(x = x)
  }
  return(data.frame(
    rows = rep.int(x = rownames(x = x), times = ncol(x = x)),
    cols = unlist(x = lapply(X = colnames(x = x), FUN = rep.int, times = nrow(x = x))),
    vals = unlist(x = x, use.names = FALSE)
  ))
}

.BlendMatrix <- function(
  # BlendMatrix => copied from Seurat for use in custom spatial blended plots
  n = 10,
  col.threshold = 0.5,
  two.colors = c("#ff0000", "#00ff00"),
  negative.color = "black"
  ) {
  if (0 > col.threshold || col.threshold > 1) {
    stop("col.threshold must be between 0 and 1")
  }
  C0 <- as.vector(col2rgb(negative.color, alpha = TRUE))
  C1 <- as.vector(col2rgb(two.colors[1], alpha = TRUE))
  C2 <- as.vector(col2rgb(two.colors[2], alpha = TRUE))
  blend_alpha <- (C1[4] + C2[4])/2
  C0 <- C0[-4]
  C1 <- C1[-4]
  C2 <- C2[-4]
  merge.weight <- min(255 / (C1 + C2 +  C0 + 0.01))
  sigmoid <- function(x) {
    return(1 / (1 + exp(-x)))
  }
  blend_color <- function(
    i,
    j,
    col.threshold,
    n,
    C0,
    C1,
    C2,
    alpha,
    merge.weight
  ) {
    c.min <- sigmoid(5 * (1 / n - col.threshold))
    c.max <- sigmoid(5 * (1 - col.threshold))
    c1_weight <- sigmoid(5 * (i / n - col.threshold))
    c2_weight <- sigmoid(5 * (j / n - col.threshold))
    c0_weight <-  sigmoid(5 * ((i + j) / (2 * n) - col.threshold))
    c1_weight <- (c1_weight - c.min) / (c.max - c.min)
    c2_weight <- (c2_weight - c.min) / (c.max - c.min)
    c0_weight <- (c0_weight - c.min) / (c.max - c.min)
    C1_length <- sqrt(sum((C1 - C0) ** 2))
    C2_length <- sqrt(sum((C2 - C0) ** 2))
    C1_unit <- (C1 - C0) / C1_length
    C2_unit <- (C2 - C0) / C2_length
    C1_weight <- C1_unit * c1_weight
    C2_weight <- C2_unit * c2_weight
    C_blend <- C1_weight * (i - 1) * C1_length / (n - 1) + C2_weight * (j - 1) * C2_length / (n - 1) + (i - 1) * (j - 1) * c0_weight * C0 / (n - 1) ** 2 + C0
    C_blend[C_blend > 255] <- 255
    C_blend[C_blend < 0] <- 0
    return(rgb(
      red = C_blend[1],
      green = C_blend[2],
      blue = C_blend[3],
      alpha = alpha,
      maxColorValue = 255
    ))
  }
  blend_matrix <- matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      blend_matrix[i, j] <- blend_color(
        i = i,
        j = j,
        col.threshold = col.threshold,
        n = n,
        C0 = C0,
        C1 = C1,
        C2 = C2,
        alpha = blend_alpha,
        merge.weight = merge.weight
      )
    }
  }
  return(blend_matrix)
}


spatial_blended_plot <- function(
  # Custom blended spatial plot
    object,
    features,
    blend.threshold = 0.5,
    pt.size = 1,
    cols = c('lightgrey', '#ff0000', '#0400ff')
    ) {

    # cowplot dependency

    ncol <- 4
    color.matrix <- .BlendMatrix(
      two.colors = cols[2:3],
      col.threshold = blend.threshold,
      negative.color = cols[1]
    )
    cols <- cols[2:3]
    colors <- list(
      color.matrix[, 1],
      color.matrix[1, ],
      as.vector(x = color.matrix)
    )

    sample_id <- names(object@images)

    data.plot <- GetAssayData(object, slot = "data", assay = "Spatial")

    data.plot <- as.data.frame(t(as.matrix((data.plot[features, ]))))
    
    data.plot <- as.data.frame(x = apply(
        X = data.plot,
        MARGIN = 2,
        FUN = function(x) {
          return(round(x = 9 * (x - min(x)) / (max(x) - min(x))))
        }
    ))
    
    data.plot[, 3] <- data.plot[, 1] + data.plot[, 2] * 10
    
    c(features, paste(features, collapse = '_'))
    
    colnames(x = data.plot) <- c(features, paste(features, collapse = '_'))

    for (i in 1:ncol(x = data.plot)) {
        data.plot[, i] <- factor(x = data.plot[, i])
    }
    data.plot <- cbind(object@images[[sample_id]]@coordinates[, c("row", "col")], data.plot)

    features <- colnames(x = data.plot)[3:ncol(x = data.plot)]
    gene1 <- features[1]
    gene2 <- features[2]
    
    plot_list <- lapply(1:length(features), FUN = function(j) {
        feature <- features[j]
        # print(feature)
    
        cols.use <- as.numeric(x = as.character(x = data.plot[, feature])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
        p <- ggplot(data.plot, aes(x = row, y = col, color = .data[[feature]])) +
            geom_point(size = pt.size) +
            scale_color_manual(values = cols.use) +
            theme_cowplot() +
            theme(axis.line=element_blank(),
            	axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                legend.position="none",
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank()
		) +
            ggtitle(feature)
        }
    )
    
    blend.legend <- .BlendMap(color.matrix = color.matrix, xlab = gene1, ylab = gene2)
    plot_list[[length(plot_list) + 1]] <- blend.legend

    final_plot <- wrap_plots(plot_list, nrow = 1)
    return(final_plot)
}

.extract_image <- function (x){
  # Function extracts the first image from a Seurat object and converts it to a matrix.
  #' @param x A Seurat object from which the image is to be extracted.
  #' @return A matrix representing the extracted image.

  # Get the first/only slice from the Seurat object
  slice <- SeuratObject::Images(x)[1]
  
  # Extract the image corresponding to the slice and convert it to a raster
  x <- SeuratObject::GetImage(x, slice=slice, mode = "raster")
  
  # Convert the raster to a matrix
  x <- as.matrix(x)
  
  # Return the matrix
  return(x)
}

.plot_image <- function(x, alpha=1) {
  # Function generates a ggplot of an image represented as a matrix of hexadecimal color values.
  #' @param x A matrix of hexadecimal color values representing an image.
  #' @param alpha The alpha transparency level to be applied to the colors. Default is 1 (no transparency).
  #' @return A ggplot object representing the image.

  # Convert the hexadecimal color values to RGB and apply the alpha transparency level
  x_rgb_alpha <- t(apply(x, c(1, 2), function(y) {
    rgb_val <- col2rgb(y) / 255
    rgb(rgb_val[1], rgb_val[2], rgb_val[3], alpha = alpha)
  }))

  # Convert the matrix of RGB color values to a raster grob
  x <- grid::rasterGrob(t(x_rgb_alpha),
    interpolate = FALSE,
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc"))
  
  # Generate a ggplot of the raster grob
  p <- ggplot() +
    annotation_custom(
      grob = x,
      xmin = 0,
      xmax = ncol(x$raster),
      ymin = 0,
      ymax = nrow(x$raster)) + 
    coord_fixed(
      xlim = c(0, ncol(x$raster)),
      ylim = c(0, nrow(x$raster))) + 
    theme_void()
  
  # Return the ggplot
  return(p)
}

.extract_coord <- function(x, img) {
    # Function extracts the spatial coordinates from a Seurat object.
    #' @param x A Seurat object from which the spatial coordinates are to be extracted.
    #' @param img An optional parameter, not used in the current function but can be used for future enhancements.
    #' @return A matrix containing the spatial coordinates of the input Seurat object.

    # Check if the Seurat object has any images, if not, stop the function
    stopifnot(
        # Stop if there are no images
        !is.null(SeuratObject::Images(x))
    )
    
    # Get the first image slice from the Seurat object
    slice <- SeuratObject::Images(x)[1]

    # Extract spatial coordinates from the Seurat object for the given image slice
    x <- as.matrix(SeuratObject::GetTissueCoordinates(x, image = slice))

    # Define the column names for the coordinates
    cnames <- c("coord_y", "coord_x")
    
    # If the column names of the extracted coordinates do not match the defined names, rename them
    if (!all(colnames(x) %in% cnames)) {
        colnames(x) <- cnames
    }

    # Return the coordinates as a matrix
    return(as.matrix(x))
}


plotSpatialScatterpie <- function (x, abund_df, img_alpha = 0.5, crop=FALSE, scatterpie_alpha = 1, pie_scale = 0.4, ...){
    # This function plots a scatterpie plot using spatial coordinates from a Seurat object.
    #' @param x A Seurat object from which the spatial coordinates are to be extracted.
    #' @param abund_df A dataframe with the deconvolution proportion. Rownames True
    #' @param img_alpha A numeric value specifying the transparency of the image. Default is 0.5.
    #' @param crop A boolean value, specifying if return also the zoom in image together with the normal. Default is FALSE
    #' @param scatterpie_alpha A numeric value specifying the transparency of the scatterpie. Default is 1.
    #' @param pie_scale A numeric value specifying the size of the pies in the scatterpie plot. Default is 0.4.
    #' @return A ggplot object representing the scatterpie plot.

    #extract cell types
    cell_types <- colnames(abund_df)

    # Extract the cell type data from the Seurat object
    y <- x@meta.data[, cell_types] %>% as.matrix()
    
    # Extract the spatial coordinates from the Seurat object
    coords <- .extract_coord(x = x)
    
    # Extract the image from the Seurat object
    img <- .extract_image(x)

    # Plot the image with specified transparency
    p <- .plot_image(x = img, alpha = img_alpha) 
    ymax <- max(p$coordinates$limits$y)
  
    # Merge the coordinates and cell type data
    df <- merge(coords, y, by = 0, all = TRUE)
    
    # Convert coord_y to absolute value
    df$coord_y <- abs(df$coord_y - ymax)

    # Add the scatterpie plot to the image plot
    p1 <- p + scatterpie::geom_scatterpie(data = df, aes(x = coord_x, 
        y = coord_y), cols = cell_types, color = NA, 
        alpha = scatterpie_alpha, pie_scale = pie_scale, ...) + 
        theme(legend.key.size = unit(1.5, "lines"),
              plot.background = element_rect(fill = "white", colour="white"),
              panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank()) +
        labs(fill = "Cell Type")

    #add also the zoom in image
    if (crop){
        library(patchwork)
        xmin <- min(df$coord_x) 
        xmax <- max(df$coord_x)
        ymin <- min(df$coord_y)
        ymax <- max(df$coord_y)

        p1 <- p1 + theme(legend.position="none")
        p2 <- p + scatterpie::geom_scatterpie(data = df, aes(x = coord_x, 
        y = coord_y), cols = cell_types, color = NA, 
        alpha = scatterpie_alpha, pie_scale = pie_scale) + 
        theme(legend.key.size = unit(1.5, "lines"),
              plot.background = element_rect(fill = "white", colour="white"),
              panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank()) +
        labs(fill = "Cell Type") +
        coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax))
        # return(wrap_plots(p1, p2, ncol=2))
        return(p2)
    } else{
        return(p1)
    }
}

plotUMAPScatterpie <- function (sample, abund_df, reduction= "umap", scatterpie_alpha = 0.8, pie_scale = 0.3, ...){
    # This function plots a scatterpie plot using UMAP embeddings from a Seurat object.
    #' @param sample A Seurat object from which the UMAP embeddings are to be extracted.
    #' @param abund_df A dataframe with the deconvolution proportion. Rownames True
    #' @param scatterpie_alpha A numeric value specifying the transparency of the scatterpie. Default is 0.8.
    #' @param pie_scale A numeric value specifying the size of the pies in the scatterpie plot. Default is 0.3.
    #' @return A ggplot object representing the scatterpie plot.

    # Check if the input is a Seurat object and if the alpha values and pie scale are numeric
    stopifnot(is(sample, "Seurat"),
              is.numeric(scatterpie_alpha), 
              is.numeric(pie_scale))

    #extract cell types
    cell_types <- colnames(abund_df)

    # Extract the cell type data from the Seurat object
    y <- sample@meta.data[, cell_types] %>% as.matrix()

    # Extract the UMAP embeddings from the Seurat object
    embeddings <- Embeddings(sample, reduction = reduction)[, c(1,2)]

    colnames(embeddings) <- c("coord_x", "coord_y" )  

    # Merge the UMAP embeddings and cell type data
    df <- merge(embeddings, y, by = 0, all = TRUE)

    # Plot the scatterpie
    p <- ggplot() +
        scatterpie::geom_scatterpie(data = df, aes(x = coord_x, 
        y = coord_y), cols = cell_types, color = NA, 
        alpha = scatterpie_alpha, pie_scale = pie_scale, ...) + 
        theme_classic() +
        theme(legend.key.size = unit(1.5, "lines"),
              plot.background = element_rect(fill = "white", colour="white"),
              panel.grid.major = element_blank(), # Remove major grid
              panel.grid.minor = element_blank(), # Remove minor grid
              panel.border = element_blank(),
              axis.line = element_line(color = "black"),
              axis.title = element_text(size = 14), # Increase axis title size
              axis.text = element_text(size = 12) # Increase axis label size
              ) +
        labs(fill = "Cell Type", x = paste0(toupper(reduction), "_1"), y = paste0(toupper(reduction), "_2"))
    
    return(p)
}

.write_seurat_to_disk <- function(seurat_obj, path, feature_df_path = NULL, assay = "Spatial", slot = "counts") {
  # Most of the code from:https://rdrr.io/github/daskelly/earlycross/src/R/Write10X.R
  # This function saves a barcodes, features and matrix of a seurat object, sparse matrix.
  #' @param seurat_obj A Seurat object that contains the data to be saved.
  #' @param path A character string that specifies the directory where the files will be saved.
  #' @param feature_df_path Path to tsv file with feature ids (ex. gene symbol + ensembl ids)
  #' @assay assay to take the matrix. Default Spatial
  #' @slot slot to take the matrix. Default counts
  #' @return No return value. This function saves files to the specified paths.
  
  # Check if the compressed versions of files already exist and delete
  print(path)
  compressed_mtx <- paste0(path, "/matrix.mtx.gz")
  compressed_features <- paste0(path, "/features.tsv.gz")
  compressed_barcodes <- paste0(path, "/barcodes.tsv.gz")
  if (file.exists(compressed_mtx)) {
    file.remove(compressed_mtx)
  }

  if (file.exists(compressed_features)) {
    file.remove(compressed_features)
  }

  if (file.exists(compressed_barcodes)) {
    file.remove(compressed_barcodes)
  }

  # Write barcodes
  gz1 <- gzfile(compressed_barcodes, "w")
  cat(colnames(seurat_obj), file = gz1, sep = "\n")
  close(gz1)

  # Write features
  if (is.null(feature_df_path)) {
    df <- data.frame(ID = "ENS_ID", gene_symbol = rownames(seurat_obj))
  }
  else {
    # add ensembl ids if provided
    # feature_df is csv file with gene_id (ENSG--- ensembl id) and gene_symbol column
    feature_df <- read.csv(feature_df_path, header = TRUE)

    # left join ensembl ids to seurat object gene symbols
    seurat_df <- data.frame(gene_symbol = rownames(seurat_obj))
    df <- seurat_df %>%
      dplyr::left_join(feature_df, by = "gene_symbol")
    
    df <- df[, c("gene_id", "gene_symbol")]
  }

  gz2 <- gzfile(compressed_features, "w")
  write.table(df, row.names = F, col.names = F, sep = "\t", quote = F, file = gz2)
  close(gz2)

  # Write matrix
  mat <- GetAssayData(seurat_obj, assay = assay, slot = slot)[, colnames(seurat_obj)]
  writeMM(mat, file = paste0(path, "/matrix.mtx")) #writeMM from Matrix
  gzip(filename = paste0(path, "/matrix.mtx")) #gzip from R.utils
}


write_seurat_spatial_coord <- function(seurat_obj, path, integrated = FALSE){
  if (integrated) {
    print("here")
    scale_factor_list <- lapply(names(seurat_obj@images), FUN = function(x) {
      # write tissue spot coordinates
      coordinates <- seurat_obj@images[[x]]@coordinates
      # remove Sample from x
      sample_alias <- gsub("sample", "", x)
      sample_path <- paste(path, sample_alias, sep = "/")
      print(sample_path)
      if (!dir.exists(sample_path)) {dir.create(sample_path, recursive = TRUE)}
      write.csv(coordinates, file = paste(sample_path, "tissue_positions_list.csv", sep = "/"))

      # write scale factors
      scale_factors <- seurat_obj@images[[x]]@scale.factors
      spot_radius <- seurat_obj@images[[x]]@spot.radius

      df <- data.frame(
        fiducial_diameter_fullres = scale_factors$fiducial,
        tissue_hires_scalef = scale_factors$hires,
        tissue_lowres_scalef = scale_factors$lowres
      )
      json_scale_factors <- toJSON(df, row.names = FALSE)
      write(json_scale_factors, file = paste(path, sample_alias, "scale_factors_json.json", sep = "/"))
      return(df)
    }
  )

  combined_df <- do.call("rbind", scale_factor_list)
  rownames(combined_df) <- names(seurat_obj@images)
  write.csv(combined_df, file = paste0(path, "combined_scale_factors.csv"))
  }
  else {
    # write tissue spot coordinates for single sample
    sample_name <- names(seurat_obj@images)
    sample_alias <- gsub("sample", "", sample_name)
    sample_path <- paste(path, sample_alias, sep = "/")
    if (!dir.exists(sample_path)) {dir.create(sample_path, recursive = TRUE)}

    coordinates <- seurat_obj@images[[sample_name]]@coordinates

    write.csv(coordinates, file = paste(path, sample_alias, "tissue_positions_list.csv", sep = "/"))
    # write scale factors
    scale_factors <- seurat_obj@images[[sample_name]]@scale.factors
    spot_radius <- seurat_obj@images[[sample_name]]@spot.radius

    df <- data.frame(
      fiducial_diameter_fullres = scale_factors$fiducial,
      tissue_hires_scalef = scale_factors$hires,
      tissue_lowres_scalef = scale_factors$lowres
    )
    json_scale_factors <- toJSON(df, row.names = FALSE)
    write(json_scale_factors, file = paste(path, sample_alias, "scale_factors_json.json", sep = "/"))
  }
}


save_counts_norm_meta <- function(seurat_obj, sample_alias, path, feature_df_path = NULL) {
    # This function saves a Seurat object, counts data, normalized data, and metadata to specified paths.
    #' @param seurat_obj A Seurat object that contains the data to be saved.
    #' @param sample_alias A character string that will be used as the filename for the saved files.
    #' @param path A character string that specifies the directory where the files will be saved.
    #' @param feature_df_path Path to tsv file with feature ids (ex. gene symbol + ensembl ids)
    #' @return No return value. This function saves files to the specified paths.

  # Save the Seurat object as an RDS file in the specified path
  sample_rds_path <- paste0(path, "/rds/")
  if (!file.exists(sample_rds_path)) {dir.create(sample_rds_path, recursive = TRUE)}
  saveRDS(seurat_obj, paste0(sample_rds_path, sample_alias, ".rds"))

  # Get the counts from the Seurat object 
  sample_counts_path <- paste0(path, "/counts/", sample_alias)
  if (!file.exists(sample_counts_path)) {dir.create(sample_counts_path, recursive = TRUE)}
  .write_seurat_to_disk(seurat_obj, sample_counts_path, feature_df_path, assay = "Spatial", slot = "counts")

  # Get the logCP10K data from the Seurat object 
  sample_lognorm_path <- paste0(path, "/logCP10K/", sample_alias)
  if (!file.exists(sample_lognorm_path)) {dir.create(sample_lognorm_path, recursive = TRUE)}
  .write_seurat_to_disk(seurat_obj, sample_lognorm_path, feature_df_path, assay = "Spatial", slot = "data")

  # Get the SCT data from the Seurat object 
  sample_SCT_path <- paste0(path, "/SCT_counts/", sample_alias)
  if (!file.exists(sample_SCT_path)) {dir.create(sample_SCT_path, recursive = TRUE)}
  .write_seurat_to_disk(seurat_obj, sample_SCT_path, feature_df_path, assay = "SCT", slot = "counts")

  sample_meta_path <- paste0(path, "/meta/")
  if (!file.exists(sample_meta_path)) {dir.create(sample_meta_path, recursive = TRUE)}

  # Save the metadata as a CSV file in the specified path
  write.csv(seurat_obj@meta.data, 
              paste0(sample_meta_path, sample_alias, ".csv"), 
              quote = TRUE
  )
}

write_dim_reduction <- function(seurat_obj, path, dim_reduction = "pca") {
  # Write dimensionality reduction to file 
  #' @param seurat_obj A Seurat object that contains the data to be saved.
  #' @param path A character string that specifies the directory  where the files will be saved.
  #' @param dim_reduction Name of dimensionality reduction to save (ex. PCA, UMAP); default is pca 
  #' @return No return value. This function saves files to the specified paths.
  
  # Check if dimensionality reduction exists in seurat object
  if (!file.exists(path)) {dir.create(path, recursive = TRUE)}

  avail_reductions <- names(seurat_obj@reductions)
  if (dim_reduction %in% avail_reductions) {
    embeddings_df <- as.data.frame(Embeddings(seurat_obj, reduction = dim_reduction))
    write.csv(embeddings_df, paste0(path, "/", dim_reduction, ".csv"))
  }
  else {
    stop(paste0("Error: Dimensionality reduction ", dim_reduction, " has not been calculated for this dataset"))
  }
}


.euclidean_distance <- function(spot1, spot2) {
  # This function calculates euclidean distance for two spots given their coords
  #' @param spot1 x and y coordinates of a spot 
  #' @param spot2 x and y coordinates of a spot 
  #' @return return a db with the euclidean distance value.
  x1 <- spot1[1]
  y1 <- spot1[2]
  x2 <- spot2[1]
  y2 <- spot2[2]
  return(sqrt((x1 - x2)^2 + (y1 - y2)^2))
}

calculateContactSpots <- function(sample, cond_subsetting, distance_threshold = 10) {
    #' This function finds the spots that are in contact with a selected group of spots.
    #' @param sample A Seurat object.
    #' @param ident Which condition/subcluster to filter within the idents 
    #' @param distance_threshold The maximum distance for a spot to be considered in contact. Default is 10.
    #' 55 µm spot's diameter, distance between centers of each spot is 100 µm.It works with 10, probably for scale factors. 
    #' @return A list with selected spots according the cond_subsetting and the contact spots
    stopifnot( cond_subsetting %in% levels(sample))
    
    # Get coordinates of all spots
    all_spots_coords <- .extract_coord(x = sample)
    
    # Get coordinates of selected spots
    selected_spots_coords <- all_spots_coords[WhichCells(sample, ident = cond_subsetting),]
    
    # Initialize matrix to store distances
    distances <- matrix(nrow = nrow(selected_spots_coords), ncol = nrow(all_spots_coords))
    
    # Calculate distances
    for (i in seq_len(nrow(selected_spots_coords))) {
    for (j in seq_len(nrow(all_spots_coords))) {
      distances[i, j] <- .euclidean_distance(selected_spots_coords[i, ], all_spots_coords[j, ])
    }
    }
    colnames(distances) <- rownames(all_spots_coords)
    rownames(distances) <- rownames(selected_spots_coords)
    
    distances <- distances[, !colnames(distances) %in% rownames(distances)]
    
    # Find columns where the minimum value in each column is less than or equal to the distance threshold
    contact_spots <- colnames(distances)[apply(distances, 2, function(x) min(x) <= distance_threshold)]
    
    return(list(contact_spots = contact_spots, selected_spots_coords = selected_spots_coords))
}


metrics_coexpression <- function(seuratObj, genes, group_var) {
    #' This function calculate the metrics for the coexpression of two genes in the group_var given
    #' @param seuratObj A Seurat object.
    #' @param genes Pair of genes to check the coexpression
    #' @param group_var Column from the metadata to separate the metrics, usually diagnosis
    #' @return A dataframe indicating the number of spots expressing one gene, both or none
    gene1 <- genes[1]
    gene2 <- genes[2]
    # Check if genes and metadata column are in the dataset
    if (!(gene1 %in% rownames(seuratObj)) | !(gene2 %in% rownames(seuratObj))) {
        stop("One or both genes are not in the dataset")
    }
    
    if (!(group_var %in% colnames(seuratObj@meta.data))) {
        stop("The metadata column is not in the dataset")
    }
    
    # Get the expression matrix
    exprMatrix <- seuratObj@assays$Spatial@counts
    
    # Get the cells expressing each gene
    gene1Cells <- which(exprMatrix[gene1, ] > 0)
    gene2Cells <- which(exprMatrix[gene2, ] > 0)
    
    # Get the metadata categories
    metadataCategories <- unique(seuratObj@meta.data[[group_var]])
    
    # Initialize the result table
    resultTable <- data.frame()
    
    # Loop over each metadata category
    for (category in metadataCategories) {
    # Get the cells in this category
    categoryCells <- which(seuratObj@meta.data[[group_var]] == category)
    
    # Calculate the number of cells expressing each gene and both
    gene1Only <- length(intersect(setdiff(gene1Cells, gene2Cells), categoryCells))
    gene2Only <- length(intersect(setdiff(gene2Cells, gene1Cells), categoryCells))
    bothGenes <- length(intersect(intersect(gene1Cells, gene2Cells), categoryCells))
    neitherGene <- length(setdiff(categoryCells, union(gene1Cells, gene2Cells)))
    
    # Calculate the percentages
    totalCells <- length(categoryCells)
    gene1OnlyPerc <- gene1Only / totalCells * 100
    gene2OnlyPerc <- gene2Only / totalCells * 100
    bothGenesPerc <- bothGenes / totalCells * 100
    neitherGenePerc <- neitherGene / totalCells * 100
    
    # Create the table for this category
    categoryTable <- data.frame(
      Expression = c(gene1, gene2, "Both", "Neither"),
      group_var = rep(category, 4),
      N_spots = c(gene1Only, gene2Only, bothGenes, neitherGene),
      Percentage = c(gene1OnlyPerc, gene2OnlyPerc, bothGenesPerc, neitherGenePerc)      
    )
    
    # Append to the result table
    resultTable <- rbind(resultTable, categoryTable)
    }
    
    return(resultTable)
}
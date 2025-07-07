#JXMM 20240522

stopifnot(
  require(tidyverse),
  require(clusterProfiler),
  require(org.Mm.eg.db),
  require(org.Hs.eg.db),
  require(ReactomePA),
  require(fgsea),
  require(msigdbr),
  require(patchwork)
)

#always same results, some github discussions claim that does not work 
set.seed(123)

plotting_gsea <- function(gsea_result, title){
  gsea_result %>% 
  # mutate(pathway = str_replace_all(pathway, "_", " ")) %>% 
  # mutate(pathway = ifelse(startsWith(pathway, "GOBP ADAPTIVE IMMUNE RESPONSE BASED"), "GOBP ADAPTIVE IMMUNE RESPONSE (SomRecomb)", pathway)) %>%
  #filter(NES > 0) %>% 
  ggplot(aes(abs(NES), reorder(pathway, abs(NES)), size = -log10(p.adjust), color = NES)) + 
  geom_point() +
  scale_color_gradient2(low = "navyblue", high = "firebrick") + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  theme_classic() +
  labs(x = "Normalized Enrichment Score", y = "", title = title, size = "-log10 padj") +
  theme(axis.text.y = element_text(size = 16)) +
  guides(size = guide_legend(order = 1)) 
  # +  scale_size(range = c(1, 9), name = "-log10(p.adjust)", limits = c(1,12), breaks = c(2,4,6,8,10))
}

# Helper function to process GSEA results
processGSEAresult <- function(gsea_result, orgDb, dbType) {
  if (nrow(gsea_result@result) == 0) {
    message(paste0("No pathways enriched at ", dbType, " database."))
    return(NULL)
  }
  gsea_result <- setReadable(gsea_result, orgDb, keyType = "ENTREZID")
  gsea_result@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", gsea_result@result$Description)
  gsea_result@result$dbType <- dbType
  names(gsea_result@result)[names(gsea_result@result) == 'Description'] <- 'pathway'
  return(gsea_result@result)
}


prepareGeneListForGSEA <- function(df, orgDb) {
  names(df) <- c('genes', 'avg_logFC')
  
  df <- df %>% arrange(desc(avg_logFC)) %>% mutate(avg_logFC = as.numeric(avg_logFC))
  # Map gene symbols to Entrez IDs
  df$entrez <- mapIds(orgDb,
                      keys=df$genes,
                      column="ENTREZID",
                      keytype="SYMBOL",
                      multiVals="first")
  
  # Prepare the gene list for GSEA
  geneList <- setNames(df$avg_logFC, df$entrez)
  geneList <- na.omit(geneList)
  return(geneList)
}


performGSEA <- function(df, dbType, species, minGSSize= 8, maxGSSize = 300 ) {
  #' Perform Gene Set Enrichment Analysis (GSEA) using either the KEGG, GO, or Reactome databases
#' for human or mouse gene lists.
#' @param df A dataframe with two columns: the first column contains gene !symbols!
#'        and the second column contains log fold changes. 
#' @param dbType A character string specifying the database for GSEA. Options are
#'        "KEGG", "GO", or "Reactome".
#' @param species A character string specifying the species. Options are "human"
#'        or "mouse".
#'
#' @return A dataframe containing the GSEA results, including enriched pathways,
#'         their descriptions, statistics
#' @examples
#' # Perform KEGG GSEA for a mouse gene list
#' result <- performGSEA(geneList, "KEGG", "mouse")
#' 
  #each method has different way of calling the species
  if(species == "human"){
    organism_kegg = "hsa"
    orgDB = org.Hs.eg.db
    organism_reactome = species
  } else if (species == "mouse"){
    organism_kegg = "mmu"
    orgDB = org.Mm.eg.db
    organism_reactome =species
  } else {stop("Please provide either 'human' or 'mouse'")}

  #need same format for all of them
  geneList <- prepareGeneListForGSEA(df, orgDB)

  if (dbType == "KEGG") {
    kegg_gsea_result <- gseKEGG(gene= geneList,
                                organism = organism_kegg,
                                pAdjustMethod="BH",
                                pvalueCutoff = 0.05,
                                minGSSize = minGSSize,
                                maxGSSize = maxGSSize,
                                verbose=FALSE,
                                seed=TRUE)
    return(processGSEAresult(kegg_gsea_result, orgDB, dbType))

  } else if (dbType == "GO") {
    go_gsea_results <- lapply(c("BP", "MF"), function(subont){
      go_gsea_result <- gseGO(gene = geneList,
                              OrgDb = orgDB,
                              pAdjustMethod="BH",
                              pvalueCutoff= 0.05,
                              minGSSize= minGSSize,
                              maxGSSize= maxGSSize,
                              verbose=FALSE,
                              seed=TRUE,
                              ont = subont)
      go_gsea_result@result <- go_gsea_result@result %>% 
        mutate(ont = subont) %>%
        mutate(Description =  paste0(ont, ": ", Description)) 
      return(go_gsea_result)
    })
    go_gsea_results <- lapply(go_gsea_results, processGSEAresult, orgDB, dbType)
    return(do.call(rbind, go_gsea_results))

  } else if (dbType == "Reactome") {
    reactome_gsea_result <- gsePathway(geneList,
                                       organism= organism_reactome,
                                       pAdjustMethod=  "BH",
                                       pvalueCutoff= 0.05,
                                       minGSSize= minGSSize,
                                       maxGSSize= maxGSSize,
                                       verbose=FALSE,
                                       seed=TRUE)
    return(processGSEAresult(reactome_gsea_result, orgDB, dbType))
    
  } else if (dbType == "Hallmarks") {
    genesets = msigdbr(species=species, category = "H")
    HALLMARK_GS = split(genesets$entrez_gene, genesets$gs_name)
    set.seed(13) #again just above fgsea call, double sure
    Hallmarks_gsea_result <- fgsea(HALLMARK_GS,
                                      geneList,
                                      minSize= minGSSize,
                                      maxSize= maxGSSize
                                  )
    Hallmarks_gsea_result[, leadingEdge := mapIdsList(
                                      x=orgDB, 
                                      keys=leadingEdge,
                                      keytype="ENTREZID", 
                                      column="SYMBOL")]
    #filter pvalue
    Hallmarks_gsea_result = Hallmarks_gsea_result[Hallmarks_gsea_result$padj < 0.05,]
    #change names to be consistent with other dbs
    names(Hallmarks_gsea_result)[names(Hallmarks_gsea_result) == 'padj'] <- 'p.adjust'

    #process manually since result is not consistent with other dbs
    if (nrow(Hallmarks_gsea_result) == 0) {
      message(paste0("No pathways enriched at ", dbType, " database."))
      return(NULL)
    }
    Hallmarks_gsea_result <- Hallmarks_gsea_result %>%
            mutate(pathway = str_replace_all(pathway, "_", " ")) %>%
            mutate(pathway = gsub("HALLMARK",
                                          "",
                                          pathway)) %>%
            mutate(dbType = dbType)
    return(Hallmarks_gsea_result)

    #possible errors
  } else {
    stop("Invalid database type. Please choose 'GO', 'KEGG', 'Reactome' or 'Hallmarks'.")
  }
}

#######   EXAMPLE   #######
#example for one db
#result <- performGSEA(df, "KEGG", "human")
# plotting_gsea(result, "Obese vs Healthy in adipocytes :)")

#example for all dbs
# dbs <- c("KEGG", "Reactome", "GO", "Hallmarks")
# results <- lapply(dbs, function(x){
#     result <- performGSEA(df, x, "human")
# })
# names(results) <- dbs
# results <- Filter(Negate(is.null), results)
# plot_results <- lapply(1:length(results), function(x){
#     plotting_gsea(results[[x]], names(results)[x])
# })
# patchwork::wrap_plots(plot_results, ncol=length(plot_results)) + plot_annotation("Title") + plot_layout(axis_titles = "collect")

# read_and_combine_excel_sheets <- function(path, object=NULL) {
#   library(readxl)
#   library(stringr)
#   sheet_names <- excel_sheets(path)
  
#   excel_data <- list()
  
#   for (sheet_name in sheet_names) {
#     if (sheet_name == "Hallmarks") {
#       excel_data[[sheet_name]] <- read_excel(path, sheet = sheet_name) %>%
#          mutate(pathway_all = str_to_title(pathway)) %>%
#         mutate(pathway_all = paste0(sheet_name, "-", pathway)) %>%
#         mutate(setSize = size) %>%
#         # dplyr::select(pathway, setSize, NES, p.adjust, core_enrichment, dbType)
#         dplyr::select(pathway, setSize, NES, p.adjust, dbType, pathway_all)
      
#     } else {
#     excel_data[[sheet_name]] <- read_excel(path, sheet = sheet_name) %>%
#       mutate(pathway_all = paste0(sheet_name, "-", pathway)) %>%
#       # dplyr::select(pathway, setSize, NES, p.adjust, core_enrichment, dbType, pathway_all)
#       dplyr::select(pathway, setSize, NES, p.adjust, dbType, pathway_all)
#     }
#   }
  
#   gsea_data <- bind_rows(excel_data)
#   gsea_data$pathway <- gsub("BP: ", "", gsea_data$pathway)
#   gsea_data$pathway <- gsub("MF: ", "", gsea_data$pathway)
#   gsea_data$pathway <- str_to_title(gsea_data$pathway)
#   return(gsea_data)
# }



run_gsea_jxmm <- function(markers, species = "mouse", minGSSize = 8, maxGSSize = 300) {
  #' Run GSEA for a given set of markers
  #' @param markers A dataframe with two columns: the first column contains gene symbols
  #'        and the second column contains log fold changes.
  #' @param species A character string specifying the species. Options are "human" or "mouse".
  #' @param minGSSize Minimum size of gene sets to consider for GSEA.
  #' @param maxGSSize Maximum size of gene sets to consider for GSEA.
  #'
  #' @return A list of GSEA results for KEGG, GO, Reactome, and Hallmarks databases.
  results <- lapply(c("KEGG", "Reactome", "GO", "Hallmarks"), function(x){

      result <- performGSEA(markers, x, species, minGSSize, maxGSSize)
      if (is.null(result)) {
        return(NULL)
      } else {
        if(x == "Hallmarks"){
          result <- result %>%
            mutate(pathway_all = str_to_title(pathway)) %>%
            mutate(core_enrichment = leadingEdge) %>%
            mutate(core_enrichment = str_replace_all(core_enrichment, ",", "/")) %>%
            mutate(setSize = size,
              dbSubtype = NULL) %>%
            dplyr::select(pathway, setSize, NES, p.adjust, dbType, pathway_all,core_enrichment)
        } else if(x == "GO"){
          result <- result %>% mutate(
                pathway_all = str_to_title(pathway),
                split_pathway = str_split(pathway, ":", simplify = TRUE),
                dbSubtype = split_pathway[, 1],
                pathway = split_pathway[, 2]
              ) %>% 
              dplyr::select(pathway, setSize, NES, p.adjust, dbType, pathway_all,core_enrichment)
        } else {
          result <- result %>% 
            mutate(pathway_all = str_to_title(pathway)) %>%
            mutate(dbSubtype = NULL) %>%
            dplyr::select(pathway, setSize, NES, p.adjust, dbType, pathway_all,core_enrichment)
        }
      }
      return(result)
  })

  results <- Filter(Negate(is.null), results)
  results <- bind_rows(results) %>%
      arrange(desc(NES)) %>%
      rownames_to_column("ID")

    return(results)
}


compute_gsea_overlap <- function(file_path1, file_path2, p_thresh = 0.05) {
  # Read the CSV files
  df1 <- read.csv(file_path1)
  df2 <- read.csv(file_path2)
  
  # enforce same ordering of pathways
  rownames(df1) <- df1$ID
  rownames(df2) <- df2$ID
  pathway_intersect <- intersect(rownames(df1), rownames(df2))
  df1 <- df1[pathway_intersect,]
  df2 <- df2[pathway_intersect,]
  
  # only focus on pathways significant in one of the two csvs
  df1$p.adjust <- stats::p.adjust(df1$pvalue, method = "BH")
  df2$p.adjust <- stats::p.adjust(df2$pvalue, method = "BH")
  pathways_of_interest <- intersect(
    rownames(df1)[which(df1$p.adjust <= p_thresh)],
    rownames(df2)[which(df2$p.adjust <= p_thresh)]
  )
  
  # Exit early if no significant overlapping pathways
  if (length(pathways_of_interest) == 0) {
    message("No pathways were found to be significant in both datasets.")
    return(NULL)
  }
  
  # Subset data frames to significant pathways
  df1 <- df1[pathways_of_interest,]
  df2 <- df2[pathways_of_interest,]
  
  # Calculate overlap metrics for each pathway
  result <- data.frame(
    ID = character(),
    Description = character(),
    pvalue1 = numeric(),
    pvalue2 = numeric(),
    padj1 = numeric(),
    padj2 = numeric(),
    num_geneset1 = integer(),
    num_geneset2 = integer(),
    num_intersect = integer(),
    overlap_perc = numeric(),
    gene_intersect = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:length(pathways_of_interest)) {
    pathway <- pathways_of_interest[i]
    
    # Extract gene sets for each pathway
    geneset1 <- strsplit(df1[pathway, "core_enrichment"], split = "/")[[1]]
    geneset2 <- strsplit(df2[pathway, "core_enrichment"], split = "/")[[1]]
    
    # Find the intersection of genes
    gene_intersect <- intersect(geneset1, geneset2)
    
    # Calculate overlap percentage (Jaccard index)
    union_size <- length(unique(c(geneset1, geneset2)))
    overlap_perc <- length(gene_intersect) / union_size
    
    # Compile result for this pathway
    pathway_result <- data.frame(
      ID = df1[pathway, "ID"],
      Description = df1[pathway, "Description"],
      pvalue1 = df1[pathway, "pvalue"],
      pvalue2 = df2[pathway, "pvalue"],
      padj1 = df1[pathway, "p.adjust"],
      padj2 = df2[pathway, "p.adjust"],
      num_geneset1 = length(geneset1),
      num_geneset2 = length(geneset2),
      num_intersect = length(gene_intersect),
      overlap_perc = overlap_perc,
      gene_intersect = paste0(gene_intersect, collapse = "/"),
      stringsAsFactors = FALSE
    )
    
    result <- rbind(result, pathway_result)
  }
  
  # Sort by overlap percentage (descending)
  result <- result[order(-result$overlap_perc),]
  
  return(result)
}
#' Hounkpe et al. housekeeping genes
#'
#' From https://academic.oup.com/nar/article/49/D1/D947/5871367#supplementary-data,
#' gkaa609_Supplemental_Files, specifically, \code{Supplementary_Table1.xlsx})
#'
#' @format ## `housekeeping_hounkpe_df`
#' A data frame with 1003 rows (genes) and 1 column:
#' \describe{
#'   \item{Gene}{Gene name}
#' }
#' @source <https://academic.oup.com/nar/article/49/D1/D947/5871367>
#' @examples
#' \dontrun{
#' # How the data was loaded
#' df <- openxlsx::read.xlsx("Supplementary_Table1.xlsx",
#' sheet = "Gene Model")
#' colnames(df) <- df[1,1]
#' housekeeping_hounkpe_df <- df[-1,,drop = FALSE]
#' usethis::use_data(housekeeping_hounkpe_df)
#' 
#' # To load
#' housekeeping_hounkpe <- housekeeping_hounkpe_df$Gene
#' }
"housekeeping_hounkpe_df"

#' Lin et al. housekeeping genes
#'
#' From https://academic.oup.com/gigascience/article/8/9/giz106/5570567, 
#' giz106_supplemental_table.xlsx.
#' 
#' Using the criterion that \code{Stability.index} should be above 0.8, there are 759 housekeeping genes
#'
#' @format ## `housekeeping_lin`
#' A data frame with 12986 rows (genes) and 3 columns:
#' \describe{
#'   \item{Gene}{Gene name}
#'   \item{Stability.index}{Stability index, the calculation the authors made to derive if a gene is "stable" or not}
#'   \item{Housekeeping}{Internal criteria defined as Stability.index larger that 0.8}
#' }
#' @source <https://academic.oup.com/gigascience/article/8/9/giz106/5570567>
#' @examples
#' \dontrun{
#' # How the data was loaded
#' df <- openxlsx::read.xlsx("giz106_supplemental_table.xlsx",
#' sheet = "human development")
#' colnames(df)[1] <- "Gene"
#' df$Housekeeping <- (df$Stability.index > 0.8)
#' df <- df[,c("Gene", "Stability.index", "Housekeeping")]
#' housekeeping_lin_df <- df
#' usethis::use_data(housekeeping_lin_df)
#' 
#' # To load
#' housekeeping_lin <- housekeeping_lin_df$Gene[housekeeping_lin_df$Housekeeping]
#' }
"housekeeping_lin_df"

#' Sun et al. microglia AD genes
#'
#' From https://www.sciencedirect.com/science/article/pii/S0092867423009716, 
#' 1-s2.0-S0092867423009716-mmc1.xlsx.
#' 
#' Using the criterion that a gene is deemed DEG if it was found to be DE in 5 or more (among the 24 microglia subtype and AD-progression combinations).
#' 
#' @format ## `microglia_sun_df`
#' A data frame with 3689 rows (genes) and 8 columns:
#' \describe{
#'   \item{Gene}{Gene name}
#'   \item{Pr..Chisq.}{P-value from MAST}
#'   \item{coef}{Coefficient for the gene in the MAST model}
#'   \item{ci.hi}{Upper confidence value for the coef}
#'   \item{ci.lo}{Lower confidence value for the coef}
#'   \item{fdr}{False discovery rate (i.e., adjusted p-value) based on the p-value}
#'   \item{groupID}{The microglia subtype the DE analysis was performed on, where "early" refers to control-vs-earlyAD and "late" refers to earlyAD-vs-lateAD}
#'   \item{DEG}{Internal criteria to deem a gene as DEG when looking at microglia as a whole}
#' }
#' @source <https://www.sciencedirect.com/science/article/pii/S0092867423009716>
#' @examples
#' \dontrun{
#' # How the data was loaded
#' df <- openxlsx::read.xlsx("1-s2.0-S0092867423009716-mmc1.xlsx",
#' sheet = "Page 10.DEGs_AD")
#' colnames(df)[1] <- "Gene"
#' microglia_sun_df <- df
#' table_vec <- table(microglia_sun_df$Gene)
#' enriched_genes <- unique(names(table_vec)[table_vec >= 5])
#' microglia_sun_df$DEG <- (microglia_sun_df$Gene %in% enriched_genes)
#' usethis::use_data(microglia_sun_df)
#' 
#' # To load
#' microglia_sun <- unique(microglia_sun_df$Gene[microglia_sun_df$DEG])
#' }
"microglia_sun_df"

#' Sun et al. microglia AD genes
#'
#' From https://www.nature.com/articles/s43587-023-00424-y, 
#' 43587_2023_424_MOESM3_ESM.xlsx.
#' 
#' Using the criterion that a gene is deemed DEG if its adjusted p-value is smaller than 0.05 in any microglia cluster
#' 
#' @format ## `microglia_prater_df`
#' A data frame with 942 rows (genes) and 9 columns:
#' \describe{
#'   \item{Gene}{Gene name}
#'   \item{baseMean}{Average mean across all donors}
#'   \item{log2FoldChange}{Log2 fold change between AD and non-AD}
#'   \item{lfcSE}{Standard error of the log fold change}
#'   \item{stat}{Test statistic from MAST}
#'   \item{pvalue}{P-value from MAST}
#'   \item{padj}{Adjust p-value from MAST}
#'   \item{cluster}{The microglia cluster that the DEG analysis was performed on}
#'   \item{DEG}{Internal criteria to deem a gene as DEG when looking at microglia as a whole}
#' }
#' @source <https://www.sciencedirect.com/science/article/pii/S0092867423009716>
#' @examples
#' \dontrun{
#' # How the data was loaded
#' df <- openxlsx::read.xlsx("43587_2023_424_MOESM3_ESM.xlsx",
#' sheet = "AD_vs_Ctrl_DEGs",
#' startRow = 3)
#' microglia_prater_df <- df
#' enriched_genes <- unique(microglia_prater_df$Gene[microglia_prater_df$padj <= 0.05])
#' microglia_prater_df$DEG <- (microglia_prater_df$Gene %in% enriched_genes)
#' usethis::use_data(microglia_prater_df)
#' 
#' # To load
#' microglia_prater <- unique(microglia_prater_df$Gene[microglia_prater_df$DEG])
#' }
"microglia_prater_df"
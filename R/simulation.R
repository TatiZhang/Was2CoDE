rm(list=ls())
set.seed(10)
library(SeuratObject)
library(Seurat)
library(IdeasCustom)
library(devtools)

n <- 10

# Mean and standard deviation vectors for cases and controls
case_mean_vec <- 2 + runif(n, min = -0.5, 0.5)  # Shifted means for case group
case_sd_vec <- rep(1, n)  # Standard deviation for case group
case_size_vec <- rep(1, n)  # Size multiplier for scaling

control_mean_vec <- 0 + runif(n, min = -0.5, 0.5)  # Shifted means for control group
control_sd_vec <- rep(1, n)  # Standard deviation for control group
control_size_vec <- rep(1, n)  # Size multiplier for scaling

# Set x and y axis limits for the plot
xlim <- c(-2, 4)
ylim <- c(0, 10)

# Generate Gaussian distributions for the case group
case_gaussian_list <- lapply(1:n, function(i) {
  mean_val <- case_mean_vec[i]
  sd_val <- case_sd_vec[i]
  xseq <- seq(xlim[1], xlim[2], length.out = 1000)
  yseq <- stats::dnorm(xseq, mean = mean_val, sd = sd_val)
  yseq <- yseq - min(yseq)
  yseq <- yseq / max(yseq) * 4 * case_size_vec[i]  # Normalize and scale
  cbind(xseq, yseq)  # Return x and y values as a matrix
})

# Generate Gaussian distributions for the control group
control_gaussian_list <- lapply(1:n, function(i) {
  mean_val <- control_mean_vec[i]
  sd_val <- control_sd_vec[i]
  xseq <- seq(xlim[1], xlim[2], length.out = 1000)
  yseq <- stats::dnorm(xseq, mean = mean_val, sd = sd_val)
  yseq <- yseq - min(yseq)
  yseq <- yseq / max(yseq) * 4 * control_size_vec[i]  # Normalize and scale
  cbind(xseq, yseq)  # Return x and y values as a matrix
})

############################################
# ---- Simulate data (same as before) ----
n <- 10
set.seed(10)

# Generate the case and control vectors
case_mean_vec <- 2 + runif(n, min = -0.5, 0.5)
case_sd_vec <- rep(1, n)
control_mean_vec <- 0 + runif(n, min = -0.5, 0.5)
control_sd_vec <- rep(1, n)

# Create meta information for individuals and cells
meta_ind <- data.frame(
  individual = 1:(2*n),
  ADNC = rep(c(1, 0), each = n)  # Binary variable for case/control status (1 = case, 0 = control)
)

# Simulate count matrix (Gaussian distributed data as placeholder for RNA counts)
count_matrix <- matrix(rnorm(200 * (2 * n), mean = c(rep(case_mean_vec, each = 200), 
                                                     rep(control_mean_vec, each = 200)),
                             sd = 1), 
                       nrow = 200 * (2 * n))

# Simulate metadata for cells
meta_cell <- data.frame(
  individual = rep(1:(2 * n), each = 200),
  donor_id = rep(1:(2 * n), each = 200),
  nCount_RNA = rowSums(count_matrix)  # Cell-specific variable based on sum of counts
)

# Simulate Seurat object for demonstration
seurat_obj <- CreateSeuratObject(counts = t(count_matrix))
seurat_obj$donor_id <- meta_cell$donor_id
seurat_obj$ADNC <- rep(c(1, 0), each = n)

# Extract the count matrix and metadata
count_matrix <- SeuratObject::LayerData(seurat_obj, features = Seurat::VariableFeatures(seurat_obj), layer = "data", assay = "RNA")

meta_cell$cell_id <- row.names(meta_cell)

meta_ind <- unique(data.frame(
  "individual" = seurat_obj$donor_id,
  "ADNC" = seurat_obj$ADNC,
  row.names = NULL
))

# Handle missing data (not applicable in this simulation but included for consistency)
for(j in 1:ncol(meta_ind)){
  if(!is.numeric(meta_ind[,j])) next()
  meta_ind[which(is.na(meta_ind[,j])), j] <- stats::median(meta_ind[,j], na.rm = TRUE)
}

# Variables for testing
var2test <- "ADNC"
var2adjust <- setdiff(colnames(meta_ind), c("individual", "ADNC"))
var2test_type <- "binary"  # ADNC is a binary variable
var_per_cell <- "nCount_RNA"  # Per cell variable

# ---- Apply was2 method (IDEAS) ----
dist_list <- IdeasCustom::ideas_dist_custom(
  count_input = count_matrix, 
  meta_cell = meta_cell, 
  meta_ind = meta_ind, 
  var_per_cell = var_per_cell, 
  var2test = var2test, 
  var2test_type = var2test_type, 
  verbose = 3
)


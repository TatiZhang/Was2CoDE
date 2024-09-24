rm(list=ls())
set.seed(10)

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

n_cases <- 10
n_controls <- 10
cells_per_donor <- 200
mean_case <- 2    
sd_case <- 1 
mean_control <- 0  
sd_control <- 1 

sim_data <- data.frame(
  donor_id = rep(1:(n_cases + n_controls), each = cells_per_donor),
  pathology = rep(c(rep("case", n_cases), rep("control", n_controls)), each = cells_per_donor),
  cell_id = 1:(cells_per_donor * (n_cases + n_controls)),
  value = c(
    rnorm(n_cases * cells_per_donor, mean = mean_case, sd = sd_case), 
    rnorm(n_controls * cells_per_donor, mean = mean_control, sd = sd_control)
  )
)

head(sim_data)

rm(list=ls())
set.seed(10)

n <- 10

# Color palettes for case and control groups
case_color_palette <- grDevices::colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                                    rgb(244, 84, 84, maxColorValue = 255)))(n)
control_color_palette <- grDevices::colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                       rgb(27, 198, 245, maxColorValue = 255)))(n)
transparent_gray <- rgb(0.5, 0.5, 0.5, 0.3)
two_letters <- substr(transparent_gray, start = 8, stop = 9)

# Add transparency to the palettes
case_color_trans_palette <- paste0(case_color_palette, two_letters)
control_color_trans_palette <- paste0(control_color_palette, two_letters)

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

# Create the plot
par(mar = c(3, 0.25, 0, 0.25))
plot(NA,
     xlim = xlim,
     ylim = ylim,
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "")

# Add control group polygons to the plot
for(i in 1:n) {
  graphics::polygon(x = c(control_gaussian_list[[i]][,1], rev(control_gaussian_list[[i]][,1])),
                    y = 2 + c(control_gaussian_list[[i]][,2], rep(0, nrow(control_gaussian_list[[i]]))),
                    col = control_color_trans_palette[i])
}

# Add case group polygons to the plot
for(i in 1:n) {
  graphics::polygon(x = c(case_gaussian_list[[i]][,1], rev(case_gaussian_list[[i]][,1])),
                    y = c(case_gaussian_list[[i]][,2], rep(0, nrow(case_gaussian_list[[i]]))),
                    col = case_color_trans_palette[i])
}

# Add rug plots to indicate the means
for(i in 1:n) {
  graphics::rug(control_mean_vec[i],
                col = control_color_palette[i],
                lwd = 2)
  graphics::rug(case_mean_vec[i],
                col = case_color_palette[i],
                lwd = 2)
}

# Add the x-axis
axis(side = 1, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)

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

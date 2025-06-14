# ================================================================
# 2D MOTILITY-INDUCED PHASE SEPARATION (MIPS) SIMULATION
# ================================================================
# Simulates clustering dynamics in active matter systems
# where particle motility decreases with local density
# ================================================================

# Install and load necessary packages
required_packages <- c("ggplot2", "gganimate", "gifski", "ggforce", 
                       "gridExtra", "viridis", "grid")
new_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(new_packages)) install.packages(new_packages)

library(ggplot2)       # For plotting
library(gganimate)     # For animation (optional)
library(ggforce)       # For advanced plotting features
library(gridExtra)     # For arranging multiple plots
library(grid)          # For low-level graphics
library(viridis)       # For color palettes

# ================================================================
# UTILITY FUNCTIONS
# ================================================================

#' Apply periodic boundary conditions to coordinates
#' 
#' @param coord Coordinate vector
#' @param domain_size Size of the square domain
#' @return Coordinates wrapped within [0, domain_size]
apply_periodic <- function(coord, domain_size) {
  coord %% domain_size
}

#' Compute displacement considering periodic boundaries
#' 
#' @param obj Current position
#' @param targ Target position
#' @param domain_size Size of the square domain
#' @return Shortest displacement vector in periodic space
wrapped_displacement <- function(obj, targ, domain_size) {
  d <- targ - obj
  d_p <- domain_size - abs(d)
  if (abs(d_p) < abs(d)) {
    d <- d_p * sign(-d)
  }
  return(d)
}

#' Compute periodic Euclidean distance between two points
#' 
#' @param x,y Current position coordinates
#' @param x_targ,y_targ Target position coordinates
#' @param domain_size Size of the square domain
#' @return Minimum distance in periodic space
periodic_distance <- function(x, y, x_targ, y_targ, domain_size) {
  delt_x <- wrapped_displacement(x, x_targ, domain_size)
  delt_y <- wrapped_displacement(y, y_targ, domain_size)
  sqrt(delt_x^2 + delt_y^2)
}

#' Fast neighbor counting using spatial binning
#' 
#' @param df_walkers Dataframe of particle positions
#' @param R Interaction radius
#' @param domain_size Size of the square domain
#' @return Vector of neighbor counts for each particle
count_all_neighbors_fast <- function(df_walkers, R, domain_size) {
  n_walkers <- nrow(df_walkers)
  neighbor_counts <- integer(n_walkers)
  
  # Spatial binning setup
  cell_size <- R
  df_walkers$cell_x <- floor(df_walkers$x / cell_size)
  df_walkers$cell_y <- floor(df_walkers$y / cell_size)
  
  # Create spatial dictionary
  cell_dict <- split(seq_len(n_walkers), 
                     paste(df_walkers$cell_x, df_walkers$cell_y, sep = "_"))
  
  # Check neighbors in 3x3 neighborhood
  for (j in seq_len(n_walkers)) {
    x_j <- df_walkers$x[j]
    y_j <- df_walkers$y[j]
    neighbor_count <- 0
    
    for (dx in -1:1) {
      for (dy in -1:1) {
        neighbor_cell <- paste(
          (df_walkers$cell_x[j] + dx) %% (domain_size / R),
          (df_walkers$cell_y[j] + dy) %% (domain_size / R),
          sep = "_"
        )
        
        if (neighbor_cell %in% names(cell_dict)) {
          neighbors <- cell_dict[[neighbor_cell]]
          neighbor_count <- neighbor_count + sum(
            sapply(neighbors, function(k) {
              periodic_distance(x_j, y_j, df_walkers$x[k], df_walkers$y[k], domain_size) < R
            })
          )
        }
      }
    }
    neighbor_counts[j] <- neighbor_count - 1  # Exclude self
  }
  return(neighbor_counts)
}

# ================================================================
# SIMULATION PARAMETERS
# ================================================================
set.seed(1288)  # Seed for reproducibility
domain_size <- 1        # Square domain [0, domain_size]

# Particle parameters
n_walkers <- 500        # Number of particles
l_min <- 1e-3           # Minimum step length
mu <- 1.2               # Lévy flight exponent (1 < μ ≤ 3)
n_substeps <- 30        # Substeps per main step

# Interaction parameters
R <- 0.05 * domain_size  # Interaction radius
gamma <- 4              # Motility reduction strength

# Stopping conditions
step_threshold <- l_min * exp(-gamma * 1)  # Mobility threshold
max_steps <- 200                           # Maximum simulation steps
consecutive_number <- 30                   # Steps below threshold to stop

# Visualization parameters
grid_res <- 32          # Resolution for density grid
plot_every_n <- 30      # Plot frequency (substeps)

# ================================================================
# SIMULATION INITIALIZATION
# ================================================================
steps <- 0
avg_step_sizes <- numeric()          # Stores average step sizes
avg_step_sizes_trimmed <- numeric()  # Stores trimmed averages
consecutive_below_threshold <- 0     # Convergence counter

# Initialize particles randomly
df_walkers <- data.frame(
  x = runif(n_walkers, 0, domain_size),
  y = runif(n_walkers, 0, domain_size),
  dx = 0,
  dy = 0
)

# ================================================================
# MAIN SIMULATION LOOP
# ================================================================
while (steps < max_steps && consecutive_below_threshold < consecutive_number) {
  steps <- steps + 1
  n_walkers <- nrow(df_walkers)
  
  # Generate Lévy-distributed step sizes and directions
  step_sizes <- l_min * (1 - runif(n_walkers))^(1 / (1 - mu))
  angles <- runif(n_walkers, 0, 2 * pi)
  
  total_step_magnitude <- numeric(n_walkers)  # Accumulator for step magnitudes
  
  for (i in 1:n_substeps) {
    # Count neighbors for motility adjustment
    neighbor_counts <- count_all_neighbors_fast(df_walkers, R, domain_size)
    
    # Adjust step size based on local density
    adjusted_step_size <- step_sizes * exp(-gamma * neighbor_counts)
    total_step_magnitude <- total_step_magnitude + adjusted_step_size / n_substeps
    
    # Update displacements
    df_walkers$dx <- adjusted_step_size * cos(angles) / n_substeps
    df_walkers$dy <- adjusted_step_size * sin(angles) / n_substeps
    
    # Update positions with periodic boundaries
    df_walkers$x <- apply_periodic(df_walkers$x + df_walkers$dx, domain_size)
    df_walkers$y <- apply_periodic(df_walkers$y + df_walkers$dy, domain_size)
    
    # Visualization (less frequent for performance)
    if (i %% plot_every_n == 0) {
      # Prepare density data
      df_walkers$grid_x <- floor(df_walkers$x * grid_res)
      df_walkers$grid_y <- floor(df_walkers$y * grid_res)
      
      df_density <- as.data.frame(table(df_walkers$grid_x, df_walkers$grid_y))
      names(df_density) <- c("grid_x", "grid_y", "count")
      df_density$grid_x <- as.numeric(as.character(df_density$grid_x))
      df_density$grid_y <- as.numeric(as.character(df_density$grid_y))
      
      # Create full grid for heatmap
      full_grid <- expand.grid(grid_x = 0:(grid_res-1), grid_y = 0:(grid_res-1))
      df_density_full <- merge(full_grid, df_density, by = c("grid_x", "grid_y"), all.x = TRUE)
      df_density_full$count[is.na(df_density_full$count)] <- 0
      df_density_full$x <- (df_density_full$grid_x + 0.5) / grid_res
      df_density_full$y <- (df_density_full$grid_y + 0.5) / grid_res
      
      # Plot 1: Particle positions
      plot1 <- ggplot() +
        geom_point(data = df_walkers, aes(x, y), size = 1) +
        labs(title = "Particle Positions",
             subtitle = paste("Step:", steps, "| Particles:", n_walkers, 
                              "\nμ:", mu, "| γ:", gamma, "| R:", R)) +
        coord_fixed(xlim = c(0, domain_size), ylim = c(0, domain_size)) +
        theme_minimal()
      
      # Plot 2: Density heatmap
      plot2 <- ggplot(df_density_full, aes(x, y, fill = count)) +
        geom_tile(width = 1/grid_res, height = 1/grid_res) +
        scale_fill_viridis_c(option = "plasma", name = "Density") +
        labs(title = "Local Density Heatmap") +
        coord_fixed() +
        theme_minimal()
      
      # Combine and display
      grid.arrange(plot1, plot2, ncol = 2)
    }
  }
  
  # Calculate trimmed average (exclude top 1%)
  sorted_steps <- sort(total_step_magnitude)
  cutoff_index <- floor(0.99 * length(sorted_steps))
  avg_trimmed <- mean(sorted_steps[1:cutoff_index])
  
  # Update convergence tracking
  avg_step_sizes_trimmed <- c(avg_step_sizes_trimmed, avg_trimmed)
  if (avg_trimmed <= step_threshold) {
    consecutive_below_threshold <- consecutive_below_threshold + 1
  } else {
    consecutive_below_threshold <- 0
  }
  
  # Progress reporting
  cat(sprintf("Step %d | Trimmed avg step: %.6f | Below threshold: %d/%d\n",
              steps, avg_trimmed, consecutive_below_threshold, consecutive_number))
}

# ================================================================
# POST-SIMULATION VISUALIZATION
# ================================================================

# Plot mobility convergence
par(mfrow = c(1, 2))

# Raw average step sizes
plot(seq_along(avg_step_sizes), log10(avg_step_sizes), type = "l", col = "blue",
     xlab = "Step", ylab = "Log10(Average Step Size)",
     main = "Mobility Decay")
abline(h = log10(step_threshold), col = "red", lty = 2)
legend("topright", legend = c("Mobility", "Threshold"), 
       col = c("blue", "red"), lty = 1:2)

# Trimmed average step sizes
plot(seq_along(avg_step_sizes_trimmed), log10(avg_step_sizes_trimmed), 
     type = "l", col = "purple",
     xlab = "Step", ylab = "Log10(Trimmed Step Size)",
     main = "Mobility Decay (99% Trimmed)")
abline(h = log10(step_threshold), col = "red", lty = 2)
legend("topright", legend = c("Trimmed Mobility", "Threshold"), 
       col = c("purple", "red"), lty = 1:2)
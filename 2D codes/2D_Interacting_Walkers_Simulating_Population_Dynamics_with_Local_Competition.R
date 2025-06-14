# ================================================================
# 2D INTERACTING WALKERS SIMULATION (BUGS MODEL)
# ================================================================
# Simulates population dynamics with movement, reproduction, and death
# where local density affects birth and death rates
# ================================================================

# Install and load necessary packages
required_packages <- c("ggplot2", "gifski", "ggforce")
new_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(new_packages)) install.packages(new_packages)

library(ggplot2)   # For plotting
library(ggforce)   # For advanced plotting features

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
  if (d > domain_size / 2) {
    d <- d - domain_size
  } else if (d < -domain_size / 2) {
    d <- d + domain_size
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
domain_size <- 1          # Square domain [0, domain_size]
initial_walkers <- 500    # Initial number of particles
l_min <- 1e-3             # Minimum step length
mu <- 2                   # Lévy flight exponent (1 < μ ≤ 3)
max_steps <- 150          # Maximum simulation steps
n_substeps <- 1           # Substeps per main step
set.seed(1288)            # Seed for reproducibility

# Interaction parameters
R <- 0.1 * domain_size    # Interaction radius
r_d0 <- 0.1               # Base death rate
r_b0 <- 1                 # Base birth rate
beta <- 0.02              # Death rate density dependence
alpha <- 0                # Birth rate density dependence

# Equilibrium detection parameters
starting_steps <- 45      # Warm-up period before checking equilibrium
equilibrium_steps <- 15   # Consecutive stable steps to stop

# ================================================================
# SIMULATION INITIALIZATION
# ================================================================
# Initialize walkers with random positions
df_walkers <- data.frame(
  x = runif(initial_walkers, 0, domain_size),
  y = runif(initial_walkers, 0, domain_size),
  dx = 0,
  dy = 0
)

# Initialize population tracking
pop_history <- numeric(max_steps + 1)
pop_history[1] <- initial_walkers

# Equilibrium tracking
equilibrium_reached <- FALSE
steps <- 0
stable_counter <- 0

# ================================================================
# MAIN SIMULATION LOOP
# ================================================================
while (steps < max_steps && !equilibrium_reached) {
  steps <- steps + 1
  n_walkers <- nrow(df_walkers)
  
  # Generate Lévy-distributed step sizes and directions
  step_sizes <- l_min * (1 - runif(n_walkers))^(1 / (1 - mu))
  angles <- runif(n_walkers, 0, 2 * pi)
  
  # Compute displacement per substep
  df_walkers$dx <- step_sizes * cos(angles) / n_substeps
  df_walkers$dy <- step_sizes * sin(angles) / n_substeps
  
  # Movement phase
  for (i in 1:n_substeps) {
    # Update positions with periodic boundaries
    df_walkers$x <- apply_periodic(df_walkers$x + df_walkers$dx, domain_size)
    df_walkers$y <- apply_periodic(df_walkers$y + df_walkers$dy, domain_size)
  }
  
  # Count neighbors for each walker
  N <- count_all_neighbors_fast(df_walkers, R, domain_size)
  
  # Compute density-dependent rates
  r_d <- r_d0 + beta * N           # Death rate
  r_b <- pmax(r_b0 - alpha * N, 0) # Birth rate (non-negative)
  
  # Convert rates to probabilities
  p_death <- 1 - exp(-r_d)
  p_birth <- 1 - exp(-r_b)
  
  # Apply death/birth events
  death_events <- runif(n_walkers) < p_death
  birth_events <- runif(n_walkers) < p_birth
  
  # Update population
  survivors <- df_walkers[!death_events, ]
  newborns <- df_walkers[birth_events, ]
  df_walkers <- rbind(survivors, newborns)
  
  # Record population size
  current_pop <- nrow(df_walkers)
  pop_history[steps + 1] <- current_pop
  
  # Equilibrium detection
  if (steps >= starting_steps) {
    # Calculate equilibrium estimate and fluctuation threshold
    recent_pop <- pop_history[(steps - equilibrium_steps + 1):steps]
    N_eq_estimate <- mean(recent_pop)
    pop_threshold <- max(1, ceiling(sqrt(2 * N_eq_estimate)))
    
    # Check stability
    if (abs(current_pop - pop_history[steps]) < pop_threshold) {
      stable_counter <- stable_counter + 1
    } else {
      stable_counter <- 0
    }
    
    # Check if equilibrium reached
    if (stable_counter >= equilibrium_steps) {
      equilibrium_reached <- TRUE
      message("\nEquilibrium reached at step ", steps, 
              ": Population stable at ", round(N_eq_estimate),
              " ±", pop_threshold)
    }
  }
  
  # Progress reporting
  cat(sprintf("Step %d: Population = %d | Stable steps: %d/%d\n", 
              steps, current_pop, stable_counter, equilibrium_steps))
}

# ================================================================
# VISUALIZATION AND RESULTS SAVING
# ================================================================
# Create population history data frame
pop_df <- data.frame(
  Step = 0:steps,
  Population = pop_history[1:(steps + 1)]
)

# Final particle positions plot
final_plot <- ggplot(df_walkers, aes(x, y)) +
  geom_point(size = 1.5, alpha = 0.6) +
  labs(title = "Final Particle Distribution",
       subtitle = paste("Step:", steps, "| Particles:", nrow(df_walkers),
                        "| μ:", mu, "| α:", alpha, "| β:", beta),
       x = "X", y = "Y") +
  coord_fixed(xlim = c(0, domain_size), ylim = c(0, domain_size)) +
  theme_minimal()

# Population dynamics plot
pop_plot <- ggplot(pop_df, aes(Step, Population)) +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  labs(title = "Population Dynamics",
       subtitle = paste("Final population:", nrow(df_walkers),
                        "| Equilibrium threshold: ±", pop_threshold),
       x = "Time Step", y = "Population Size") +
  theme_minimal()

# Display plots
print(final_plot)
print(pop_plot)
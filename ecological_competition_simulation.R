# Install necessary packages if not installed
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(ggplot2)

############################
### UTILITY FUNCTIONS ###
############################
# Apply periodic boundary conditions to coordinates
apply_periodic <- function(coord, domain_size) {
  coord %% domain_size
}

# Compute displacement considering periodic boundaries
wrapped_displacement <- function(obj, targ, domain_size) {
  d <- targ - obj
  if (d > domain_size / 2) {
    d <- d - domain_size
  } else if (d < -domain_size / 2) {
    d <- d + domain_size
  }
  return(d)
}

# Compute periodic Euclidean distance between two points
periodic_distance <- function(x, y, x_targ, y_targ, domain_size) {
  delt_x <- wrapped_displacement(x, x_targ, domain_size)
  delt_y <- wrapped_displacement(y, y_targ, domain_size)
  sqrt(delt_x^2 + delt_y^2)
}

# Count neighbors by species using spatial binning
count_neighbors_by_species <- function(df_walkers, R, domain_size) {
  n_walkers <- nrow(df_walkers)
  
  # Handle case with no interaction radius
  if (R == 0) {
    return(list(same = integer(n_walkers), other = integer(n_walkers)))
  }
  
  # Initialize neighbor counters
  neighbors_same <- integer(n_walkers)
  neighbors_other <- integer(n_walkers)
  
  # Create spatial grid for efficient neighbor search
  cell_size <- R
  df_walkers$cell_x <- floor(df_walkers$x / cell_size)
  df_walkers$cell_y <- floor(df_walkers$y / cell_size)
  cell_dict <- split(seq_len(n_walkers), 
                     paste(df_walkers$cell_x, df_walkers$cell_y, sep = "_"))
  
  n_cells <- floor(domain_size / cell_size)  # Cells per axis for wrap-around
  
  # Find neighbors for each individual
  for (j in seq_len(n_walkers)) {
    x_j <- df_walkers$x[j]
    y_j <- df_walkers$y[j]
    sp_j <- df_walkers$species[j]
    
    count_same <- 0
    count_other <- 0
    
    # Check 3x3 neighborhood of cells
    for (dx in -1:1) {
      for (dy in -1:1) {
        cx <- (df_walkers$cell_x[j] + dx) %% n_cells
        cy <- (df_walkers$cell_y[j] + dy) %% n_cells
        neighbor_cell <- paste(cx, cy, sep = "_")
        
        if (neighbor_cell %in% names(cell_dict)) {
          neighbors <- cell_dict[[neighbor_cell]]
          for (k in neighbors) {
            if (k != j) {
              dist <- periodic_distance(x_j, y_j, df_walkers$x[k], 
                                        df_walkers$y[k], domain_size)
              if (dist < R) {
                if (df_walkers$species[k] == sp_j) {
                  count_same <- count_same + 1
                } else {
                  count_other <- count_other + 1
                }
              }
            }
          }
        }
      }
    }
    neighbors_same[j] <- count_same
    neighbors_other[j] <- count_other
  }
  
  return(list(same = neighbors_same, other = neighbors_other))
}

# Create visualization of current state
plot_walkers <- function(df_walkers, step_num, pop_A, pop_B, domain_size,
                         R, r_b0, r_d0, alfa, beta, kappa_b, kappa_d,
                         mu_A, mu_B, l_min,
                         gamma_fast, gamma_slow, R_MIPS) {
  
  subtitle_text <- sprintf(
    "STEP: %d | A: %d | B: %d | R: %.2f\nr_b0: %.2f | r_d0: %.2f | α: %.2f | β: %.2f\nκ_b: %.2f | κ_d: %.2f | μ_A: %.2f | μ_B: %.2f",
    step_num, pop_A, pop_B, R, r_b0, r_d0, alfa, beta, kappa_b, kappa_d, mu_A, mu_B
  )
  
  ggplot() +
    annotate("rect", xmin = 0, xmax = domain_size, ymin = 0, ymax = domain_size,
             fill = NA, color = "blue", linetype = "dashed", linewidth = 1) +
    geom_point(data = df_walkers, aes(x = x, y = y, color = species),
               size = 1.5, alpha = 0.7) +
    scale_color_manual(values = c("A" = "dodgerblue", "B" = "firebrick")) +
    labs(
      title = "Ecological Competition Simulation",
      subtitle = subtitle_text,
      x = "X",
      y = "Y",
      color = "Species"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom"
    ) +
    coord_fixed(xlim = c(0, domain_size), ylim = c(0, domain_size))
}


############################
### SIMULATION PARAMETERS ###
############################

# Set seed for reproducibility
set.seed(1287)

### DOMAIN AND POPULATION ###
domain_size <- 1              # Square domain [0, domain_size]
initial_walkers_A <- 250      # Initial population of species A
initial_walkers_B <- 250      # Initial population of species B
max_steps <- 200              # Maximum simulation steps

### MOVEMENT PARAMETERS ###
l_min <- 1e-3                 # Minimum step length
mu_A <- 2.0                   # Lévy exponent for species A (higher = longer jumps)
mu_B <- 1.5                   # Lévy exponent for species B

### INTERACTION PARAMETERS ###
R <- 0.1 * domain_size        # Interaction radius
r_d0 <- 0.1                   # Base death rate
r_b0 <- 0.5                   # Base birth rate
beta <- 0.02                  # Density-dependent death coefficient
kappa_b <- 0.03               # Competition effect on birth
kappa_d <- 0.03               # Competition effect on death
alfa <- 0                     # Density-dependent birth coefficient

### MOVEMENT ADJUSTMENT PARAMETERS ###
gamma_slow <- 0.4             # Crowding slow-down factor
gamma_fast <- 0.04            # Isolation speed-up factor
R_MIPS <- R * 0.5             # Radius for movement adjustment

### SIMULATION CONTROL ###
n_substeps <- 30              # Movement substeps per main step
plot_every_n <- 1             # Plot every N steps
plot_substeps <- FALSE        # Plot substeps (not recommended)

### EQUILIBRIUM DETECTION ###
equilibrium_steps <- 15       # Steps for population stability check
warming_steps <- 45           # Initial steps before equilibrium check
step_threshold <- l_min * exp(-gamma_slow)  # Step size threshold for MIPS
consecutive_number <- 15       # Consecutive steps below threshold for MIPS

############################
### SIMULATION EXECUTION ###
############################

# Initialize walkers for both species
df_walkers <- rbind(
  data.frame(x = runif(initial_walkers_A, 0, domain_size),
             y = runif(initial_walkers_A, 0, domain_size),
             dx = 0, dy = 0,
             species = factor(rep("A", initial_walkers_A))),
  data.frame(x = runif(initial_walkers_B, 0, domain_size),
             y = runif(initial_walkers_B, 0, domain_size),
             dx = 0, dy = 0,
             species = factor(rep("B", initial_walkers_B)))
)

# Prepare population history per species
pop_history_A <- numeric(max_steps + 1)
pop_history_B <- numeric(max_steps + 1)
pop_history_A[1] <- sum(df_walkers$species == "A")
pop_history_B[1] <- sum(df_walkers$species == "B")

# Simulation control variables
stable_counter <- 0
consecutive_below_threshold <- 0        # To count how many times the average step size is below the threshold
equilibrium_reached <- FALSE
steps <- 0
avg_step_sizes_trimmed <- numeric()      # To store trimmed step size averages

# Initialize plots:
plot.new()
dev.new()

# Main simulation loop
while (steps < max_steps && !equilibrium_reached) {
  steps <- steps + 1
  
  n_walkers <- nrow(df_walkers)
  
  # Generate species-specific mu values
  mu_vec <- ifelse(df_walkers$species == "A", mu_A, mu_B)
  
  step_sizes <- l_min * (1 - runif(n_walkers))^(1 / (1 - mu_vec))
  angles <- runif(n_walkers, 0, 2 * pi)
  
  # This will accumulate the step magnitude for all substeps
  total_step_magnitude <- numeric(n_walkers)
  
  for (i in 1:n_substeps) {
    # Count neighbors for each walker
    neighbors <- count_neighbors_by_species(df_walkers, R_MIPS, domain_size)
    N_same <- neighbors$same
    N_other <- neighbors$other
    
    # Adjust step size to slow down if there are many same-species neighbors
    adjusted_step_size <- step_sizes * exp(-gamma_slow * N_same + gamma_fast * N_other)
    total_step_magnitude <- total_step_magnitude + adjusted_step_size / n_substeps
    
    # Update the displacement based on adjusted step size
    df_walkers$dx <- adjusted_step_size * cos(angles) / n_substeps
    df_walkers$dy <- adjusted_step_size * sin(angles) / n_substeps
    
    # Update positions with periodic boundary conditions
    df_walkers$x <- apply_periodic(df_walkers$x + df_walkers$dx, domain_size)
    df_walkers$y <- apply_periodic(df_walkers$y + df_walkers$dy, domain_size)
    
    if (plot_substeps) {
      plot_sub <- ggplot(df_walkers, aes(x = x, y = y, color = species)) +
        geom_point(size = 1, alpha = 0.65) +
        scale_color_manual(values = c("A" = "blue", "B" = "red")) +
        labs(
          title = paste("Main Step:", steps, "| Substep:", i),
          subtitle = paste("A:", sum(df_walkers$species == "A"), " B:", sum(df_walkers$species == "B")),
          x = "X", y = "Y"
        ) +
        theme_minimal() +
        coord_fixed(xlim = c(0, domain_size), ylim = c(0, domain_size))
      
      dev.hold()
      print(plot_sub)
      Sys.sleep(2)  # Pause for 2 seconds
      dev.flush()
    }
    
  }
  
  if(n_substeps >1 && steps >= warming_steps){
    # Equilibrium MIPS
    # Remove top 1% of steps (potential outliers)
    sorted_steps <- sort(total_step_magnitude)
    cutoff_index <- floor(0.99 * length(sorted_steps))
    filtered_steps <- sorted_steps[1:cutoff_index]
    avg_trimmed <- mean(filtered_steps)
    
    avg_step_sizes_trimmed <- c(avg_step_sizes_trimmed, avg_trimmed)
    # Check if trimmed average step size is below threshold
    if (avg_trimmed <= step_threshold) {
      consecutive_below_threshold <- consecutive_below_threshold + 1
    } else {
      consecutive_below_threshold <- 0
    }
    
    cat(sprintf("Step %d | Avg. trimmed step size = %.8f vs threshold = %.8f\n", steps, avg_trimmed, step_threshold))
    # Print stopping condition progress
    cat(sprintf("Consecutive steps below threshold: %d / %d\n", consecutive_below_threshold, consecutive_number))
    if(consecutive_below_threshold >= consecutive_number){
      equilibrium_reached <- TRUE
    }
  }
  
  #Equilibrium population
  # Count neighbors by species
  neighbors <- count_neighbors_by_species(df_walkers, R, domain_size)
  N_same <- neighbors$same
  N_other <- neighbors$other
  
  # Death rate: own base + beta * same species neighbors + kappa * other species neighbors
  r_d <- r_d0 + beta * N_same +  kappa_d * N_other
  
  # Birth rate: adjusted by alfa and same species neighbors (could be customized)
  r_b <- pmax(r_b0 - alfa * N_same - kappa_b * N_other, 0)
  
  p_death <- 1 - exp(-r_d)
  p_birth <- 1 - exp(-r_b)
  
  random_vals_death <- runif(n_walkers)
  random_vals_birth <- runif(n_walkers)
  
  survivor_indices <- which(random_vals_death >= p_death)
  df_survivors <- df_walkers[survivor_indices, ]
  
  birth_indices <- which(random_vals_birth < p_birth)
  df_birth <- df_walkers[birth_indices, ]
  
  # New walkers born at parent's location, with same species
  df_walkers <- rbind(df_survivors, df_birth)
  
  pop_history_A[steps + 1] <- sum(df_walkers$species == "A")
  pop_history_B[steps + 1] <- sum(df_walkers$species == "B")
  
  # Stop simulation if one species goes extinct
  if (pop_history_A[steps + 1] == 0 || pop_history_B[steps + 1] == 0) {
    cat("Extinction event: Species A or B reached 0 individuals at step", steps, "\n")
    equilibrium_reached <- TRUE
    break
  }
  
  # Check equilibrium on total population (can be refined)
  total_pop <- pop_history_A[steps + 1] + pop_history_B[steps + 1]
  if (steps >= max(warming_steps, equilibrium_steps)) {
    recent_pops <- pop_history_A[(steps - equilibrium_steps + 1):steps] + pop_history_B[(steps - equilibrium_steps + 1):steps]
    N_eq_estimate <- mean(recent_pops)
    pop_threshold <- max(1, ceiling(sqrt(2 * N_eq_estimate )))
    
    if (abs(total_pop - (pop_history_A[steps] + pop_history_B[steps])) < pop_threshold) {
      stable_counter <- stable_counter + 1
    } else {
      stable_counter <- 0
    }
    cat(sprintf("Stable steps in a row: %d/%d | Population A: %d, B: %d\n\n", stable_counter, equilibrium_steps, pop_history_A[steps + 1], pop_history_B[steps + 1]))
    
    if (stable_counter >= equilibrium_steps && steps >= warming_steps) {
      equilibrium_reached <- TRUE
      message("Equilibrium reached at step ", steps, ": population stable within threshold ", pop_threshold)
    }
  }
  
  # Plot every plot_every_n steps
  if (steps %% plot_every_n == 0) {
    plot <- plot_walkers(
      df_walkers = df_walkers,
      step_num = steps,
      pop_A = pop_history_A[steps + 1],
      pop_B = pop_history_B[steps + 1],
      domain_size = domain_size,
      R = R,
      r_b0 = r_b0,
      r_d0 = r_d0,
      alfa = alfa,
      beta = beta,
      kappa_b = kappa_b,
      kappa_d = kappa_d,
      mu_A = mu_A,
      mu_B = mu_B,
      l_min = l_min,
      gamma_fast = gamma_fast,
      gamma_slow = gamma_slow,
      R_MIPS = R_MIPS
    )
    
    dev.hold()
    print(plot)
    dev.flush()
  }
}

# Final plot
plot <- plot_walkers(
  df_walkers = df_walkers,
  step_num = steps,
  pop_A = pop_history_A[steps + 1],
  pop_B = pop_history_B[steps + 1],
  domain_size = domain_size,
  R = R,
  r_b0 = r_b0,
  r_d0 = r_d0,
  alfa = alfa,
  beta = beta,
  kappa_b = kappa_b,
  kappa_d = kappa_d,
  mu_A = mu_A,
  mu_B = mu_B,
  l_min = l_min,
  gamma_fast = gamma_fast,
  gamma_slow = gamma_slow,
  R_MIPS = R_MIPS
)

print(plot)


# Plot population history per species over time
df_pop <- data.frame(
  Step = 0:steps,
  Population_A = pop_history_A[1:(steps + 1)],
  Population_B = pop_history_B[1:(steps + 1)]
)

plot_pop <- ggplot(df_pop, aes(x = Step)) +
  geom_line(aes(y = Population_A, color = "Species A")) +
  geom_point(aes(y = Population_A, color = "Species A")) +
  geom_line(aes(y = Population_B, color = "Species B")) +
  geom_point(aes(y = Population_B, color = "Species B")) +
  geom_hline(yintercept = N_eq_estimate, linetype = "dashed", color = "red") +
  annotate("text", x = steps/2, y = N_eq_estimate + 10, 
           label = sprintf("Est. Equilibrium: %.0f", N_eq_estimate), color = "red")+
  scale_color_manual(values = c("Species A" = "blue", "Species B" = "red")) +
  labs(title = "Population Size vs Steps",
       subtitle = paste("STEP:", steps,
                        "| Population A:", pop_history_A[steps + 1],
                        "| Population B:", pop_history_B[steps + 1],
                        "| Interaction radius R:", R,"| R_MIPS:", R_MIPS,
                        "\nr_b0:", r_b0, "| r_d0:", r_d0, "| α:", alfa,
                        "| β:", beta, "| κ_b:", kappa_b, "| κ_d:", kappa_d,
                        "| mu_A:", mu_A, "| mu_B:", mu_B, "| l_min:", l_min),
       x = "Step",
       y = "Population Size",
       color = "Species") +
  theme_minimal()

print(plot_pop)

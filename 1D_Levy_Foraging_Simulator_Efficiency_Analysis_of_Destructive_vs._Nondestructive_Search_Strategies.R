# Install and load necessary packages
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

# ================================================================
# SIMULATION PARAMETERS
# ================================================================
set.seed(200)  # Seed for reproducibility
n_simulations <- 10000000    # Number of simulations per (mu, lambda) pair, 
                           #Nondestructive: 10000000, Destructive: 100000
n_partitions <- 40          # Number of mu values to test in interval (1,3]
max_steps <- 40000000       # Max steps per walker (prevents infinite loops)

# Random walk parameters
l_min <- 1                  # Minimum step length (r_v in the model)
r_v <- l_min                # Direct vision range = minimum step length
mu_max <- 3                 # Maximum Lévy exponent
mu_min <- 1                 # Minimum Lévy exponent
mu_values <- seq(           # Create sequence of mu values to test
  from = mu_min + (mu_max - mu_min)/n_partitions,
  to = mu_max,
  length.out = n_partitions
)
warning_count <- 0          # Counter for simulations hitting max_steps
lambda_values <- c(10, 100, 1000, 10000)  # Target spacings to test
results <- data.frame()     # Storage for all results

# Vectors for correlation and covariance metrics
cor_N_l <- numeric(n_partitions)   # Pearson correlation between N and avg step length
cov_N_l <- numeric(n_partitions)   # Covariance between N and avg step length

# ================================================================
# MAIN SIMULATION LOOP
# ================================================================
for (lambda in lambda_values) {
  print(paste("Current lambda:", lambda))
  
  # Domain setup - ADJUST FOR FORAGING TYPE:
  # Destructive foraging: L = 2*lambda, start at s = lambda
  # Nondestructive foraging: L = lambda, start at s = r_v (as shown)
  L <- lambda              # Domain size (0 to L)
  s <- r_v                 # Starting position (r_v from left boundary)
  
  # Initialize efficiency vectors for current lambda
  efficiency1 <- numeric(n_partitions)
  efficiency2 <- numeric(n_partitions)
  efficiency3 <- numeric(n_partitions)
  
  for (j in seq_along(mu_values)) {
    mu <- mu_values[j]
    print(paste("Current mu:", mu))
    
    # Storage for simulation results
    absorption_times <- numeric(n_simulations)  # Steps until absorption (N)
    distance <- numeric(n_simulations)          # Total distance traveled
    
    # Run individual simulations
    for (i in 1:n_simulations) {
      position <- s        # Start at initial position
      steps <- 0           # Step counter
      total_dist <- 0      # Distance accumulator
      
      # Random walk until absorption or max steps
      while (position > 0 && position < L && steps < max_steps) {
        steps <- steps + 1
        
        # Random direction: 50% left, 50% right
        sign <- ifelse(runif(1) < 0.5, -1, 1)
        
        # Generate Lévy-distributed step using inverse transform sampling
        # Correct implementation based on Eq. 2.4: ℓ = ℓ_min * (1 - u)^{1/(1-μ)}
        step_sizes <- l_min * (1 - runif(1))^(1/(1 - mu))
        
        # Calculate movement and new position
        move <- step_sizes * sign
        next_position <- position + move
        
        # Handle boundary absorption:
        # - Partial steps when crossing boundaries
        # - Full steps when staying within domain
        if (next_position <= 0) {
          # Absorbed at left boundary (0)
          total_dist <- total_dist + abs(position)  # Actual distance to boundary
        } else if (next_position >= L) {
          # Absorbed at right boundary (L)
          total_dist <- total_dist + (L - position)  # Actual distance to boundary
        } else {
          # Still in domain - full step taken
          total_dist <- total_dist + abs(move)
        }
        
        position <- next_position  # Update position
      }
      
      # Handle max steps reached
      if (steps == max_steps) {
        warning_count <- warning_count + 1
        warning_msg <- paste("Simulation", i, "mu =", mu, "lambda =", lambda, 
                             ": Max steps reached")
        warning(warning_msg)
        print(warning_msg)
      }
      
      # Store results for this simulation
      absorption_times[i] <- steps
      distance[i] <- total_dist
    }
    
    # Identify successful simulations (absorbed before max_steps)
    successful <- absorption_times < max_steps
    
    # ================================================================
    # EFFICIENCY CALCULATIONS
    # ================================================================
    # Three different efficiency estimators:
    # η₁ = 1/<d>         (inverse of mean distance)
    efficiency1[j] <- 1 / mean(distance[successful])
    
    # η₂ = 1/(<N><ℓ>)    (based on mean steps and mean step length)
    mean_N <- mean(absorption_times[successful])
    mean_l <- mean(distance[successful] / absorption_times[successful])
    efficiency2[j] <- 1 / (mean_N * mean_l)
    
    # η₃ = <1/d>         (mean of inverse distances)
    efficiency3[j] <- mean(1 / distance[successful])
    
    # ================================================================
    # CORRELATION ANALYSIS
    # ================================================================
    N_vals <- absorption_times[successful]          # Number of steps (N)
    l_vals <- distance[successful] / N_vals         # Average step length (ℓ)
    
    # Pearson correlation between N and ℓ
    cor_N_l[j] <- cor(N_vals, l_vals)
    
    # Covariance between N and ℓ
    cov_N_l[j] <- cov(N_vals, l_vals)
  }
  
  # Store results for current lambda
  results <- rbind(
    results,
    data.frame(
      mu = mu_values,
      efficiency1 = efficiency1,
      efficiency2 = efficiency2,
      efficiency3 = efficiency3,
      lambda = lambda,
      cor_N_l = cor_N_l,
      cov_N_l = cov_N_l  # Covariance metric
    )
  )
}

# ================================================================
# VISUALIZATION SECTION
# ================================================================
# Create labels for lambda values (with mathematical notation)
lambda_labels <- c(
  `10` = expression(10),
  `100` = expression(10^2),
  `1000` = expression(10^3),
  `10000` = expression(10^4)
)

# Plot 1: Efficiency vs. mu (η₁ estimator)
ggplot(results, aes(x = mu, y = lambda * efficiency1, color = factor(lambda))) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 3)) +
  scale_color_discrete(
    name = expression(lambda),
    labels = lambda_labels
  ) +
  labs(
    title = expression("Efficiency " * lambda * eta[1] * " vs " * mu),
    x = expression(mu),
    y = expression(lambda * eta[1](mu))
  ) +
  theme_minimal(base_size = 14)

# Plot 2: Pearson correlation between N and average step length
ggplot(results, aes(x = mu, y = cor_N_l, color = factor(lambda))) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 3)) +
  scale_color_discrete(
    name = expression(lambda),
    labels = lambda_labels
  ) +
  labs(
    title = expression("Correlation: " * rho(N, bar(italic(l)))),
    x = expression(mu),
    y = expression(rho(N, bar(italic(l))))
  ) +
  theme_minimal(base_size = 14)

# Plot 3: Covariance between N and average step length
ggplot(results, aes(x = mu, y = cov_N_l, color = factor(lambda))) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 3)) +
  scale_color_discrete(
    name = expression(lambda),
    labels = lambda_labels
  ) +
  labs(
    title = expression("Covariance: " * Cov(N, bar(italic(l)))),
    x = expression(mu),
    y = expression(Cov(N, bar(italic(l))))
  ) +
  theme_minimal(base_size = 14)

# ================================================================
# ERROR ANALYSIS SECTION
# ================================================================
# Calculate corrected efficiency using covariance
results$predicted_eta2_with_cov <- 1 / (1 / results$efficiency2 + results$cov_N_l)

# Relative error between η₁ and covariance-corrected η₂
results$eta1_relative_error <- abs(
  (results$efficiency1 - results$predicted_eta2_with_cov) / 
    results$efficiency1
)

# Plot relative error with covariance correction
ggplot(results, aes(x = mu, y = eta1_relative_error, color = factor(lambda))) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 3)) +
  scale_color_discrete(
    name = expression(lambda),
    labels = lambda_labels
  ) +
  labs(
    title = expression("Relative Error: " * eta[1] * " vs " * eta[2]^"corr"),
    x = expression(mu),
    y = expression(epsilon)
  ) +
  theme_minimal(base_size = 14)

# Calculate relative error without covariance correction
results$eta2_relative_error_naive <- abs(
  (results$efficiency1 - results$efficiency2) / 
    results$efficiency1
)

# Plot naive relative error
ggplot(results, aes(x = mu, y = eta2_relative_error_naive, color = factor(lambda))) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 3)) +
  scale_color_discrete(
    name = expression(lambda),
    labels = lambda_labels
  ) +
  labs(
    title = expression("Relative Error: " * eta[1] * " vs " * eta[2]),
    x = expression(mu),
    y = expression(epsilon)
  ) +
  theme_minimal(base_size = 14)

# Plot covariance-corrected efficiency
ggplot(results, aes(x = mu, y = lambda * predicted_eta2_with_cov, color = factor(lambda))) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 3)) +
  scale_color_discrete(
    name = expression(lambda),
    labels = lambda_labels
  ) +
  labs(
    title = expression("Covariance-Corrected " * lambda * eta[2]),
    x = expression(mu),
    y = expression(lambda * eta[2]^"corr"(mu))
  ) +
  theme_minimal(base_size = 14)
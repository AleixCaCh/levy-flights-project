# Lévy Flights and Competition Simulations Repository

This repository contains R scripts for simulating Lévy flights and competition dynamics in 1D and 2D environments. The simulations explore search efficiency, motility-induced phase separation (MIPS), and competitive interactions between species.

## Table of Contents
1. [1D Efficiency Simulation](#1d-efficiency-simulation)
2. [2D Interacting Walkers Simulation](#2d-interacting-walkers-simulation)
3. [2D Motility-Induced Phase Separation](#2d-motility-induced-phase-separation)
4. [2D Species Competition Simulation](#2d-species-competition-simulation)
5. [Dependencies](#dependencies)
6. [Usage](#usage)
7. [Contributing](#contributing)

---

## 1D Efficiency Simulation
`1D_Levy_Foraging_Simulator_Efficiency_Analysis_of_Destructive_and_Nondestructive_Search_Strategies.R`

### Purpose
Simulates Lévy flights in 1D environments to analyze search efficiency and boundary absorption dynamics for both destructive and non-destructive foraging strategies.

### Key Features
- **Lévy flight generation** with exponent μ ∈ (1,3]
- Three efficiency estimators: 
  - η₁ = 1/〈d〉 
  - η₂ = 1/(〈N〉〈ℓ〉)
  - η₃ = 〈1/d〉
- Correlation analysis between step count (N) and step length (ℓ)
- Error analysis with covariance correction
- Domain size adaptation for different foraging types

### Parameters
```r
n_simulations <- 10000000   # Simulations per (μ, λ)
n_partitions <- 40           # μ values in (1,3]
lambda_values <- c(10, 100, 1000, 10000)  # Target spacings
max_steps <- 40000000        # Max steps per walker
```

### Outputs
1. Efficiency vs. μ for different λ values
2. Correlation ρ(N, ℓ) vs. μ
3. Covariance Cov(N, ℓ) vs. μ
4. Relative error analysis between estimators

---

## 2D Interacting Walkers Simulation
`2D_Interacting_Walkers_Simulating_Population_Dynamics_with_Local_Competition.R`

### Purpose
Simulates population dynamics with density-dependent birth/death rates and Lévy flight movement in 2D environments.

### Key Features
- Periodic boundary conditions
- Spatial binning for efficient neighbor counting
- Density-dependent rates
- Equilibrium detection algorithm based on population stability
- Real-time population tracking

### Parameters
```r
domain_size <- 1             # Simulation area
initial_walkers <- 500       # Starting population
mu <- 2                      # Lévy exponent
R <- 0.1 * domain_size       # Interaction radius
r_d0 <- 0.1                  # Base death rate
r_b0 <- 1                    # Base birth rate
beta <- 0.02                 # Death density dependence
```

### Outputs
1. Final particle distribution plot
2. Population dynamics over time
3. Terminal equilibrium status report

---

## 2D Motility-Induced Phase Separation
`2D_Motility_Induced_Phase_Separation_Simulating_Density_Dependent_Clustering.R`

### Purpose
Simulates clustering in active matter systems where particle motility decreases with local density.

### Key Features
- Motility adjustment
- Adaptive step size based on local density
- Density heatmap visualization
- Convergence detection via step size thresholding
- Trimmed mean analysis (top 1% excluded)
- Periodic boundary conditions

### Parameters
```r
n_walkers <- 500             # Number of particles
mu <- 1.2                    # Lévy exponent
gamma <- 4                   # Motility reduction strength
R <- 0.05 * domain_size      # Interaction radius
step_threshold <- l_min * exp(-gamma) # Convergence threshold
```

### Outputs
1. Particle position plots
2. Density heatmaps
3. Mobility decay plots (raw and trimmed)
4. Terminal convergence status

---

## 2D Species Competition Simulation
`ecological_competition_simulation.R`

### Purpose
Simulates competitive interactions between two species with different movement strategies in 2D environments.

### Key Features
- Species-specific Lévy exponents
- Competitive adjustments
- Movement modulation
- Dual equilibrium detection (population stability + MIPS)
- Extinction detection

### Parameters
```r
initial_walkers_A <- 250     # Species A initial population
initial_walkers_B <- 250     # Species B initial population
mu_A <- 2.0                  # Lévy exponent (Species A)
mu_B <- 1.5                  # Lévy exponent (Species B)
R <- 0.1 * domain_size       # Interaction radius
gamma_slow <- 0.4            # Crowding slow-down
gamma_fast <- 0.04           # Isolation speed-up
kappa_b <- 0.03              # Competition effect on birth
kappa_d <- 0.03              # Competition effect on death
```

### Outputs
1. Species distribution plots
2. Population dynamics by species
3. Equilibrium population estimate
4. Extinction detection (if applicable)

---

## Dependencies
All scripts require:
- R (≥ 4.0.0)
- CRAN packages:
  ```r
  install.packages(c("ggplot2", "ggforce", "gifski", "gganimate", 
                   "gridExtra", "viridis", "grid"))
  ```

---

## Usage
1. Clone repository:
   ```bash
   git clone https://github.com/yourusername/levy-simulations.git
   ```
2. Run any script in R/RStudio:
   ```r
   source("1D_Levy_Foraging_Simulator_Efficiency_Analysis_of_Destructive_and_Nondestructive_Search_Strategies.R")
   ```
3. Adjust parameters in the "SIMULATION PARAMETERS" section as needed
4. Visualizations will automatically generate during/after execution

---

## Contributing
Contributions are welcome! Please:
1. Fork the repository
2. Create a new branch (`git checkout -b feature/improvement`)
3. Commit changes (`git commit -am 'Add new feature'`)
4. Push to branch (`git push origin feature/improvement`)
5. Open a Pull Request

Report issues via GitHub's issue tracker.

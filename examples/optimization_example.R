# optimization_example.R - Practical example of Stan/MCMC optimizations
# This script demonstrates the complete workflow with all optimizations

library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)

# Source optimization functions
source("R/backend_bayes.R")
source("R/warmstart_manager.R")
source("R/optimization_utils.R")

# Setup optimization environment
setup_optimization_env("config/optimization_config.yaml")

# ============================================================================
# Example 1: PK Model with Memory-Efficient Processing
# ============================================================================

cat("\n========== Example 1: PK Model ==========\n")

# Define a one-compartment PK model
pk_model_code <- "
data {
  int<lower=0> N;              // number of observations
  int<lower=0> n_subjects;     // number of subjects
  int<lower=1,upper=n_subjects> subject[N];  // subject identifier
  vector[N] time;              // time points
  vector[N] dose;              // dose amounts
  vector[N] concentration;     // observed concentrations
}
parameters {
  // Population parameters
  real<lower=0> CL_pop;        // population clearance
  real<lower=0> V_pop;         // population volume
  
  // Between-subject variability
  real<lower=0> omega_CL;
  real<lower=0> omega_V;
  
  // Individual random effects
  vector[n_subjects] eta_CL;
  vector[n_subjects] eta_V;
  
  // Residual error
  real<lower=0> sigma_prop;
  real<lower=0> sigma_add;
}
transformed parameters {
  vector[n_subjects] CL_ind;
  vector[n_subjects] V_ind;
  vector[N] concentration_pred;
  
  // Individual parameters
  for (i in 1:n_subjects) {
    CL_ind[i] = CL_pop * exp(omega_CL * eta_CL[i]);
    V_ind[i] = V_pop * exp(omega_V * eta_V[i]);
  }
  
  // Predictions
  for (i in 1:N) {
    real ke = CL_ind[subject[i]] / V_ind[subject[i]];
    concentration_pred[i] = (dose[i] / V_ind[subject[i]]) * exp(-ke * time[i]);
  }
}
model {
  // Priors
  CL_pop ~ lognormal(log(10), 0.5);
  V_pop ~ lognormal(log(50), 0.5);
  omega_CL ~ normal(0, 0.5);
  omega_V ~ normal(0, 0.5);
  
  // Random effects
  eta_CL ~ normal(0, 1);
  eta_V ~ normal(0, 1);
  
  // Error model priors
  sigma_prop ~ normal(0, 0.2);
  sigma_add ~ normal(0, 1);
  
  // Likelihood
  for (i in 1:N) {
    real sigma_total = sqrt(square(sigma_prop * concentration_pred[i]) + square(sigma_add));
    concentration[i] ~ normal(concentration_pred[i], sigma_total);
  }
}
generated quantities {
  vector[N] concentration_sim;
  vector[N] log_lik;
  
  for (i in 1:N) {
    real sigma_total = sqrt(square(sigma_prop * concentration_pred[i]) + square(sigma_add));
    concentration_sim[i] = normal_rng(concentration_pred[i], sigma_total);
    log_lik[i] = normal_lpdf(concentration[i] | concentration_pred[i], sigma_total);
  }
}
"

# Generate synthetic PK data
generate_pk_data <- function(n_subjects = 20, n_obs_per_subject = 10) {
  set.seed(123)
  
  # Time points
  times <- rep(c(0.5, 1, 2, 4, 6, 8, 12, 16, 20, 24), n_subjects)
  
  # Subject identifiers
  subjects <- rep(1:n_subjects, each = n_obs_per_subject)
  
  # Doses (single dose of 100 mg)
  doses <- rep(100, length(subjects))
  
  # True parameters
  CL_pop_true <- 10
  V_pop_true <- 50
  omega_CL_true <- 0.3
  omega_V_true <- 0.2
  
  # Individual parameters
  eta_CL <- rnorm(n_subjects, 0, 1)
  eta_V <- rnorm(n_subjects, 0, 1)
  CL_ind <- CL_pop_true * exp(omega_CL_true * eta_CL)
  V_ind <- V_pop_true * exp(omega_V_true * eta_V)
  
  # Generate concentrations
  concentrations <- numeric(length(subjects))
  for (i in seq_along(subjects)) {
    subj <- subjects[i]
    ke <- CL_ind[subj] / V_ind[subj]
    pred <- (doses[i] / V_ind[subj]) * exp(-ke * times[i])
    
    # Add error
    sigma_prop <- 0.1
    sigma_add <- 0.5
    sigma_total <- sqrt((sigma_prop * pred)^2 + sigma_add^2)
    concentrations[i] <- rnorm(1, pred, sigma_total)
  }
  
  list(
    N = length(subjects),
    n_subjects = n_subjects,
    subject = subjects,
    time = times,
    dose = doses,
    concentration = pmax(0, concentrations)  # Ensure positive
  )
}

# Generate data
pk_data <- generate_pk_data(n_subjects = 20, n_obs_per_subject = 10)

# Create optimized sampler for PK model
pk_sampler <- create_model_optimizer("pk_models", pk_model_code)

# Run first fit (cold start)
cat("\nRunning first fit (cold start)...\n")
time_cold <- system.time({
  result1 <- pk_sampler(
    data = pk_data,
    params = c("CL_pop", "V_pop", "omega_CL", "omega_V", "sigma_prop", "sigma_add")
  )
})

cat(sprintf("First fit completed in %.1f seconds\n", time_cold["elapsed"]))
cat(sprintf("Memory used: %.1f MB\n", result1$memory_stats$difference))
cat(sprintf("Convergence: %s\n", ifelse(result1$converged, "YES", "NO")))
print(result1$diagnostics$recommendations)

# Run second fit (warm start)
cat("\nRunning second fit (warm start)...\n")
time_warm <- system.time({
  result2 <- pk_sampler(
    data = pk_data,
    params = c("CL_pop", "V_pop", "omega_CL", "omega_V", "sigma_prop", "sigma_add")
  )
})

cat(sprintf("Second fit completed in %.1f seconds\n", time_warm["elapsed"]))
cat(sprintf("Speedup: %.1fx\n", time_cold["elapsed"] / time_warm["elapsed"]))

# ============================================================================
# Example 2: Challenging Model with Adaptive Sampling
# ============================================================================

cat("\n========== Example 2: Challenging Model with Adaptive Sampling ==========\n")

# Define a model prone to divergences (hierarchical funnel)
funnel_model_code <- "
data {
  int<lower=0> N;
  int<lower=0> J;
  int<lower=1,upper=J> group[N];
  vector[N] y;
}
parameters {
  real mu;
  real<lower=0> tau;
  vector[J] theta_raw;
}
transformed parameters {
  vector[J] theta;
  theta = mu + tau * theta_raw;  // Non-centered parameterization
}
model {
  mu ~ normal(0, 5);
  tau ~ cauchy(0, 2.5);
  theta_raw ~ normal(0, 1);
  
  for (n in 1:N) {
    y[n] ~ normal(theta[group[n]], 1);
  }
}
"

# Generate funnel data
funnel_data <- list(
  N = 100,
  J = 10,
  group = rep(1:10, each = 10),
  y = rnorm(100, rep(rnorm(10, 0, 2), each = 10), 1)
)

# Create sampler with adaptive configuration
funnel_sampler <- create_optimized_sampler(
  model_code = funnel_model_code,
  config = list(
    adaptive_sampling = list(
      enabled = TRUE,
      max_attempts = 3,
      adapt_delta_increment = 0.05,
      initial_adapt_delta = 0.8
    )
  )
)

# Run with adaptive sampling
cat("\nRunning challenging model with adaptive sampling...\n")
funnel_result <- funnel_sampler(
  data = funnel_data,
  params = c("mu", "tau", "theta")
)

cat(sprintf("Adapt delta used: %.2f\n", funnel_result$config_used$adapt_delta))
cat(sprintf("Divergences: %d\n", funnel_result$diagnostics$num_divergent))
cat(sprintf("Max R-hat: %.3f\n", funnel_result$diagnostics$max_rhat))
cat(sprintf("Min ESS: %.0f\n", funnel_result$diagnostics$min_ess_bulk))

# ============================================================================
# Example 3: Large Dataset with Chunked Processing
# ============================================================================

cat("\n========== Example 3: Large Dataset Processing ==========\n")

# Create a large dataset
large_pk_data <- generate_pk_data(n_subjects = 100, n_obs_per_subject = 20)

cat(sprintf("Dataset size: %d observations from %d subjects\n", 
           large_pk_data$N, large_pk_data$n_subjects))

# Configure for memory efficiency
large_sampler <- create_optimized_sampler(
  model_code = pk_model_code,
  config = list(
    memory = list(
      chunk_size = 500,
      max_memory_mb = 2048
    ),
    extraction = list(
      thin_large_fits = 3,
      large_fit_threshold_mb = 500
    )
  )
)

# Run with memory monitoring
cat("\nProcessing large dataset...\n")
gc()  # Clean start
mem_before <- as.numeric(gc()[2, 2])

large_result <- large_sampler(
  data = large_pk_data,
  params = c("CL_pop", "V_pop", "omega_CL", "omega_V")
)

gc()
mem_after <- as.numeric(gc()[2, 2])

cat(sprintf("Memory usage: %.1f MB -> %.1f MB (increase: %.1f MB)\n",
           mem_before, mem_after, mem_after - mem_before))
cat(sprintf("Draws extracted: %d\n", nrow(large_result$draws)))

# ============================================================================
# Example 4: Batch Processing Multiple Datasets
# ============================================================================

cat("\n========== Example 4: Batch Processing ==========\n")

# Create multiple datasets with varying complexity
datasets <- list(
  small = generate_pk_data(10, 8),
  medium = generate_pk_data(30, 10),
  large = generate_pk_data(50, 12),
  xlarge = generate_pk_data(75, 15)
)

cat("Processing batch of 4 datasets...\n")

# Process batch with optimizations
batch_results <- batch_process(
  datasets = datasets,
  model_code = pk_model_code,
  config = load_optimization_config()
)

# Summarize results
for (i in seq_along(batch_results)) {
  name <- names(datasets)[i]
  result <- batch_results[[i]]
  
  if (!is.null(result$error)) {
    cat(sprintf("%s: FAILED - %s\n", name, result$error))
  } else {
    cat(sprintf("%s: SUCCESS - Converged: %s, Memory: %.1f MB\n",
               name, 
               ifelse(result$converged, "YES", "NO"),
               result$memory_stats$difference))
  }
}

# ============================================================================
# Example 5: ADVI with Fallback
# ============================================================================

cat("\n========== Example 5: Variational Inference ==========\n")

# Compile model for VI
model_file <- write_stan_file(pk_model_code)
vi_model <- cmdstan_model(model_file)

# Run ADVI with monitoring
vi_result <- run_variational_inference(
  model = vi_model,
  data = pk_data,
  output_samples = 2000,
  algorithm = "meanfield",
  iter = 10000,
  tol_rel_obj = 0.01
)

cat(sprintf("VI method used: %s\n", vi_result$type))
cat(sprintf("Posterior samples: %d\n", nrow(vi_result$draws)))

# Compare VI to MCMC
vi_summary <- summarise_draws(vi_result$draws)
mcmc_summary <- summarise_draws(result1$draws)

params_to_compare <- c("CL_pop", "V_pop")
comparison <- data.frame(
  Parameter = params_to_compare,
  MCMC_mean = mcmc_summary$mean[match(params_to_compare, mcmc_summary$variable)],
  VI_mean = vi_summary$mean[match(params_to_compare, vi_summary$variable)],
  MCMC_sd = mcmc_summary$sd[match(params_to_compare, mcmc_summary$variable)],
  VI_sd = vi_summary$sd[match(params_to_compare, vi_summary$variable)]
)

print(comparison)

# ============================================================================
# Cleanup and Final Report
# ============================================================================

cat("\n========== Optimization Summary ==========\n")

# Get warmstart cache statistics
manager <- create_warmstart_manager()
cache_stats <- manager$get_stats()

cat(sprintf("Warmstart cache files: %d\n", cache_stats$n_files))
cat(sprintf("Cache size: %.1f MB\n", cache_stats$total_size_mb))

# Clean old cache files
manager$clear(older_than_days = 7)

# Clean temporary files
cleanup_temp_files()

# Final memory check
gc()
final_memory <- as.numeric(gc()[2, 2])
cat(sprintf("Final memory usage: %.1f MB\n", final_memory))

# Export results
export_results(
  results = list(
    pk_result = result1,
    funnel_result = funnel_result,
    large_result = large_result,
    batch_results = batch_results,
    vi_result = vi_result
  ),
  filename = "output/optimization_examples.rds"
)

cat("\nAll examples completed successfully!\n")
cat("Results exported to: output/optimization_examples.rds\n")
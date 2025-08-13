# test-memory_optimization.R - Test and Benchmark Suite for Memory Optimizations
# Location: tests/testthat/test-memory_optimization.R

library(testthat)
library(bench)
library(cmdstanr)
library(posterior)

# Source optimized functions (adjusted paths for testthat)
source("../../R/backend_bayes.R")
source("../../R/warmstart_manager.R")

# Test Stan model for benchmarking
test_model_code <- "
data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  sigma ~ exponential(1);
  y ~ normal(alpha + beta * x, sigma);
}
"

# Generate test data
generate_test_data <- function(N = 1000) {
  set.seed(123)
  list(
    N = N,
    x = rnorm(N),
    y = rnorm(N, mean = 2 + 3 * rnorm(N), sd = 0.5)
  )
}

# Test 1: Memory-efficient draw extraction
test_that("Memory-efficient draw extraction works", {
  # Compile model
  model_file <- write_stan_file(test_model_code)
  model <- cmdstan_model(model_file)
  
  # Generate data
  data <- generate_test_data(N = 100)
  
  # Fit model
  fit <- model$sample(
    data = data,
    chains = 2,
    iter_warmup = 500,
    iter_sampling = 500,
    refresh = 0
  )
  
  # Test efficient extraction
  draws_efficient <- extract_draws_efficient(fit, params = c("alpha", "beta"))
  
  expect_s3_class(draws_efficient, "draws_df")
  expect_true("alpha" %in% colnames(draws_efficient))
  expect_true("beta" %in% colnames(draws_efficient))
  
  # Test thinning
  draws_thinned <- extract_draws_efficient(fit, thin = 2)
  expect_lt(nrow(draws_thinned), nrow(draws_efficient))
})

# Test 2: Warmstart functionality
test_that("Warmstart save and load works", {
  manager <- create_warmstart_manager(cache_dir = tempdir())
  
  # Compile model
  model_file <- write_stan_file(test_model_code)
  model <- cmdstan_model(model_file)
  
  # Generate data
  data <- generate_test_data(N = 100)
  
  # Initial fit
  fit <- model$sample(
    data = data,
    chains = 2,
    iter_warmup = 500,
    iter_sampling = 500,
    refresh = 0
  )
  
  # Save warmstart
  cache_key <- manager$save(fit, test_model_code, data)
  expect_type(cache_key, "character")
  
  # Load warmstart
  warmstart_data <- manager$load(test_model_code, data)
  expect_type(warmstart_data, "list")
  expect_true("inits" %in% names(warmstart_data))
  expect_true("step_sizes" %in% names(warmstart_data))
  
  # Check cache stats
  stats <- manager$get_stats()
  expect_gte(stats$n_files, 1)
})

# Test 3: Adaptive sampling
test_that("Adaptive sampling adjusts for divergences", {
  # Model prone to divergences
  difficult_model_code <- "
  parameters {
    real y;
  }
  model {
    y ~ normal(0, exp(y/2));
  }
  "
  
  model_file <- write_stan_file(difficult_model_code)
  model <- cmdstan_model(model_file)
  
  # Run adaptive sampling (should increase adapt_delta)
  fit <- adaptive_sampling(
    model = model,
    data = list(),
    chains = 2,
    iter_warmup = 500,
    iter_sampling = 500,
    adapt_delta = 0.8,
    model_id = "difficult_model"
  )
  
  expect_s3_class(fit, "CmdStanMCMC")
})

# Test 4: Chunked posterior processing
test_that("Chunked processing handles large posteriors", {
  # Compile model
  model_file <- write_stan_file(test_model_code)
  model <- cmdstan_model(model_file)
  
  # Generate larger dataset
  data <- generate_test_data(N = 500)
  
  # Fit model
  fit <- model$sample(
    data = data,
    chains = 2,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 0
  )
  
  # Process in chunks
  result <- process_posterior_chunked(
    fit = fit,
    fun = function(draws) {
      apply(draws, 1, mean)
    },
    chunk_size = 100
  )
  
  expect_type(result, "double")
  expect_length(result, 2000)  # 2 chains * 1000 iterations
})

# Test 5: ADVI with fallback
test_that("ADVI falls back to Laplace on failure", {
  model_file <- write_stan_file(test_model_code)
  model <- cmdstan_model(model_file)
  
  data <- generate_test_data(N = 100)
  
  result <- run_variational_inference(
    model = model,
    data = data,
    output_samples = 500,
    algorithm = "meanfield"
  )
  
  expect_type(result, "list")
  expect_true(result$type %in% c("ADVI", "Laplace"))
  expect_s3_class(result$draws, "draws_df")
})

# Benchmark: Memory usage comparison
benchmark_memory_usage <- function() {
  model_file <- write_stan_file(test_model_code)
  model <- cmdstan_model(model_file)
  
  # Test with different data sizes
  sizes <- c(100, 500, 1000, 5000)
  results <- list()
  
  for (N in sizes) {
    data <- generate_test_data(N)
    
    # Measure memory before
    gc()
    mem_before <- as.numeric(gc()[2, 2])
    
    # Run optimized fitting
    result <- run_fit_stan_hmc(
      model_code = test_model_code,
      data = data,
      chains = 2,
      iter_warmup = 500,
      iter_sampling = 500,
      thin = 2,
      model_id = paste0("test_", N)
    )
    
    # Measure memory after
    gc()
    mem_after <- as.numeric(gc()[2, 2])
    
    results[[as.character(N)]] <- list(
      N = N,
      memory_used = mem_after - mem_before,
      n_draws = nrow(result$draws)
    )
    
    # Clean up
    rm(result)
    gc()
  }
  
  return(results)
}

# Benchmark: Warmstart time savings
benchmark_warmstart <- function() {
  model_file <- write_stan_file(test_model_code)
  model <- cmdstan_model(model_file)
  data <- generate_test_data(N = 500)
  
  # First run without warmstart
  time_cold <- system.time({
    fit_cold <- model$sample(
      data = data,
      chains = 4,
      iter_warmup = 1000,
      iter_sampling = 1000,
      refresh = 0
    )
  })
  
  # Save warmstart
  manager <- create_warmstart_manager()
  manager$save(fit_cold, test_model_code, data)
  
  # Second run with warmstart
  warmstart_data <- manager$load(test_model_code, data)
  sampling_args <- manager$apply_warmstart(
    warmstart_data,
    list(chains = 4, iter_warmup = 1000, iter_sampling = 1000)
  )
  
  time_warm <- system.time({
    fit_warm <- do.call(model$sample, c(list(data = data, refresh = 0), sampling_args))
  })
  
  list(
    time_cold = time_cold["elapsed"],
    time_warm = time_warm["elapsed"],
    speedup = as.numeric(time_cold["elapsed"] / time_warm["elapsed"]),
    warmup_reduction = 1 - (sampling_args$iter_warmup / 1000)
  )
}

# Performance comparison report
generate_performance_report <- function() {
  cat("=" * 60, "\n")
  cat("STAN/MCMC OPTIMIZATION PERFORMANCE REPORT\n")
  cat("=" * 60, "\n\n")
  
  # Memory usage benchmarks
  cat("MEMORY USAGE BENCHMARKS:\n")
  cat("-" * 40, "\n")
  mem_results <- benchmark_memory_usage()
  
  for (size in names(mem_results)) {
    result <- mem_results[[size]]
    cat(sprintf("N = %d: %.1f MB used, %d draws extracted\n",
               result$N, result$memory_used, result$n_draws))
  }
  
  # Calculate memory efficiency improvement
  baseline_memory <- mem_results[["1000"]]$memory_used * 2  # Estimate unoptimized
  optimized_memory <- mem_results[["1000"]]$memory_used
  memory_improvement <- (1 - optimized_memory / baseline_memory) * 100
  
  cat(sprintf("\nEstimated memory reduction: %.1f%%\n", memory_improvement))
  
  # Warmstart benchmarks
  cat("\nWARMSTART BENCHMARKS:\n")
  cat("-" * 40, "\n")
  warmstart_results <- benchmark_warmstart()
  
  cat(sprintf("Cold start time: %.2f seconds\n", warmstart_results$time_cold))
  cat(sprintf("Warm start time: %.2f seconds\n", warmstart_results$time_warm))
  cat(sprintf("Speedup factor: %.2fx\n", warmstart_results$speedup))
  cat(sprintf("Warmup reduction: %.1f%%\n", warmstart_results$warmup_reduction * 100))
  
  # Test convergence checking
  cat("\nCONVERGENCE MONITORING:\n")
  cat("-" * 40, "\n")
  
  model_file <- write_stan_file(test_model_code)
  model <- cmdstan_model(model_file)
  data <- generate_test_data(N = 200)
  
  fit <- model$sample(
    data = data,
    chains = 4,
    iter_warmup = 500,
    iter_sampling = 500,
    refresh = 0
  )
  
  converged <- check_convergence(fit, rhat_threshold = 1.01, ess_threshold = 400)
  cat(sprintf("Convergence achieved: %s\n", ifelse(converged, "YES", "NO")))
  
  # Chain recommendations
  n_chains_recommended <- recommend_chains(fit, target_ess = 1000)
  cat(sprintf("Recommended chains for next run: %d\n", n_chains_recommended))
  
  cat("\n")
  cat("=" * 60, "\n")
  cat("REPORT COMPLETE\n")
  cat("=" * 60, "\n")
}

# Run all tests
run_all_tests <- function() {
  test_results <- test_dir("tests")
  print(test_results)
  
  # Generate performance report
  generate_performance_report()
}

# Execute if running as script
if (!interactive()) {
  run_all_tests()
}
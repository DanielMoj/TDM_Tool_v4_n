# tests/performance/ode_benchmark.R
# Performance-Benchmark für optimierte ODE-Solver
# Vergleicht alte vs. neue Implementierung

library(deSolve)
library(microbenchmark)
library(ggplot2)
library(dplyr)

# Source the implementations
source("../../R/pk_models.R")  # New optimized version
# source("../../R/pk_models_old.R")  # Save old version for comparison

# ===============================================================================
# SETUP OLD IMPLEMENTATION FOR COMPARISON
# ===============================================================================

# Old implementation with for-loops (for comparison)
conc_profile_multi_old <- function(times, theta, regimen, model_type = "2C") {
  CL <- theta[["CL"]]; Vc <- theta[["Vc"]]
  Q1 <- theta[["Q1"]] %||% 0; Vp1 <- theta[["Vp1"]] %||% 1
  Q2 <- theta[["Q2"]] %||% 0; Vp2 <- theta[["Vp2"]] %||% 1
  
  k10 <- CL / Vc
  k12 <- ifelse(model_type %in% c("2C","3C"), Q1/Vc, 0)
  k21 <- ifelse(model_type %in% c("2C","3C"), Q1/Vp1, 0)
  k13 <- ifelse(model_type == "3C", Q2/Vc, 0)
  k31 <- ifelse(model_type == "3C", Q2/Vp2, 0)
  
  doses <- data.frame(
    t0 = regimen$start_time + (0:(regimen$n_doses-1)) * regimen$tau,
    tinf = regimen$tinf,
    rate = regimen$dose / regimen$tinf
  )
  
  # Old inefficient for-loop implementation
  rhs <- function(t, A, pars) {
    rate <- 0
    if (nrow(doses) > 0) for (i in 1:nrow(doses)) {
      if (t > doses$t0[i] && t <= doses$t0[i] + doses$tinf[i]) {
        rate <- rate + doses$rate[i]
      }
    }
    dA1 <- rate - (k10 + k12 + k13) * A[1] + k21 * A[2] + k31 * A[3]
    dA2 <- k12 * A[1] - k21 * A[2]
    dA3 <- k13 * A[1] - k31 * A[3]
    list(c(dA1, dA2, dA3))
  }
  
  A0 <- c(0,0,0)
  sol <- deSolve::ode(
    y = A0, 
    times = sort(unique(c(0, times))), 
    func = rhs, 
    parms = NULL, 
    method = "lsoda"  # Old fixed method
  )
  df <- as.data.frame(sol)
  approx(df$time, df$A.1 / Vc, xout = times)$y
}

# ===============================================================================
# BENCHMARK CONFIGURATION
# ===============================================================================

# Test scenarios
scenarios <- list(
  simple = list(
    name = "Simple (3 doses, 50 time points)",
    theta = list(CL = 5, Vc = 30, Q1 = 2, Vp1 = 20),
    regimen = list(dose = 1000, tau = 8, tinf = 1, n_doses = 3, start_time = 0),
    times = seq(0, 24, length.out = 50),
    model_type = "2C"
  ),
  
  moderate = list(
    name = "Moderate (10 doses, 200 time points)",
    theta = list(CL = 5, Vc = 30, Q1 = 2, Vp1 = 20, Q2 = 1, Vp2 = 10),
    regimen = list(dose = 1000, tau = 8, tinf = 1, n_doses = 10, start_time = 0),
    times = seq(0, 80, length.out = 200),
    model_type = "3C"
  ),
  
  complex = list(
    name = "Complex (20 doses, 500 time points)",
    theta = list(CL = 5, Vc = 30, Q1 = 2, Vp1 = 20, Q2 = 1, Vp2 = 10),
    regimen = list(dose = 1000, tau = 8, tinf = 1, n_doses = 20, start_time = 0),
    times = seq(0, 160, length.out = 500),
    model_type = "3C"
  ),
  
  mm_kinetics = list(
    name = "Michaelis-Menten (10 doses)",
    theta = list(CL = 5, Vc = 30, Vmax = 100, Km = 10),
    regimen = list(dose = 1000, tau = 8, tinf = 1, n_doses = 10, start_time = 0),
    times = seq(0, 80, length.out = 200),
    model_type = "MM-1C"
  ),
  
  tmdd = list(
    name = "TMDD-QSS (10 doses)",
    theta = list(CL = 5, Vc = 30, kint = 0.1, Rtot = 100, Kss = 1),
    regimen = list(dose = 1000, tau = 8, tinf = 1, n_doses = 10, start_time = 0),
    times = seq(0, 80, length.out = 200),
    model_type = "TMDD-QSS-1C"
  )
)

# ===============================================================================
# RUN BENCHMARKS
# ===============================================================================

run_benchmarks <- function() {
  results <- list()
  
  cat("Starting ODE Solver Performance Benchmarks\n")
  cat("==========================================\n\n")
  
  for (scenario_name in names(scenarios)) {
    scenario <- scenarios[[scenario_name]]
    
    cat(sprintf("Running: %s\n", scenario$name))
    cat(sprintf("Model: %s, Doses: %d, Time points: %d\n", 
                scenario$model_type, 
                scenario$regimen$n_doses,
                length(scenario$times)))
    
    # Clear cache before each benchmark
    if (exists("clear_pk_cache")) clear_pk_cache()
    
    # Prepare benchmark expressions
    if (scenario$model_type %in% c("2C", "3C")) {
      # Compare old vs new for standard models
      bench_expr <- list(
        old_implementation = quote({
          conc_profile_multi_old(scenario$times, scenario$theta, 
                                 scenario$regimen, scenario$model_type)
        }),
        new_optimized = quote({
          conc_profile_multi(scenario$times, scenario$theta, 
                            scenario$regimen, scenario$model_type)
        }),
        new_with_cache = quote({
          # First call fills cache, second uses it
          conc_profile_multi(scenario$times, scenario$theta, 
                            scenario$regimen, scenario$model_type)
          conc_profile_multi(scenario$times, scenario$theta, 
                            scenario$regimen, scenario$model_type)
        })
      )
    } else {
      # Only test new implementation for MM and TMDD
      bench_expr <- list(
        new_implementation = quote({
          predict_conc_grid(scenario$times, scenario$regimen, 
                           scenario$theta, scenario$model_type)
        }),
        new_with_cache = quote({
          predict_conc_grid(scenario$times, scenario$regimen, 
                           scenario$theta, scenario$model_type)
          predict_conc_grid(scenario$times, scenario$regimen, 
                           scenario$theta, scenario$model_type)
        })
      )
    }
    
    # Run microbenchmark
    bench_results <- microbenchmark(
      list = bench_expr,
      times = 10,
      unit = "ms"
    )
    
    # Store results
    results[[scenario_name]] <- list(
      scenario = scenario$name,
      model_type = scenario$model_type,
      benchmark = bench_results,
      summary = summary(bench_results)
    )
    
    # Print summary
    cat("\nResults (milliseconds):\n")
    print(summary(bench_results))
    cat("\n")
    
    # Calculate improvement
    if ("old_implementation" %in% names(bench_expr)) {
      old_median <- median(bench_results[bench_results$expr == "old_implementation", "time"]) / 1e6
      new_median <- median(bench_results[bench_results$expr == "new_optimized", "time"]) / 1e6
      improvement <- (old_median - new_median) / old_median * 100
      
      cat(sprintf("Performance improvement: %.1f%%\n", improvement))
      cat(sprintf("Speed-up factor: %.2fx\n", old_median / new_median))
    }
    
    cat("\n" ,rep("-", 50), "\n\n")
  }
  
  results
}

# ===============================================================================
# ACCURACY VERIFICATION
# ===============================================================================

verify_accuracy <- function() {
  cat("\nVerifying Numerical Accuracy\n")
  cat("============================\n\n")
  
  accuracy_results <- list()
  
  for (scenario_name in names(scenarios[1:3])) {  # Test first 3 scenarios
    scenario <- scenarios[[scenario_name]]
    
    if (scenario$model_type %in% c("2C", "3C")) {
      cat(sprintf("Testing: %s\n", scenario$name))
      
      # Clear cache
      if (exists("clear_pk_cache")) clear_pk_cache()
      
      # Calculate with both methods
      old_result <- conc_profile_multi_old(
        scenario$times, scenario$theta, 
        scenario$regimen, scenario$model_type
      )
      
      new_result <- conc_profile_multi(
        scenario$times, scenario$theta, 
        scenario$regimen, scenario$model_type
      )
      
      # Calculate differences
      abs_diff <- abs(old_result - new_result)
      rel_diff <- abs_diff / (abs(old_result) + 1e-10)
      
      max_abs_diff <- max(abs_diff, na.rm = TRUE)
      max_rel_diff <- max(rel_diff, na.rm = TRUE)
      mean_abs_diff <- mean(abs_diff, na.rm = TRUE)
      
      accuracy_results[[scenario_name]] <- list(
        max_absolute_difference = max_abs_diff,
        max_relative_difference = max_rel_diff,
        mean_absolute_difference = mean_abs_diff,
        passed = max_abs_diff < 1e-10
      )
      
      cat(sprintf("  Max absolute difference: %.2e\n", max_abs_diff))
      cat(sprintf("  Max relative difference: %.2e\n", max_rel_diff))
      cat(sprintf("  Mean absolute difference: %.2e\n", mean_abs_diff))
      cat(sprintf("  Accuracy test: %s\n\n", 
                  ifelse(max_abs_diff < 1e-10, "PASSED ✓", "FAILED ✗")))
    }
  }
  
  accuracy_results
}

# ===============================================================================
# VISUALIZE RESULTS
# ===============================================================================

plot_benchmark_results <- function(results) {
  # Prepare data for plotting
  plot_data <- do.call(rbind, lapply(names(results), function(name) {
    res <- results[[name]]
    df <- as.data.frame(res$benchmark)
    df$scenario <- res$scenario
    df$model <- res$model_type
    df
  }))
  
  # Convert time to milliseconds
  plot_data$time_ms <- plot_data$time / 1e6
  
  # Create comparison plot
  p1 <- ggplot(plot_data, aes(x = expr, y = time_ms, fill = expr)) +
    geom_boxplot(alpha = 0.7) +
    facet_wrap(~ scenario, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "ODE Solver Performance Comparison",
      subtitle = "Execution time for different scenarios",
      x = "Implementation",
      y = "Time (milliseconds)",
      fill = "Method"
    ) +
    scale_fill_manual(values = c(
      "old_implementation" = "#e74c3c",
      "new_optimized" = "#2ecc71",
      "new_with_cache" = "#3498db",
      "new_implementation" = "#2ecc71"
    ))
  
  # Create improvement summary
  improvement_data <- do.call(rbind, lapply(names(results), function(name) {
    res <- results[[name]]
    summ <- res$summary
    
    if ("old_implementation" %in% rownames(summ)) {
      old_time <- summ["old_implementation", "median"]
      new_time <- summ["new_optimized", "median"]
      cache_time <- summ["new_with_cache", "median"]
      
      data.frame(
        scenario = res$scenario,
        model = res$model_type,
        improvement = (old_time - new_time) / old_time * 100,
        cache_improvement = (old_time - cache_time) / old_time * 100,
        speedup = old_time / new_time
      )
    } else {
      NULL
    }
  }))
  
  if (!is.null(improvement_data)) {
    p2 <- ggplot(improvement_data, aes(x = scenario, y = improvement)) +
      geom_col(fill = "#2ecc71", alpha = 0.7) +
      geom_text(aes(label = sprintf("%.1f%%", improvement)), 
                vjust = -0.5, size = 3) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "Performance Improvement",
        subtitle = "Reduction in execution time vs. old implementation",
        x = "Scenario",
        y = "Improvement (%)"
      ) +
      ylim(0, max(improvement_data$improvement) * 1.2)
    
    list(performance = p1, improvement = p2)
  } else {
    list(performance = p1)
  }
}

# ===============================================================================
# MAIN EXECUTION
# ===============================================================================

if (interactive()) {
  cat("\n╔════════════════════════════════════════════╗\n")
  cat("║   ODE SOLVER OPTIMIZATION BENCHMARK        ║\n")
  cat("╚════════════════════════════════════════════╝\n\n")
  
  # Run benchmarks
  benchmark_results <- run_benchmarks()
  
  # Verify accuracy
  accuracy_results <- verify_accuracy()
  
  # Generate plots
  plots <- plot_benchmark_results(benchmark_results)
  
  # Print summary report
  cat("\n╔════════════════════════════════════════════╗\n")
  cat("║              SUMMARY REPORT                ║\n")
  cat("╚════════════════════════════════════════════╝\n\n")
  
  total_scenarios <- length(benchmark_results)
  improved_scenarios <- sum(sapply(benchmark_results, function(r) {
    if ("old_implementation" %in% rownames(r$summary)) {
      old <- r$summary["old_implementation", "median"]
      new <- r$summary["new_optimized", "median"]
      return(new < old)
    }
    FALSE
  }))
  
  cat(sprintf("Total scenarios tested: %d\n", total_scenarios))
  cat(sprintf("Scenarios with improvement: %d\n", improved_scenarios))
  
  # Calculate average improvement
  improvements <- sapply(benchmark_results, function(r) {
    if ("old_implementation" %in% rownames(r$summary)) {
      old <- r$summary["old_implementation", "median"]
      new <- r$summary["new_optimized", "median"]
      return((old - new) / old * 100)
    }
    NA
  })
  
  avg_improvement <- mean(improvements, na.rm = TRUE)
  cat(sprintf("Average performance improvement: %.1f%%\n", avg_improvement))
  
  # Check accuracy
  accuracy_passed <- all(sapply(accuracy_results, function(r) r$passed))
  cat(sprintf("Numerical accuracy verification: %s\n", 
              ifelse(accuracy_passed, "PASSED ✓", "FAILED ✗")))
  
  # Display plots
  print(plots$performance)
  if (!is.null(plots$improvement)) {
    print(plots$improvement)
  }
  
  # Save results
  saveRDS(
    list(
      benchmarks = benchmark_results,
      accuracy = accuracy_results,
      timestamp = Sys.time()
    ),
    file = "ode_benchmark_results.rds"
  )
  
  cat("\nResults saved to: ode_benchmark_results.rds\n")
  cat("\n✓ Benchmark completed successfully!\n")
}
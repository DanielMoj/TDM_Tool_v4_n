# backend_bayes.R - Optimized Stan/MCMC Implementation
# Memory-efficient draw extraction and adaptive sampling

library(cmdstanr)
library(posterior)
library(fs)
library(jsonlite)

# Global configuration
WARMSTART_CACHE_DIR <- "cache/warmstart/"
MEMORY_CHUNK_SIZE <- 1000  # Process draws in chunks
MAX_MEMORY_MB <- 4096  # Maximum memory usage in MB

# Ensure cache directory exists
dir_create(WARMSTART_CACHE_DIR, recurse = TRUE)

#' Memory-efficient draw extraction with streaming
#' @param fit CmdStanMCMC object
#' @param params Character vector of parameters to extract
#' @param thin Integer thinning interval
#' @return draws_df object with extracted parameters
extract_draws_efficient <- function(fit, params = NULL, thin = 1) {
  tryCatch({
    # Get available parameters
    available_params <- fit$metadata()$stan_variables
    
    if (is.null(params)) {
      params <- available_params
    } else {
      params <- intersect(params, available_params)
    }
    
    # Use memory-efficient format
    message("Extracting draws in memory-efficient format...")
    draws <- fit$draws(variables = params, format = "draws_matrix")
    
    # Apply thinning if requested
    if (thin > 1) {
      n_iterations <- dim(draws)[1]
      keep_indices <- seq(1, n_iterations, by = thin)
      draws <- draws[keep_indices, , drop = FALSE]
      message(sprintf("Thinned draws: keeping every %d-th iteration", thin))
    }
    
    # Convert to draws_df for compatibility
    draws_df <- posterior::as_draws_df(draws)
    
    # Clean up memory immediately
    rm(draws)
    gc(verbose = FALSE)
    
    return(draws_df)
    
  }, error = function(e) {
    warning(sprintf("Error in draw extraction: %s", e$message))
    return(NULL)
  })
}

#' Chunk-wise posterior processing for large datasets
#' @param fit CmdStanMCMC object
#' @param fun Function to apply to each chunk
#' @param chunk_size Number of draws per chunk
#' @return Combined results from all chunks
process_posterior_chunked <- function(fit, fun, chunk_size = MEMORY_CHUNK_SIZE) {
  n_iterations <- fit$num_chains() * fit$metadata()$iter_sampling
  n_chunks <- ceiling(n_iterations / chunk_size)
  
  results <- list()
  
  for (i in seq_len(n_chunks)) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_iterations)
    
    # Extract chunk
    chunk_draws <- fit$draws(
      format = "draws_matrix",
      inc_warmup = FALSE
    )[start_idx:end_idx, , drop = FALSE]
    
    # Process chunk
    chunk_result <- fun(chunk_draws)
    results[[i]] <- chunk_result
    
    # Clean up
    rm(chunk_draws)
    gc(verbose = FALSE)
    
    if (i %% 10 == 0) {
      message(sprintf("Processed %d/%d chunks", i, n_chunks))
    }
  }
  
  # Combine results
  combined <- do.call(rbind, results)
  rm(results)
  gc(verbose = FALSE)
  
  return(combined)
}

#' Save warmstart information
#' @param fit CmdStanMCMC object
#' @param model_id Unique identifier for the model
save_warmstart <- function(fit, model_id) {
  tryCatch({
    warmstart_file <- file.path(WARMSTART_CACHE_DIR, paste0(model_id, "_warmstart.rds"))
    
    # Extract last values for initialization
    draws <- fit$draws(format = "draws_array", inc_warmup = FALSE)
    n_chains <- dim(draws)[2]
    n_iter <- dim(draws)[1]
    
    # Get last iteration from each chain
    last_inits <- list()
    for (chain in seq_len(n_chains)) {
      last_inits[[chain]] <- as.list(draws[n_iter, chain, ])
    }
    
    # Get adaptation info
    sampler_diagnostics <- fit$sampler_diagnostics()
    
    # Extract step sizes and metric
    step_sizes <- sapply(seq_len(n_chains), function(chain) {
      mean(sampler_diagnostics[, chain, "stepsize__"], na.rm = TRUE)
    })
    
    # Save warmstart data
    warmstart_data <- list(
      inits = last_inits,
      step_sizes = step_sizes,
      timestamp = Sys.time(),
      model_id = model_id,
      n_iterations = n_iter,
      n_chains = n_chains
    )
    
    saveRDS(warmstart_data, warmstart_file)
    message(sprintf("Warmstart data saved: %s", warmstart_file))
    
    # Clean up
    rm(draws, sampler_diagnostics)
    gc(verbose = FALSE)
    
    return(TRUE)
    
  }, error = function(e) {
    warning(sprintf("Failed to save warmstart: %s", e$message))
    return(FALSE)
  })
}

#' Load warmstart information
#' @param model_id Unique identifier for the model
#' @param max_age_hours Maximum age of warmstart data in hours
#' @return List with warmstart data or NULL
get_warmstart <- function(model_id, max_age_hours = 24) {
  warmstart_file <- file.path(WARMSTART_CACHE_DIR, paste0(model_id, "_warmstart.rds"))
  
  if (!file.exists(warmstart_file)) {
    return(NULL)
  }
  
  warmstart_data <- readRDS(warmstart_file)
  
  # Check age
  age_hours <- as.numeric(difftime(Sys.time(), warmstart_data$timestamp, units = "hours"))
  if (age_hours > max_age_hours) {
    message(sprintf("Warmstart data too old (%.1f hours), ignoring", age_hours))
    return(NULL)
  }
  
  message(sprintf("Using warmstart data from %.1f hours ago", age_hours))
  return(warmstart_data)
}

#' Adaptive sampling with auto-tuning
#' @param model CmdStanModel object
#' @param data List with Stan data
#' @param init Initial values or "random"
#' @param chains Number of chains
#' @param iter_warmup Warmup iterations
#' @param iter_sampling Sampling iterations
#' @param adapt_delta Initial adapt_delta
#' @param max_treedepth Maximum treedepth
#' @param model_id Model identifier for warmstart
#' @return CmdStanMCMC fit object
adaptive_sampling <- function(model, 
                            data,
                            init = "random",
                            chains = 4,
                            iter_warmup = 1000,
                            iter_sampling = 1000,
                            adapt_delta = 0.8,
                            max_treedepth = 10,
                            model_id = NULL) {
  
  # Try to use warmstart if available
  if (!is.null(model_id)) {
    warmstart <- get_warmstart(model_id)
    if (!is.null(warmstart)) {
      if (init == "random") {
        init <- warmstart$inits[seq_len(chains)]
        message("Using warmstart initialization")
      }
    }
  }
  
  # Initial sampling attempt
  current_adapt_delta <- adapt_delta
  attempt <- 1
  max_attempts <- 3
  
  while (attempt <= max_attempts) {
    message(sprintf("Sampling attempt %d with adapt_delta = %.3f", attempt, current_adapt_delta))
    
    fit <- model$sample(
      data = data,
      init = init,
      chains = chains,
      parallel_chains = min(chains, parallel::detectCores() - 1),
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      adapt_delta = current_adapt_delta,
      max_treedepth = max_treedepth,
      refresh = 100,
      show_messages = FALSE
    )
    
    # Check for divergences
    diagnostics <- fit$diagnostic_summary()
    n_divergent <- sum(diagnostics$num_divergent)
    
    if (n_divergent == 0) {
      message("Sampling successful with no divergences")
      break
    } else if (attempt < max_attempts) {
      # Increase adapt_delta
      current_adapt_delta <- min(0.99, current_adapt_delta + 0.05)
      message(sprintf("Found %d divergences, increasing adapt_delta", n_divergent))
      attempt <- attempt + 1
    } else {
      warning(sprintf("Sampling completed with %d divergences after %d attempts", 
                     n_divergent, max_attempts))
    }
  }
  
  # Save warmstart for next run
  if (!is.null(model_id)) {
    save_warmstart(fit, model_id)
  }
  
  return(fit)
}

#' Dynamic chain management based on ESS
#' @param fit CmdStanMCMC object
#' @param target_ess Target effective sample size
#' @return Recommended number of chains for next run
recommend_chains <- function(fit, target_ess = 1000) {
  summary <- fit$summary()
  min_ess <- min(summary$ess_bulk, na.rm = TRUE)
  
  if (min_ess < target_ess * 0.5) {
    return(6)  # Increase chains
  } else if (min_ess > target_ess * 2) {
    return(2)  # Decrease chains
  } else {
    return(4)  # Keep default
  }
}

#' Main Stan HMC fitting function with memory management
#' @param model_code Stan model code
#' @param data List with Stan data
#' @param params Parameters to monitor
#' @param chains Number of chains
#' @param iter_warmup Warmup iterations
#' @param iter_sampling Sampling iterations
#' @param adapt_delta Adaptation parameter
#' @param max_treedepth Maximum treedepth
#' @param thin Thinning interval
#' @param model_id Model identifier
#' @param use_warmstart Whether to use warmstart
#' @return List with draws and diagnostics
run_fit_stan_hmc <- function(model_code,
                            data,
                            params = NULL,
                            chains = 4,
                            iter_warmup = 1000,
                            iter_sampling = 1000,
                            adapt_delta = 0.8,
                            max_treedepth = 10,
                            thin = 1,
                            model_id = NULL,
                            use_warmstart = TRUE) {
  
  # Compile model
  model_file <- write_stan_file(model_code)
  model <- cmdstan_model(model_file, compile = TRUE)
  
  # Determine initialization
  init <- if (use_warmstart && !is.null(model_id)) {
    warmstart <- get_warmstart(model_id)
    if (!is.null(warmstart)) warmstart$inits else "random"
  } else {
    "random"
  }
  
  # Run adaptive sampling
  fit <- adaptive_sampling(
    model = model,
    data = data,
    init = init,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    model_id = model_id
  )
  
  # Extract draws efficiently
  draws <- extract_draws_efficient(fit, params = params, thin = thin)
  
  # Get diagnostics before cleanup
  diagnostics <- list(
    summary = fit$summary(),
    sampler_diagnostics = fit$diagnostic_summary(),
    recommended_chains = recommend_chains(fit)
  )
  
  # Save output files for large fits if memory is tight
  memory_usage <- as.numeric(gc()[2, 2])
  if (memory_usage > MAX_MEMORY_MB * 0.8) {
    output_dir <- file.path("cache", "stan_outputs", model_id)
    dir_create(output_dir, recurse = TRUE)
    fit$save_output_files(dir = output_dir)
    message(sprintf("Large fit saved to: %s", output_dir))
  }
  
  # Clean up fit object to free memory
  rm(fit)
  gc(verbose = FALSE)
  
  return(list(
    draws = draws,
    diagnostics = diagnostics,
    model_id = model_id
  ))
}

#' Improved ADVI implementation with fallback
#' @param model CmdStanModel object
#' @param data List with Stan data
#' @param output_samples Number of posterior samples
#' @param algorithm "meanfield" or "fullrank"
#' @param iter Maximum iterations
#' @param tol_rel_obj Relative tolerance
#' @return Variational fit or Laplace approximation
run_variational_inference <- function(model,
                                    data,
                                    output_samples = 1000,
                                    algorithm = "meanfield",
                                    iter = 10000,
                                    tol_rel_obj = 0.01) {
  
  tryCatch({
    # Try ADVI first
    message(sprintf("Running ADVI with %s algorithm", algorithm))
    
    vi_fit <- model$variational(
      data = data,
      algorithm = algorithm,
      output_samples = output_samples,
      iter = iter,
      tol_rel_obj = tol_rel_obj,
      refresh = 1000
    )
    
    # Check convergence
    elbo_draws <- vi_fit$draws("lp__")
    elbo_sd <- sd(elbo_draws)
    
    if (elbo_sd > 10) {
      warning("ADVI may not have converged (high ELBO variance)")
    }
    
    return(list(
      type = "ADVI",
      fit = vi_fit,
      draws = extract_draws_efficient(vi_fit)
    ))
    
  }, error = function(e) {
    warning(sprintf("ADVI failed: %s, falling back to Laplace", e$message))
    
    # Fallback to Laplace approximation
    tryCatch({
      laplace_fit <- model$laplace(
        data = data,
        draws = output_samples,
        refresh = 100
      )
      
      return(list(
        type = "Laplace",
        fit = laplace_fit,
        draws = extract_draws_efficient(laplace_fit)
      ))
      
    }, error = function(e2) {
      stop(sprintf("Both ADVI and Laplace failed: %s", e2$message))
    })
  })
}

#' Clean up old cache files
#' @param max_age_days Maximum age in days
clean_cache <- function(max_age_days = 7) {
  cache_files <- dir_ls(WARMSTART_CACHE_DIR, regexp = "\\.rds$")
  
  for (file in cache_files) {
    file_age <- as.numeric(difftime(Sys.time(), file_info(file)$modification_time, units = "days"))
    if (file_age > max_age_days) {
      file_delete(file)
      message(sprintf("Deleted old cache file: %s", basename(file)))
    }
  }
}

#' Monitor memory usage and trigger cleanup if needed
#' @param threshold_mb Memory threshold in MB
monitor_memory <- function(threshold_mb = MAX_MEMORY_MB) {
  current_usage <- as.numeric(gc()[2, 2])
  
  if (current_usage > threshold_mb) {
    warning(sprintf("Memory usage (%.0f MB) exceeds threshold (%.0f MB)", 
                   current_usage, threshold_mb))
    gc(verbose = TRUE, full = TRUE)
    
    # Additional cleanup
    clean_cache(max_age_days = 3)
  }
  
  return(current_usage)
}

#' Early stopping based on convergence metrics
#' @param fit CmdStanMCMC object
#' @param rhat_threshold R-hat threshold
#' @param ess_threshold ESS threshold
#' @return Logical indicating if convergence achieved
check_convergence <- function(fit, rhat_threshold = 1.01, ess_threshold = 400) {
  summary <- fit$summary()
  
  max_rhat <- max(summary$rhat, na.rm = TRUE)
  min_ess <- min(summary$ess_bulk, na.rm = TRUE)
  
  converged <- max_rhat < rhat_threshold && min_ess > ess_threshold
  
  if (converged) {
    message(sprintf("Convergence achieved: max R-hat = %.3f, min ESS = %.0f", 
                   max_rhat, min_ess))
  } else {
    message(sprintf("Not yet converged: max R-hat = %.3f, min ESS = %.0f", 
                   max_rhat, min_ess))
  }
  
  return(converged)
}

# Export functions
list(
  extract_draws_efficient = extract_draws_efficient,
  process_posterior_chunked = process_posterior_chunked,
  save_warmstart = save_warmstart,
  get_warmstart = get_warmstart,
  adaptive_sampling = adaptive_sampling,
  recommend_chains = recommend_chains,
  run_fit_stan_hmc = run_fit_stan_hmc,
  run_variational_inference = run_variational_inference,
  clean_cache = clean_cache,
  monitor_memory = monitor_memory,
  check_convergence = check_convergence
)
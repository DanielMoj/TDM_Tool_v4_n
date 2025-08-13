# optimization_utils.R - Utility functions for Stan/MCMC optimizations

library(yaml)
library(logger)
library(fs)

#' Load optimization configuration
#' @param config_file Path to configuration file
#' @return List with configuration settings
load_optimization_config <- function(config_file = "config/optimization_config.yaml") {
  if (!file_exists(config_file)) {
    warning(sprintf("Config file not found: %s, using defaults", config_file))
    return(get_default_config())
  }
  
  config <- read_yaml(config_file)
  log_info("Loaded optimization config from: {config_file}")
  return(config)
}

#' Get default configuration
#' @return List with default settings
get_default_config <- function() {
  list(
    memory = list(
      max_memory_mb = 4096,
      chunk_size = 1000,
      gc_threshold_mb = 3500,
      cache_max_size_mb = 500
    ),
    warmstart = list(
      enabled = TRUE,
      cache_dir = "cache/warmstart",
      max_age_hours = 24
    ),
    adaptive_sampling = list(
      enabled = TRUE,
      max_attempts = 3,
      adapt_delta_increment = 0.05
    ),
    chains = list(
      auto_adjust = TRUE,
      min_chains = 2,
      max_chains = 8,
      target_ess = 1000
    )
  )
}

#' Create optimized Stan sampler with config
#' @param model_code Stan model code
#' @param config Configuration list or path to config file
#' @return Function that runs optimized sampling
create_optimized_sampler <- function(model_code, config = NULL) {
  
  # Load config
  if (is.null(config)) {
    config <- load_optimization_config()
  } else if (is.character(config)) {
    config <- load_optimization_config(config)
  }
  
  # Compile model once
  model_file <- write_stan_file(model_code)
  model <- cmdstan_model(model_file)
  
  # Initialize warmstart manager if enabled
  warmstart_manager <- NULL
  if (config$warmstart$enabled) {
    warmstart_manager <- create_warmstart_manager(
      cache_dir = config$warmstart$cache_dir,
      max_cache_size_mb = config$warmstart$cache_max_size_mb
    )
  }
  
  # Return sampler function
  function(data, 
           params = NULL,
           model_type = "default",
           override_config = list()) {
    
    # Merge override config
    sampling_config <- modifyList(config, override_config)
    
    # Get model-specific settings if available
    if (model_type %in% names(sampling_config$models)) {
      model_settings <- sampling_config$models[[model_type]]
    } else {
      model_settings <- list(
        adapt_delta = sampling_config$adaptive_sampling$initial_adapt_delta,
        max_treedepth = 10,
        thin = 1
      )
    }
    
    # Check memory before starting
    current_memory <- monitor_memory(sampling_config$memory$gc_threshold_mb)
    log_info("Current memory usage: {current_memory} MB")
    
    # Load warmstart if available
    warmstart_data <- NULL
    if (!is.null(warmstart_manager)) {
      warmstart_data <- warmstart_manager$load(
        model_code = model_code,
        data = data,
        max_age_hours = sampling_config$warmstart$max_age_hours
      )
    }
    
    # Prepare sampling arguments
    sampling_args <- list(
      data = data,
      chains = sampling_config$chains$min_chains,
      iter_warmup = 1000,
      iter_sampling = 1000,
      adapt_delta = model_settings$adapt_delta,
      max_treedepth = model_settings$max_treedepth
    )
    
    # Apply warmstart if available
    if (!is.null(warmstart_data) && !is.null(warmstart_manager)) {
      sampling_args <- warmstart_manager$apply_warmstart(warmstart_data, sampling_args)
      log_info("Applied warmstart configuration")
    }
    
    # Run adaptive sampling
    if (sampling_config$adaptive_sampling$enabled) {
      fit <- adaptive_sampling(
        model = model,
        data = data,
        init = sampling_args$init %||% "random",
        chains = sampling_args$chains,
        iter_warmup = sampling_args$iter_warmup,
        iter_sampling = sampling_args$iter_sampling,
        adapt_delta = sampling_args$adapt_delta,
        max_treedepth = sampling_args$max_treedepth,
        model_id = digest::digest(list(model_code, data))
      )
    } else {
      fit <- do.call(model$sample, sampling_args)
    }
    
    # Check convergence
    converged <- check_convergence(
      fit,
      rhat_threshold = sampling_config$diagnostics$max_rhat,
      ess_threshold = sampling_config$diagnostics$min_ess_bulk
    )
    
    # Extract draws efficiently
    thin <- if (object.size(fit) > sampling_config$extraction$large_fit_threshold_mb * 1024^2) {
      sampling_config$extraction$thin_large_fits
    } else {
      model_settings$thin
    }
    
    draws <- extract_draws_efficient(fit, params = params, thin = thin)
    
    # Save warmstart for next run
    if (!is.null(warmstart_manager)) {
      warmstart_manager$save(fit, model_code, data)
    }
    
    # Get diagnostics
    diagnostics <- extract_diagnostics(fit, config = sampling_config)
    
    # Clean up
    if (sampling_config$memory$temp_file_cleanup) {
      cleanup_temp_files()
    }
    
    # Monitor final memory
    final_memory <- monitor_memory(sampling_config$memory$gc_threshold_mb)
    log_info("Final memory usage: {final_memory} MB")
    
    # Return results
    list(
      draws = draws,
      diagnostics = diagnostics,
      converged = converged,
      config_used = model_settings,
      memory_stats = list(
        initial = current_memory,
        final = final_memory,
        difference = final_memory - current_memory
      )
    )
  }
}

#' Extract comprehensive diagnostics
#' @param fit Stan fit object
#' @param config Configuration list
#' @return List with diagnostics
extract_diagnostics <- function(fit, config) {
  summary <- fit$summary()
  sampler_diag <- fit$diagnostic_summary()
  
  # Basic diagnostics
  diagnostics <- list(
    max_rhat = max(summary$rhat, na.rm = TRUE),
    min_ess_bulk = min(summary$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(summary$ess_tail, na.rm = TRUE),
    num_divergent = sum(sampler_diag$num_divergent),
    num_max_treedepth = sum(sampler_diag$num_max_treedepth),
    ebfmi = sampler_diag$ebfmi
  )
  
  # Check against thresholds
  diagnostics$passed <- all(
    diagnostics$max_rhat <= config$diagnostics$max_rhat,
    diagnostics$min_ess_bulk >= config$diagnostics$min_ess_bulk,
    diagnostics$min_ess_tail >= config$diagnostics$min_ess_tail,
    diagnostics$num_divergent <= config$diagnostics$max_divergences,
    diagnostics$num_max_treedepth <= config$diagnostics$max_treedepth_hits,
    all(diagnostics$ebfmi >= config$diagnostics$energy_threshold, na.rm = TRUE)
  )
  
  # Add recommendations
  diagnostics$recommendations <- generate_recommendations(diagnostics, config)
  
  return(diagnostics)
}

#' Generate recommendations based on diagnostics
#' @param diagnostics Diagnostic results
#' @param config Configuration list
#' @return Character vector with recommendations
generate_recommendations <- function(diagnostics, config) {
  recommendations <- character()
  
  if (diagnostics$max_rhat > config$diagnostics$max_rhat) {
    recommendations <- c(recommendations, 
                        sprintf("Increase iterations (R-hat = %.3f)", diagnostics$max_rhat))
  }
  
  if (diagnostics$min_ess_bulk < config$diagnostics$min_ess_bulk) {
    recommendations <- c(recommendations,
                        sprintf("Increase chains or iterations (ESS = %.0f)", diagnostics$min_ess_bulk))
  }
  
  if (diagnostics$num_divergent > 0) {
    recommendations <- c(recommendations,
                        sprintf("Increase adapt_delta (%d divergences)", diagnostics$num_divergent))
  }
  
  if (diagnostics$num_max_treedepth > 0) {
    recommendations <- c(recommendations,
                        sprintf("Increase max_treedepth (%d hits)", diagnostics$num_max_treedepth))
  }
  
  if (length(recommendations) == 0) {
    recommendations <- "Sampling successful, no adjustments needed"
  }
  
  return(recommendations)
}

#' Clean up temporary files
#' @param pattern File pattern to match
#' @param max_age_hours Maximum age in hours
cleanup_temp_files <- function(pattern = "stan_model.*", max_age_hours = 24) {
  temp_dir <- tempdir()
  temp_files <- dir_ls(temp_dir, regexp = pattern)
  
  for (file in temp_files) {
    file_age <- as.numeric(difftime(Sys.time(), file_info(file)$modification_time, units = "hours"))
    if (file_age > max_age_hours) {
      file_delete(file)
      log_debug("Deleted temp file: {basename(file)}")
    }
  }
}

#' Setup optimization environment
#' @param config_file Path to configuration file
#' @return Invisible NULL
setup_optimization_env <- function(config_file = "config/optimization_config.yaml") {
  config <- load_optimization_config(config_file)
  
  # Create directories
  dir_create(config$warmstart$cache_dir, recurse = TRUE)
  dir_create(config$extraction$output_dir, recurse = TRUE)
  dir_create(dirname(config$monitoring$log_file), recurse = TRUE)
  
  # Setup logging
  if (config$monitoring$enabled) {
    log_appender(appender_file(config$monitoring$log_file))
    log_threshold(if (config$monitoring$verbose) DEBUG else INFO)
    log_info("Optimization environment initialized")
  }
  
  # Set memory limits
  options(
    stan.max_memory = config$memory$max_memory_mb,
    mc.cores = config$chains$parallel_chains
  )
  
  invisible(NULL)
}

#' Create model-specific optimizer
#' @param model_type Type of model (pk_models, pd_models, pkpd_models)
#' @param model_code Stan model code
#' @return Optimized sampler function
create_model_optimizer <- function(model_type, model_code) {
  config <- load_optimization_config()
  
  if (!model_type %in% names(config$models)) {
    stop(sprintf("Unknown model type: %s", model_type))
  }
  
  # Create specialized sampler
  sampler <- create_optimized_sampler(model_code, config)
  
  # Return wrapped function with model type
  function(data, params = NULL, ...) {
    sampler(data, params = params, model_type = model_type, ...)
  }
}

#' Batch process multiple datasets
#' @param datasets List of datasets
#' @param model_code Stan model code
#' @param config Configuration list
#' @return List of results
batch_process <- function(datasets, model_code, config = NULL) {
  sampler <- create_optimized_sampler(model_code, config)
  
  results <- list()
  total <- length(datasets)
  
  for (i in seq_along(datasets)) {
    log_info("Processing dataset {i}/{total}")
    
    tryCatch({
      results[[i]] <- sampler(datasets[[i]])
      
      # Clean up between runs
      gc()
      
      # Save intermediate results
      if (i %% 5 == 0) {
        saveRDS(results, "cache/batch_results_temp.rds")
      }
      
    }, error = function(e) {
      log_error("Failed on dataset {i}: {e$message}")
      results[[i]] <- list(error = e$message)
    })
  }
  
  return(results)
}

#' Export results with compression
#' @param results Results object
#' @param filename Output filename
#' @param config Configuration list
export_results <- function(results, filename, config = NULL) {
  if (is.null(config)) {
    config <- load_optimization_config()
  }
  
  export_config <- config$export
  
  # Prepare export object
  export_obj <- list(
    results = results,
    timestamp = Sys.time(),
    config = if (export_config$include_metadata) config else NULL,
    diagnostics = if (export_config$include_diagnostics) results$diagnostics else NULL
  )
  
  # Save with compression
  if (export_config$compress) {
    saveRDS(export_obj, filename, compress = "xz")
    log_info("Results exported with compression: {filename}")
  } else {
    saveRDS(export_obj, filename, compress = FALSE)
    log_info("Results exported: {filename}")
  }
}

# Initialize environment on load
if (interactive()) {
  setup_optimization_env()
}
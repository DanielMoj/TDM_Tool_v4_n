# warmstart_manager.R - Advanced Warmstart Management System
# Handles initialization caching and parameter persistence

library(digest)
library(jsonlite)
library(fs)

#' WarmstartManager R6 Class
#' @description Manages warmstart caching with advanced features
#' @export
WarmstartManager <- R6::R6Class(
  "WarmstartManager",
  
  private = list(
    cache_dir = NULL,
    max_cache_size_mb = 500,
    cache_index = list(),
    
    #' Generate cache key from model and data
    get_cache_key = function(model_code, data) {
      # Create unique hash from model code and data structure
      model_hash <- digest(model_code, algo = "md5")
      data_hash <- digest(names(data), algo = "md5")
      paste0(model_hash, "_", data_hash)
    },
    
    #' Check cache size and clean if needed
    manage_cache_size = function() {
      cache_files <- dir_ls(private$cache_dir, regexp = "\\.rds$")
      
      if (length(cache_files) == 0) return()
      
      # Get file sizes
      file_info <- file_info(cache_files)
      total_size_mb <- sum(file_info$size) / 1024^2
      
      if (total_size_mb > private$max_cache_size_mb) {
        # Remove oldest files
        file_info <- file_info[order(file_info$modification_time), ]
        
        while (total_size_mb > private$max_cache_size_mb * 0.8 && nrow(file_info) > 0) {
          oldest_file <- rownames(file_info)[1]
          file_delete(oldest_file)
          message(sprintf("Removed old cache file: %s", basename(oldest_file)))
          
          file_info <- file_info[-1, ]
          total_size_mb <- sum(file_info$size) / 1024^2
        }
      }
    },
    
    #' Load cache index
    load_index = function() {
      index_file <- file.path(private$cache_dir, "cache_index.json")
      if (file_exists(index_file)) {
        private$cache_index <- fromJSON(index_file)
      }
    },
    
    #' Save cache index
    save_index = function() {
      index_file <- file.path(private$cache_dir, "cache_index.json")
      write_json(private$cache_index, index_file, pretty = TRUE)
    }
  ),
  
  public = list(
    #' Initialize WarmstartManager
    #' @param cache_dir Cache directory path
    #' @param max_cache_size_mb Maximum cache size in MB
    initialize = function(cache_dir = "cache/warmstart", max_cache_size_mb = 500) {
      private$cache_dir <- cache_dir
      private$max_cache_size_mb <- max_cache_size_mb
      
      # Create cache directory
      dir_create(private$cache_dir, recurse = TRUE)
      
      # Load existing index
      private$load_index()
      
      message(sprintf("WarmstartManager initialized with cache at: %s", private$cache_dir))
    },
    
    #' Save warmstart data with metadata
    #' @param fit Stan fit object
    #' @param model_code Stan model code
    #' @param data Stan data
    #' @param metadata Additional metadata
    save = function(fit, model_code, data, metadata = list()) {
      cache_key <- private$get_cache_key(model_code, data)
      cache_file <- file.path(private$cache_dir, paste0(cache_key, ".rds"))
      
      # Extract warmstart data
      warmstart_data <- self$extract_warmstart_data(fit)
      warmstart_data$metadata <- metadata
      warmstart_data$timestamp <- Sys.time()
      warmstart_data$cache_key <- cache_key
      
      # Save to cache
      saveRDS(warmstart_data, cache_file)
      
      # Update index
      private$cache_index[[cache_key]] <- list(
        file = basename(cache_file),
        timestamp = as.character(warmstart_data$timestamp),
        metadata = metadata
      )
      private$save_index()
      
      # Manage cache size
      private$manage_cache_size()
      
      message(sprintf("Warmstart saved with key: %s", cache_key))
      return(cache_key)
    },
    
    #' Load warmstart data
    #' @param model_code Stan model code
    #' @param data Stan data
    #' @param max_age_hours Maximum age in hours
    #' @return Warmstart data or NULL
    load = function(model_code, data, max_age_hours = 24) {
      cache_key <- private$get_cache_key(model_code, data)
      cache_file <- file.path(private$cache_dir, paste0(cache_key, ".rds"))
      
      if (!file_exists(cache_file)) {
        return(NULL)
      }
      
      warmstart_data <- readRDS(cache_file)
      
      # Check age
      age_hours <- as.numeric(difftime(Sys.time(), warmstart_data$timestamp, units = "hours"))
      if (age_hours > max_age_hours) {
        message(sprintf("Warmstart cache expired (%.1f hours old)", age_hours))
        return(NULL)
      }
      
      message(sprintf("Loaded warmstart (%.1f hours old)", age_hours))
      return(warmstart_data)
    },
    
    #' Extract comprehensive warmstart data from fit
    #' @param fit Stan fit object
    #' @return List with warmstart components
    extract_warmstart_data = function(fit) {
      # Get draws
      draws <- fit$draws(format = "draws_array", inc_warmup = FALSE)
      n_chains <- dim(draws)[2]
      n_iter <- dim(draws)[1]
      
      # Extract last iteration for each chain
      inits <- list()
      for (chain in seq_len(n_chains)) {
        inits[[chain]] <- as.list(draws[n_iter, chain, ])
      }
      
      # Get sampler parameters
      sampler_params <- fit$sampler_diagnostics()
      
      # Extract adaptation parameters
      step_sizes <- sapply(seq_len(n_chains), function(chain) {
        steps <- sampler_params[, chain, "stepsize__"]
        mean(steps[length(steps) - 100:length(steps)], na.rm = TRUE)  # Last 100 iterations
      })
      
      # Extract metric (mass matrix) if available
      metric <- tryCatch({
        fit$inv_metric()
      }, error = function(e) NULL)
      
      # Get convergence statistics
      summary <- fit$summary()
      
      warmstart_data <- list(
        inits = inits,
        step_sizes = step_sizes,
        metric = metric,
        n_chains = n_chains,
        n_iterations = n_iter,
        convergence = list(
          max_rhat = max(summary$rhat, na.rm = TRUE),
          min_ess = min(summary$ess_bulk, na.rm = TRUE),
          n_divergent = sum(fit$diagnostic_summary()$num_divergent)
        )
      )
      
      return(warmstart_data)
    },
    
    #' Apply warmstart to sampling arguments
    #' @param warmstart_data Warmstart data from load()
    #' @param sampling_args List of sampling arguments
    #' @return Modified sampling arguments
    apply_warmstart = function(warmstart_data, sampling_args = list()) {
      if (is.null(warmstart_data)) {
        return(sampling_args)
      }
      
      # Set initial values
      if (is.null(sampling_args$init)) {
        n_chains <- sampling_args$chains %||% 4
        sampling_args$init <- warmstart_data$inits[seq_len(n_chains)]
      }
      
      # Set step sizes if available
      if (!is.null(warmstart_data$step_sizes)) {
        sampling_args$step_size <- warmstart_data$step_sizes
      }
      
      # Set metric if available
      if (!is.null(warmstart_data$metric)) {
        sampling_args$inv_metric <- warmstart_data$metric
      }
      
      # Adjust warmup based on previous convergence
      if (!is.null(warmstart_data$convergence)) {
        if (warmstart_data$convergence$max_rhat < 1.01 && 
            warmstart_data$convergence$min_ess > 1000) {
          # Good convergence, reduce warmup
          sampling_args$iter_warmup <- sampling_args$iter_warmup %||% 1000
          sampling_args$iter_warmup <- round(sampling_args$iter_warmup * 0.5)
          message(sprintf("Reducing warmup to %d iterations based on previous convergence", 
                         sampling_args$iter_warmup))
        }
      }
      
      return(sampling_args)
    },
    
    #' Get cache statistics
    #' @return List with cache statistics
    get_stats = function() {
      cache_files <- dir_ls(private$cache_dir, regexp = "\\.rds$")
      
      if (length(cache_files) == 0) {
        return(list(
          n_files = 0,
          total_size_mb = 0,
          oldest_file = NULL,
          newest_file = NULL
        ))
      }
      
      file_info <- file_info(cache_files)
      
      list(
        n_files = length(cache_files),
        total_size_mb = sum(file_info$size) / 1024^2,
        oldest_file = basename(cache_files[which.min(file_info$modification_time)]),
        newest_file = basename(cache_files[which.max(file_info$modification_time)]),
        cache_index = names(private$cache_index)
      )
    },
    
    #' Clear cache
    #' @param older_than_days Remove files older than this many days
    clear = function(older_than_days = NULL) {
      cache_files <- dir_ls(private$cache_dir, regexp = "\\.rds$")
      
      if (is.null(older_than_days)) {
        # Remove all cache files
        file_delete(cache_files)
        private$cache_index <- list()
        private$save_index()
        message("Cache cleared completely")
      } else {
        # Remove old files
        file_info <- file_info(cache_files)
        cutoff_time <- Sys.time() - older_than_days * 24 * 3600
        old_files <- cache_files[file_info$modification_time < cutoff_time]
        
        if (length(old_files) > 0) {
          file_delete(old_files)
          
          # Update index
          for (file in old_files) {
            key <- gsub("\\.rds$", "", basename(file))
            private$cache_index[[key]] <- NULL
          }
          private$save_index()
          
          message(sprintf("Removed %d old cache files", length(old_files)))
        }
      }
    },
    
    #' Export cache to archive
    #' @param output_file Archive file path
    export_cache = function(output_file = "warmstart_cache.tar.gz") {
      temp_dir <- tempdir()
      export_dir <- file.path(temp_dir, "warmstart_export")
      dir_create(export_dir)
      
      # Copy cache files
      cache_files <- dir_ls(private$cache_dir)
      file_copy(cache_files, export_dir)
      
      # Create archive
      tar(output_file, files = export_dir, compression = "gzip")
      
      # Clean up
      dir_delete(export_dir)
      
      message(sprintf("Cache exported to: %s", output_file))
    },
    
    #' Import cache from archive
    #' @param archive_file Archive file path
    import_cache = function(archive_file) {
      if (!file_exists(archive_file)) {
        stop("Archive file does not exist")
      }
      
      # Extract to cache directory
      untar(archive_file, exdir = dirname(private$cache_dir))
      
      # Reload index
      private$load_index()
      
      message(sprintf("Cache imported from: %s", archive_file))
    }
  )
)

#' Create global warmstart manager instance
#' @export
create_warmstart_manager <- function(cache_dir = "cache/warmstart", 
                                    max_cache_size_mb = 500) {
  WarmstartManager$new(cache_dir = cache_dir, max_cache_size_mb = max_cache_size_mb)
}

# Convenience functions for backward compatibility

#' Simple warmstart save function
#' @export
save_warmstart_simple <- function(fit, model_id) {
  manager <- create_warmstart_manager()
  manager$save(fit, model_code = model_id, data = list(), 
               metadata = list(model_id = model_id))
}

#' Simple warmstart load function
#' @export
load_warmstart_simple <- function(model_id, max_age_hours = 24) {
  manager <- create_warmstart_manager()
  manager$load(model_code = model_id, data = list(), 
               max_age_hours = max_age_hours)
}

#' Use warmstart in Stan sampling
#' @export
use_warmstart <- function(model, data, sampling_args = list(), 
                         model_code = NULL, max_age_hours = 24) {
  manager <- create_warmstart_manager()
  
  # Try to load warmstart
  warmstart_data <- manager$load(
    model_code = model_code %||% deparse(substitute(model)),
    data = data,
    max_age_hours = max_age_hours
  )
  
  # Apply warmstart to sampling arguments
  sampling_args <- manager$apply_warmstart(warmstart_data, sampling_args)
  
  # Run sampling
  fit <- do.call(model$sample, c(list(data = data), sampling_args))
  
  # Save new warmstart
  manager$save(fit, model_code = model_code, data = data)
  
  return(fit)
}
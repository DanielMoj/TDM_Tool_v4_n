# R/cache_manager.R
# LRU (Least Recently Used) Cache Manager for Stan Models
# Implements memory-efficient caching with automatic eviction

#' @title LRU Cache for Stan Models
#' @description R6 class implementing an LRU cache with memory monitoring
#' @export
LRUCache <- R6::R6Class(
  "LRUCache",
  
  private = list(
    .cache = NULL,           # Environment to store models
    .access_times = NULL,    # Environment to track access times
    .access_counts = NULL,   # Environment to track access frequency
    .model_sizes = NULL,     # Environment to track memory usage
    .max_size = 10,          # Maximum number of cached models
    .total_memory = 0,       # Total memory usage in bytes
    
    #' Calculate object size in memory
    get_object_size = function(obj) {
      tryCatch({
        as.numeric(object.size(obj))
      }, error = function(e) {
        0
      })
    },
    
    #' Evict least recently used model
    evict_lru = function() {
      if (length(ls(private$.cache)) == 0) return(invisible(NULL))
      
      # Get all keys with their access times
      keys <- ls(private$.cache)
      access_times <- sapply(keys, function(k) {
        if (exists(k, envir = private$.access_times)) {
          get(k, envir = private$.access_times)
        } else {
          0
        }
      })
      
      # Find least recently used
      lru_key <- keys[which.min(access_times)]
      
      # Log eviction
      message(sprintf(
        "Evicting cached model: %s (last accessed: %s, size: %.2f MB)", 
        lru_key,
        format(as.POSIXct(access_times[lru_key], origin = "1970-01-01"), "%Y-%m-%d %H:%M:%S"),
        private$.model_sizes[[lru_key]] / 1024^2
      ))
      
      # Update total memory before eviction
      private$.total_memory <- private$.total_memory - private$.model_sizes[[lru_key]]
      
      # Remove from all environments
      rm(list = lru_key, envir = private$.cache)
      rm(list = lru_key, envir = private$.access_times)
      rm(list = lru_key, envir = private$.access_counts)
      rm(list = lru_key, envir = private$.model_sizes)
      
      # Force garbage collection after eviction
      gc(verbose = FALSE)
      
      invisible(NULL)
    }
  ),
  
  public = list(
    #' @description Initialize LRU Cache
    #' @param max_size Maximum number of models to cache (default: 10)
    initialize = function(max_size = 10) {
      private$.max_size <- max_size
      private$.cache <- new.env(parent = emptyenv())
      private$.access_times <- new.env(parent = emptyenv())
      private$.access_counts <- new.env(parent = emptyenv())
      private$.model_sizes <- new.env(parent = emptyenv())
      private$.total_memory <- 0
      
      message(sprintf("LRU Cache initialized with max size: %d models", max_size))
    },
    
    #' @description Get or compile a Stan model
    #' @param stan_file Path to Stan file
    #' @param force_recompile Force recompilation even if cached
    #' @return Compiled Stan model
    get = function(stan_file, force_recompile = FALSE) {
      if (!file.exists(stan_file)) {
        stop("Stan file not found: ", stan_file)
      }
      
      # Create cache key from file name and modification time
      key <- digest::digest(c(stan_file, file.mtime(stan_file)))
      
      # Check if we should force recompile
      if (force_recompile && exists(key, envir = private$.cache)) {
        message("Force recompiling model: ", basename(stan_file))
        rm(list = key, envir = private$.cache)
        if (exists(key, envir = private$.access_times)) {
          rm(list = key, envir = private$.access_times)
        }
        if (exists(key, envir = private$.access_counts)) {
          rm(list = key, envir = private$.access_counts)
        }
        if (exists(key, envir = private$.model_sizes)) {
          private$.total_memory <- private$.total_memory - private$.model_sizes[[key]]
          rm(list = key, envir = private$.model_sizes)
        }
      }
      
      # Check cache
      if (exists(key, envir = private$.cache)) {
        # Update access time and count
        assign(key, as.numeric(Sys.time()), envir = private$.access_times)
        current_count <- get(key, envir = private$.access_counts)
        assign(key, current_count + 1, envir = private$.access_counts)
        
        message(sprintf(
          "Using cached Stan model: %s (access count: %d)", 
          basename(stan_file),
          current_count + 1
        ))
        
        return(get(key, envir = private$.cache))
      }
      
      # Check if we need to evict before adding new model
      current_size <- length(ls(private$.cache))
      if (current_size >= private$.max_size) {
        message(sprintf(
          "Cache full (%d/%d), evicting least recently used model",
          current_size, private$.max_size
        ))
        private$evict_lru()
      }
      
      # Compile new model
      message("Compiling Stan model: ", basename(stan_file))
      
      model <- tryCatch({
        if (requireNamespace("cmdstanr", quietly = TRUE)) {
          cmdstanr::cmdstan_model(stan_file, compile = TRUE)
        } else {
          stop("cmdstanr package not available")
        }
      }, error = function(e) {
        stop(sprintf("Failed to compile Stan model %s: %s", basename(stan_file), e$message))
      })
      
      # Calculate model size
      model_size <- private$get_object_size(model)
      
      # Store in cache with metadata
      assign(key, model, envir = private$.cache)
      assign(key, as.numeric(Sys.time()), envir = private$.access_times)
      assign(key, 1, envir = private$.access_counts)
      assign(key, model_size, envir = private$.model_sizes)
      
      # Update total memory
      private$.total_memory <- private$.total_memory + model_size
      
      message(sprintf(
        "Model cached successfully (size: %.2f MB, total cache: %.2f MB)",
        model_size / 1024^2,
        private$.total_memory / 1024^2
      ))
      
      return(model)
    },
    
    #' @description Clear the entire cache
    clear = function() {
      keys <- ls(private$.cache)
      if (length(keys) > 0) {
        message(sprintf("Clearing %d cached models", length(keys)))
        rm(list = keys, envir = private$.cache)
        rm(list = keys, envir = private$.access_times)
        rm(list = keys, envir = private$.access_counts)
        rm(list = keys, envir = private$.model_sizes)
      }
      private$.total_memory <- 0
      gc(verbose = FALSE)
      message("Cache cleared and memory freed")
    },
    
    #' @description Get cache statistics
    #' @return List with cache statistics
    get_stats = function() {
      keys <- ls(private$.cache)
      n_models <- length(keys)
      
      if (n_models == 0) {
        return(list(
          n_models = 0,
          total_size_mb = 0,
          models = data.frame(),
          max_size = private$.max_size
        ))
      }
      
      # Collect model information
      model_info <- data.frame(
        key = keys,
        size_mb = sapply(keys, function(k) private$.model_sizes[[k]] / 1024^2),
        access_count = sapply(keys, function(k) get(k, envir = private$.access_counts)),
        last_access = sapply(keys, function(k) {
          as.POSIXct(get(k, envir = private$.access_times), origin = "1970-01-01")
        }),
        stringsAsFactors = FALSE
      )
      
      # Sort by last access time (most recent first)
      model_info <- model_info[order(model_info$last_access, decreasing = TRUE), ]
      rownames(model_info) <- NULL
      
      list(
        n_models = n_models,
        total_size_mb = private$.total_memory / 1024^2,
        models = model_info,
        max_size = private$.max_size
      )
    },
    
    #' @description Print cache status
    print = function() {
      stats <- self$get_stats()
      cat("LRU Stan Model Cache\n")
      cat(sprintf("  Models: %d / %d\n", stats$n_models, stats$max_size))
      cat(sprintf("  Total size: %.2f MB\n", stats$total_size_mb))
      
      if (stats$n_models > 0) {
        cat("\nCached models:\n")
        for (i in 1:nrow(stats$models)) {
          cat(sprintf("  %d. Size: %.2f MB, Accessed: %d times, Last: %s\n",
                     i,
                     stats$models$size_mb[i],
                     stats$models$access_count[i],
                     format(stats$models$last_access[i], "%Y-%m-%d %H:%M:%S")))
        }
      }
      invisible(self)
    }
  )
)

# Global cache instance
.stan_lru_cache <- NULL

#' @title Get or create global LRU cache instance
#' @param max_size Maximum number of models to cache
#' @return LRUCache instance
#' @export
get_lru_cache <- function(max_size = 10) {
  if (is.null(.stan_lru_cache)) {
    .stan_lru_cache <<- LRUCache$new(max_size = max_size)
  }
  .stan_lru_cache
}

#' @title Memory monitoring function
#' @description Get comprehensive memory statistics
#' @return List with memory information
#' @export
get_memory_stats <- function() {
  # Get R memory usage
  mem_info <- gc(verbose = FALSE)
  
  # Get cache statistics
  cache <- get_lru_cache()
  cache_stats <- cache$get_stats()
  
  # System memory (if available)
  system_mem <- tryCatch({
    if (Sys.info()["sysname"] == "Linux") {
      mem_output <- system("free -m", intern = TRUE)
      total_line <- grep("^Mem:", mem_output, value = TRUE)
      if (length(total_line) > 0) {
        values <- as.numeric(strsplit(gsub("\\s+", " ", total_line), " ")[[1]][2:4])
        list(
          total_mb = values[1],
          used_mb = values[2],
          free_mb = values[3]
        )
      } else NULL
    } else if (Sys.info()["sysname"] == "Darwin") {
      # macOS
      vm_output <- system("vm_stat", intern = TRUE)
      page_size <- as.numeric(gsub("[^0-9]", "", grep("page size", vm_output, value = TRUE)[1]))
      free_pages <- as.numeric(gsub("[^0-9]", "", grep("Pages free", vm_output, value = TRUE)[1]))
      active_pages <- as.numeric(gsub("[^0-9]", "", grep("Pages active", vm_output, value = TRUE)[1]))
      
      list(
        free_mb = (free_pages * page_size) / 1024^2,
        active_mb = (active_pages * page_size) / 1024^2
      )
    } else {
      NULL
    }
  }, error = function(e) NULL)
  
  # Compile results
  stats <- list(
    timestamp = Sys.time(),
    r_memory = list(
      used_mb = sum(mem_info[, "used"]) / 1024,  # Convert from Kb to Mb
      gc_trigger_mb = sum(mem_info[, "max used"]) / 1024
    ),
    cache = list(
      n_models = cache_stats$n_models,
      max_models = cache_stats$max_size,
      total_size_mb = cache_stats$total_size_mb,
      models = cache_stats$models
    ),
    system = system_mem
  )
  
  class(stats) <- c("memory_stats", "list")
  stats
}

#' @title Print memory statistics
#' @param x Memory stats object
#' @param ... Additional arguments (ignored)
#' @export
print.memory_stats <- function(x, ...) {
  cat("=== Memory Statistics ===\n")
  cat("Timestamp:", format(x$timestamp, "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  cat("R Memory Usage:\n")
  cat(sprintf("  Used: %.2f MB\n", x$r_memory$used_mb))
  cat(sprintf("  GC Trigger: %.2f MB\n", x$r_memory$gc_trigger_mb))
  cat("\n")
  
  cat("Stan Model Cache:\n")
  cat(sprintf("  Models: %d / %d\n", x$cache$n_models, x$cache$max_models))
  cat(sprintf("  Total Size: %.2f MB\n", x$cache$total_size_mb))
  
  if (x$cache$n_models > 0) {
    cat("\n  Top 3 models by size:\n")
    top_models <- head(x$cache$models[order(x$cache$models$size_mb, decreasing = TRUE), ], 3)
    for (i in 1:nrow(top_models)) {
      cat(sprintf("    - %.2f MB (accessed %d times)\n", 
                 top_models$size_mb[i], 
                 top_models$access_count[i]))
    }
  }
  
  if (!is.null(x$system)) {
    cat("\nSystem Memory:\n")
    if (!is.null(x$system$total_mb)) {
      cat(sprintf("  Total: %.0f MB\n", x$system$total_mb))
      cat(sprintf("  Used: %.0f MB\n", x$system$used_mb))
      cat(sprintf("  Free: %.0f MB\n", x$system$free_mb))
    } else if (!is.null(x$system$free_mb)) {
      cat(sprintf("  Free: %.0f MB\n", x$system$free_mb))
      cat(sprintf("  Active: %.0f MB\n", x$system$active_mb))
    }
  }
  
  invisible(x)
}

#' @title Force garbage collection with logging
#' @param verbose Print collection results
#' @export
force_gc <- function(verbose = TRUE) {
  before <- gc(verbose = FALSE)
  gc_result <- gc(verbose = FALSE, full = TRUE)
  after <- gc(verbose = FALSE)
  
  if (verbose) {
    freed_mb <- (sum(before[, "used"]) - sum(after[, "used"])) / 1024
    message(sprintf("Garbage collection completed. Freed: %.2f MB", max(0, freed_mb)))
  }
  
  invisible(gc_result)
}
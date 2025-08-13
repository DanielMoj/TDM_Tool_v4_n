# R/load_all.R
# Central loading function for all dependencies and source files
# Modernized to use dependency management instead of library() calls

#' Load all required source files and dependencies
#' @param check_dependencies Whether to check package dependencies
#' @param install_missing Whether to install missing packages
#' @return invisible TRUE
#' @export
load_all_sources <- function(check_dependencies = TRUE, install_missing = FALSE) {
  
  # Load dependency management system
  dep_mgr_path <- file.path("R", "dependencies.R")
  if (file.exists(dep_mgr_path)) {
    source(dep_mgr_path, local = FALSE)
    message("✓ Loaded dependency management system")
  } else {
    stop("CRITICAL: dependencies.R not found - cannot manage package dependencies")
  }
  
  # Check and load dependencies
  if (check_dependencies) {
    tryCatch({
      load_dependencies(install_missing = install_missing)
      message("✓ All package dependencies satisfied")
    }, error = function(e) {
      stop(sprintf("Dependency check failed: %s", e$message))
    })
  }
  
  # CRITICAL: Load utils.R first - contains %||% operator and other utilities
  utils_path <- file.path("R", "utils.R")
  if (file.exists(utils_path)) {
    source(utils_path, local = FALSE)
    message("✓ Loaded core utilities: utils.R")
  } else {
    stop("CRITICAL: utils.R not found - this file contains essential utilities including %||% operator")
  }
  
  # Define source files in loading order
  # Dependencies first, then modules
  source_files <- c(
    # Core utilities (utils.R already loaded above)
    "utils/config",
    "utils/helpers",
    "utils/validators",
    
    # Safe I/O operations (requires utils.R for %||%)
    "safe_io",
    
    # Database and authentication
    "core/db_connection",
    "core/auth_functions",
    "core/audit_logger",
    
    # Error monitoring (now uses :: notation)
    "error_monitor",
    
    # Health checks
    "health_checks",
    
    # FHIR modules (require utils.R for %||%)
    "fhir_auth",
    "fhir_cache",
    "fhir_circuit",
    "fhir",
    
    # PK/PD specific
    "models/pk_models",
    "models/pk_parameters",
    "models/dosing_regimens",
    "pk_calculations",
    
    # Analysis functions
    "run_fit_jags",  # Main fitting function
    "diagnostics",
    "optimization",
    "sensitivity",
    
    # Data processing
    "data/data_import",
    "data/data_validation",
    "data/data_transformation",
    
    # Plotting functions
    "plots/plot_concentrations",
    "plots/plot_diagnostics",
    "plots/plot_parameters",
    "plotting_functions",
    
    # Report generation
    "reports/report_generator",
    "reports/export_functions"
  )
  
  # Track loading results
  loading_results <- list(
    loaded = character(),
    skipped = character(),
    failed = character()
  )
  
  # Source all files with error handling
  for (file in source_files) {
    file_path <- file.path("R", paste0(file, ".R"))
    if (file.exists(file_path)) {
      tryCatch({
        source(file_path, local = FALSE)
        loading_results$loaded <- c(loading_results$loaded, file_path)
        message(paste("✓ Loaded:", file_path))
      }, error = function(e) {
        loading_results$failed <- c(loading_results$failed, file_path)
        warning(paste("✗ Failed to load:", file_path, "-", e$message))
      })
    } else {
      loading_results$skipped <- c(loading_results$skipped, file_path)
      # Only show info for optional files
      if (!grepl("^(utils/|core/|plots/|reports/|data/)", file)) {
        message(paste("ℹ File not found (skipping):", file_path))
      }
    }
  }
  
  # Load all modules
  module_dir <- file.path("R", "modules")
  if (dir.exists(module_dir)) {
    module_files <- list.files(
      path = module_dir,
      pattern = "^mod_.*\\.R$",
      full.names = TRUE
    )
    
    for (module_file in module_files) {
      tryCatch({
        source(module_file, local = FALSE)
        loading_results$loaded <- c(loading_results$loaded, module_file)
        message(paste("✓ Loaded module:", basename(module_file)))
      }, error = function(e) {
        loading_results$failed <- c(loading_results$failed, module_file)
        warning(paste("✗ Failed to load module:", module_file, "-", e$message))
      })
    }
  }
  
  # Set global options
  options(
    shiny.maxRequestSize = 50*1024^2,  # 50MB upload limit
    shiny.sanitize.errors = TRUE,
    digits = 4,
    scipen = 999
  )
  
  # Initialize database pool if config exists
  if (file.exists("config/database.yml")) {
    message("Checking database configuration...")
    # Check if DBI is available
    if (requireNamespace("DBI", quietly = TRUE) && 
        requireNamespace("pool", quietly = TRUE)) {
      message("  Database packages available")
      # initialize_db_pool() would be defined in db_connection.R
    } else {
      warning("  Database packages not available - database features disabled")
    }
  }
  
  # Verify critical functions are available
  if (!exists("%||%")) {
    stop("CRITICAL: %||% operator not loaded - check utils.R")
  }
  
  # Summary
  message("\n========================================")
  message("Source Loading Summary:")
  message(sprintf("  ✓ Loaded: %d files", length(loading_results$loaded)))
  if (length(loading_results$skipped) > 0) {
    message(sprintf("  ℹ Skipped: %d files (not found)", length(loading_results$skipped)))
  }
  if (length(loading_results$failed) > 0) {
    message(sprintf("  ✗ Failed: %d files", length(loading_results$failed)))
    warning("Some files failed to load. Check warnings above for details.")
  }
  message("========================================")
  message("✓ Application sources loaded successfully")
  message("✓ %||% operator is available globally")
  
  invisible(TRUE)
}

#' Clean up resources on app shutdown
#' @export
cleanup_resources <- function() {
  # Close database connections if DBI is available
  if (requireNamespace("DBI", quietly = TRUE) && 
      exists("db_pool") && !is.null(db_pool)) {
    if (requireNamespace("pool", quietly = TRUE)) {
      pool::poolClose(db_pool)
      message("Database pool closed")
    }
  }
  
  # Clean up temporary files
  temp_files <- list.files(tempdir(), pattern = "^pk_", full.names = TRUE)
  if (length(temp_files) > 0) {
    unlink(temp_files)
    message(paste("Cleaned up", length(temp_files), "temporary files"))
  }
  
  invisible(TRUE)
}

#' Check system health before starting
#' @export
system_health_check <- function() {
  
  checks <- list()
  
  # Check R version
  r_version <- R.version
  checks$r_version <- list(
    status = r_version$major >= "4",
    message = sprintf("R version: %s.%s", r_version$major, r_version$minor)
  )
  
  # Check critical packages
  critical_packages <- c("shiny", "ggplot2", "dplyr")
  for (pkg in critical_packages) {
    checks[[paste0("package_", pkg)]] <- list(
      status = requireNamespace(pkg, quietly = TRUE),
      message = sprintf("Package %s: %s", 
                       pkg, 
                       ifelse(requireNamespace(pkg, quietly = TRUE), "OK", "MISSING"))
    )
  }
  
  # Check critical directories
  critical_dirs <- c("R", "config", "data")
  for (dir in critical_dirs) {
    checks[[paste0("dir_", dir)]] <- list(
      status = dir.exists(dir),
      message = sprintf("Directory %s: %s", 
                       dir, 
                       ifelse(dir.exists(dir), "OK", "MISSING"))
    )
  }
  
  # Check write permissions
  temp_test <- tempfile()
  write_test <- tryCatch({
    writeLines("test", temp_test)
    unlink(temp_test)
    TRUE
  }, error = function(e) FALSE)
  
  checks$write_permission <- list(
    status = write_test,
    message = sprintf("Write permission: %s", ifelse(write_test, "OK", "FAILED"))
  )
  
  # Print results
  message("\n========================================")
  message("System Health Check Results:")
  message("========================================")
  
  all_passed <- TRUE
  for (check_name in names(checks)) {
    check <- checks[[check_name]]
    status_symbol <- ifelse(check$status, "✓", "✗")
    message(sprintf("  %s %s", status_symbol, check$message))
    if (!check$status) all_passed <- FALSE
  }
  
  message("========================================")
  
  if (!all_passed) {
    warning("Some health checks failed. The application may not function correctly.")
  } else {
    message("✓ All health checks passed")
  }
  
  return(all_passed)
}

# Auto-load if running in interactive mode for development
if (interactive()) {
  message("\n========================================")
  message("PK/PD Application - Development Mode")
  message("========================================")
  
  # Run health check first
  if (system_health_check()) {
    # Load sources
    message("\nLoading application sources...")
    load_all_sources(check_dependencies = TRUE, install_missing = FALSE)
  } else {
    warning("System health check failed. Please resolve issues before continuing.")
  }
}
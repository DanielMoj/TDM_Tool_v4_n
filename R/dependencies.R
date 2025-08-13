# R/dependencies.R
# Central Package Dependency Management for TDMx Application
# 
# This file manages all package dependencies for the application.
# It should be sourced at application startup (e.g., in app.R)
# 
# Benefits:
# - Single source of truth for dependencies
# - Clear version requirements
# - Graceful error handling
# - No namespace pollution from library() calls

#' Check and load required packages
#'
#' @param install_missing Logical, whether to attempt installing missing packages
#' @return Logical indicating success
#' @export
load_dependencies <- function(install_missing = FALSE) {
  
  # Define all required packages with minimum versions
  required_packages <- list(
    # Core Shiny packages
    shiny = "1.7.0",
    shinydashboard = NULL,
    shinyjs = NULL,
    shinyWidgets = NULL,
    bslib = NULL,
    
    # Data visualization
    ggplot2 = "3.4.0",
    plotly = "4.10.0",
    DT = "0.26",
    
    # Data manipulation
    dplyr = "1.1.0",
    tidyr = "1.3.0",
    data.table = NULL,
    
    # Database
    DBI = "1.1.0",
    RPostgres = NULL,
    RSQLite = NULL,
    pool = NULL,
    
    # Bayesian inference
    rjags = NULL,
    coda = NULL,
    
    # Utilities
    jsonlite = "1.8.0",
    yaml = "2.3.0",
    digest = "0.6.30",
    lubridate = "1.9.0",
    glue = NULL,
    readr = NULL,
    tibble = NULL,
    
    # PK/PD specific
    deSolve = NULL,
    numDeriv = NULL,
    MASS = NULL,
    
    # Reporting
    rmarkdown = NULL,
    knitr = NULL,
    
    # Async/parallel
    promises = NULL,
    future = NULL,
    
    # Additional analysis packages
    zoo = NULL,
    posterior = NULL,
    bayesplot = NULL,
    matrixStats = NULL
  )
  
  # Optional packages (warn if missing but don't fail)
  optional_packages <- list(
    cmdstanr = NULL,  # For Stan models
    rstan = NULL,     # Alternative Stan interface
    httr = NULL,      # For API calls
    httr2 = NULL,     # Modern HTTP client
    sendmailR = NULL, # For email alerts
    later = NULL      # For scheduling
  )
  
  # Check required packages
  missing_required <- character()
  version_mismatches <- character()
  
  for (pkg_name in names(required_packages)) {
    min_version <- required_packages[[pkg_name]]
    
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      missing_required <- c(missing_required, pkg_name)
    } else if (!is.null(min_version)) {
      # Check version requirement
      installed_version <- as.character(packageVersion(pkg_name))
      if (package_version(installed_version) < package_version(min_version)) {
        version_mismatches <- c(
          version_mismatches,
          sprintf("%s (installed: %s, required: >= %s)", 
                  pkg_name, installed_version, min_version)
        )
      }
    }
  }
  
  # Check optional packages
  missing_optional <- character()
  
  for (pkg_name in names(optional_packages)) {
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      missing_optional <- c(missing_optional, pkg_name)
    }
  }
  
  # Handle missing required packages
  if (length(missing_required) > 0) {
    if (install_missing) {
      message("Installing missing required packages...")
      for (pkg in missing_required) {
        tryCatch({
          install.packages(pkg, quiet = TRUE)
          message(sprintf("  ✓ Installed %s", pkg))
        }, error = function(e) {
          warning(sprintf("  ✗ Failed to install %s: %s", pkg, e$message))
        })
      }
      
      # Re-check after installation
      still_missing <- character()
      for (pkg in missing_required) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
          still_missing <- c(still_missing, pkg)
        }
      }
      
      if (length(still_missing) > 0) {
        stop(sprintf(
          "Required packages still missing after installation attempt:\n  %s\n\nPlease install manually with:\n  install.packages(c(%s))",
          paste(still_missing, collapse = ", "),
          paste(sprintf('"%s"', still_missing), collapse = ", ")
        ))
      }
    } else {
      stop(sprintf(
        "Required packages missing:\n  %s\n\nInstall with:\n  install.packages(c(%s))",
        paste(missing_required, collapse = ", "),
        paste(sprintf('"%s"', missing_required), collapse = ", ")
      ))
    }
  }
  
  # Handle version mismatches
  if (length(version_mismatches) > 0) {
    warning(sprintf(
      "Package version requirements not met:\n  %s\n\nConsider updating with:\n  update.packages()",
      paste(version_mismatches, collapse = "\n  ")
    ))
  }
  
  # Report optional packages
  if (length(missing_optional) > 0) {
    message(sprintf(
      "Optional packages not installed (some features may be unavailable):\n  %s",
      paste(missing_optional, collapse = ", ")
    ))
  }
  
  # Success message
  message(sprintf(
    "✓ Package dependency check complete: %d required packages available",
    length(required_packages)
  ))
  
  return(TRUE)
}

#' Get package dependency report
#'
#' @return Data frame with package status
#' @export
get_dependency_report <- function() {
  
  all_packages <- c(
    names(formals(load_dependencies)$required_packages),
    names(formals(load_dependencies)$optional_packages)
  )
  
  report <- data.frame(
    package = character(),
    installed = logical(),
    version = character(),
    required = logical(),
    status = character(),
    stringsAsFactors = FALSE
  )
  
  for (pkg in all_packages) {
    installed <- requireNamespace(pkg, quietly = TRUE)
    version <- if (installed) as.character(packageVersion(pkg)) else NA
    required <- pkg %in% names(formals(load_dependencies)$required_packages)
    status <- if (!installed) {
      if (required) "MISSING (REQUIRED)" else "MISSING (OPTIONAL)"
    } else {
      "OK"
    }
    
    report <- rbind(report, data.frame(
      package = pkg,
      installed = installed,
      version = version,
      required = required,
      status = status,
      stringsAsFactors = FALSE
    ))
  }
  
  report <- report[order(report$required, report$package, decreasing = c(TRUE, FALSE)), ]
  rownames(report) <- NULL
  
  return(report)
}

#' Check if specific package is available
#'
#' @param package_name Name of package to check
#' @param min_version Minimum required version (optional)
#' @return Logical indicating if package meets requirements
#' @export
check_package <- function(package_name, min_version = NULL) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    return(FALSE)
  }
  
  if (!is.null(min_version)) {
    installed_version <- packageVersion(package_name)
    return(installed_version >= package_version(min_version))
  }
  
  return(TRUE)
}

#' Conditional package loading helper
#'
#' This function provides a safe way to use package functions
#' without library() calls
#'
#' @param package_name Name of package
#' @param error_if_missing Whether to error or return NULL if missing
#' @return Package namespace or NULL
#' @export
use_package <- function(package_name, error_if_missing = TRUE) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    if (error_if_missing) {
      stop(sprintf(
        "Package '%s' is required but not installed. Install with: install.packages('%s')",
        package_name, package_name
      ))
    } else {
      return(NULL)
    }
  }
  
  # Return the package namespace for use with ::
  return(getNamespace(package_name))
}

#' Create package manifest for reproducibility
#'
#' @param file Path to save manifest (NULL for return only)
#' @return List with package information
#' @export
create_package_manifest <- function(file = NULL) {
  
  # Get all loaded namespaces
  loaded_packages <- loadedNamespaces()
  
  manifest <- list(
    created_at = Sys.time(),
    r_version = R.version.string,
    platform = R.version$platform,
    packages = list()
  )
  
  for (pkg in loaded_packages) {
    # Skip base packages
    if (pkg %in% c("base", "methods", "utils", "grDevices", 
                   "graphics", "stats", "datasets")) {
      next
    }
    
    pkg_info <- list(
      version = as.character(packageVersion(pkg)),
      loaded = pkg %in% search()
    )
    
    # Try to get package source
    desc <- tryCatch(
      packageDescription(pkg),
      error = function(e) NULL
    )
    
    if (!is.null(desc)) {
      pkg_info$source <- desc$Repository %||% desc$RemoteType %||% "unknown"
      pkg_info$built <- desc$Built %||% NA
    }
    
    manifest$packages[[pkg]] <- pkg_info
  }
  
  # Save if file specified
  if (!is.null(file)) {
    jsonlite::write_json(manifest, file, pretty = TRUE, auto_unbox = TRUE)
    message(sprintf("Package manifest saved to: %s", file))
  }
  
  return(manifest)
}

#' Install packages from manifest
#'
#' @param manifest_file Path to manifest file
#' @param upgrade Whether to upgrade existing packages
#' @export
install_from_manifest <- function(manifest_file, upgrade = FALSE) {
  
  if (!file.exists(manifest_file)) {
    stop("Manifest file not found: ", manifest_file)
  }
  
  manifest <- jsonlite::read_json(manifest_file)
  
  message(sprintf(
    "Installing packages from manifest created on %s",
    manifest$created_at
  ))
  
  packages_to_install <- character()
  packages_to_upgrade <- character()
  
  for (pkg_name in names(manifest$packages)) {
    pkg_info <- manifest$packages[[pkg_name]]
    required_version <- package_version(pkg_info$version)
    
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      packages_to_install <- c(packages_to_install, pkg_name)
    } else if (upgrade) {
      current_version <- packageVersion(pkg_name)
      if (current_version < required_version) {
        packages_to_upgrade <- c(packages_to_upgrade, pkg_name)
      }
    }
  }
  
  # Install missing packages
  if (length(packages_to_install) > 0) {
    message(sprintf("Installing %d missing packages...", length(packages_to_install)))
    install.packages(packages_to_install)
  }
  
  # Upgrade packages if requested
  if (upgrade && length(packages_to_upgrade) > 0) {
    message(sprintf("Upgrading %d packages...", length(packages_to_upgrade)))
    install.packages(packages_to_upgrade)
  }
  
  message("Package installation from manifest complete")
}

# Package groups for specific functionality
.package_groups <- list(
  core = c("shiny", "shinydashboard", "DT", "ggplot2", "dplyr"),
  database = c("DBI", "RPostgres", "RSQLite", "pool"),
  bayesian = c("rjags", "coda"),
  stan = c("cmdstanr", "rstan", "posterior", "bayesplot"),
  pk_pd = c("deSolve", "numDeriv", "MASS"),
  reporting = c("rmarkdown", "knitr"),
  api = c("httr", "httr2", "jsonlite")
)

#' Check if package group is available
#'
#' @param group Name of package group
#' @return Logical vector indicating package availability
#' @export
check_package_group <- function(group) {
  if (!group %in% names(.package_groups)) {
    stop(sprintf("Unknown package group: %s. Available groups: %s",
                 group, paste(names(.package_groups), collapse = ", ")))
  }
  
  packages <- .package_groups[[group]]
  result <- sapply(packages, function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
  })
  
  if (!all(result)) {
    missing <- packages[!result]
    message(sprintf("Missing packages in group '%s': %s",
                    group, paste(missing, collapse = ", ")))
  }
  
  return(result)
}
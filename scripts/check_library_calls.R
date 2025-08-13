#!/usr/bin/env Rscript
# scripts/check_library_calls.R
# Script to identify remaining library() and require() calls in R files
# Usage: Rscript scripts/check_library_calls.R

cat("=====================================\n")
cat("Package Management Migration Checker\n")
cat("=====================================\n\n")

# Function to check for library/require calls
check_library_calls <- function(dir = "R", exclude_files = c("dependencies.R")) {
  
  # Find all R files
  r_files <- list.files(
    path = dir,
    pattern = "\\.R$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  # Exclude specified files
  for (exclude in exclude_files) {
    r_files <- r_files[!grepl(exclude, r_files)]
  }
  
  cat(sprintf("Checking %d R files in %s/\n\n", length(r_files), dir))
  
  # Results storage
  results <- list(
    library_calls = list(),
    require_calls = list(),
    clean_files = character(),
    summary = list()
  )
  
  # Check each file
  for (file in r_files) {
    content <- readLines(file, warn = FALSE)
    
    # Find library() calls
    library_lines <- grep("^\\s*library\\s*\\(", content)
    if (length(library_lines) > 0) {
      results$library_calls[[file]] <- list(
        lines = library_lines,
        calls = content[library_lines]
      )
    }
    
    # Find require() calls
    require_lines <- grep("^\\s*require\\s*\\(", content)
    if (length(require_lines) > 0) {
      results$require_calls[[file]] <- list(
        lines = require_lines,
        calls = content[require_lines]
      )
    }
    
    # Track clean files
    if (length(library_lines) == 0 && length(require_lines) == 0) {
      results$clean_files <- c(results$clean_files, file)
    }
  }
  
  # Generate summary
  results$summary <- list(
    total_files = length(r_files),
    files_with_library = length(results$library_calls),
    files_with_require = length(results$require_calls),
    clean_files = length(results$clean_files),
    total_library_calls = sum(sapply(results$library_calls, function(x) length(x$lines))),
    total_require_calls = sum(sapply(results$require_calls, function(x) length(x$lines)))
  )
  
  return(results)
}

# Function to print results
print_results <- function(results) {
  
  # Print summary
  cat("SUMMARY\n")
  cat("-------\n")
  cat(sprintf("Total files checked: %d\n", results$summary$total_files))
  cat(sprintf("Files with library(): %d\n", results$summary$files_with_library))
  cat(sprintf("Files with require(): %d\n", results$summary$files_with_require))
  cat(sprintf("Clean files: %d\n", results$summary$clean_files))
  cat(sprintf("Total library() calls: %d\n", results$summary$total_library_calls))
  cat(sprintf("Total require() calls: %d\n", results$summary$total_require_calls))
  cat("\n")
  
  # Print library() calls
  if (length(results$library_calls) > 0) {
    cat("FILES WITH library() CALLS\n")
    cat("--------------------------\n")
    for (file in names(results$library_calls)) {
      cat(sprintf("\n📁 %s\n", file))
      info <- results$library_calls[[file]]
      for (i in seq_along(info$lines)) {
        cat(sprintf("  Line %d: %s\n", info$lines[i], info$calls[i]))
      }
    }
    cat("\n")
  } else {
    cat("✅ No library() calls found!\n\n")
  }
  
  # Print require() calls
  if (length(results$require_calls) > 0) {
    cat("FILES WITH require() CALLS\n")
    cat("--------------------------\n")
    for (file in names(results$require_calls)) {
      cat(sprintf("\n📁 %s\n", file))
      info <- results$require_calls[[file]]
      for (i in seq_along(info$lines)) {
        cat(sprintf("  Line %d: %s\n", info$lines[i], info$calls[i]))
      }
    }
    cat("\n")
  } else {
    cat("✅ No require() calls found!\n\n")
  }
  
  # Migration status
  if (results$summary$files_with_library == 0 && results$summary$files_with_require == 0) {
    cat("🎉 MIGRATION COMPLETE!\n")
    cat("All files are using namespace notation (::)\n")
  } else {
    cat("⚠️  MIGRATION NEEDED\n")
    cat(sprintf("%d files still need migration\n", 
                results$summary$files_with_library + results$summary$files_with_require))
  }
}

# Function to generate migration script
generate_migration_script <- function(results, output_file = "migrate_packages.R") {
  
  if (results$summary$files_with_library == 0 && results$summary$files_with_require == 0) {
    cat("\nNo migration script needed - all files are clean!\n")
    return(invisible(NULL))
  }
  
  cat(sprintf("\nGenerating migration script: %s\n", output_file))
  
  script_lines <- c(
    "# Auto-generated migration script",
    sprintf("# Generated: %s", Sys.Date()),
    "# This script shows the changes needed for each file",
    "",
    "# Files to migrate:"
  )
  
  # Add files with library() calls
  for (file in names(results$library_calls)) {
    script_lines <- c(script_lines,
                     "",
                     sprintf("# File: %s", file),
                     "# Replace these library() calls with :: notation:")
    
    info <- results$library_calls[[file]]
    for (call in info$calls) {
      # Extract package name
      pkg_match <- regmatches(call, regexpr("library\\s*\\(\\s*['\"]?([^)'\"]+)", call))
      if (length(pkg_match) > 0) {
        pkg_name <- gsub("library\\s*\\(\\s*['\"]?", "", pkg_match)
        script_lines <- c(script_lines,
                         sprintf("#   %s", call),
                         sprintf("#   -> Use %s:: for functions from this package", pkg_name))
      }
    }
  }
  
  # Add files with require() calls
  for (file in names(results$require_calls)) {
    script_lines <- c(script_lines,
                     "",
                     sprintf("# File: %s", file),
                     "# Replace these require() calls with requireNamespace():")
    
    info <- results$require_calls[[file]]
    for (call in info$calls) {
      # Extract package name
      pkg_match <- regmatches(call, regexpr("require\\s*\\(\\s*['\"]?([^)'\"]+)", call))
      if (length(pkg_match) > 0) {
        pkg_name <- gsub("require\\s*\\(\\s*['\"]?", "", pkg_match)
        script_lines <- c(script_lines,
                         sprintf("#   %s", call),
                         sprintf('#   -> if (!requireNamespace("%s", quietly = TRUE)) stop("Package %s needed")', 
                                pkg_name, pkg_name))
      }
    }
  }
  
  # Add migration template
  script_lines <- c(script_lines,
                   "",
                   "# MIGRATION TEMPLATE:",
                   "# For each file:",
                   "# 1. Remove library() calls",
                   "# 2. Add dependency check at the beginning:",
                   "#    if (!requireNamespace('package', quietly = TRUE)) {",
                   "#      stop('Package needed. Install with: install.packages(\"package\")')",
                   "#    }",
                   "# 3. Use package::function() notation throughout the file",
                   "# 4. Test the file: source('path/to/file.R')")
  
  writeLines(script_lines, output_file)
  cat(sprintf("Migration script written to: %s\n", output_file))
}

# Main execution
main <- function() {
  # Check R directory
  results <- check_library_calls("R", exclude_files = c("dependencies.R"))
  
  # Print results
  print_results(results)
  
  # Generate migration script if needed
  if (results$summary$files_with_library > 0 || results$summary$files_with_require > 0) {
    generate_migration_script(results)
    
    cat("\nNEXT STEPS:\n")
    cat("-----------\n")
    cat("1. Review the migration script: migrate_packages.R\n")
    cat("2. Migrate each file manually\n")
    cat("3. Test after each migration\n")
    cat("4. Run this script again to verify\n")
  }
  
  # Return results for programmatic use
  invisible(results)
}

# Run if executed as script
if (!interactive()) {
  results <- main()
} else {
  cat("Run main() to check for library/require calls\n")
}
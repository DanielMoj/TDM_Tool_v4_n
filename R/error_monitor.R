#' Error Monitoring System for TDMx
#' 
#' Central error tracking, aggregation, and alerting system
#' 
#' @description
#' This module provides comprehensive error monitoring with:
#' - Real-time error tracking
#' - Error aggregation and classification
#' - Automatic alerting for critical errors
#' - Performance metrics tracking
#' - Error pattern detection
#'
#' @note This module requires the following packages to be installed:
#' DBI, jsonlite, digest, lubridate
#' These are loaded via R/dependencies.R or checked with requireNamespace()

# Check for required packages at module load time
.check_error_monitor_dependencies <- function() {
  required_packages <- c("DBI", "jsonlite", "digest", "lubridate")
  missing_packages <- character()
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
    }
  }
  
  if (length(missing_packages) > 0) {
    stop(sprintf(
      "Error monitor requires the following packages: %s\nPlease install them with: install.packages(c(%s))",
      paste(missing_packages, collapse = ", "),
      paste(sprintf('"%s"', missing_packages), collapse = ", ")
    ))
  }
}

# Check dependencies when module is sourced
.check_error_monitor_dependencies()

# Configuration
.error_monitor_config <- new.env(parent = emptyenv())
.error_monitor_config$max_errors_in_memory <- 1000
.error_monitor_config$alert_threshold_critical <- 5
.error_monitor_config$alert_threshold_warning <- 20
.error_monitor_config$alert_window_minutes <- 5
.error_monitor_config$db_enabled <- FALSE
.error_monitor_config$email_alerts <- FALSE
.error_monitor_config$slack_webhook <- NULL

# Error storage
.error_store <- new.env(parent = emptyenv())
.error_store$errors <- list()
.error_store$metrics <- list()
.error_store$alerts_sent <- list()

#' Initialize error monitoring system
#' 
#' @param db_connection Optional database connection for persistence
#' @param config List with configuration overrides
#' @export
init_error_monitor <- function(db_connection = NULL, config = list()) {
  
  # Apply configuration
  for (name in names(config)) {
    if (name %in% ls(.error_monitor_config)) {
      .error_monitor_config[[name]] <- config[[name]]
    }
  }
  
  # Set up database if provided
  if (!is.null(db_connection)) {
    .error_monitor_config$db_connection <- db_connection
    .error_monitor_config$db_enabled <- TRUE
    
    # Create tables if they don't exist
    create_error_tables(db_connection)
  }
  
  # Initialize metrics
  .error_store$metrics <- list(
    total_errors = 0,
    errors_by_severity = list(CRITICAL = 0, ERROR = 0, WARNING = 0, INFO = 0),
    errors_by_module = list(),
    error_rate_per_minute = numeric(),
    last_error_time = NULL,
    uptime_start = Sys.time()
  )
  
  # Set up cleanup scheduler
  schedule_cleanup()
  
  message("Error monitoring system initialized")
  return(TRUE)
}

#' Log an error event
#' 
#' @param error_type Type/category of error
#' @param message Error message
#' @param severity One of: CRITICAL, ERROR, WARNING, INFO
#' @param module Module where error occurred
#' @param context Additional context (list)
#' @param stack_trace Optional stack trace
#' @return Error ID
#' @export
log_error <- function(error_type, message, severity = "ERROR", 
                     module = NA, context = list(), stack_trace = NULL) {
  
  # Validate severity
  valid_severities <- c("CRITICAL", "ERROR", "WARNING", "INFO")
  if (!severity %in% valid_severities) {
    severity <- "ERROR"
  }
  
  # Create error record using digest:: for ID generation
  error_id <- generate_error_id()
  
  error_record <- list(
    id = error_id,
    timestamp = Sys.time(),
    error_type = error_type,
    message = message,
    severity = severity,
    module = module,
    context = context,
    stack_trace = stack_trace %||% capture_stack_trace(),
    session_id = get_session_id(),
    user = get_current_user(),
    environment = Sys.getenv("R_ENV", "production")
  )
  
  # Store in memory (limited)
  if (length(.error_store$errors) >= .error_monitor_config$max_errors_in_memory) {
    # Remove oldest errors
    .error_store$errors <- tail(.error_store$errors, 
                               .error_monitor_config$max_errors_in_memory - 1)
  }
  .error_store$errors[[error_id]] <- error_record
  
  # Update metrics
  update_error_metrics(error_record)
  
  # Check alert conditions
  check_alert_conditions(error_record)
  
  # Persist to database if enabled
  if (.error_monitor_config$db_enabled) {
    store_error_in_db(error_record)
  }
  
  # Write to log file
  write_error_log(error_record)
  
  return(error_id)
}

#' Create database tables for error storage
create_error_tables <- function(con) {
  # Using DBI:: namespace
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS error_log (
      id TEXT PRIMARY KEY,
      timestamp TIMESTAMP,
      error_type TEXT,
      message TEXT,
      severity TEXT,
      module TEXT,
      context TEXT,
      stack_trace TEXT,
      session_id TEXT,
      user_name TEXT,
      environment TEXT
    )
  ")
  
  DBI::dbExecute(con, "
    CREATE INDEX IF NOT EXISTS idx_error_timestamp ON error_log(timestamp);
    CREATE INDEX IF NOT EXISTS idx_error_severity ON error_log(severity);
    CREATE INDEX IF NOT EXISTS idx_error_type ON error_log(error_type);
  ")
}

#' Store error in database
store_error_in_db <- function(error_record) {
  con <- .error_monitor_config$db_connection
  if (is.null(con)) return(FALSE)
  
  tryCatch({
    # Convert complex fields to JSON using jsonlite::
    context_json <- jsonlite::toJSON(error_record$context, auto_unbox = TRUE)
    stack_json <- jsonlite::toJSON(error_record$stack_trace, auto_unbox = TRUE)
    
    # Use DBI:: for database operations
    DBI::dbExecute(con, "
      INSERT INTO error_log (id, timestamp, error_type, message, severity, 
                            module, context, stack_trace, session_id, user_name, environment)
      VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ", list(
      error_record$id,
      error_record$timestamp,
      error_record$error_type,
      error_record$message,
      error_record$severity,
      error_record$module,
      context_json,
      stack_json,
      error_record$session_id,
      error_record$user,
      error_record$environment
    ))
    
    return(TRUE)
    
  }, error = function(e) {
    warning(sprintf("Failed to store error in database: %s", e$message))
    return(FALSE)
  })
}

#' Update error metrics
update_error_metrics <- function(error_record) {
  metrics <- .error_store$metrics
  
  # Update totals
  metrics$total_errors <- metrics$total_errors + 1
  metrics$errors_by_severity[[error_record$severity]] <- 
    (metrics$errors_by_severity[[error_record$severity]] %||% 0) + 1
  
  # Update by module
  if (!is.na(error_record$module)) {
    metrics$errors_by_module[[error_record$module]] <- 
      (metrics$errors_by_module[[error_record$module]] %||% 0) + 1
  }
  
  # Update error rate using lubridate::
  current_minute <- lubridate::floor_date(Sys.time(), "minute")
  minute_key <- as.character(current_minute)
  
  if (length(metrics$error_rate_per_minute) == 0 || 
      names(metrics$error_rate_per_minute)[length(metrics$error_rate_per_minute)] != minute_key) {
    metrics$error_rate_per_minute[minute_key] <- 1
  } else {
    metrics$error_rate_per_minute[minute_key] <- 
      metrics$error_rate_per_minute[minute_key] + 1
  }
  
  # Keep only last hour of rate data using lubridate::
  cutoff <- Sys.time() - lubridate::hours(1)
  keep_times <- as.POSIXct(names(metrics$error_rate_per_minute)) >= cutoff
  metrics$error_rate_per_minute <- metrics$error_rate_per_minute[keep_times]
  
  metrics$last_error_time <- error_record$timestamp
  
  .error_store$metrics <- metrics
}

#' Check alert conditions and send alerts if needed
check_alert_conditions <- function(error_record) {
  
  # Check for critical errors
  if (error_record$severity == "CRITICAL") {
    send_alert(
      level = "CRITICAL",
      subject = sprintf("CRITICAL ERROR: %s", error_record$error_type),
      message = format_alert_message(error_record),
      error_id = error_record$id
    )
    return()
  }
  
  # Check error rate thresholds using lubridate::
  window_start <- Sys.time() - lubridate::minutes(.error_monitor_config$alert_window_minutes)
  recent_errors <- Filter(function(e) e$timestamp >= window_start, .error_store$errors)
  
  error_count <- length(recent_errors)
  critical_count <- sum(sapply(recent_errors, function(e) e$severity == "CRITICAL"))
  
  # Alert on high error rate
  if (error_count >= .error_monitor_config$alert_threshold_warning) {
    if (!has_recent_alert("high_error_rate")) {
      send_alert(
        level = "WARNING",
        subject = sprintf("High error rate: %d errors in %d minutes", 
                         error_count, .error_monitor_config$alert_window_minutes),
        message = format_rate_alert(recent_errors)
      )
    }
  }
  
  # Alert on multiple critical errors
  if (critical_count >= .error_monitor_config$alert_threshold_critical) {
    if (!has_recent_alert("critical_errors")) {
      send_alert(
        level = "CRITICAL",
        subject = sprintf("Multiple critical errors: %d in %d minutes", 
                         critical_count, .error_monitor_config$alert_window_minutes),
        message = format_rate_alert(Filter(function(e) e$severity == "CRITICAL", recent_errors))
      )
    }
  }
  
  # Check for error patterns
  check_error_patterns(recent_errors)
}

#' Check for error patterns
check_error_patterns <- function(errors) {
  if (length(errors) < 5) return()
  
  # Check for repeated errors
  error_types <- sapply(errors, function(e) e$error_type)
  type_counts <- table(error_types)
  
  for (error_type in names(type_counts)) {
    if (type_counts[error_type] >= 5) {
      if (!has_recent_alert(paste0("repeated_", error_type))) {
        send_alert(
          level = "WARNING",
          subject = sprintf("Repeated error pattern: %s (%d times)", 
                           error_type, type_counts[error_type]),
          message = sprintf("Error '%s' has occurred %d times in the last %d minutes",
                           error_type, type_counts[error_type], 
                           .error_monitor_config$alert_window_minutes)
        )
      }
    }
  }
  
  # Check for module-specific issues
  modules <- sapply(errors, function(e) e$module)
  modules <- modules[!is.na(modules)]
  if (length(modules) > 0) {
    module_counts <- table(modules)
    
    for (module in names(module_counts)) {
      if (module_counts[module] >= 10) {
        if (!has_recent_alert(paste0("module_", module))) {
          send_alert(
            level = "WARNING",
            subject = sprintf("Module experiencing issues: %s", module),
            message = sprintf("Module '%s' has reported %d errors in the last %d minutes",
                            module, module_counts[module], 
                            .error_monitor_config$alert_window_minutes)
          )
        }
      }
    }
  }
}

#' Send alert notification
send_alert <- function(level, subject, message, error_id = NULL) {
  alert_record <- list(
    timestamp = Sys.time(),
    level = level,
    subject = subject,
    message = message,
    error_id = error_id
  )
  
  # Store alert record
  alert_key <- paste0(level, "_", gsub("[^a-zA-Z0-9]", "_", subject))
  .error_store$alerts_sent[[alert_key]] <- alert_record
  
  # Console output
  cat(sprintf("\n[%s ALERT] %s\n%s\n\n", level, subject, message))
  
  # Email alert if configured
  if (.error_monitor_config$email_alerts) {
    send_email_alert(alert_record)
  }
  
  # Slack webhook if configured
  if (!is.null(.error_monitor_config$slack_webhook)) {
    send_slack_alert(alert_record)
  }
  
  # Log alert
  log_alert(alert_record)
}

#' Check if alert was recently sent
has_recent_alert <- function(alert_key) {
  alert <- .error_store$alerts_sent[[alert_key]]
  if (is.null(alert)) return(FALSE)
  
  # Check if alert was sent in last 30 minutes using lubridate::
  alert$timestamp > (Sys.time() - lubridate::minutes(30))
}

#' Schedule cleanup of old errors
schedule_cleanup <- function() {
  # This would use later:: package in production
  # For now, just a placeholder
  # later::later(cleanup_old_errors, delay = 3600)  # Every hour
}

#' Clean up old errors from memory
cleanup_old_errors <- function() {
  # Keep only last 24 hours of errors using lubridate::
  cutoff <- Sys.time() - lubridate::hours(24)
  
  .error_store$errors <- Filter(
    function(e) e$timestamp >= cutoff, 
    .error_store$errors
  )
  
  # Clean old alerts using lubridate::
  alert_cutoff <- Sys.time() - lubridate::hours(6)
  .error_store$alerts_sent <- Filter(
    function(a) a$timestamp >= alert_cutoff,
    .error_store$alerts_sent
  )
}

#' Get error statistics
#' @export
get_error_stats <- function(hours = 24) {
  # Calculate time window using lubridate::
  since <- Sys.time() - lubridate::hours(hours)
  
  recent_errors <- Filter(
    function(e) e$timestamp >= since,
    .error_store$errors
  )
  
  if (length(recent_errors) == 0) {
    return(list(
      total = 0,
      by_severity = list(),
      by_module = list(),
      by_type = list(),
      rate_per_hour = 0,
      recent_errors = list()
    ))
  }
  
  # Aggregate statistics
  severities <- sapply(recent_errors, function(e) e$severity)
  modules <- sapply(recent_errors, function(e) e$module)
  types <- sapply(recent_errors, function(e) e$error_type)
  
  list(
    total = length(recent_errors),
    by_severity = table(severities),
    by_module = table(modules[!is.na(modules)]),
    by_type = table(types),
    rate_per_hour = length(recent_errors) / hours,
    recent_errors = tail(recent_errors, 10)
  )
}

#' Get error trends for visualization
#' @export
get_error_trends <- function(hours = 24, interval = "hour") {
  # Calculate time window using lubridate::
  since <- Sys.time() - lubridate::hours(hours)
  
  recent_errors <- Filter(
    function(e) e$timestamp >= since,
    .error_store$errors
  )
  
  if (length(recent_errors) == 0) {
    return(data.frame(
      time = character(),
      count = integer(),
      severity = character()
    ))
  }
  
  # Create time bins using lubridate::
  timestamps <- sapply(recent_errors, function(e) e$timestamp)
  severities <- sapply(recent_errors, function(e) e$severity)
  
  # Round to interval using lubridate::
  binned_times <- lubridate::floor_date(
    as.POSIXct(timestamps, origin = "1970-01-01"),
    interval
  )
  
  # Aggregate by time and severity
  trend_data <- aggregate(
    list(count = rep(1, length(binned_times))),
    list(time = binned_times, severity = severities),
    sum
  )
  
  trend_data
}

# ===== Helper Functions =====

#' Generate unique error ID using digest::
generate_error_id <- function() {
  paste0("ERR_", 
         format(Sys.time(), "%Y%m%d%H%M%S"), 
         "_",
         substr(digest::digest(runif(1)), 1, 8))
}

#' Capture current stack trace
capture_stack_trace <- function() {
  calls <- sys.calls()
  if (length(calls) > 2) {
    # Remove this function and log_error from trace
    calls <- calls[1:(length(calls) - 2)]
  }
  lapply(calls, deparse)
}

#' Get current session ID
get_session_id <- function() {
  if (exists("session") && !is.null(session$token)) {
    return(session$token)
  }
  return(NA)
}

#' Get current user
get_current_user <- function() {
  if (exists("session") && !is.null(session$user)) {
    return(session$user)
  }
  return(Sys.getenv("USER", "unknown"))
}

#' Write error to log file
write_error_log <- function(error_record) {
  log_dir <- "log/errors"
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  
  log_file <- file.path(log_dir, sprintf("errors_%s.log", format(Sys.Date(), "%Y%m%d")))
  
  log_entry <- sprintf(
    "[%s] %s | %s | %s | %s | %s",
    format(error_record$timestamp, "%Y-%m-%d %H:%M:%S"),
    error_record$severity,
    error_record$error_type,
    error_record$module,
    error_record$message,
    error_record$user
  )
  
  cat(log_entry, "\n", file = log_file, append = TRUE)
}

#' Format alert message
format_alert_message <- function(error_record) {
  # Using jsonlite:: for JSON formatting
  sprintf(
    "Error Details:\n\nType: %s\nSeverity: %s\nModule: %s\nMessage: %s\nTime: %s\nUser: %s\nSession: %s\nEnvironment: %s\n\nContext: %s",
    error_record$error_type,
    error_record$severity,
    error_record$module %||% "Unknown",
    error_record$message,
    error_record$timestamp,
    error_record$user,
    error_record$session_id %||% "N/A",
    error_record$environment,
    jsonlite::toJSON(error_record$context, pretty = TRUE, auto_unbox = TRUE)
  )
}

#' Format rate alert message
format_rate_alert <- function(errors) {
  summary <- table(sapply(errors, function(e) e$error_type))
  
  paste0(
    "Error Summary:\n\n",
    "Total errors: ", length(errors), "\n",
    "Time window: Last ", .error_monitor_config$alert_window_minutes, " minutes\n\n",
    "Error types:\n",
    paste(sprintf("  - %s: %d", names(summary), summary), collapse = "\n"),
    "\n\nMost recent errors:\n",
    paste(sapply(tail(errors, 5), function(e) {
      sprintf("  [%s] %s: %s", 
              format(e$timestamp, "%H:%M:%S"), 
              e$error_type, 
              substr(e$message, 1, 50))
    }), collapse = "\n")
  )
}

#' Send email alert (placeholder)
send_email_alert <- function(alert_record) {
  # Implementation would depend on mail configuration
  # Using sendmailR, mailR, or blastula package
  message("Email alert would be sent: ", alert_record$subject)
}

#' Send Slack alert (placeholder)
send_slack_alert <- function(alert_record) {
  # Implementation would use httr to post to webhook
  # httr::POST(
  #   url = .error_monitor_config$slack_webhook,
  #   body = jsonlite::toJSON(list(text = alert_record$message)),
  #   encode = "json"
  # )
  message("Slack alert would be sent: ", alert_record$subject)
}

#' Log alert to file
log_alert <- function(alert_record) {
  log_dir <- "log/alerts"
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  
  log_file <- file.path(log_dir, sprintf("alerts_%s.log", format(Sys.Date(), "%Y%m%d")))
  
  log_entry <- sprintf(
    "[%s] %s | %s | %s",
    format(alert_record$timestamp, "%Y-%m-%d %H:%M:%S"),
    alert_record$level,
    alert_record$subject,
    substr(alert_record$message, 1, 100)
  )
  
  cat(log_entry, "\n", file = log_file, append = TRUE)
}
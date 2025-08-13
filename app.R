# app.R - Modularized PK/PD Analysis Platform
# Version 2.1 - Modernized Package Management (no library() calls)
# 
# This application uses namespace references (::) instead of library()
# to avoid namespace pollution and conflicts

# Load all dependencies and modules with single source call
# This loads R/dependencies.R first, then all other modules
source("R/load_all.R")

# Check system and load all sources
if (!interactive()) {
  # Production mode - strict checking
  system_health_check()
  load_all_sources(check_dependencies = TRUE, install_missing = FALSE)
} else {
  # Development mode - already loaded in load_all.R
  if (!exists("load_all_sources")) {
    source("R/load_all.R")
    load_all_sources()
  }
}

# Define UI using shiny:: namespace
ui <- shinydashboard::dashboardPage(
  shinydashboard::dashboardHeader(
    title = "PK/PD Analysis Platform",
    tags$li(
      class = "dropdown",
      tags$a(
        href = "#",
        class = "dropdown-toggle",
        shiny::icon("user"),
        tags$span(class = "hidden-xs", shiny::textOutput("username_display"))
      )
    ),
    tags$li(
      class = "dropdown",
      shiny::actionButton("logout_btn", "Logout", icon = shiny::icon("sign-out-alt"))
    )
  ),
  
  shinydashboard::dashboardSidebar(
    shinydashboard::sidebarMenu(
      id = "sidebar",
      shinydashboard::menuItem("Dashboard", tabName = "dashboard", icon = shiny::icon("dashboard")),
      shinydashboard::menuItem("Data Import", tabName = "data", icon = shiny::icon("upload")),
      shinydashboard::menuItem("Model Fitting", tabName = "fitting", icon = shiny::icon("chart-line")),
      shinydashboard::menuItem("Diagnostics", tabName = "diagnostics", icon = shiny::icon("stethoscope")),
      shinydashboard::menuItem("Dose Optimization", tabName = "optimization", icon = shiny::icon("calculator")),
      shinydashboard::menuItem("Reports", tabName = "reports", icon = shiny::icon("file-pdf")),
      shinydashboard::menuItem("Administration", tabName = "admin", icon = shiny::icon("cogs"))
    )
  ),
  
  shinydashboard::dashboardBody(
    # Custom CSS
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side {
          background-color: #f4f4f4;
        }
        .box {
          border-radius: 5px;
          box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }
      "))
    ),
    
    # Enable shinyjs if available
    if (requireNamespace("shinyjs", quietly = TRUE)) {
      shinyjs::useShinyjs()
    } else {
      tags$div()
    },
    
    # Authentication UI (shown initially)
    shiny::div(id = "auth_ui", 
               if (exists("mod_auth_ui")) mod_auth_ui("auth") else tags$div()),
    
    # Main content (hidden initially)
    shinyjs::hidden(
      shiny::div(
        id = "main_ui",
        shinydashboard::tabItems(
          # Dashboard tab
          shinydashboard::tabItem(
            tabName = "dashboard",
            shiny::h2("Dashboard"),
            shiny::fluidRow(
              shinydashboard::box(
                title = "Welcome",
                status = "primary",
                width = 12,
                shiny::p("Welcome to the PK/PD Analysis Platform."),
                shiny::p("Use the sidebar to navigate through different analysis modules.")
              )
            ),
            shiny::fluidRow(
              shinydashboard::valueBox(
                value = "0",
                subtitle = "Analyses Completed",
                icon = shiny::icon("check-circle"),
                color = "green"
              ),
              shinydashboard::valueBox(
                value = "0",
                subtitle = "Models Fitted",
                icon = shiny::icon("chart-line"),
                color = "blue"
              ),
              shinydashboard::valueBox(
                value = "Ready",
                subtitle = "System Status",
                icon = shiny::icon("server"),
                color = "yellow"
              )
            )
          ),
          
          # Data Import tab
          shinydashboard::tabItem(
            tabName = "data",
            shiny::h2("Data Import"),
            shinydashboard::box(
              title = "Upload Data",
              status = "info",
              width = 12,
              shiny::fileInput("data_file", "Choose CSV File",
                              accept = c("text/csv", ".csv")),
              shiny::actionButton("load_example", "Load Example Data",
                                 class = "btn-warning"),
              shiny::hr(),
              shiny::verbatimTextOutput("upload_status"),
              DT::DTOutput("data_preview")
            )
          ),
          
          # Model Fitting tab
          shinydashboard::tabItem(
            tabName = "fitting",
            if (exists("mod_fitting_ui")) {
              mod_fitting_ui("fitting")
            } else {
              shiny::h2("Model Fitting")
            }
          ),
          
          # Diagnostics tab
          shinydashboard::tabItem(
            tabName = "diagnostics",
            if (exists("mod_diagnostics_ui")) {
              mod_diagnostics_ui("diagnostics")
            } else {
              shiny::h2("Diagnostics")
            }
          ),
          
          # Dose Optimization tab
          shinydashboard::tabItem(
            tabName = "optimization",
            if (exists("mod_optimization_ui")) {
              mod_optimization_ui("optimization")
            } else {
              shiny::h2("Dose Optimization")
            }
          ),
          
          # Reports tab
          shinydashboard::tabItem(
            tabName = "reports",
            shiny::h2("Reports"),
            shinydashboard::box(
              title = "Generate Report",
              status = "warning",
              width = 12,
              shiny::selectInput("report_type", "Report Type",
                                choices = c("Full Analysis", "Summary", "Technical")),
              shiny::downloadButton("download_report", "Download Report",
                                   class = "btn-primary")
            )
          ),
          
          # Admin tab
          shinydashboard::tabItem(
            tabName = "admin",
            if (exists("mod_admin_ui")) {
              mod_admin_ui("admin")
            } else {
              shiny::h2("Administration")
            }
          )
        )
      )
    )
  )
)

# Define Server using shiny:: namespace
server <- function(input, output, session) {
  
  # Initialize reactive values
  app_data <- shiny::reactiveVal(NULL)
  
  # Authentication module (if available)
  if (exists("mod_auth_server")) {
    auth <- mod_auth_server("auth")
  } else {
    # Mock auth for development
    auth <- list(
      is_logged_in = function() TRUE,
      get_user = function() "dev_user",
      logout = function() NULL
    )
  }
  
  # Show/hide UI based on authentication
  shiny::observe({
    if (auth$is_logged_in()) {
      if (requireNamespace("shinyjs", quietly = TRUE)) {
        shinyjs::hide("auth_ui")
        shinyjs::show("main_ui")
      }
      output$username_display <- shiny::renderText(auth$get_user())
    } else {
      if (requireNamespace("shinyjs", quietly = TRUE)) {
        shinyjs::show("auth_ui")
        shinyjs::hide("main_ui")
      }
    }
  })
  
  # Logout handler
  shiny::observeEvent(input$logout_btn, {
    auth$logout()
  })
  
  # Data upload handler
  shiny::observeEvent(input$data_file, {
    shiny::req(input$data_file)
    
    # Use readr:: if available, otherwise base R
    df <- if (requireNamespace("readr", quietly = TRUE)) {
      readr::read_csv(input$data_file$datapath, show_col_types = FALSE)
    } else {
      utils::read.csv(input$data_file$datapath)
    }
    
    app_data(df)
    
    output$upload_status <- shiny::renderPrint({
      cat("File uploaded successfully!\n")
      cat("Rows:", nrow(df), "\n")
      cat("Columns:", ncol(df))
    })
    
    output$data_preview <- DT::renderDT({
      DT::datatable(df, options = list(pageLength = 10))
    })
  })
  
  # Load example data
  shiny::observeEvent(input$load_example, {
    # Create example PK data
    set.seed(123)
    n_subjects <- 10
    times <- c(0, 0.5, 1, 2, 4, 6, 8, 12, 24)
    
    example_data <- expand.grid(
      ID = 1:n_subjects,
      TIME = times
    )
    
    # Simulate concentrations (simple 1-compartment model)
    example_data$DV <- with(example_data, {
      dose <- 1000
      cl <- 5 * exp(rnorm(n_subjects, 0, 0.3))[ID]
      v <- 30 * exp(rnorm(n_subjects, 0, 0.2))[ID]
      ke <- cl/v
      
      conc <- (dose/v) * exp(-ke * TIME)
      conc * exp(rnorm(length(conc), 0, 0.1))  # Add residual error
    })
    
    example_data$DOSE <- ifelse(example_data$TIME == 0, 1000, 0)
    example_data$MDV <- ifelse(example_data$TIME == 0, 1, 0)
    
    app_data(example_data)
    
    output$upload_status <- shiny::renderPrint({
      cat("Example data loaded!\n")
      cat("Subjects:", n_subjects, "\n")
      cat("Observations per subject:", length(times))
    })
    
    output$data_preview <- DT::renderDT({
      DT::datatable(example_data, options = list(pageLength = 10))
    })
  })
  
  # Module servers (if available)
  if (exists("mod_fitting_server")) {
    fitting_results <- mod_fitting_server("fitting", app_data)
  }
  
  if (exists("mod_diagnostics_server")) {
    mod_diagnostics_server("diagnostics", fitting_results)
  }
  
  if (exists("mod_optimization_server")) {
    mod_optimization_server("optimization", fitting_results)
  }
  
  if (exists("mod_admin_server")) {
    mod_admin_server("admin")
  }
  
  # Report generation
  output$download_report <- shiny::downloadHandler(
    filename = function() {
      paste0("pkpd_report_", format(Sys.Date(), "%Y%m%d"), ".pdf")
    },
    content = function(file) {
      # Check if rmarkdown is available
      if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        warning("rmarkdown package required for report generation")
        return(NULL)
      }
      
      # Create temporary report file
      temp_report <- tempfile(fileext = ".Rmd")
      
      # Simple report template
      report_content <- '
---
title: "PK/PD Analysis Report"
date: "`r Sys.Date()`"
output: pdf_document
---

## Summary

This report was generated on `r Sys.Date()`.

## Data Summary

```{r, echo=FALSE}
if (exists("app_data") && !is.null(app_data())) {
  summary(app_data())
} else {
  print("No data loaded")
}
```
'
      
      writeLines(report_content, temp_report)
      
      # Render the report
      rmarkdown::render(
        temp_report,
        output_file = file,
        quiet = TRUE
      )
    }
  )
  
  # Cleanup on session end
  session$onSessionEnded(function() {
    if (exists("cleanup_resources")) {
      cleanup_resources()
    }
  })
}

# Run the application using shiny::
shiny::shinyApp(ui = ui, server = server)
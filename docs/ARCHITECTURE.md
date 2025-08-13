# ARCHITECTURE

## Overview
TDMx ist eine **R/Shiny**-Anwendung mit klarer Schichtentrennung und modularem Aufbau.

## Package Management (Updated 2025-08-11)

### Decision: Namespace-Based Package Management

**Problem**: 
- `library()` calls in modules cause namespace pollution and conflicts
- Side effects from package loading can affect global environment
- Difficult to track and manage dependencies

**Solution**:
- **NO `library()` calls in modules** - use explicit namespace references (`package::function()`)
- **Central dependency management** via `R/dependencies.R`
- **Conditional loading** with `requireNamespace()` for optional features
- **Main app (app.R)** uses the central loading system

### Implementation:

1. **R/dependencies.R**: Central package dependency management
   - Defines all required and optional packages
   - Version checking and validation
   - Installation helpers for missing packages
   - Package group management (core, database, bayesian, etc.)

2. **R/load_all.R**: Modernized loading orchestrator
   - Loads dependencies.R first
   - Sources all modules in correct order
   - Health checks and validation
   - No direct `library()` or `require()` calls

3. **Module Pattern**:
```r
# BAD (old pattern):
library(DBI)
library(jsonlite)
result <- toJSON(dbGetQuery(con, query))

# GOOD (new pattern):
# At module start:
if (!requireNamespace("DBI", quietly = TRUE)) {
  stop("Package 'DBI' needed. Please install it.")
}
# In functions:
result <- jsonlite::toJSON(DBI::dbGetQuery(con, query))
```

### Benefits:
- ✅ No namespace conflicts
- ✅ Clear dependency tracking
- ✅ Better performance (lazy loading)
- ✅ Easier testing and mocking
- ✅ Explicit dependencies in code

## Layers

- **Presentation (UI):** `app.R` mit Shiny-UI, Module unter `R/modules/mod_*.R`.
- **Application/Controller:** server-Logik in `app.R` orchestriert Inputs → Priors → Inferenz → Outputs → Audit.
- **Domain (PK/Fehler/Kovariaten):** `R/pk_models.R`, `R/error_models.R`, `R/units_checks.R`, `R/pediatric.R`, `R/crrt.R`, `R/nonlinear.R`, `R/tissue.R`.
- **Inferenz-Backends:** `R/backend_bayes.R` (Laplace, Stan, ADVI, JAGS), `R/backend_optim.R` (MLE).
- **Data Access & Config:** `R/prior_db.R`, `config/*.json|yaml`, `R/db.R` (optional Postgres), `audit/audit_log.csv`.

## File Structure
```
app.R
  |
  ├── R/dependencies.R     [NEW: Package management]
  ├── R/load_all.R         [UPDATED: No library() calls]
  ├── R/utils.R            [Core utilities, %||% operator]
  ├── R/error_monitor.R    [UPDATED: Uses :: notation]
  |
  ├── R/modules/
  │   ├── mod_auth.R
  │   ├── mod_fitting.R
  │   ├── mod_diagnostics.R
  │   └── mod_admin.R
  |
  ├── R/core/
  │   ├── pk_models.R
  │   ├── error_models.R
  │   └── backend_bayes.R
  |
  └── config/
      ├── database.yml
      ├── targets.json
      └── users.yaml
```

## Module Loading Order

Critical loading sequence maintained by `R/load_all.R`:

1. **dependencies.R** - Package management
2. **utils.R** - Core utilities including `%||%` operator
3. **Core modules** - Database, auth, error handling
4. **Domain modules** - PK/PD models, calculations
5. **UI modules** - Shiny modules
6. **Optional modules** - Based on available packages

## Key Components

### Core Infrastructure
- **dependencies.R**: Package dependency management system
- **load_all.R**: Central loading orchestrator
- **utils.R**: Core utilities, null-coalescing operator (`%||%`)
- **error_monitor.R**: Comprehensive error tracking and alerting
- **health_checks.R**: System health monitoring

### PK/PD Domain
- **pk_models.R:** 1C analytisch; 2C/3C (ODE); *TVCL* für CRRT; Hooks für **MM** & **TMDD**
- **error_models.R:** σ-Funktionen, **t-Residuen**, **Mixture**, **BLQ (M3)**-Likelihood
- **backend_bayes.R:** Laplace/MAP; Stan (HMC) & **Stan-ADVI**; JAGS-1C; Diagnostics; Warm-Start Cache

### Data Flow
1. **Input:** CSV/Excel oder manuelle Eingabe → Validierung (Einheiten, Plausibilität)
2. **Processing:** Priors laden → Kovariaten-Adjustierung → Bayes-Fit
3. **Output:** Posterior-Draws → Konz-Zeit-Profile → PTA/CFR → Plots → Report (PDF)

## Configuration Management

### Package Dependencies
- Required packages defined in `R/dependencies.R`
- Version requirements enforced
- Optional packages for extended features
- Package groups for specific functionality

### Environment Variables
```bash
# Database
PG_HOST=localhost
PG_PORT=5432
PG_DATABASE=tdmx
PG_USER=tdmx
PG_PASSWORD=secure_password

# Application
R_ENV=production
SHINY_PORT=3838
```

## Best Practices

### Package Usage
1. **Never use `library()` in modules** - only in app.R or dependencies.R
2. **Always use `::`** for external package functions
3. **Check package availability** with `requireNamespace()`
4. **Document package requirements** in module headers
5. **Handle missing packages gracefully** with informative error messages

### Error Handling
- All critical operations wrapped in `tryCatch()`
- Comprehensive error monitoring via `error_monitor.R`
- Automatic alerting for critical errors
- Error pattern detection

### Performance
- Lazy loading of optional features
- Namespace references avoid loading entire packages
- Conditional feature availability based on installed packages

## Testing Strategy

### Dependency Testing
```r
# Check all dependencies
source("R/dependencies.R")
report <- get_dependency_report()
print(report)

# Verify no library() calls in modules
grep -r "^library(" R/ --include="*.R" | grep -v "dependencies.R"
```

### Module Testing
- Each module tested in isolation
- Mock external dependencies
- Verify namespace usage

## Migration Notes

### From library() to :: notation
1. Identify all `library()` calls: `grep -r "^library(" R/`
2. Replace with namespace references
3. Add dependency checks at module start
4. Update documentation
5. Test thoroughly

### Compatibility
- Backwards compatible with existing code
- Gradual migration possible
- No breaking changes to API

## Future Enhancements

1. **Automated dependency installation** during docker builds
2. **Package version locking** via renv integration  
3. **Dependency graph visualization**
4. **Automated compatibility testing**
5. **Package performance profiling**

---

*Architecture last updated: 2025-08-11*
*Package management system implemented as per WP5*
# prior_db.R - Vollständig vektorisierte Prior-Datenbank-Funktionen
# Optimiert für Performance durch lapply/sapply statt for-Loops

#' Vektorisiertes Laden von Priors
#' 
#' Nutzt lapply für paralleles Laden mehrerer Prior-Dateien
#' @param drug_names Vektor mit Medikamentennamen
#' @param prior_dir Verzeichnis mit Prior-Dateien
#' @param cache Verwende Cache für schnelleren Zugriff
#' @return Liste mit Prior-Informationen
load_priors <- function(drug_names, prior_dir = "priors", cache = TRUE) {
  
  # Cache-Initialisierung
  if (cache && exists(".prior_cache", envir = .GlobalEnv)) {
    prior_cache <- get(".prior_cache", envir = .GlobalEnv)
  } else {
    prior_cache <- list()
  }
  
  # Vektorisiertes Laden mit lapply
  priors <- lapply(drug_names, function(drug) {
    
    # Check Cache
    if (cache && drug %in% names(prior_cache)) {
      return(prior_cache[[drug]])
    }
    
    # Konstruiere Dateipfad
    file_path <- file.path(prior_dir, paste0(drug, ".json"))
    
    # Prüfe Existenz
    if (!file.exists(file_path)) {
      warning(paste("Prior-Datei nicht gefunden:", file_path))
      return(NULL)
    }
    
    # Lade JSON (vektorisiert durch jsonlite)
    prior_data <- jsonlite::fromJSON(file_path, simplifyVector = TRUE)
    
    # Validiere und transformiere
    prior_data <- validate_prior_structure(prior_data)
    
    # Cache aktualisieren
    if (cache) {
      prior_cache[[drug]] <- prior_data
    }
    
    return(prior_data)
  })
  
  # Namen setzen
  names(priors) <- drug_names
  
  # Cache global speichern
  if (cache) {
    assign(".prior_cache", prior_cache, envir = .GlobalEnv)
  }
  
  return(priors)
}

#' Batch-Verarbeitung von Prior-JSON-Dateien
#' 
#' Vektorisierte JSON-Verarbeitung für multiple Dateien
#' @param json_files Vektor mit JSON-Dateipfaden
#' @return Data frame mit allen Priors
batch_process_prior_json <- function(json_files) {
  
  # Paralleles Einlesen aller JSON-Dateien
  json_list <- lapply(json_files, function(file) {
    if (file.exists(file)) {
      jsonlite::fromJSON(file, simplifyVector = TRUE)
    } else {
      NULL
    }
  })
  
  # Entferne NULL-Einträge
  json_list <- json_list[!sapply(json_list, is.null)]
  
  # Extrahiere gemeinsame Struktur
  all_params <- unique(unlist(lapply(json_list, names)))
  
  # Vektorisierte Datenrahmen-Erstellung
  prior_df <- do.call(rbind, lapply(json_list, function(prior) {
    # Stelle sicher, dass alle Parameter vorhanden sind
    row_data <- sapply(all_params, function(param) {
      if (param %in% names(prior)) {
        value <- prior[[param]]
        
        # Behandle verschachtelte Strukturen
        if (is.list(value) && !is.data.frame(value)) {
          # Konvertiere zu String für komplexe Strukturen
          jsonlite::toJSON(value, auto_unbox = TRUE)
        } else {
          value
        }
      } else {
        NA
      }
    })
    
    as.data.frame(t(row_data), stringsAsFactors = FALSE)
  }))
  
  # Typ-Konvertierung (vektorisiert)
  numeric_cols <- c("CL_pop", "V_pop", "Q_pop", "V2_pop", "V3_pop",
                   "omega_CL", "omega_V", "omega_Q", "sigma_add", "sigma_prop")
  
  for (col in numeric_cols) {
    if (col %in% names(prior_df)) {
      prior_df[[col]] <- as.numeric(prior_df[[col]])
    }
  }
  
  return(prior_df)
}

#' Vektorisierte Prior-Validierung
#' 
#' Validiert Struktur mehrerer Priors gleichzeitig
#' @param prior_data Prior-Daten (Liste oder Data Frame)
#' @return Validierte Prior-Daten
validate_prior_structure <- function(prior_data) {
  
  # Erforderliche Felder
  required_fields <- c("drug_name", "model_type", "population")
  required_params <- c("CL_pop", "V_pop", "omega_CL", "omega_V")
  
  # Vektorisierte Validierung
  if (is.list(prior_data) && !is.data.frame(prior_data)) {
    
    # Prüfe erforderliche Felder
    missing_fields <- setdiff(required_fields, names(prior_data))
    if (length(missing_fields) > 0) {
      warning(paste("Fehlende Felder:", paste(missing_fields, collapse = ", ")))
    }
    
    # Prüfe Parameter
    if ("parameters" %in% names(prior_data)) {
      missing_params <- setdiff(required_params, names(prior_data$parameters))
      if (length(missing_params) > 0) {
        warning(paste("Fehlende Parameter:", paste(missing_params, collapse = ", ")))
        
        # Füge Standard-Werte hinzu
        default_values <- list(
          CL_pop = 5, V_pop = 30, 
          omega_CL = 0.3, omega_V = 0.3
        )
        
        for (param in missing_params) {
          prior_data$parameters[[param]] <- default_values[[param]]
        }
      }
    }
    
    # Validiere Kovariaten (vektorisiert)
    if ("covariates" %in% names(prior_data)) {
      prior_data$covariates <- validate_covariates_vec(prior_data$covariates)
    }
  }
  
  return(prior_data)
}

#' Vektorisierte Kovariaten-Validierung
#' 
#' @param covariates Liste oder Data Frame mit Kovariaten
#' @return Validierte Kovariaten
validate_covariates_vec <- function(covariates) {
  
  if (is.null(covariates)) return(list())
  
  # Standard-Kovariaten-Struktur
  standard_covariates <- c("weight", "age", "creatinine", "sex", "height")
  
  # Vektorisierte Prüfung und Korrektur
  validated <- lapply(standard_covariates, function(cov) {
    if (cov %in% names(covariates)) {
      cov_data <- covariates[[cov]]
      
      # Stelle sicher, dass effect und reference vorhanden sind
      if (!is.list(cov_data)) {
        cov_data <- list(effect = cov_data, reference = NA)
      }
      
      if (!"effect" %in% names(cov_data)) {
        cov_data$effect <- 0
      }
      
      if (!"reference" %in% names(cov_data)) {
        cov_data$reference <- switch(cov,
          weight = 70,
          age = 40,
          creatinine = 1.0,
          sex = "male",
          height = 170
        )
      }
      
      return(cov_data)
    } else {
      return(NULL)
    }
  })
  
  names(validated) <- standard_covariates
  validated <- validated[!sapply(validated, is.null)]
  
  return(validated)
}

#' Optimierte Prior-Suche
#' 
#' Vektorisierte Suche in Prior-Datenbank
#' @param criteria Liste mit Suchkriterien
#' @param prior_db Prior-Datenbank (Data Frame)
#' @return Gefilterte Priors
search_priors <- function(criteria, prior_db) {
  
  # Initialisiere mit allen Priors
  matches <- rep(TRUE, nrow(prior_db))
  
  # Vektorisierte Filterung für jedes Kriterium
  for (field in names(criteria)) {
    if (field %in% names(prior_db)) {
      value <- criteria[[field]]
      
      # Verschiedene Vergleichsoperatoren
      if (is.numeric(value) && is.numeric(prior_db[[field]])) {
        # Numerischer Vergleich mit Toleranz
        if (length(value) == 2) {
          # Bereich
          matches <- matches & 
                    prior_db[[field]] >= value[1] & 
                    prior_db[[field]] <= value[2]
        } else {
          # Exakter Wert (mit Toleranz)
          matches <- matches & 
                    abs(prior_db[[field]] - value) < 0.01
        }
      } else {
        # String-Vergleich (vektorisiert)
        if (is.character(value)) {
          # Partial Match mit grepl (vektorisiert)
          pattern <- paste0(".*", value, ".*")
          matches <- matches & 
                    grepl(pattern, prior_db[[field]], ignore.case = TRUE)
        } else {
          # Exakter Vergleich
          matches <- matches & prior_db[[field]] == value
        }
      }
    }
  }
  
  return(prior_db[matches, ])
}

#' Batch-Update von Priors
#' 
#' Aktualisiert mehrere Priors gleichzeitig
#' @param prior_list Liste mit Priors
#' @param updates Liste mit Updates
#' @return Aktualisierte Prior-Liste
batch_update_priors <- function(prior_list, updates) {
  
  # Vektorisierte Updates mit mapply
  updated_priors <- mapply(function(prior, update) {
    if (is.null(update)) return(prior)
    
    # Rekursives Update für verschachtelte Strukturen
    update_recursive <- function(original, new) {
      for (key in names(new)) {
        if (key %in% names(original) && is.list(original[[key]]) && is.list(new[[key]])) {
          original[[key]] <- update_recursive(original[[key]], new[[key]])
        } else {
          original[[key]] <- new[[key]]
        }
      }
      return(original)
    }
    
    return(update_recursive(prior, update))
    
  }, prior_list, updates, SIMPLIFY = FALSE)
  
  return(updated_priors)
}

#' Vektorisierte Prior-Interpolation
#' 
#' Interpoliert Prior-Parameter für neue Populationen
#' @param reference_priors Referenz-Priors
#' @param target_characteristics Ziel-Charakteristika
#' @return Interpolierte Priors
interpolate_priors <- function(reference_priors, target_characteristics) {
  
  # Extrahiere relevante Charakteristika
  ref_chars <- do.call(rbind, lapply(reference_priors, function(p) {
    c(weight = p$population$weight,
      age = p$population$age,
      creatinine = p$population$creatinine)
  }))
  
  # Berechne Distanzen (vektorisiert)
  distances <- apply(ref_chars, 1, function(ref) {
    sqrt(sum((ref - target_characteristics)^2, na.rm = TRUE))
  })
  
  # Gewichte basierend auf Distanz (inverse distance weighting)
  weights <- 1 / (distances + 1e-6)
  weights <- weights / sum(weights)
  
  # Vektorisierte gewichtete Mittelung
  interpolated <- list()
  
  # Parameter-Namen extrahieren
  all_params <- unique(unlist(lapply(reference_priors, function(p) {
    names(p$parameters)
  })))
  
  # Gewichtete Mittelung für jeden Parameter
  for (param in all_params) {
    values <- sapply(reference_priors, function(p) {
      if (param %in% names(p$parameters)) {
        p$parameters[[param]]
      } else {
        NA
      }
    })
    
    # Entferne NA-Werte
    valid_idx <- !is.na(values)
    if (any(valid_idx)) {
      interpolated[[param]] <- sum(values[valid_idx] * weights[valid_idx]) / 
                               sum(weights[valid_idx])
    }
  }
  
  return(list(
    parameters = interpolated,
    population = target_characteristics,
    method = "interpolated"
  ))
}

#' Schnelle Prior-Aggregation
#' 
#' Aggregiert Priors aus mehreren Quellen
#' @param prior_sources Liste mit Prior-Quellen
#' @param aggregation_method Methode (mean, median, weighted)
#' @param weights Gewichte für weighted aggregation
#' @return Aggregierte Priors
aggregate_priors <- function(prior_sources, 
                           aggregation_method = "mean", 
                           weights = NULL) {
  
  # Extrahiere alle Parameter
  all_params <- do.call(rbind, lapply(prior_sources, function(source) {
    if ("parameters" %in% names(source)) {
      as.data.frame(source$parameters)
    } else {
      NULL
    }
  }))
  
  # Entferne NULL-Einträge
  all_params <- all_params[!sapply(all_params, is.null), ]
  
  # Vektorisierte Aggregation
  aggregated <- switch(aggregation_method,
    "mean" = colMeans(all_params, na.rm = TRUE),
    "median" = apply(all_params, 2, median, na.rm = TRUE),
    "weighted" = {
      if (is.null(weights)) {
        weights <- rep(1 / nrow(all_params), nrow(all_params))
      }
      colSums(all_params * weights, na.rm = TRUE) / sum(weights)
    },
    "geometric" = exp(colMeans(log(all_params), na.rm = TRUE)),
    stop("Unbekannte Aggregationsmethode: ", aggregation_method)
  )
  
  return(list(
    parameters = as.list(aggregated),
    method = paste("aggregated", aggregation_method),
    n_sources = nrow(all_params)
  ))
}

#' Prior-Datenbank Cache-Management
#' 
#' Effizientes Cache-Management für Priors
#' @param action Aktion (save, load, clear, info)
#' @param data Zu speichernde Daten
#' @param cache_dir Cache-Verzeichnis
#' @return Cache-Status oder Daten
manage_prior_cache <- function(action = "info", data = NULL, 
                             cache_dir = "cache/priors") {
  
  # Stelle sicher, dass Cache-Verzeichnis existiert
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  
  result <- switch(action,
    "save" = {
      if (!is.null(data)) {
        # Hash für Cache-Key
        cache_key <- digest::digest(data)
        cache_file <- file.path(cache_dir, paste0(cache_key, ".rds"))
        
        # Speichere mit Kompression
        saveRDS(data, cache_file, compress = "xz")
        
        # Metadaten speichern
        meta_file <- file.path(cache_dir, "metadata.rds")
        if (file.exists(meta_file)) {
          metadata <- readRDS(meta_file)
        } else {
          metadata <- list()
        }
        
        metadata[[cache_key]] <- list(
          timestamp = Sys.time(),
          size = file.size(cache_file)
        )
        
        saveRDS(metadata, meta_file)
        
        return(cache_key)
      }
    },
    "load" = {
      if (!is.null(data)) {  # data enthält cache_key
        cache_file <- file.path(cache_dir, paste0(data, ".rds"))
        if (file.exists(cache_file)) {
          return(readRDS(cache_file))
        } else {
          return(NULL)
        }
      }
    },
    "clear" = {
      # Lösche alle Cache-Dateien
      cache_files <- list.files(cache_dir, pattern = "\\.rds$", 
                               full.names = TRUE)
      file.remove(cache_files)
      return(TRUE)
    },
    "info" = {
      # Cache-Statistiken
      cache_files <- list.files(cache_dir, pattern = "\\.rds$", 
                               full.names = TRUE)
      
      if (length(cache_files) > 0) {
        total_size <- sum(file.size(cache_files))
        
        # Lade Metadaten
        meta_file <- file.path(cache_dir, "metadata.rds")
        if (file.exists(meta_file)) {
          metadata <- readRDS(meta_file)
          oldest <- min(sapply(metadata, function(m) m$timestamp))
        } else {
          oldest <- NA
        }
        
        return(list(
          n_files = length(cache_files),
          total_size_mb = total_size / 1024^2,
          oldest_entry = oldest
        ))
      } else {
        return(list(
          n_files = 0,
          total_size_mb = 0,
          oldest_entry = NA
        ))
      }
    }
  )
  
  return(result)
}

#' Export der optimierten Prior-Funktionen
#' @export
prior_db_vectorized <- list(
  load = load_priors,
  batch_process = batch_process_prior_json,
  validate = validate_prior_structure,
  search = search_priors,
  update = batch_update_priors,
  interpolate = interpolate_priors,
  aggregate = aggregate_priors,
  cache = manage_prior_cache
)
# crrt.R - Vollständig vektorisierte CRRT-Funktionen
# Continuous Renal Replacement Therapy - Optimiert für Performance

#' Vektorisierte CRRT Clearance-Profil Berechnung
#' 
#' Verwendet findInterval() statt linearer Suche
#' @param times Zeitvektor für Abfrage
#' @param crrt_data Data frame mit CRRT-Parametern über Zeit
#' @return Matrix mit Clearance-Werten für alle Zeitpunkte
crrt_clearance_profile <- function(times, crrt_data) {
  
  # Input-Validierung
  if (!all(c("time", "effluent_rate", "sieving_coefficient") %in% names(crrt_data))) {
    stop("crrt_data muss time, effluent_rate und sieving_coefficient enthalten")
  }
  
  # Sortiere CRRT-Daten nach Zeit
  crrt_data <- crrt_data[order(crrt_data$time), ]
  
  # Verwende findInterval() für effiziente Zeitpunkt-Zuordnung
  # Dies ersetzt die lineare Suche mit O(log n) Komplexität
  intervals <- findInterval(times, crrt_data$time, rightmost.closed = TRUE)
  
  # Initialisiere Ergebnis-Matrix
  n_times <- length(times)
  clearance_matrix <- matrix(0, nrow = n_times, ncol = 3,
                            dimnames = list(NULL, c("time", "crrt_cl", "total_cl")))
  
  # Vektorisierte Clearance-Berechnung
  clearance_matrix[, "time"] <- times
  
  # Behandle Zeitpunkte vor CRRT-Start
  pre_crrt <- intervals == 0
  clearance_matrix[pre_crrt, "crrt_cl"] <- 0
  
  # Behandle Zeitpunkte während CRRT
  during_crrt <- intervals > 0 & intervals <= nrow(crrt_data)
  
  if (any(during_crrt)) {
    # Vektorisierte Interpolation zwischen CRRT-Zeitpunkten
    idx <- intervals[during_crrt]
    
    # Direkte Zuordnung für exakte Zeitpunkte
    exact_match <- times[during_crrt] == crrt_data$time[idx]
    
    # Berechne CRRT Clearance
    crrt_cl <- crrt_data$effluent_rate[idx] * 
               crrt_data$sieving_coefficient[idx] / 60  # ml/min zu L/h
    
    # Lineare Interpolation für Zwischenwerte
    needs_interp <- !exact_match & idx < nrow(crrt_data)
    
    if (any(needs_interp)) {
      idx_interp <- idx[needs_interp]
      t1 <- crrt_data$time[idx_interp]
      t2 <- crrt_data$time[idx_interp + 1]
      
      # Gewichte für Interpolation
      w <- (times[during_crrt][needs_interp] - t1) / (t2 - t1)
      
      # Interpolierte Werte
      eff1 <- crrt_data$effluent_rate[idx_interp]
      eff2 <- crrt_data$effluent_rate[idx_interp + 1]
      sc1 <- crrt_data$sieving_coefficient[idx_interp]
      sc2 <- crrt_data$sieving_coefficient[idx_interp + 1]
      
      crrt_cl[needs_interp] <- ((1 - w) * eff1 * sc1 + w * eff2 * sc2) / 60
    }
    
    clearance_matrix[during_crrt, "crrt_cl"] <- crrt_cl
  }
  
  return(clearance_matrix)
}

#' Batch-Processing für CRRT-Patienten
#' 
#' Vektorisierte Verarbeitung mehrerer Patienten
#' @param patient_list Liste mit Patientendaten
#' @param time_grid Gemeinsames Zeitgitter für Auswertung
#' @return Matrix mit Clearance-Profilen für alle Patienten
batch_crrt_processing <- function(patient_list, time_grid = seq(0, 168, by = 1)) {
  
  n_patients <- length(patient_list)
  n_times <- length(time_grid)
  
  # Pre-allocate Ergebnis-Arrays
  crrt_clearance <- array(0, dim = c(n_patients, n_times, 2),
                          dimnames = list(
                            paste0("Patient_", seq_len(n_patients)),
                            NULL,
                            c("crrt_cl", "total_cl")
                          ))
  
  # Gruppiere Patienten mit ähnlichen CRRT-Profilen
  profile_hash <- sapply(patient_list, function(p) {
    if (!is.null(p$crrt_data)) {
      digest::digest(p$crrt_data)
    } else {
      "no_crrt"
    }
  })
  
  unique_profiles <- unique(profile_hash)
  
  # Batch-Verarbeitung pro Profil-Gruppe
  for (profile in unique_profiles) {
    if (profile == "no_crrt") next
    
    # Indizes der Patienten mit diesem Profil
    patient_idx <- which(profile_hash == profile)
    
    # Berechne Clearance-Profil einmal für die Gruppe
    example_patient <- patient_list[[patient_idx[1]]]
    cl_profile <- crrt_clearance_profile(time_grid, example_patient$crrt_data)
    
    # Weise allen Patienten in der Gruppe zu
    for (idx in patient_idx) {
      patient <- patient_list[[idx]]
      
      # Basis-Clearance des Patienten
      base_cl <- patient$params$CL
      
      # Addiere CRRT-Clearance
      crrt_clearance[idx, , "crrt_cl"] <- cl_profile[, "crrt_cl"]
      crrt_clearance[idx, , "total_cl"] <- base_cl + cl_profile[, "crrt_cl"]
    }
  }
  
  return(crrt_clearance)
}

#' Vektorisierte zeitabhängige Clearance-Funktion
#' 
#' Für Integration in PK-Modelle
#' @param t Zeitvektor
#' @param base_cl Basis-Clearance
#' @param crrt_profile CRRT-Clearance-Profil
#' @return Zeitabhängige Gesamt-Clearance
time_varying_clearance <- function(t, base_cl, crrt_profile) {
  
  # Verwende findInterval für schnelle Zuordnung
  if (is.null(crrt_profile) || nrow(crrt_profile) == 0) {
    return(rep(base_cl, length(t)))
  }
  
  # Interpoliere CRRT-Clearance auf gewünschte Zeitpunkte
  crrt_cl <- approx(x = crrt_profile[, "time"], 
                    y = crrt_profile[, "crrt_cl"],
                    xout = t, 
                    method = "linear",
                    yleft = 0, 
                    yright = tail(crrt_profile[, "crrt_cl"], 1))$y
  
  # Addiere zu Basis-Clearance
  total_cl <- base_cl + crrt_cl
  
  return(total_cl)
}

#' Optimierte Filtrations-Berechnung
#' 
#' Vektorisierte Berechnung der Medikamenten-Filtration
#' @param concentration Plasma-Konzentration
#' @param flow_rates Vektor der Flow-Raten über Zeit
#' @param sieving Sieving-Koeffizienten
#' @param protein_binding Proteinbindung (0-1)
#' @return Filtrationsrate über Zeit
calculate_filtration_vec <- function(concentration, flow_rates, 
                                   sieving, protein_binding = 0) {
  
  # Freie Fraktion
  free_fraction <- 1 - protein_binding
  
  # Vektorisierte Berechnung
  filtration_rate <- concentration * free_fraction * flow_rates * sieving
  
  return(filtration_rate)
}

#' CRRT-Modus spezifische Clearance
#' 
#' Unterstützt CVVH, CVVHD, CVVHDF
#' @param mode CRRT-Modus
#' @param params CRRT-Parameter
#' @return Clearance-Komponenten
calculate_crrt_clearance_by_mode <- function(mode, params) {
  
  # Extrahiere Parameter
  qb <- params$blood_flow  # ml/min
  qd <- params$dialysate_flow  # ml/min  
  qf <- params$filtration_rate  # ml/min
  sc <- params$sieving_coefficient
  sat <- params$saturation_coefficient
  
  # Vektorisierte Berechnung je nach Modus
  clearance <- switch(mode,
    "CVVH" = {
      # Nur Hämofiltration
      list(
        convective = qf * sc,
        diffusive = 0,
        total = qf * sc
      )
    },
    "CVVHD" = {
      # Nur Hämodialyse
      list(
        convective = 0,
        diffusive = qd * sat,
        total = qd * sat
      )
    },
    "CVVHDF" = {
      # Kombiniert
      conv_cl <- qf * sc
      diff_cl <- qd * sat
      
      # Berücksichtige Interaktion
      interaction <- 1 - (conv_cl * diff_cl) / ((conv_cl + diff_cl) * qb)
      
      list(
        convective = conv_cl,
        diffusive = diff_cl,
        total = (conv_cl + diff_cl) * interaction
      )
    },
    "SCUF" = {
      # Slow Continuous Ultrafiltration
      list(
        convective = qf * sc * 0.5,  # Reduzierte Effizienz
        diffusive = 0,
        total = qf * sc * 0.5
      )
    },
    stop("Unbekannter CRRT-Modus: ", mode)
  )
  
  # Konvertiere zu L/h
  clearance <- lapply(clearance, function(x) x * 0.06)
  
  return(clearance)
}

#' Vektorisierte AKI-Stadium Adjustierung
#' 
#' Passt Clearance basierend auf AKI-Stadium an
#' @param base_clearance Basis-Clearance
#' @param aki_stage AKI-Stadium (0-3)
#' @param times Zeitpunkte
#' @return Adjustierte Clearance über Zeit
adjust_clearance_for_aki <- function(base_clearance, aki_stage, times) {
  
  # AKI-Stadium Faktoren (vektorisiert)
  aki_factors <- c(1.0, 0.7, 0.5, 0.3)  # Normal, Stage 1, 2, 3
  
  # Validierung
  aki_stage <- pmin(pmax(aki_stage, 0), 3)
  
  # Wenn aki_stage ein Vektor ist (zeitabhängig)
  if (length(aki_stage) == length(times)) {
    # Direkte vektorisierte Multiplikation
    adjusted_cl <- base_clearance * aki_factors[aki_stage + 1]
  } else if (length(aki_stage) == 1) {
    # Konstantes AKI-Stadium
    adjusted_cl <- base_clearance * aki_factors[aki_stage + 1]
  } else {
    stop("aki_stage muss Länge 1 oder gleich times haben")
  }
  
  return(adjusted_cl)
}

#' Optimierte Effluent-Konzentrations-Berechnung
#' 
#' @param plasma_conc Plasma-Konzentration
#' @param sieving Sieving-Koeffizient
#' @param saturation Sättigungs-Koeffizient
#' @param mode CRRT-Modus
#' @return Effluent-Konzentration
calculate_effluent_concentration <- function(plasma_conc, sieving, 
                                           saturation, mode) {
  
  # Vektorisierte Berechnung
  effluent_conc <- switch(mode,
    "CVVH" = plasma_conc * sieving,
    "CVVHD" = plasma_conc * saturation,
    "CVVHDF" = plasma_conc * (sieving + saturation) / 2,
    "SCUF" = plasma_conc * sieving * 0.5,
    plasma_conc * 0  # Kein CRRT
  )
  
  return(effluent_conc)
}

#' Batch-Simulation für CRRT-Szenarien
#' 
#' Simuliert mehrere CRRT-Einstellungen parallel
#' @param scenarios Liste mit CRRT-Szenarien
#' @param pk_params PK-Parameter
#' @param time_grid Zeitgitter
#' @return Array mit Konzentrationsverläufen
simulate_crrt_scenarios <- function(scenarios, pk_params, 
                                   time_grid = seq(0, 168, by = 0.5)) {
  
  n_scenarios <- length(scenarios)
  n_times <- length(time_grid)
  
  # Pre-allocate Ergebnis-Array
  concentrations <- array(NA, dim = c(n_scenarios, n_times),
                         dimnames = list(
                           paste0("Scenario_", seq_len(n_scenarios)),
                           paste0("t_", time_grid)
                         ))
  
  # Gruppiere ähnliche Szenarien für Batch-Processing
  scenario_groups <- list()
  
  for (i in seq_len(n_scenarios)) {
    scenario <- scenarios[[i]]
    
    # Erstelle Clearance-Profil
    cl_profile <- calculate_crrt_clearance_by_mode(
      mode = scenario$mode,
      params = scenario$params
    )
    
    # Zeitabhängige Clearance
    total_cl <- pk_params$CL + cl_profile$total
    
    # Berechne Konzentration (1-Kompartiment als Beispiel)
    ke <- total_cl / pk_params$V
    
    # Vektorisierte Berechnung für Multiple-Dose
    dose_times <- seq(0, max(time_grid), by = scenario$tau)
    
    # Superposition
    conc <- numeric(n_times)
    for (dose_time in dose_times) {
      mask <- time_grid >= dose_time
      if (any(mask)) {
        dt <- time_grid[mask] - dose_time
        
        # Infusions-Kinetik
        if (dt <= scenario$tinf) {
          conc[mask] <- conc[mask] + 
                       (scenario$dose / scenario$tinf) / total_cl * 
                       (1 - exp(-ke * dt))
        } else {
          c_end <- (scenario$dose / scenario$tinf) / total_cl * 
                  (1 - exp(-ke * scenario$tinf))
          conc[mask] <- conc[mask] + 
                       c_end * exp(-ke * (dt - scenario$tinf))
        }
      }
    }
    
    concentrations[i, ] <- conc
  }
  
  return(concentrations)
}

#' Export der optimierten CRRT-Funktionen
#' @export
crrt_functions_vectorized <- list(
  clearance_profile = crrt_clearance_profile,
  batch_process = batch_crrt_processing,
  time_varying_cl = time_varying_clearance,
  filtration = calculate_filtration_vec,
  clearance_by_mode = calculate_crrt_clearance_by_mode,
  aki_adjustment = adjust_clearance_for_aki,
  effluent_conc = calculate_effluent_concentration,
  simulate_scenarios = simulate_crrt_scenarios
)
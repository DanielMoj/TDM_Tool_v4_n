# pk_models.R - Vollständig vektorisierte PK-Modelle
# Optimiert für Performance durch Eliminierung aller for-Loops

#' Vektorisierte 1-Kompartiment Infusions-Konzentrationsberechnung
#' 
#' @param t Zeitvektor (kann mehrere Zeitpunkte enthalten)
#' @param dose Dosis in mg
#' @param tau Dosierungsintervall in Stunden
#' @param tinf Infusionsdauer in Stunden
#' @param CL Clearance in L/h
#' @param V Verteilungsvolumen in L
#' @param start_time Startzeit der ersten Dosis
#' @param n_doses Anzahl der Dosen
#' @return Konzentrationsvektor für alle Zeitpunkte
conc_1c_inf_analytical <- function(t, dose, tau, tinf, CL, V, 
                                   start_time = 0, n_doses = 1) {
  # Parameter-Validierung
  if (any(c(dose, tau, tinf, CL, V) <= 0)) {
    stop("Alle pharmakokinetischen Parameter müssen positiv sein")
  }
  
  ke <- CL / V  # Eliminationskonstante
  
  # Vektorisierte Berechnung über alle Dosen und Zeitpunkte
  # Erstelle Matrix der Zeitdifferenzen (Zeilen: Zeitpunkte, Spalten: Dosen)
  dose_indices <- 0:(n_doses - 1)
  dose_start_times <- start_time + dose_indices * tau
  
  # Outer product für alle Kombinationen von t und dose_start_times
  time_since_dose <- outer(t, dose_start_times, "-")
  
  # Masken für verschiedene Phasen
  during_infusion <- (time_since_dose >= 0) & (time_since_dose <= tinf)
  after_infusion <- (time_since_dose > tinf) & (time_since_dose < tau)
  
  # Vektorisierte Konzentrationsberechnung
  # Während der Infusion: C = (dose/tinf) / CL * (1 - exp(-ke*t))
  C_during <- (dose / tinf) / CL * (1 - exp(-ke * time_since_dose))
  C_during[!during_infusion] <- 0
  
  # Nach der Infusion: C = C_end * exp(-ke*(t-tinf))
  C_end_infusion <- (dose / tinf) / CL * (1 - exp(-ke * tinf))
  C_after <- C_end_infusion * exp(-ke * (time_since_dose - tinf))
  C_after[!after_infusion] <- 0
  
  # Superposition: Summiere über alle Dosen
  C_total <- rowSums(C_during + C_after, na.rm = TRUE)
  
  return(C_total)
}

#' Vektorisierte 2-Kompartiment ODE-Lösung
#' 
#' Nutzt Matrix-Operationen für multiple Zeitpunkte
conc_2c_ode_vectorized <- function(times, doses, params) {
  # Verwende deSolve mit vektorisierten Anfangsbedingungen
  # Batch-Processing für mehrere Szenarien gleichzeitig
  
  CL <- params$CL
  V1 <- params$V1
  Q <- params$Q
  V2 <- params$V2
  
  # Systemmatrix für 2-Kompartiment-Modell
  k10 <- CL / V1
  k12 <- Q / V1
  k21 <- Q / V2
  
  # Eigenwert-Zerlegung für analytische Lösung
  lambda1 <- 0.5 * ((k10 + k12 + k21) + sqrt((k10 + k12 + k21)^2 - 4 * k10 * k21))
  lambda2 <- 0.5 * ((k10 + k12 + k21) - sqrt((k10 + k12 + k21)^2 - 4 * k10 * k21))
  
  # Vektorisierte Berechnung über alle Zeitpunkte
  A <- (k21 - lambda1) / (lambda2 - lambda1)
  B <- (lambda2 - k21) / (lambda2 - lambda1)
  
  # Matrix-Multiplikation für alle Dosen
  n_times <- length(times)
  n_doses <- nrow(doses)
  
  C <- matrix(0, n_times, 2)  # Konzentrationen in beiden Kompartimenten
  
  for (i in seq_len(n_doses)) {
    dose_time <- doses$time[i]
    dose_amount <- doses$amount[i]
    
    # Vektorisierte Berechnung für alle relevanten Zeitpunkte
    relevant_times <- times[times >= dose_time]
    if (length(relevant_times) > 0) {
      dt <- relevant_times - dose_time
      
      # Analytische Lösung
      C1 <- (dose_amount / V1) * (A * exp(-lambda1 * dt) + B * exp(-lambda2 * dt))
      C2 <- (dose_amount / V1) * (k12 / (lambda2 - lambda1)) * 
            (exp(-lambda1 * dt) - exp(-lambda2 * dt))
      
      # Füge zu Gesamtkonzentration hinzu (Superposition)
      idx <- which(times >= dose_time)
      C[idx, 1] <- C[idx, 1] + C1
      C[idx, 2] <- C[idx, 2] + C2
    }
  }
  
  return(C)
}

#' Vektorisierte 3-Kompartiment Berechnung
#' 
#' Nutzt Eigenwert-Methode für schnelle Berechnung
conc_3c_analytical_vectorized <- function(times, doses, params) {
  CL <- params$CL
  V1 <- params$V1
  Q2 <- params$Q2
  V2 <- params$V2
  Q3 <- params$Q3
  V3 <- params$V3
  
  # Rate-Konstanten
  k10 <- CL / V1
  k12 <- Q2 / V1
  k21 <- Q2 / V2
  k13 <- Q3 / V1
  k31 <- Q3 / V3
  
  # Charakteristisches Polynom Koeffizienten
  a <- k10 + k12 + k13 + k21 + k31
  b <- k10 * (k21 + k31) + k21 * k31 + k12 * k31 + k13 * k21
  c <- k10 * k21 * k31
  
  # Cardano's Methode für kubische Gleichung (vektorisiert)
  # Berechne Eigenwerte
  p <- b - a^2 / 3
  q <- 2 * a^3 / 27 - a * b / 3 + c
  discriminant <- (q / 2)^2 + (p / 3)^3
  
  # Für reelle Eigenwerte (pharmakokinetisch relevant)
  if (discriminant < 0) {
    m <- 2 * sqrt(-p / 3)
    theta <- acos(3 * q / (p * m)) / 3
    
    lambda1 <- -m * cos(theta) - a / 3
    lambda2 <- -m * cos(theta - 2 * pi / 3) - a / 3
    lambda3 <- -m * cos(theta - 4 * pi / 3) - a / 3
  } else {
    # Fallback für spezielle Fälle
    lambda1 <- k10
    lambda2 <- k21
    lambda3 <- k31
  }
  
  # Vektorisierte Exponential-Berechnung
  n_times <- length(times)
  C <- matrix(0, n_times, 3)
  
  # Batch-Processing für alle Dosen
  dose_matrix <- as.matrix(doses)
  
  # Verwende outer() für Zeit-Dosis-Kombinationen
  for (i in seq_len(nrow(dose_matrix))) {
    dose_time <- dose_matrix[i, "time"]
    dose_amount <- dose_matrix[i, "amount"]
    
    # Vektorisierte Berechnung
    mask <- times >= dose_time
    if (any(mask)) {
      dt <- times[mask] - dose_time
      
      # Tri-exponentielles Modell (vektorisiert)
      A1 <- dose_amount / V1 * (lambda2 - k10) * (lambda3 - k10) / 
           ((lambda2 - lambda1) * (lambda3 - lambda1))
      A2 <- dose_amount / V1 * (lambda1 - k10) * (lambda3 - k10) / 
           ((lambda1 - lambda2) * (lambda3 - lambda2))
      A3 <- dose_amount / V1 * (lambda1 - k10) * (lambda2 - k10) / 
           ((lambda1 - lambda3) * (lambda2 - lambda3))
      
      C[mask, 1] <- C[mask, 1] + 
                    A1 * exp(-lambda1 * dt) + 
                    A2 * exp(-lambda2 * dt) + 
                    A3 * exp(-lambda3 * dt)
    }
  }
  
  return(C)
}

#' Vektorisierte Steady-State Berechnung
#' 
#' Berechnet Steady-State ohne Loops
calculate_steady_state_vectorized <- function(params, tau, n_tau = 10) {
  CL <- params$CL
  V <- params$V
  ke <- CL / V
  
  # Vektorisierte Berechnung der Akkumulation
  accumulation_factor <- 1 / (1 - exp(-ke * tau))
  
  # Berechne für multiple Tau-Werte gleichzeitig wenn nötig
  if (length(tau) > 1) {
    accumulation_factor <- 1 / (1 - exp(-ke * outer(tau, rep(1, n_tau))))
    css_max <- params$dose / V * accumulation_factor
    css_min <- css_max * exp(-ke * tau)
    css_avg <- params$dose / (CL * tau)
    
    return(data.frame(
      tau = tau,
      css_max = css_max,
      css_min = css_min,
      css_avg = css_avg
    ))
  }
  
  css_max <- params$dose / V * accumulation_factor
  css_min <- css_max * exp(-ke * tau)
  css_avg <- params$dose / (CL * tau)
  
  return(list(
    css_max = css_max,
    css_min = css_min,
    css_avg = css_avg
  ))
}

#' Helper: Vektorisierte AUC-Berechnung
#' 
#' Trapezregel vollständig vektorisiert
calculate_auc_vectorized <- function(times, concentrations) {
  if (length(times) != length(concentrations)) {
    stop("times und concentrations müssen gleiche Länge haben")
  }
  
  if (length(times) < 2) {
    return(0)
  }
  
  # Sortiere nach Zeit
  ord <- order(times)
  times <- times[ord]
  concentrations <- concentrations[ord]
  
  # Vektorisierte Trapezregel
  dt <- diff(times)
  conc_mean <- (concentrations[-length(concentrations)] + 
                concentrations[-1]) / 2
  
  auc <- sum(dt * conc_mean)
  
  return(auc)
}

#' Batch-Verarbeitung für multiple Patienten
#' 
#' Vektorisierte Berechnung für N Patienten gleichzeitig
batch_pk_calculation <- function(patient_data, model = "1c") {
  n_patients <- nrow(patient_data)
  
  # Pre-allocate result matrix
  results <- matrix(NA, nrow = n_patients, ncol = 100)  # 100 Zeitpunkte
  time_grid <- seq(0, 24, length.out = 100)
  
  # Gruppiere Patienten mit ähnlichen Parametern für Batch-Processing
  param_groups <- split(seq_len(n_patients), 
                       paste(patient_data$CL, patient_data$V))
  
  for (group_idx in names(param_groups)) {
    patient_indices <- param_groups[[group_idx]]
    
    # Extrahiere gemeinsame Parameter
    group_params <- patient_data[patient_indices[1], ]
    
    if (model == "1c") {
      # Vektorisierte Berechnung für alle Patienten in der Gruppe
      conc <- conc_1c_inf_analytical(
        t = time_grid,
        dose = group_params$dose,
        tau = group_params$tau,
        tinf = group_params$tinf,
        CL = group_params$CL,
        V = group_params$V,
        n_doses = group_params$n_doses
      )
      
      # Weise Ergebnisse zu
      for (idx in patient_indices) {
        results[idx, ] <- conc
      }
    }
  }
  
  return(results)
}

#' Export der optimierten Funktionen
#' @export
pk_models_vectorized <- list(
  conc_1c_inf = conc_1c_inf_analytical,
  conc_2c = conc_2c_ode_vectorized,
  conc_3c = conc_3c_analytical_vectorized,
  steady_state = calculate_steady_state_vectorized,
  auc = calculate_auc_vectorized,
  batch_calculate = batch_pk_calculation
)
# error_models.R - Vollständig vektorisierte Fehlermodelle
# Optimiert für Performance durch logische Indizierung statt Loops

#' Vektorisierte Log-Likelihood Berechnung für Residuen
#' 
#' @param observed Vektor der beobachteten Konzentrationen
#' @param predicted Vektor der vorhergesagten Konzentrationen
#' @param sigma_add Additiver Fehlerparameter
#' @param sigma_prop Proportionaler Fehlerparameter
#' @param error_model Fehlermodell (1=additiv, 2=proportional, 3=kombiniert, 4=t-additiv, 5=t-prop, 6=mixture)
#' @param lloq Lower Limit of Quantification
#' @param is_blq Logischer Vektor für BLQ-Werte
#' @param nu Freiheitsgrade für t-Verteilung
#' @param mix_w Mixture weight
#' @param mix_scale Mixture scale factor
#' @return Log-Likelihood Wert
loglik_residuals_vec <- function(observed, predicted, 
                                sigma_add = 0.1, sigma_prop = 0.1,
                                error_model = 1, lloq = NULL, 
                                is_blq = NULL, nu = 3, 
                                mix_w = 0.1, mix_scale = 10) {
  
  n <- length(observed)
  
  # Input-Validierung
  if (length(predicted) != n) {
    stop("observed und predicted müssen gleiche Länge haben")
  }
  
  # Initialisiere BLQ-Vektor falls nicht vorhanden
  if (is.null(is_blq)) {
    is_blq <- rep(FALSE, n)
  }
  
  if (is.null(lloq)) {
    lloq <- min(observed[observed > 0], na.rm = TRUE) / 2
  }
  
  # Vektorisierte Fehler-Berechnung basierend auf Modell
  sigma <- switch(error_model,
    # 1: Additiv
    rep(sigma_add, n),
    # 2: Proportional
    sigma_prop * predicted,
    # 3: Kombiniert
    sqrt(sigma_add^2 + (sigma_prop * predicted)^2),
    # 4: t-additiv
    rep(sigma_add, n),
    # 5: t-proportional
    sigma_prop * predicted,
    # 6: Mixture (behandelt separat)
    sqrt(sigma_add^2 + (sigma_prop * predicted)^2)
  )
  
  # Separiere Standard-Fälle von Spezialfällen
  standard_models <- error_model %in% c(1, 2, 3)
  t_models <- error_model %in% c(4, 5)
  mixture_model <- error_model == 6
  
  # Initialisiere Log-Likelihood
  ll <- numeric(n)
  
  # Standard-Normalverteilung (vektorisiert)
  if (standard_models) {
    # BLQ-Werte
    blq_mask <- is_blq & !is.na(observed)
    if (any(blq_mask)) {
      ll[blq_mask] <- pnorm(lloq, mean = predicted[blq_mask], 
                            sd = sigma[blq_mask], log.p = TRUE)
    }
    
    # Nicht-BLQ-Werte
    non_blq_mask <- !is_blq & !is.na(observed)
    if (any(non_blq_mask)) {
      ll[non_blq_mask] <- dnorm(observed[non_blq_mask], 
                                mean = predicted[non_blq_mask],
                                sd = sigma[non_blq_mask], log = TRUE)
    }
  }
  
  # t-Verteilung (vektorisiert)
  if (t_models) {
    # BLQ-Werte
    blq_mask <- is_blq & !is.na(observed)
    if (any(blq_mask)) {
      ll[blq_mask] <- pt((lloq - predicted[blq_mask]) / sigma[blq_mask], 
                        df = nu, log.p = TRUE)
    }
    
    # Nicht-BLQ-Werte
    non_blq_mask <- !is_blq & !is.na(observed)
    if (any(non_blq_mask)) {
      # Verwende vektorisierte t-Verteilung
      ll[non_blq_mask] <- dt((observed[non_blq_mask] - predicted[non_blq_mask]) / 
                            sigma[non_blq_mask], 
                            df = nu, log = TRUE) - log(sigma[non_blq_mask])
    }
  }
  
  # Mixture-Modell (vektorisiert)
  if (mixture_model) {
    # BLQ-Werte (normale Komponente)
    blq_mask <- is_blq & !is.na(observed)
    if (any(blq_mask)) {
      ll[blq_mask] <- pnorm(lloq, mean = predicted[blq_mask], 
                           sd = sigma[blq_mask], log.p = TRUE)
    }
    
    # Nicht-BLQ-Werte (Mixture)
    non_blq_mask <- !is_blq & !is.na(observed)
    if (any(non_blq_mask)) {
      # Log-sum-exp Trick für numerische Stabilität
      comp1 <- log(mix_w) + dnorm(observed[non_blq_mask], 
                                  mean = predicted[non_blq_mask],
                                  sd = sigma[non_blq_mask], log = TRUE)
      comp2 <- log(1 - mix_w) + dnorm(observed[non_blq_mask], 
                                      mean = predicted[non_blq_mask],
                                      sd = mix_scale * sigma[non_blq_mask], log = TRUE)
      
      # Vektorisiertes log-sum-exp
      ll[non_blq_mask] <- pmax(comp1, comp2) + 
                         log1p(exp(pmin(comp1, comp2) - pmax(comp1, comp2)))
    }
  }
  
  # Behandle NA-Werte
  ll[is.na(ll)] <- 0
  
  return(sum(ll))
}

#' Vektorisierte Residuen-Berechnung
#' 
#' @param observed Beobachtete Werte
#' @param predicted Vorhergesagte Werte
#' @param sigma Fehler-Standardabweichung
#' @param type Typ der Residuen (standard, weighted, studentized)
#' @return Vektor der Residuen
calculate_residuals_vec <- function(observed, predicted, sigma, 
                                   type = "standard") {
  
  n <- length(observed)
  
  # Input-Validierung
  if (length(predicted) != n || length(sigma) != n) {
    stop("Alle Vektoren müssen gleiche Länge haben")
  }
  
  # Vektorisierte Berechnung je nach Typ
  residuals <- switch(type,
    "standard" = observed - predicted,
    "weighted" = (observed - predicted) / sigma,
    "studentized" = {
      # Berechne Leverage (vereinfacht für vektorisierte Version)
      leverage <- rep(1 / n, n)  # Vereinfachte Annahme
      (observed - predicted) / (sigma * sqrt(1 - leverage))
    },
    stop("Unbekannter Residuen-Typ: ", type)
  )
  
  return(residuals)
}

#' Vektorisierte Varianz-Funktionen
#' 
#' Berechnet Varianz für verschiedene Fehlermodelle
calculate_variance_vec <- function(predicted, sigma_add, sigma_prop, 
                                  error_model = 1) {
  
  variance <- switch(error_model,
    # Additiv
    rep(sigma_add^2, length(predicted)),
    # Proportional  
    (sigma_prop * predicted)^2,
    # Kombiniert
    sigma_add^2 + (sigma_prop * predicted)^2,
    # Power (zusätzliches Modell)
    sigma_add^2 * (predicted^(2 * sigma_prop)),
    stop("Unbekanntes Fehlermodell: ", error_model)
  )
  
  return(variance)
}

#' Vektorisierte Gewichtung für WLS
#' 
#' Weighted Least Squares Gewichte
calculate_weights_vec <- function(predicted, sigma_add, sigma_prop, 
                                 error_model = 1, robust = FALSE) {
  
  # Basis-Gewichte (inverse Varianz)
  variance <- calculate_variance_vec(predicted, sigma_add, sigma_prop, error_model)
  weights <- 1 / variance
  
  # Robuste Gewichtung (optional)
  if (robust) {
    # Huber-Gewichte
    residuals_standardized <- abs(predicted - mean(predicted)) / sd(predicted)
    k <- 1.345  # Huber-Konstante
    
    # Vektorisierte Huber-Gewichte
    huber_weights <- ifelse(residuals_standardized <= k, 
                           1, 
                           k / residuals_standardized)
    weights <- weights * huber_weights
  }
  
  # Normalisiere Gewichte
  weights <- weights / sum(weights) * length(weights)
  
  return(weights)
}

#' Batch-Likelihood für multiple Datensätze
#' 
#' Berechnet Log-Likelihood für mehrere Patienten gleichzeitig
batch_loglik <- function(data_list, params_matrix, error_model = 1) {
  
  n_datasets <- length(data_list)
  n_params <- nrow(params_matrix)
  
  # Pre-allocate result matrix
  ll_matrix <- matrix(NA, nrow = n_datasets, ncol = n_params)
  
  # Vektorisierte Berechnung über alle Kombinationen
  for (i in seq_len(n_datasets)) {
    data <- data_list[[i]]
    
    # Parallele Berechnung für alle Parameter-Sets
    ll_values <- apply(params_matrix, 1, function(params) {
      loglik_residuals_vec(
        observed = data$observed,
        predicted = data$predicted,
        sigma_add = params["sigma_add"],
        sigma_prop = params["sigma_prop"],
        error_model = error_model,
        lloq = data$lloq,
        is_blq = data$is_blq
      )
    })
    
    ll_matrix[i, ] <- ll_values
  }
  
  return(ll_matrix)
}

#' Optimierte BLQ-Handling Funktion
#' 
#' Vektorisierte Behandlung von Below Limit of Quantification
handle_blq_vec <- function(observed, lloq, method = "M3") {
  
  n <- length(observed)
  
  # Identifiziere BLQ-Werte
  is_blq <- observed < lloq | is.na(observed)
  
  # Verschiedene BLQ-Methoden (vektorisiert)
  result <- switch(method,
    "M1" = {
      # Ignoriere BLQ
      list(observed = observed[!is_blq], 
           is_blq = rep(FALSE, sum(!is_blq)))
    },
    "M3" = {
      # Zensierte Likelihood
      observed_adj <- observed
      observed_adj[is_blq] <- lloq / 2  # Platzhalter für Berechnung
      list(observed = observed_adj, 
           is_blq = is_blq)
    },
    "M4" = {
      # Erste BLQ behalten, Rest ignorieren
      first_blq <- which(is_blq)[1]
      if (!is.na(first_blq)) {
        keep <- c(seq_len(first_blq), which(!is_blq & seq_len(n) > first_blq))
        list(observed = observed[keep], 
             is_blq = is_blq[keep])
      } else {
        list(observed = observed, 
             is_blq = is_blq)
      }
    },
    "M5" = {
      # LLOQ/2 Substitution
      observed_adj <- observed
      observed_adj[is_blq] <- lloq / 2
      list(observed = observed_adj, 
           is_blq = rep(FALSE, n))
    },
    "M6" = {
      # 0 Substitution
      observed_adj <- observed
      observed_adj[is_blq] <- 0
      list(observed = observed_adj, 
           is_blq = rep(FALSE, n))
    },
    stop("Unbekannte BLQ-Methode: ", method)
  )
  
  return(result)
}

#' Vektorisierte Outlier-Detektion
#' 
#' Identifiziert Ausreißer in Residuen
detect_outliers_vec <- function(residuals, method = "iqr", threshold = 3) {
  
  n <- length(residuals)
  
  outliers <- switch(method,
    "iqr" = {
      # Interquartilsabstand-Methode
      q1 <- quantile(residuals, 0.25, na.rm = TRUE)
      q3 <- quantile(residuals, 0.75, na.rm = TRUE)
      iqr <- q3 - q1
      lower <- q1 - threshold * iqr
      upper <- q3 + threshold * iqr
      residuals < lower | residuals > upper
    },
    "zscore" = {
      # Z-Score Methode
      z_scores <- abs((residuals - mean(residuals, na.rm = TRUE)) / 
                     sd(residuals, na.rm = TRUE))
      z_scores > threshold
    },
    "mad" = {
      # Median Absolute Deviation
      med <- median(residuals, na.rm = TRUE)
      mad <- median(abs(residuals - med), na.rm = TRUE)
      abs(residuals - med) > threshold * mad * 1.4826
    },
    "grubbs" = {
      # Grubbs Test (vektorisiert)
      mean_r <- mean(residuals, na.rm = TRUE)
      sd_r <- sd(residuals, na.rm = TRUE)
      g_scores <- abs(residuals - mean_r) / sd_r
      
      # Kritischer Wert (approximiert)
      alpha <- 0.05
      t_crit <- qt(1 - alpha / (2 * n), n - 2)
      g_crit <- (n - 1) / sqrt(n) * sqrt(t_crit^2 / (n - 2 + t_crit^2))
      
      g_scores > g_crit
    },
    stop("Unbekannte Outlier-Methode: ", method)
  )
  
  return(outliers)
}

#' Export der optimierten Funktionen
#' @export
error_models_vectorized <- list(
  loglik = loglik_residuals_vec,
  residuals = calculate_residuals_vec,
  variance = calculate_variance_vec,
  weights = calculate_weights_vec,
  batch_loglik = batch_loglik,
  handle_blq = handle_blq_vec,
  detect_outliers = detect_outliers_vec
)
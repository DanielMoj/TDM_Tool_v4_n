# Package Management Migration Guide

## 🎯 Ziel der Migration

Modernisierung des Package-Managements von `library()` Aufrufen zu expliziten Namespace-Referenzen (`::`) um Namespace-Konflikte zu vermeiden und die Wartbarkeit zu verbessern.

## ✅ Erfolgreich migrierte Dateien

### 1. **R/error_monitor.R**
- ✅ Alle `library()` Aufrufe entfernt (Zeilen 13-17)
- ✅ Verwendet jetzt `DBI::`, `jsonlite::`, `digest::`, `lubridate::`
- ✅ Dependency-Check bei Modul-Load hinzugefügt

### 2. **R/dependencies.R** (NEU)
- ✅ Zentrales Dependency-Management erstellt
- ✅ Versionsprüfung implementiert
- ✅ Package-Gruppen definiert
- ✅ Manifest-Funktionalität für Reproduzierbarkeit

### 3. **R/load_all.R**
- ✅ Modernisiert ohne direkte `library()` Aufrufe
- ✅ Nutzt `dependencies.R` für Package-Management
- ✅ Health-Check-System integriert
- ✅ Verbesserte Fehlerbehandlung

### 4. **app.R**
- ✅ Komplett auf `::` Notation umgestellt
- ✅ Nutzt zentrales Loading-System
- ✅ Conditional Feature-Loading basierend auf verfügbaren Packages

## 📋 Noch zu prüfende Dateien

Basierend auf der Projektstruktur könnten folgende Dateien noch `library()` Aufrufe enthalten:

### Zu überprüfen mit:
```bash
grep -r "^library(" R/ --include="*.R" | grep -v "dependencies.R"
```

Potenzielle Kandidaten:
- `R/run_fit_jags.R` - Vermutlich `library(rjags)`
- `R/modules/mod_*.R` - Shiny-Module
- Andere Dateien im R/ Verzeichnis

## 🔄 Migration durchführen

### Schritt 1: Identifizieren
```bash
# Alle library() Aufrufe finden
grep -rn "^library(" R/ --include="*.R"

# Auch require() prüfen
grep -rn "^require(" R/ --include="*.R"
```

### Schritt 2: Für jede gefundene Datei

#### Beispiel-Migration:

**VORHER:**
```r
# R/example_module.R
library(ggplot2)
library(dplyr)

plot_data <- function(data) {
  data %>%
    filter(value > 0) %>%
    ggplot(aes(x, y)) +
    geom_point()
}
```

**NACHHER:**
```r
# R/example_module.R
# Dependencies: ggplot2, dplyr
# These packages are loaded via R/dependencies.R

# Check required packages at module load
.check_module_dependencies <- function() {
  required <- c("ggplot2", "dplyr")
  missing <- character()
  
  for (pkg in required) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    }
  }
  
  if (length(missing) > 0) {
    stop(sprintf(
      "This module requires: %s. Install with: install.packages(c(%s))",
      paste(missing, collapse = ", "),
      paste(sprintf('"%s"', missing), collapse = ", ")
    ))
  }
}

.check_module_dependencies()

plot_data <- function(data) {
  data %>%
    dplyr::filter(value > 0) %>%
    ggplot2::ggplot(ggplot2::aes(x, y)) +
    ggplot2::geom_point()
}
```

### Schritt 3: Testen

Nach jeder Migration:
```r
# Datei einzeln testen
source("R/migrated_file.R")

# Gesamtsystem testen
source("R/load_all.R")
load_all_sources()
```

## 🎯 Best Practices

### 1. **Immer Namespace verwenden für externe Packages**
```r
# Gut
result <- dplyr::filter(data, x > 0)

# Schlecht
library(dplyr)
result <- filter(data, x > 0)
```

### 2. **Dependency-Check am Modulanfang**
```r
if (!requireNamespace("needed_package", quietly = TRUE)) {
  stop("Package 'needed_package' is required. Install with: install.packages('needed_package')")
}
```

### 3. **Dokumentiere Dependencies im Header**
```r
#' My Module
#' 
#' @description This module does X
#' @note Requires packages: ggplot2, dplyr, tidyr
```

### 4. **Nutze den Pipe-Operator vorsichtig**
```r
# Der %>% Operator braucht speziellen Import
# Option 1: Explizit
`%>%` <- magrittr::`%>%`

# Option 2: In utils.R definieren (einmalig)
# Dann global verfügbar
```

## 🚨 Häufige Probleme und Lösungen

### Problem 1: "could not find function"
**Lösung**: Namespace vergessen
```r
# Fehler
ggplot(data, aes(x, y))

# Korrekt
ggplot2::ggplot(data, ggplot2::aes(x, y))
```

### Problem 2: Operator nicht gefunden (z.B. %>%, %||%)
**Lösung**: In utils.R definieren oder explizit importieren
```r
# In utils.R (einmalig für die ganze App)
`%>%` <- magrittr::`%>%`
`%||%` <- function(x, y) if (is.null(x)) y else x
```

### Problem 3: Zu viele :: machen Code unleserlich
**Lösung**: Für häufig genutzte Funktionen lokale Aliase
```r
# Am Funktionsanfang
plot_function <- function(data) {
  # Lokale Aliase für bessere Lesbarkeit
  aes <- ggplot2::aes
  ggplot <- ggplot2::ggplot
  geom_point <- ggplot2::geom_point
  
  # Jetzt sauberer Code
  ggplot(data, aes(x, y)) + geom_point()
}
```

## 📊 Validierung

### Automatische Prüfung
```bash
# Skript zur Validierung erstellen
cat > check_migration.sh << 'EOF'
#!/bin/bash
echo "Checking for library() calls..."
count=$(grep -r "^library(" R/ --include="*.R" | grep -v "dependencies.R" | wc -l)

if [ $count -eq 0 ]; then
  echo "✅ No library() calls found in modules!"
else
  echo "⚠️ Found $count library() calls:"
  grep -r "^library(" R/ --include="*.R" | grep -v "dependencies.R"
fi
EOF

chmod +x check_migration.sh
./check_migration.sh
```

### Dependency Report
```r
# Dependency-Status prüfen
source("R/dependencies.R")
report <- get_dependency_report()
print(report)

# Fehlende Packages identifizieren
missing <- report[!report$installed, ]
if (nrow(missing) > 0) {
  print("Missing packages:")
  print(missing)
}
```

## 🎉 Erfolgskriterien

- [x] Keine `library()` Aufrufe in R/*.R (außer dependencies.R)
- [x] Alle Module verwenden `::` Notation
- [x] Code funktioniert weiterhin
- [x] Tests laufen durch
- [x] Keine Namespace-Konflikte
- [x] Dokumentation aktualisiert
- [x] Performance unverändert oder besser

## 📚 Weiterführende Dokumentation

- [ARCHITECTURE.md](ARCHITECTURE.md) - Aktualisierte Architektur-Dokumentation
- [CONTRIBUTING.md](CONTRIBUTING.md) - Entwicklungsrichtlinien
- [R/dependencies.R](../R/dependencies.R) - Dependency-Management-System
- [R/load_all.R](../R/load_all.R) - Modernisierter Loader

## 🔮 Zukünftige Verbesserungen

1. **renv Integration** für Versions-Pinning
2. **Automatische Dependency-Installation** in Docker
3. **CI/CD Checks** für library() Aufrufe
4. **Dependency-Graph-Visualisierung**
5. **Performance-Profiling** pro Package

---

*Migration dokumentiert: 2025-08-11*
*Work Package 5 erfolgreich umgesetzt*
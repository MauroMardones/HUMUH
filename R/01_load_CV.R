# =============================================================================
# 01_load_CV.R
# Lectura y recorte espacial de datos de Cajas Verdes 2017-2024
# Golfo de Cádiz - Demarcación Sudatlántica
# =============================================================================

library(data.table)
library(tidyverse)
library(fst)

# --- Crear carpetas si no existen --------------------------------------------
dir.create("outputs",      showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/shp",  showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/CSV",  showWarnings = FALSE, recursive = TRUE)
dir.create("figs",         showWarnings = FALSE, recursive = TRUE)
# --- Zona de estudio (bounding box) ------------------------------------------
lat_min <- 34; lat_max <- 38
lon_min <- -8; lon_max <- -5

# --- Rutas base --------------------------------------------------------------
path_001 <- "data/Junta Andalucía/Cajas verdes/CAJAS_NUEVAS_2017-102025_001"
path_002 <- "data/Junta Andalucía/Cajas verdes/CAJAS_NUEVAS_2017-102025_002/Resultados_EQ_SLFP_PET314554"

# Años con dos particiones en _001 (2017-2021) y en _002 (2023-2024)
# 2022 solo tiene partición _1 en _001
years_001 <- 2017:2022
years_002 <- 2023:2024

# --- Función de lectura con filtro espacial inmediato -----------------------
read_cv <- function(path, year) {

  f1 <- file.path(path, paste0("Resultado_2_CAJAS_NUEVAS_", year, ".csv"))
  f2 <- file.path(path, paste0("Resultado_2_CAJAS_NUEVAS_", year, "_2.csv"))

  # Leer partición 1
  dt1 <- fread(f1,
               sep = ";",
               dec = ".",
               encoding = "UTF-8",
               showProgress = FALSE)

  # Leer partición 2 si existe
  if (file.exists(f2)) {
    dt2 <- fread(f2,
                 sep = ";",
                 dec = ".",
                 encoding = "UTF-8",
                 showProgress = FALSE)
    dt <- rbindlist(list(dt1, dt2), use.names = TRUE)
    rm(dt1, dt2)
  } else {
    dt <- dt1
    rm(dt1)
  }

  # Renombrar columna con caracteres especiales
  setnames(dt, "\xbf\xa1EN PUERTO?",  "EN_PUERTO", skip_absent = TRUE)
  setnames(dt, "\xef\xbb\xbfFECHA",   "FECHA",     skip_absent = TRUE)
  # Alternativa si el BOM llega diferente:
  old <- grep("FECHA|fecha", names(dt), value = TRUE)[1]
  if (!is.na(old) && old != "FECHA") setnames(dt, old, "FECHA")

  # Filtro espacial inmediato para reducir memoria
  dt <- dt[LATITUD >= lat_min & LATITUD <= lat_max &
           LONGITUD >= lon_min & LONGITUD <= lon_max]

  # Añadir columna de año
  dt[, ANYO := year]

  gc()
  message(sprintf("  %d: %s filas tras filtro ZE", year, format(nrow(dt), big.mark = ".")))
  return(dt)
}

# --- Loop 2017-2024 ----------------------------------------------------------
message("Leyendo Cajas Verdes 2017-2024...")

lista_cv <- vector("list", length(2017:2024))
names(lista_cv) <- as.character(2017:2024)

for (yr in years_001) {
  lista_cv[[as.character(yr)]] <- read_cv(path_001, yr)
}

for (yr in years_002) {
  lista_cv[[as.character(yr)]] <- read_cv(path_002, yr)
}

# --- Unir todos los años -----------------------------------------------------
CV_all <- rbindlist(lista_cv, use.names = TRUE, fill = TRUE)
#elimina la lista de data.tables por año de la memoria RAM — ya no 
# la necesitas porque todo está unido en CV_all.
# fuerza al garbage collector de R a liberar esa memoria de inmediato,
# en lugar de esperar a que R lo haga solo. 
# Con 67M filas es importante hacerlo explícito 
#para no quedarte sin RAM antes de continuar con el formateo.
rm(lista_cv); gc()

message(sprintf("\nTotal filas CV 2017-2024: %s", format(nrow(CV_all), big.mark = ".")))
message(sprintf("Columnas: %s", paste(names(CV_all), collapse = ", ")))

# --- Guardar como .fst para relectura rápida ---------------------------------

write_fst(CV_all, "outputs/CV_2017_2024_ZE.fst", compress = 75)

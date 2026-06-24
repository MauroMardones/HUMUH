# =============================================================================
# 02_format_CV.R
# Formateo de datos CV 2017-2024
# - Código CFR desde CFPO
# - Timestamp
# - Velocidad recalculada (Haversine)
# - Rumbo COG
# - Columnas auxiliares estandarizadas
# =============================================================================

library(data.table)
library(geosphere)
library(lubridate)

# CV_all debe estar cargado desde 01_load_CV.R
# o bien leer desde .fst: CV_all <- fst::read_fst("outputs/CV_2017_2024_ZE.fst", as.data.table = TRUE)

# --- 1. Renombrar columna con caracteres especiales -------------------------
en_puerto_col <- grep("PUERTO|puerto", names(CV_all), value = TRUE)
if (length(en_puerto_col) > 0) setnames(CV_all, en_puerto_col, "EN_PUERTO")

# --- 2. Eliminar duplicados --------------------------------------------------
message("Eliminando duplicados...")
n_antes <- nrow(CV_all)
CV_all <- unique(CV_all)
message(sprintf("  Filas eliminadas: %s", format(n_antes - nrow(CV_all), big.mark = ".")))

# --- 3. Código CFR desde CFPO -----------------------------------------------
# Formato: "ESP" + zeros hasta 9 dígitos + CFPO → 12 caracteres
message("Creando VE_REF (CFR)...")
CV_all[, VE_REF := paste0("ESP", formatC(CFPO, width = 9, flag = "0"))]
CV_all[, CFPO := NULL]

# --- 4. Timestamp ------------------------------------------------------------
message("Creando SI_TIMESTAMP...")
# fast.time no maneja timezones; usamos as.POSIXct con format explícito
# Es más rápido que lubridate::dmy_hms sobre 67M filas
CV_all[, SI_TIMESTAMP := as.POSIXct(
  paste(FECHA, HORA),
  format = "%d/%m/%Y %H:%M:%S",
  tz = "Europe/Madrid"
)]
CV_all[, c("FECHA", "HORA") := NULL]

# Ordenar por barco y timestamp (necesario para lag correcto)
message("Ordenando por VE_REF y timestamp...")
setorder(CV_all, VE_REF, SI_TIMESTAMP)

# --- 5. Velocidad Haversine y TDIFF -----------------------------------------
# Agrupa por barco + día para no conectar registros entre días distintos
message("Calculando SI_SPCA y SI_TDIFF...")
CV_all[, fecha_dia := as.Date(SI_TIMESTAMP, tz = "Europe/Madrid")]

CV_all[, `:=`(
  SI_DISTANCECA = if (.N > 1) c(NA_real_,
    distHaversine(
      cbind(LONGITUD[-.N], LATITUD[-.N]),
      cbind(LONGITUD[-1],  LATITUD[-1])
    )
  ) else NA_real_,
  SI_TDIFF = if (.N > 1) c(NA_real_, as.numeric(diff(SI_TIMESTAMP), units = "secs")) else NA_real_
), by = .(VE_REF, fecha_dia)]

# Velocidad en nudos (m/s → nudos: × 1.94384)
CV_all[, SI_SPCA := (SI_DISTANCECA / SI_TDIFF) * 1.94384]

# --- 6. Rumbo COG -----------------------------------------------------------
message("Calculando SI_COG...")
CV_all[, SI_COG := if (.N > 1) {
  b <- bearing(
    cbind(LONGITUD[-.N], LATITUD[-.N]),
    cbind(LONGITUD[-1],  LATITUD[-1])
  )
  c(NA_real_, ifelse(b < 0, b + 360, b))
} else NA_real_, by = .(VE_REF, fecha_dia)]

CV_all[, fecha_dia := NULL]

# --- 7. Código de marea FT_REF ----------------------------------------------
CV_all[, FT_REF := paste0(
  VE_REF, "_",
  format(as.Date(SI_TIMESTAMP, tz = "Europe/Madrid"), "%Y_%m_%d")
)]

# --- 8. Eliminar NAs (primera fila de cada grupo) ---------------------------
CV_all <- na.omit(CV_all, cols = c("SI_SPCA", "SI_COG"))

message(sprintf("Filas tras na.omit: %s", format(nrow(CV_all), big.mark = ".")))

# --- 9. Renombrar columnas al estándar --------------------------------------
setnames(CV_all,
  old = c("LATITUD", "LONGITUD", "VELOCIDAD", "RUMBO"),
  new = c("SI_LATI", "SI_LONG", "SI_SP",     "SI_HE")
)

# --- 10. Columnas auxiliares vacías -----------------------------------------
CV_all[, `:=`(
  LE_MET4   = NA_character_,
  LE_MET6   = NA_character_,
  SI_FSTATUS = NA,
  SI_FOPER   = NA_character_,
  SU_ISOB   = NA,
  SI_OGT    = NA
)]

# --- 11. Selección y orden de columnas --------------------------------------
setcolorder(CV_all, c(
  "VE_REF", "FT_REF", "ANYO", "SI_TIMESTAMP",
  "SI_LATI", "SI_LONG", "SI_SP", "SI_SPCA", "SI_HE", "SI_COG",
  "SI_DISTANCECA", "SI_TDIFF",
  "EN_PUERTO",
  "LE_MET4", "LE_MET6", "SI_FSTATUS", "SI_FOPER", "SU_ISOB", "SI_OGT"
))

message("Formato completado.")
message(sprintf("Dimensiones finales: %s filas x %s columnas",
  format(nrow(CV_all), big.mark = "."), ncol(CV_all)))

# --- Guardar -----------------------------------------------------------------
fst::write_fst(CV_all, "outputs/CV_2017_2024_ZE_fmt.fst", compress = 75)
saveRDS(CV_all, "outputs/CV_2017_2024_ZE_fmt.rds")

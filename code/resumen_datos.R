# =============================================================================
# resumen_datos.R
# Resumen de los datos para la slide "Data at a glance" (WGSFD_Turku_2026.Rmd)
#   - Años analizados, n° barcos, n° registros, n° flotas (LE_MET4)
#   - Cobertura espacial (celdas 0.05° con esfuerzo)
#   - Pérdida de datos en el cruce cajas verdes × notas de venta
# Salida: outputs/resumen_datos.csv  (lo lee la presentación al hacer Knit)
# Correr desde la raíz del proyecto (HUMUH.Rproj)
# =============================================================================

library(data.table)
library(fst)
library(marmap)

res      <- 0.05          # misma resolución que 04_effort_grid.R
anios_nv <- 2018:2024     # años con notas de venta
leer_NV  <- TRUE          # TRUE = separa la pérdida en "sin nota de venta" y
                          # "arte excluido" (lee los xlsx de NV, tarda unos min)
# Caja para la cobertura sobre plataforma (0-200 m), lado español del golfo
caja <- c(xmin = -7.5, xmax = -5.0, ymin = 36.0, ymax = 37.3)

fmt_n <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE)
fmt_M <- function(x) paste0(format(round(x / 1e6, 1), decimal.mark = ","), " M")
cel   <- function(x) floor(x / res) * res + res / 2

# --- 1. Cajas verdes formateadas (antes del cruce) ---------------------------
cv <- read_fst("outputs/CV_2017_2024_ZE_fmt.fst",
               columns = c("VE_REF", "FT_REF", "ANYO", "SI_LATI", "SI_LONG"),
               as.data.table = TRUE)
n_ping_cv_all <- nrow(cv)
n_barc_cv_all <- uniqueN(cv$VE_REF)
anios_cv      <- range(cv$ANYO, na.rm = TRUE)

cv <- cv[ANYO %in% anios_nv]   # comparación justa: solo años con NV
n_ping_cv  <- nrow(cv)
n_trip_cv  <- uniqueN(cv$FT_REF)
n_barc_cv  <- uniqueN(cv$VE_REF)

# Huella de la flota: celdas con al menos un ping
huella <- unique(cv[, .(lon_cel = cel(SI_LONG), lat_cel = cel(SI_LATI))])

# --- 2. Cajas verdes × notas de venta (salida de 03) -------------------------
cvnv <- read_fst("outputs/CV_NV_2018_2024_ZE.fst",
                 columns = c("VE_REF", "FT_REF", "ANYO", "LE_MET4"),
                 as.data.table = TRUE)
n_ping_nv <- nrow(cvnv)
n_trip_nv <- uniqueN(cvnv$FT_REF)
n_barc_nv <- uniqueN(cvnv$VE_REF)
flotas    <- sort(unique(cvnv$LE_MET4))
anios_an  <- range(cvnv$ANYO, na.rm = TRUE)

# --- 3. Pérdida en el cruce ---------------------------------------------------
perd_ping <- 100 * (1 - n_ping_nv / n_ping_cv)
perd_trip <- 100 * (1 - n_trip_nv / n_trip_cv)
perd_barc <- 100 * (1 - n_barc_nv / n_barc_cv)

perd_sin_nv <- perd_excl <- NA_real_
if (leer_NV) {
  library(readxl)
  path_nv <- "data/Junta Andalucía/Notas de Venta_2018_2025/20251001 IEO NV 2018-2025"
  ft_nv <- unique(rbindlist(lapply(anios_nv, function(yr) {
    f  <- list.files(path_nv, pattern = paste0("NV ", yr, "\\.xlsx$"), full.names = TRUE)
    dt <- as.data.table(read_xlsx(f))[, .(ID_UDPROD, FECHAVENTA)]
    dt[, FT_REF := paste0(ID_UDPROD, "_", format(as.Date(FECHAVENTA), "%Y_%m_%d"))]
    dt[, .(FT_REF)]
  }))$FT_REF)
  con_nv      <- cv$FT_REF %in% ft_nv
  perd_sin_nv <- 100 * (1 - mean(con_nv))                     # pings sin nota de venta
  perd_excl   <- 100 * (sum(con_nv) - n_ping_nv) / n_ping_cv  # arrastre, cerco, sin código
}

# --- 4. Cobertura espacial ----------------------------------------------------
esf <- read_fst("outputs/esfuerzo_2018_2024.fst", as.data.table = TRUE)
cel_esf <- unique(esf[LE_MET4 != "TOTAL" & horas_pesca > 0, .(lon_cel, lat_cel)])
cel_esf[, area_km2 := (res * 111.32) * (res * 111.32 * cos(lat_cel * pi / 180))]
area_esf <- sum(cel_esf$area_km2)
pct_huella <- 100 * nrow(cel_esf) / nrow(huella)

# Plataforma 0-200 m dentro de la caja (celdas 0.05° con profundidad media en rango)
bat <- getNOAA.bathy(lon1 = unname(caja["xmin"]), lon2 = unname(caja["xmax"]),
                     lat1 = unname(caja["ymin"]), lat2 = unname(caja["ymax"]), resolution = 1)
bat <- as.data.table(fortify.bathy(bat))
bat <- bat[, .(z = mean(z)), by = .(lon_cel = cel(x), lat_cel = cel(y))]
plataforma <- bat[z < 0 & z >= -200, .(lon_cel, lat_cel)]
en_caja <- cel_esf[lon_cel >= caja["xmin"] & lon_cel <= caja["xmax"] &
                   lat_cel >= caja["ymin"] & lat_cel <= caja["ymax"]]
pct_plat <- 100 * nrow(merge(en_caja, plataforma, by = c("lon_cel", "lat_cel"))) /
            nrow(plataforma)

# --- 4b. Tamaño de la base de datos -----------------------------------------
fmt_GB <- function(b) paste0(format(round(b / 1024^3, 1), decimal.mark = ","), " GB")
csv_cv <- list.files("data/Junta Andalucía/Cajas verdes", recursive = TRUE,
                     full.names = TRUE,
                     pattern = paste0("CAJAS_NUEVAS_(", paste(anios_cv[1]:anios_cv[2], collapse = "|"),
                                      ")(_2)?\\.csv$"))
gb_csv  <- sum(file.size(csv_cv), na.rm = TRUE)
gb_fmt  <- file.size("outputs/CV_2017_2024_ZE_fmt.fst")
gb_cvnv <- file.size("outputs/CV_NV_2018_2024_ZE.fst")
message(sprintf("Raw CSV: %s (%d archivos) · fmt.fst: %s · CV_NV.fst: %s",
                fmt_GB(gb_csv), length(csv_cv), fmt_GB(gb_fmt), fmt_GB(gb_cvnv)))

# --- 5. Tabla resumen ---------------------------------------------------------
resumen <- data.table(
  id    = c("years", "vessels", "records", "fleets", "coverage", "loss"),
  label = c("Years analysed", "Vessels", "Records (pings)",
            "Fleets (LE_MET4)", "Spatial coverage", "Loss in sales-note join"),
  value = c(
    sprintf("%d–%d", anios_an[1], anios_an[2]),
    fmt_n(n_barc_nv),
    fmt_M(n_ping_nv),
    as.character(length(flotas)),
    sprintf("%.0f%%", pct_plat),
    sprintf("%.0f%%", perd_ping)
  ),
  detail = c(
    sprintf("%d years with sales notes (green boxes %d–%d)",
            diff(anios_an) + 1, anios_cv[1], anios_cv[2]),
    sprintf("of %s vessels with green box", fmt_n(n_barc_cv_all)),
    sprintf("of %s pings in study area · raw CSV %s → .fst %s → linked %s",
            fmt_M(n_ping_cv_all), fmt_GB(gb_csv), fmt_GB(gb_fmt), fmt_GB(gb_cvnv)),
    paste(flotas, collapse = " · "),
    sprintf("of 0–200 m shelf cells (0.05°) · %s km² fished · %.0f%% of fleet footprint",
            fmt_n(round(area_esf)), pct_huella),
    if (leer_NV)
      sprintf("of pings: %.0f%% no sales note + %.0f%% excluded gears · trips %.0f%% · vessels %.0f%%",
              perd_sin_nv, perd_excl, perd_trip, perd_barc)
    else
      sprintf("of pings %d–%d · trips %.0f%% · vessels %.0f%%",
              min(anios_nv), max(anios_nv), perd_trip, perd_barc)
  )
)

print(resumen[, .(label, value, detail)])
fwrite(resumen, "outputs/resumen_datos.csv")
message("Guardado: outputs/resumen_datos.csv")

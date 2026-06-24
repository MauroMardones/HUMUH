# =============================================================================
# 09_export_shp.R
# Exporta resultados del análisis como GeoPackage (.gpkg) y Shapefile (.shp)
# Capas: esfuerzo acumulado, core fishing ground, tendencia GLM
# =============================================================================
#
# OBJETOS DE ENTRADA
# ------------------
# Todos se leen desde outputs/ (no necesitan estar en memoria)
#
# Si retomas desde cero:
#   library(data.table); library(fst)
#   esfuerzo  <- fst::read_fst("outputs/esfuerzo_2018_2024.fst",        as.data.table = TRUE)
#   cfg       <- fst::read_fst("outputs/core_fishing_ground_2018_2024.fst", as.data.table = TRUE)
#   tendencia <- fst::read_fst("outputs/tendencia_glm_2018_2024.fst",   as.data.table = TRUE)
#
# =============================================================================

library(data.table)
library(fst)
library(sf)

dir.create("outputs/shp", showWarnings = FALSE)

artes_obj <- c("DRB", "HMD", "FPO", "GN", "GNS", "GTN", "GTR", "LL", "LLS", "LHP")

artes_label <- c(
  DRB = "DRB (Rastra)",
  HMD = "HMD (Draga hidráulica)",
  FPO = "FPO (Nasas/Trampas)",
  GN  = "GN (Enmalle s/e)",
  GNS = "GNS (Enmalle fondo)",
  GTN = "GTN (Enmalle-trasmallo)",
  GTR = "GTR (Trasmallo)",
  LL  = "LL (Palangre s/e)",
  LLS = "LLS (Palangre calado)",
  LHP = "LHP (Línea de mano)"
)

# --- Leer objetos ------------------------------------------------------------
message("Leyendo outputs...")
esfuerzo  <- fst::read_fst("outputs/esfuerzo_2018_2024.fst",
                             as.data.table = TRUE)
cfg       <- fst::read_fst("outputs/core_fishing_ground_2018_2024.fst",
                             as.data.table = TRUE)
tendencia <- fst::read_fst("outputs/tendencia_glm_2018_2024.fst",
                             as.data.table = TRUE)

# =============================================================================
# FUNCIÓN: celda 0.01° → polígono sf
# Recibe un data.table con lon_cel y lat_cel (centroides)
# Devuelve un sf con geometría de polígonos cuadrados de 0.01°
# =============================================================================
celdas_a_sf <- function(dt) {
  res <- 0.05 / 2  # semiancho de la celda (grilla 0.05°)
  polys <- lapply(seq_len(nrow(dt)), function(i) {
    x <- dt$lon_cel[i]
    y <- dt$lat_cel[i]
    st_polygon(list(matrix(c(
      x - res, y - res,
      x + res, y - res,
      x + res, y + res,
      x - res, y + res,
      x - res, y - res
    ), ncol = 2, byrow = TRUE)))
  })
  st_sf(dt, geometry = st_sfc(polys, crs = 4326))
}

# =============================================================================
# 1. ESFUERZO ACUMULADO 2018-2024
# Una capa por arte en el gpkg + shapefile individual
# =============================================================================
message("Exportando esfuerzo acumulado...")

acumulado <- esfuerzo[LE_MET4 %in% artes_obj,
  .(horas_pesca = sum(horas_pesca, na.rm = TRUE),
    n_barcos    = sum(n_barcos,    na.rm = TRUE),
    n_anos      = uniqueN(ANYO)),
  by = .(lon_cel, lat_cel, csq, LE_MET4)]
acumulado <- acumulado[horas_pesca > 0]
acumulado[, gear_name := artes_label[LE_MET4]]

gpkg_esf <- "outputs/shp/esfuerzo_acumulado_GoC_2018_2024.gpkg"
if (file.exists(gpkg_esf)) file.remove(gpkg_esf)

for (art in artes_obj) {
  dat <- acumulado[LE_MET4 == art]
  if (nrow(dat) == 0) next
  sf_layer <- celdas_a_sf(dat)
  st_write(sf_layer, dsn = gpkg_esf,
           layer      = paste0("effort_", art),
           driver     = "GPKG",
           append     = TRUE,
           quiet      = TRUE)
  message(sprintf("  Esfuerzo %s: %s celdas", art, format(nrow(dat), big.mark = ".")))
}
message(sprintf("Guardado: %s", gpkg_esf))


# =============================================================================
# 2. CORE FISHING GROUND (acumulado 2018-2024)
# =============================================================================
message("Exportando core fishing ground...")

cfg_export <- cfg[LE_MET4 %in% artes_obj, .(
  lon_cel, lat_cel, csq, LE_MET4,
  horas_pesca,
  core_75  = as.integer(core_75),
  core_90  = as.integer(core_90),
  categoria = fcase(
    core_75 == 1, "CFG 75%",
    core_90 == 1, "CFG 90%",
    default      = "Periphery"
  ),
  gear_name = artes_label[LE_MET4]
)]

gpkg_cfg <- "outputs/shp/core_fishing_ground_GoC_2018_2024.gpkg"
if (file.exists(gpkg_cfg)) file.remove(gpkg_cfg)

for (art in artes_obj) {
  dat <- cfg_export[LE_MET4 == art]
  if (nrow(dat) == 0) next
  sf_layer <- celdas_a_sf(dat)
  st_write(sf_layer, dsn = gpkg_cfg,
           layer  = paste0("cfg_", art),
           driver = "GPKG",
           append = TRUE,
           quiet  = TRUE)
  message(sprintf("  CFG %s: %s celdas", art, format(nrow(dat), big.mark = ".")))
}
message(sprintf("Guardado: %s", gpkg_cfg))


# =============================================================================
# 3. TENDENCIA GLM POR CELDA
# =============================================================================
message("Exportando tendencia GLM...")

tend_export <- tendencia[LE_MET4 %in% artes_obj, .(
  lon_cel, lat_cel, csq, LE_MET4,
  slope     = round(slope, 4),
  pval      = round(pval,  4),
  r2        = round(r2,    3),
  n_anos,
  sig       = as.integer(sig),
  direccion,
  gear_name = artes_label[LE_MET4]
)]

gpkg_glm <- "outputs/shp/tendencia_glm_GoC_2018_2024.gpkg"
if (file.exists(gpkg_glm)) file.remove(gpkg_glm)

for (art in artes_obj) {
  dat <- tend_export[LE_MET4 == art]
  if (nrow(dat) == 0) next
  sf_layer <- celdas_a_sf(dat)
  st_write(sf_layer, dsn = gpkg_glm,
           layer  = paste0("trend_", art),
           driver = "GPKG",
           append = TRUE,
           quiet  = TRUE)
  message(sprintf("  GLM %s: %s celdas (%s sig.)",
    art,
    format(nrow(dat), big.mark = "."),
    format(sum(dat$sig), big.mark = ".")))
}
message(sprintf("Guardado: %s", gpkg_glm))


# =============================================================================
# 4. SHAPEFILE INDIVIDUAL — esfuerzo total (todos los artes juntos)
# Útil para visualización rápida en QGIS sin abrir capas individuales
# =============================================================================
message("Exportando shapefile esfuerzo total (todos los artes)...")

esf_total <- esfuerzo[LE_MET4 == "TOTAL" & horas_pesca > 0,
  .(horas_pesca = sum(horas_pesca, na.rm = TRUE),
    n_barcos    = sum(n_barcos,    na.rm = TRUE)),
  by = .(lon_cel, lat_cel)]

sf_total <- celdas_a_sf(esf_total)
st_write(sf_total,
         dsn    = "outputs/shp/esfuerzo_TOTAL_GoC_2018_2024.shp",
         driver = "ESRI Shapefile",
         delete_dsn = TRUE,
         quiet  = TRUE)
message(sprintf("Guardado: esfuerzo_TOTAL_GoC_2018_2024.shp (%s celdas)",
  format(nrow(esf_total), big.mark = ".")))


message("\n=== Script 09 completado ===")
message("Archivos en outputs/shp/:")
message("  esfuerzo_acumulado_GoC_2018_2024.gpkg  — 10 capas (una por arte)")
message("  core_fishing_ground_GoC_2018_2024.gpkg — 10 capas (una por arte)")
message("  tendencia_glm_GoC_2018_2024.gpkg       — 10 capas (una por arte)")
message("  esfuerzo_TOTAL_GoC_2018_2024.shp       — todos los artes combinados")

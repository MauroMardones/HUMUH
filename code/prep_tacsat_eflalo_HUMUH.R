# =============================================================================
# prep_tacsat_eflalo_HUMUH.R
#
# Puente entre HUMUH y el flujo ICES VMS & Logbook Data Call.
# Convierte CV_NV (cajas verdes + notas de venta) en:
#   code/ICES_datacall/Data/tacsat_<año>.RData  (objeto tacsat_<año>)
#   code/ICES_datacall/Data/eflalo_<año>.RData  (objeto eflalo_<año>)
#
# Los scripts de ICES (code/ICES_datacall/*.R) NO se modifican.
#
# Uso (desde la raíz del proyecto HUMUH):
#   source("code/prep_tacsat_eflalo_HUMUH.R")
#   setwd("code/ICES_datacall")
#   source("0_global.R")
#   yearsToSubmit <- 2024   # en consola, sin tocar 0_global.R
#   intThres      <- 1      # en consola: los pings CV son cada ~3 min y el
#                           # valor por defecto (5 min) los borraría casi todos
#                           # como pseudo-duplicados en el script 1
#   source("1_eflalo_tacsat_preprocessing.R")
#   source("2_eflalo_tacsat_analysis.R")
#   source("3_data_submission.R")
#   setwd("../..")
#
# Formato de referencia: code/ICES_datacall/b_EFLALO & TACSAT Formats.md
# =============================================================================

library(data.table)
library(readxl)

# Rectángulo estadístico ICES (30' lat x 1° lon), sin depender de vmstools.
# Filas: 01 = 36°00'-36°30'N. Columnas: A0-A3 (44-40°W), B0..M9 (sin "I").
# Ej.: 54°N 2°E -> "37F2"; Golfo de Cádiz -> 01E3, 02E2, 03E4...
rect_ices <- function(lon, lat) { 
  fila  <- floor((lat - 36) * 2) + 1
  letra <- ifelse(lon < -40, "A",
                  c("B","C","D","E","F","G","H","J","K","L","M")[floor((lon + 40) / 10) + 1])
  num   <- ifelse(lon < -40, floor(lon + 44), floor(lon + 40) %% 10)
  r <- sprintf("%02d%s%d", fila, letra, num)
  r[fila < 1 | fila > 99 | lon < -44 | lon >= 70 | is.na(lon) | is.na(lat)] <- NA_character_
  r
}

# --- 0. Parámetros -----------------------------------------------------------
anios       <- 2024
ruta_cvnv   <- "outputs/CV_NV_2018_2024_ZE.fst"
ruta_nv     <- "data/Junta Andalucía/Notas de Venta_2018_2025/20251001 IEO NV 2018-2025"
ruta_censo  <- "data/censo_flota.csv"   # opcional: VE_REF, VE_LEN, VE_KW, VE_TON
ruta_metier <- "code/ICES_datacall/resources/RDB_ISSG_Metier_list.csv"
ruta_grupos <- "code/ICES_ref/Metier_Subgroup_Species.csv"
ruta_salida <- "code/ICES_datacall/Data"
rcg         <- "NAtl"                   # Golfo de Cádiz = RCG North Atlantic

quitar_pings_puerto      <- TRUE  # usa EN_PUERTO de CV (mejor que los puertos de vmstools)
unificar_sin_especificar <- TRUE  # GN -> GNS, LL -> LLS (GN/LL no tienen métier N6 en NAtl)

dir.create(ruta_salida, showWarnings = FALSE, recursive = TRUE)

# --- 1. Censo de flota (opcional) --------------------------------------------
# Sin censo los scripts ICES corren igual (eflalo_clean() de 0_global.R, que haría
# stop() por VE_LEN/VE_TON, no se llama en los scripts 1-3), pero quedan en NA:
#   - kwHour / kwDays (esfuerzo en kW)       -> script 2 l.644, script 3 l.45
#   - VesselLengthRange (LENGTHCAT)          -> script 3 l.156-165
#   - GEARWIDTHKM / área barrida (usa VE_LEN y VE_KW) -> script 2 l.660
# Las horas de pesca (INTV), capturas, c-squares y métiers sí se calculan.
# Fuente sugerida: EU Fleet Register (CFR = VE_REF, LOA, potencia principal, GT).
req_censo <- c("VE_REF", "VE_LEN", "VE_KW", "VE_TON")
if (file.exists(ruta_censo)) {
  censo <- fread(ruta_censo)
  faltan <- setdiff(req_censo, names(censo))
  if (length(faltan)) stop("Columnas faltantes en el censo: ", paste(faltan, collapse = ", "))
  censo <- unique(censo[, ..req_censo], by = "VE_REF")
} else {
  warning("Sin censo de flota (", ruta_censo, "): VE_LEN, VE_KW y VE_TON quedan en NA. ",
          "Horas de pesca OK; kW-horas, clase de eslora y área barrida saldrán en NA.")
  censo <- data.table(VE_REF = character(), VE_LEN = numeric(),
                      VE_KW = numeric(), VE_TON = numeric())
}

# --- 2. Tablas de referencia ICES / RCG --------------------------------------
# Métiers N6 válidos para la región
met_validos <- fread(ruta_metier)[RCG == rcg, unique(Metier_level6)]

# Especie FAO -> grupo de especies (MOL, CEP, CRU, DEF, SPF, LPF...)
# Fuente: ices-eg/RCGs/Metiers/Reference_lists/Metier Subgroup Species.xlsx
grupos <- fread(ruta_grupos)[, .(FAO, GROUP)]

# Especies de las NV de Andalucía que no están en la lista RCG
# (solo las de grupo inequívoco; el resto queda como "UNK" y se reporta)
grupos_extra <- data.table(
  FAO   = c("DXL", "GFU", "JDV", "EKW", "KDH", "SVV", "MPO", "RLI", "DHQ"),
  GROUP = c("MOL", "CRU", "CRU", "CRU", "MOL", "DEF", "DEF", "DEF", "DEF")
  # DXL coquina, GFU camarón periquito, JDV centolla, EKW camarón rayado,
  # KDH casco tirreno, SVV rascacio de altura, MPO peje obispo,
  # RLI peluda imperial, DHQ sargo
)
grupos <- unique(rbind(grupos, grupos_extra[!FAO %in% grupos$FAO]), by = "FAO")

# Segmento de flota DCF (VE_FLT) según arte
seg_dcf <- function(g) fcase(
  g %in% c("DRB", "HMD"),                     "DRB",
  g %in% c("FPO", "FIX"),                     "FPO",
  g %in% c("LHP", "LLS", "LLD", "LL"),        "HOK",
  g %in% c("GNS", "GTR", "GTN", "GN"),        "DFN",
  default = "PGP"
)

# --- 3. Cajas verdes + NV (ya unidas) ----------------------------------------
cols_cv <- c("FT_REF", "VE_REF", "ANYO", "SI_TIMESTAMP", "SI_LATI", "SI_LONG",
             "SI_SP", "SI_SPCA", "SI_HE", "EN_PUERTO", "LE_MET4", "FECHAVENTA")
CV_all <- fst::read_fst(ruta_cvnv, columns = cols_cv, as.data.table = TRUE)[ANYO %in% anios]

if (unificar_sin_especificar) {
  CV_all[, LE_MET4 := fcase(LE_MET4 == "GN", "GNS",
                            LE_MET4 == "LL", "LLS",
                            default = LE_MET4)]
}

# Velocidad: instantánea; si falta, la calculada
CV_all[, SI_SP := fifelse(is.na(SI_SP), SI_SPCA, SI_SP)]

# ICES pide fecha y hora en UTC. SI_TIMESTAMP se creó en Europe/Madrid
# (R/02_format_CV.R); format(..., tz = "UTC") convierte el instante correctamente.
CV_all[, `:=`(
  SI_DATE = format(SI_TIMESTAMP, "%d/%m/%Y", tz = "UTC"),
  SI_TIME = format(SI_TIMESTAMP, "%H:%M",    tz = "UTC")
)]

# --- 4. Bucle por año --------------------------------------------------------
for (yr in anios) {

  message("\n==== Año ", yr, " ====")
  CV <- CV_all[ANYO == yr]

  # 4.1 TACSAT ----------------------------------------------------------------
  tacsat <- CV[, .(VE_COU = "ESP", VE_REF, SI_LATI, SI_LONG,
                   SI_DATE, SI_TIME, SI_SP, SI_HE, EN_PUERTO)]
  if (quitar_pings_puerto) tacsat <- tacsat[EN_PUERTO == 0]
  tacsat[, EN_PUERTO := NULL]
  message(sprintf("tacsat: %s pings", format(nrow(tacsat), big.mark = ",")))

  # 4.2 Viajes (FT_REF = barco + día de venta) --------------------------------
  viajes <- CV[, .(
    VE_REF     = VE_REF[1],
    LE_GEAR    = LE_MET4[1],          # arte dominante (kg) de la NV, ver R/03_join_CV_NV.R
    FECHAVENTA = FECHAVENTA[1],
    t_ini      = min(SI_TIMESTAMP),   # salida  = primer ping
    t_fin      = max(SI_TIMESTAMP)    # llegada = último ping
  ), by = FT_REF]

  # Rectángulo ICES más frecuente entre los pings fuera de puerto
  pos <- CV[EN_PUERTO == 0, .(FT_REF, LE_RECT = rect_ices(SI_LONG, SI_LATI))]
  rect <- pos[!is.na(LE_RECT), .N, by = .(FT_REF, LE_RECT)][order(FT_REF, -N)][
    , .(LE_RECT = LE_RECT[1]), by = FT_REF]
  rm(pos)

  # 4.3 Captura por especie desde las NV originales ----------------------------
  f_nv <- list.files(ruta_nv, pattern = paste0("NV ", yr, "\\.xlsx$"), full.names = TRUE)
  if (length(f_nv) != 1) stop("No encuentro un único xlsx de NV para ", yr)
  nv <- as.data.table(read_xlsx(f_nv))[, .(VE_REF = ID_UDPROD, FECHAVENTA, FAO,
                                           TOTAL_KILOS, TOTAL_EUROS)]
  nv[, FECHAVENTA := as.Date(FECHAVENTA)]
  nv[, FT_REF := paste0(VE_REF, "_", format(FECHAVENTA, "%Y_%m_%d"))]
  nv <- nv[FT_REF %in% viajes$FT_REF & !is.na(FAO) & FAO != ""]

  nv_sp <- nv[, .(KG   = sum(TOTAL_KILOS, na.rm = TRUE),
                  EURO = sum(TOTAL_EUROS, na.rm = TRUE)), by = .(FT_REF, FAO)]
  rm(nv)

  # Formato ancho: LE_KG_<FAO>, LE_EURO_<FAO>
  kg <- dcast(nv_sp, FT_REF ~ FAO, value.var = "KG",   fill = 0)
  eu <- dcast(nv_sp, FT_REF ~ FAO, value.var = "EURO", fill = 0)
  sp <- setdiff(names(kg), "FT_REF")
  setnames(kg, sp, paste0("LE_KG_",   sp))
  setnames(eu, sp, paste0("LE_EURO_", sp))

  # 4.4 Métier nivel 6 ---------------------------------------------------------
  # Grupo de especies con más kg entre los que dan un métier válido para ese arte.
  # Malla desconocida -> ">0" en artes con malla; "0" en anzuelos.
  nv_sp <- merge(nv_sp, grupos, by = "FAO", all.x = TRUE)
  sin_grupo <- nv_sp[is.na(GROUP), .(KG = sum(KG)), by = FAO][order(-KG)]
  nv_sp[is.na(GROUP), GROUP := "UNK"]

  kg_grp <- nv_sp[, .(KG = sum(KG)), by = .(FT_REF, GROUP)]
  kg_grp <- merge(kg_grp, viajes[, .(FT_REF, LE_GEAR)], by = "FT_REF")
  kg_grp[, LE_MET := paste0(LE_GEAR, "_", GROUP, "_",
                            fifelse(LE_GEAR %in% c("LHP", "LLS", "LLD"), "0_0_0", ">0_0_0"))]
  met <- kg_grp[LE_MET %in% met_validos][order(FT_REF, -KG)][
    , .(LE_MET = LE_MET[1]), by = FT_REF]

  # 4.5 EFLALO -----------------------------------------------------------------
  eflalo <- viajes[, .(
    VE_REF,
    VE_FLT   = seg_dcf(LE_GEAR),
    VE_COU   = "ESP",
    FT_REF,
    FT_DCOU  = "ESP",
    FT_DHAR  = NA_character_,                          # UN/LOCODE: no usado por los scripts
    FT_DDAT  = format(t_ini, "%d/%m/%Y", tz = "UTC"),
    FT_DTIME = format(t_ini, "%H:%M",    tz = "UTC"),
    FT_LCOU  = "ESP",
    FT_LHAR  = NA_character_,
    FT_LDAT  = format(t_fin, "%d/%m/%Y", tz = "UTC"),
    FT_LTIME = format(t_fin, "%H:%M",    tz = "UTC"),
    LE_ID    = paste0(FT_REF, "_1"),                   # un evento por viaje
    LE_CDAT  = format(FECHAVENTA, "%d/%m/%Y"),
    LE_GEAR,
    LE_MSZ   = NA_real_,                               # malla desconocida
    LE_DIV   = "27.9.a"                                # Golfo de Cádiz
  )]

  eflalo <- merge(eflalo, rect,  by = "FT_REF", all.x = TRUE)
  eflalo <- merge(eflalo, met,   by = "FT_REF", all.x = TRUE)
  eflalo <- merge(eflalo, censo, by = "VE_REF", all.x = TRUE)
  eflalo <- merge(eflalo, kg,    by = "FT_REF", all.x = TRUE)
  eflalo <- merge(eflalo, eu,    by = "FT_REF", all.x = TRUE)

  # 4.6 Controles ----------------------------------------------------------------
  n0 <- nrow(eflalo)
  message(sprintf("eflalo: %s viajes", format(n0, big.mark = ",")))

  message("Viajes por arte y métier:")
  print(eflalo[, .N, by = .(LE_GEAR, LE_MET)][order(LE_GEAR, -N)])

  if (nrow(sin_grupo)) {
    message("Especies sin grupo RCG (asignadas a UNK), kg:")
    print(head(sin_grupo, 20))
  }

  sin_met <- eflalo[is.na(LE_MET), .N, by = LE_GEAR]
  if (nrow(sin_met)) {
    message("Viajes SIN métier N6 válido (el script 1 de ICES los descartará):")
    print(sin_met)
  }

  sin_censo <- eflalo[is.na(VE_LEN) | is.na(VE_KW), unique(VE_REF)]
  if (length(sin_censo)) {
    message(sprintf("Aviso: %d barcos (%s viajes) sin eslora/kW: se mantienen, con kW-horas y clase de eslora en NA.",
                    length(sin_censo), eflalo[VE_REF %in% sin_censo, .N]))
  }

  message(sprintf("Viajes en tacsat sin eflalo: %d",
                  length(setdiff(unique(CV$FT_REF), eflalo$FT_REF))))

  # 4.7 Guardar con los nombres que espera 1_eflalo_tacsat_preprocessing.R ------
  nm_t <- paste0("tacsat_", yr)
  nm_e <- paste0("eflalo_", yr)
  assign(nm_t, as.data.frame(tacsat))
  assign(nm_e, as.data.frame(eflalo))
  save(list = nm_t, file = file.path(ruta_salida, paste0(nm_t, ".RData")))
  save(list = nm_e, file = file.path(ruta_salida, paste0(nm_e, ".RData")))
  message("Guardado: ", file.path(ruta_salida, paste0(c(nm_t, nm_e), ".RData")))

  rm(list = c(nm_t, nm_e)); rm(CV, tacsat, viajes, rect, nv_sp, kg, eu, kg_grp, met, eflalo)
  gc()
}

rm(CV_all); gc()

# =============================================================================
# 03_join_CV_NV.R
# Lectura de Notas de Venta 2018-2024, join con CV_all y asignación LE_MET4
# =============================================================================

library(data.table)
library(readxl)

# CV_all debe estar en memoria desde 02_format_CV.R
CV_all <- fst::read_fst("outputs/CV_2017_2024_ZE_fmt.fst", as.data.table = TRUE)

# --- Ruta NV -----------------------------------------------------------------
path_nv <- "data/Junta Andalucía/Notas de Venta_2018_2025/20251001 IEO NV 2018-2025"

# --- 1. Lectura NV 2018-2024 -------------------------------------------------
message("Leyendo Notas de Venta...")

years_nv <- 2018:2024
lista_nv <- lapply(years_nv, function(yr) {
  f <- list.files(path_nv, pattern = paste0("NV ", yr, "\\.xlsx$"), full.names = TRUE)
  dt <- as.data.table(readxl::read_xlsx(f))
  dt[, ANYO_NV := yr]
  message(sprintf("  NV %d: %s filas", yr, format(nrow(dt), big.mark = ",")))
  return(dt)
})

NV_all <- rbindlist(lista_nv, use.names = TRUE, fill = TRUE)
rm(lista_nv); gc()

message(sprintf("Total NV 2018-2024: %s filas", format(nrow(NV_all), big.mark = ",")))

# --- 2. Estandarizar NV ------------------------------------------------------
# Renombrar ID_UDPROD → VE_REF
setnames(NV_all, "ID_UDPROD", "VE_REF")

# FECHAVENTA: readxl devuelve Date directamente desde xlsx
NV_all[, FECHAVENTA := as.Date(FECHAVENTA)]

# Crear FT_REF (clave de unión: barco + día)
NV_all[, FT_REF := paste0(VE_REF, "_", format(FECHAVENTA, "%Y_%m_%d"))]

# --- 3. Colapsar NV a una fila por FT_REF ------------------------------------
# Cada marea tiene múltiples especies vendidas → join cartesiano si no se colapsa.
# Estrategia: arte dominante = el que acumula más kilos en esa marea.
# Se suman kilos y euros totales por FT_REF.

message("Colapsando NV a una fila por FT_REF...")

# Arte dominante por FT_REF
arte_dominante <- NV_all[, .(TOTAL_KILOS_ARTE = sum(TOTAL_KILOS, na.rm = TRUE)),
                          by = .(FT_REF, ARTE)]
arte_dominante <- arte_dominante[order(FT_REF, -TOTAL_KILOS_ARTE)]
arte_dominante <- arte_dominante[, .SD[1], by = FT_REF]  # primer arte = mayor captura

# Totales por FT_REF
nv_totales <- NV_all[, .(
  VE_REF    = VE_REF[1],
  FECHAVENTA = FECHAVENTA[1],
  ANYO_NV   = ANYO_NV[1],
  TOTAL_KILOS = sum(TOTAL_KILOS, na.rm = TRUE),
  TOTAL_EUROS = sum(TOTAL_EUROS, na.rm = TRUE),
  N_ESPECIES  = .N
), by = FT_REF]

# Unir arte dominante con totales
NV_key <- merge(nv_totales, arte_dominante[, .(FT_REF, ARTE)], by = "FT_REF")

message(sprintf("NV colapsada: %s mareas únicas", format(nrow(NV_key), big.mark = ",")))
rm(NV_all, arte_dominante, nv_totales); gc()

# --- 4. Join CV × NV ---------------------------------------------------------
message("Cruzando CV con NV...")
CV_NV <- merge(CV_all, NV_key, by = "FT_REF", all.x = FALSE)
# all.x = FALSE: inner join → solo pings con nota de venta asociada
# 2017 queda fuera (no hay NV); pings en puerto sin venta también quedan fuera

# Limpiar columnas duplicadas del merge (VE_REF, ARTE)
if ("VE_REF.y" %in% names(CV_NV)) CV_NV[, VE_REF.y := NULL]
if ("VE_REF.x" %in% names(CV_NV)) setnames(CV_NV, "VE_REF.x", "VE_REF")
# ARTE.x = arte del CV (puede estar vacío); ARTE.y = arte dominante de NV (por kilos)
# Nos quedamos con ARTE de NV para la clasificación LE_MET4
if ("ARTE.x" %in% names(CV_NV)) CV_NV[, ARTE.x := NULL]
if ("ARTE.y" %in% names(CV_NV)) setnames(CV_NV, "ARTE.y", "ARTE")

message(sprintf("CV_NV: %s filas", format(nrow(CV_NV), big.mark = ",")))

# --- 5. Asignación LE_MET4 desde ARTE ----------------------------------------
# Matching por grepl sobre patrones ASCII → inmune a problemas de encoding
# (los xlsx de la Junta de Andalucía tienen tildes con encoding inconsistente)
#
# Códigos DCF (Data Collection Framework) asignados:
#   DRB  – Rastra / draga remolcada
#   HMD  – Draga hidráulica
#   LHP  – Línea de mano (pole & line): manual y mecanizada
#   LHP  – Línea de mano (pole & line): manual, mecanizada y anzuelos
#          → "ANZUELOS Y LÍNEAS (SIN ESPECIFICAR)" se unifica aquí;
#            operativamente equivalente a línea de mano
#   GNS  – Red de enmalle de fondo anclada
#   GN   – Red de enmalle sin especificar
#   GTN  – Trasmallo combinado
#   GTR  – Trasmallo
#   FPO  – Nasa / trampa para crustáceos
#   FIX  – Trampa fija
#   LLS  – Palangre de fondo calado
#   LL   – Palangre sin especificar
#   LLD  – Palangre de deriva
#   LTL  – Curricán / troll
#   NK   – Arte no especificado


# aqui junto los artes como por ejemplo FIX Y FPO Y dejo fuera los que no usaremos como 
# LLD, NK, Cerco, LTL

CV_NV[, LE_MET4 := fcase(
  grepl("RASTRA",                                    ARTE, ignore.case = TRUE), "DRB",
  grepl("DRAGAS|HIDRA|HIDR.ULIC",                   ARTE, ignore.case = TRUE), "HMD",
  grepl("MANO.*MANUAL|MANUAL.*MANO",                ARTE, ignore.case = TRUE), "LHP",
  grepl("MANO.*MECAN|MECAN.*MANO",                  ARTE, ignore.case = TRUE), "LHP",
  grepl("ANZUELOS",                                 ARTE, ignore.case = TRUE), "LHP",
  grepl("ENMALLE DE FONDO|ENMALLE.*FONDO",          ARTE, ignore.case = TRUE), "GNS",
  grepl("ENMALLE.*SIN ESPECIF",                     ARTE, ignore.case = TRUE), "GN",
  grepl("COMBINADA.*TRASMALLO|TRASMALLO.*COMBINADA", ARTE, ignore.case = TRUE), "GTN",
  grepl("^TRASMALLOS$",                             ARTE, ignore.case = TRUE), "GTR",
  grepl("NASAS",                                    ARTE, ignore.case = TRUE), "FPO",
  grepl("TRAMPAS",                                  ARTE, ignore.case = TRUE), "FPO",
  grepl("PALANGRES CALADOS",                        ARTE, ignore.case = TRUE), "LLS",
  grepl("PALANGRES.*SIN ESPECIF",                   ARTE, ignore.case = TRUE), "LL",
  default = NA_character_
)]

# Distribución de registros por metier
message("Registros por LE_MET4:")
print(table(CV_NV$LE_MET4, useNA = "ifany"))

# Artes sin código asignado
sin_codigo <- CV_NV[is.na(LE_MET4), .N, by = ARTE][order(-N)]
if (nrow(sin_codigo) > 0) {
  message("Artes sin LE_MET4 asignado:")
  print(sin_codigo)
}

# --- 6. Eliminar artes no objetivo --------------------------------------------
# grepl para evitar problemas de encoding con tildes (EMBARCACIÓN, etc.)
n_antes <- nrow(CV_NV)
CV_NV <- CV_NV[!grepl("ARRASTRE|CERCO", ARTE, ignore.case = TRUE)]
CV_NV <- CV_NV[!is.na(LE_MET4)]
message(sprintf("Filas eliminadas (arrastre/cerco/sin código): %s",
  format(n_antes - nrow(CV_NV), big.mark = ",")))
message(sprintf("CV_NV final: %s filas", format(nrow(CV_NV), big.mark = ",")))

# --- 7. Guardar --------------------------------------------------------------
fst::write_fst(CV_NV, "outputs/CV_NV_2018_2024_ZE.fst", compress = 75)
CV_NV[, .N, by = .(ARTE, LE_MET4)][order(ARTE)]


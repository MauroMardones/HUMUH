# =============================================================================
# 04_effort_grid.R
# Filtro por velocidad de pesca, agregación de esfuerzo en celdas 0.01°×0.01°
# (~1×1 km en el Golfo de Cádiz)
# Incluye código C-square (estándar ICES) por celda
# =============================================================================

library(data.table)

# CV_NV debe estar en memoria desde 03_join_CV_NV.R
CV_NV <- fst::read_fst("outputs/CV_NV_2018_2024_ZE.fst", as.data.table = TRUE)

# =============================================================================
# FUNCIÓN C-SQUARE (estándar ICES)
# Genera el código identificador de celda geográfica a resolución 0.01°
# Fuente: Eastwood et al. (2003), adaptada del Rmd original del proyecto
# =============================================================================
csquare.cod <- function(lon, lat, degrees) {
  if (length(lon) != length(lat))
    stop("length of longitude not equal to length of latitude")
  if (!degrees %in% c(10, 5, 1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001))
    stop("degrees not in range: c(10,5,1,0.5,0.1,0.05,0.01,0.005,0.001)")

  dims <- length(lon)
  quadrants <- array(NA,
    dim      = c(5, 6, dims),
    dimnames = list(
      c("globalQuadrant","intmQuadrant1","intmQuadrant2","intmQuadrant3","intmQuadrant4"),
      c("quadrantDigit","latDigit","lonDigit","latRemain","lonRemain","code"),
      seq_len(dims)
    )
  )

  quadrants["globalQuadrant","quadrantDigit",] <-
    4 - (((2 * floor(1 + (lon/200))) - 1) * ((2 * floor(1 + (lat/200))) + 1))
  quadrants["globalQuadrant","latDigit",] <- floor(abs(lat)/10)
  quadrants["globalQuadrant","lonDigit",] <- floor(abs(lon)/10)
  quadrants["globalQuadrant","latRemain",] <-
    round(abs(lat) - (quadrants["globalQuadrant","latDigit",] * 10), 7)
  quadrants["globalQuadrant","lonRemain",] <-
    round(abs(lon) - (quadrants["globalQuadrant","lonDigit",] * 10), 7)
  quadrants["globalQuadrant","code",] <-
    quadrants["globalQuadrant","quadrantDigit",] * 1000 +
    quadrants["globalQuadrant","latDigit",]      * 100  +
    quadrants["globalQuadrant","lonDigit",]

  for (q in 1:4) {
    prev <- if (q == 1) "globalQuadrant" else paste0("intmQuadrant", q-1)
    curr <- paste0("intmQuadrant", q)
    quadrants[curr,"quadrantDigit",] <-
      (2 * floor(quadrants[prev,"latRemain",] * 0.2)) +
      floor(quadrants[prev,"lonRemain",] * 0.2) + 1
    quadrants[curr,"latDigit",] <- floor(quadrants[prev,"latRemain",])
    quadrants[curr,"lonDigit",] <- floor(quadrants[prev,"lonRemain",])
    quadrants[curr,"latRemain",] <-
      round((quadrants[prev,"latRemain",] - quadrants[curr,"latDigit",]) * 10, 7)
    quadrants[curr,"lonRemain",] <-
      round((quadrants[prev,"lonRemain",] - quadrants[curr,"lonDigit",]) * 10, 7)
    quadrants[curr,"code",] <-
      quadrants[curr,"quadrantDigit",] * 100 +
      quadrants[curr,"latDigit",]      * 10  +
      quadrants[curr,"lonDigit",]
  }

  if (degrees == 10)
    return(quadrants["globalQuadrant","code",])
  if (degrees == 5)
    return(paste0(quadrants["globalQuadrant","code",],":",
                  quadrants["intmQuadrant1","quadrantDigit",]))
  if (degrees == 1)
    return(paste0(quadrants["globalQuadrant","code",],":",
                  quadrants["intmQuadrant1","code",]))
  if (degrees == 0.5)
    return(paste0(quadrants["globalQuadrant","code",],":",
                  quadrants["intmQuadrant1","code",],":",
                  quadrants["intmQuadrant2","quadrantDigit",]))
  if (degrees == 0.1)
    return(paste0(quadrants["globalQuadrant","code",],":",
                  quadrants["intmQuadrant1","code",],":",
                  quadrants["intmQuadrant2","code",]))
  if (degrees == 0.05)
    return(paste0(quadrants["globalQuadrant","code",],":",
                  quadrants["intmQuadrant1","code",],":",
                  quadrants["intmQuadrant2","code",],":",
                  quadrants["intmQuadrant3","quadrantDigit",]))
  if (degrees == 0.01)
    return(paste0(quadrants["globalQuadrant","code",],":",
                  quadrants["intmQuadrant1","code",],":",
                  quadrants["intmQuadrant2","code",],":",
                  quadrants["intmQuadrant3","code",]))
  if (degrees == 0.005)
    return(paste0(quadrants["globalQuadrant","code",],":",
                  quadrants["intmQuadrant1","code",],":",
                  quadrants["intmQuadrant2","code",],":",
                  quadrants["intmQuadrant3","code",],":",
                  quadrants["intmQuadrant4","quadrantDigit",]))
  if (degrees == 0.001)
    return(paste0(quadrants["globalQuadrant","code",],":",
                  quadrants["intmQuadrant1","code",],":",
                  quadrants["intmQuadrant2","code",],":",
                  quadrants["intmQuadrant3","code",],":",
                  quadrants["intmQuadrant4","code",]))
}

# =============================================================================
# 1. FILTRO POR VELOCIDAD DE PESCA POR ARTE
# =============================================================================
message("Filtrando por velocidad de pesca...")

CV_NV[, SI_SPCA := as.numeric(SI_SPCA)]

CV_fil <- CV_NV[
  (LE_MET4 == "DRB" & SI_SPCA >= 0.1 & SI_SPCA <= 1.0) |
    (LE_MET4 == "HMD" & SI_SPCA >= 1.0 & SI_SPCA <= 3.0) |
    (LE_MET4 == "FPO" & SI_SPCA >= 0.1 & SI_SPCA <= 3.0) |
    (LE_MET4 == "GN"  & SI_SPCA >= 0.1 & SI_SPCA <= 4.5) |
    (LE_MET4 == "GNS" & SI_SPCA >= 0.1 & SI_SPCA <= 2.0) |
    (LE_MET4 == "GTN" & SI_SPCA >= 0.1 & SI_SPCA <= 2.0) |
    (LE_MET4 == "GTR" & SI_SPCA >= 0.1 & SI_SPCA <= 2.0) |
    (LE_MET4 == "LL"  & SI_SPCA >= 0.1 & SI_SPCA <= 2.5) |
    (LE_MET4 == "LLS" & SI_SPCA >= 0.1 & SI_SPCA <= 2.5) |
    (LE_MET4 == "LHP" & SI_SPCA >= 0.1 & SI_SPCA <= 1.5)
]

message(sprintf("Pings en estado de pesca: %s", format(nrow(CV_fil), big.mark = ".")))

# =============================================================================
# 2. ASIGNAR CELDA 0.01° Y CÓDIGO C-SQUARE
# =============================================================================
message("Asignando celdas y códigos C-square...")

# --- Resolución de la grilla — descomenta UNA opción -------------------------
res <- 0.05   # ~5 km (≈ 4.5 km E-O × 5.6 km N-S a 37°N)
# res <- 0.01   # ~1 km (≈ 0.9 km E-O × 1.1 km N-S a 37°N)

CV_fil[, lon_cel := floor(SI_LONG / res) * res + res / 2]
CV_fil[, lat_cel := floor(SI_LATI / res) * res + res / 2]

centroides <- unique(CV_fil[, .(lon_cel, lat_cel)])
centroides[, csq      := csquare.cod(lon_cel, lat_cel, degrees = 0.05)]
centroides[, area_km2 := (res * 111.32) * (res * 111.32 * cos(lat_cel * pi / 180))]

CV_fil[centroides, `:=`(csq = i.csq, area_km2 = i.area_km2), on = .(lon_cel, lat_cel)]

message(sprintf("Celdas únicas: %s", format(nrow(centroides), big.mark = ".")))
message(sprintf("Ejemplo C-square: %s", centroides$csq[1]))

# =============================================================================
# 3. AGREGAR ESFUERZO POR CELDA × AÑO × ARTE
# =============================================================================
message("Agregando esfuerzo...")

# Por arte
esfuerzo_arte <- CV_fil[, .(
  horas_pesca = sum(SI_TDIFF, na.rm = TRUE) / 3600,
  n_barcos    = uniqueN(VE_REF)
), by = .(lon_cel, lat_cel, csq, ANYO, LE_MET4)]

# Total (todos los artes)
esfuerzo_total <- CV_fil[, .(
  horas_pesca = sum(SI_TDIFF, na.rm = TRUE) / 3600,
  n_barcos    = uniqueN(VE_REF)
), by = .(lon_cel, lat_cel, csq, ANYO)]
esfuerzo_total[, LE_MET4 := "TOTAL"]

esfuerzo <- rbindlist(list(esfuerzo_arte, esfuerzo_total), use.names = TRUE)
rm(CV_fil, esfuerzo_arte, esfuerzo_total, centroides); gc()

message(sprintf("Celdas × año × arte: %s", format(nrow(esfuerzo), big.mark = ".")))
message(sprintf("Artes presentes: %s",
  paste(sort(unique(esfuerzo$LE_MET4)), collapse = ", ")))

# =============================================================================
# 4. GUARDAR
# =============================================================================
fst::write_fst(esfuerzo, "outputs/esfuerzo_2018_2024.fst", compress = 75)
message("Script 04 completado. Campo 'csq' incluido en esfuerzo.")

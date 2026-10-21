# =============================================================================
# hintzen2019_agregacion_persistencia.R
#
# Agregación, zonas inaccesibles/desfavorables y persistencia espacial del
# esfuerzo pesquero a escala fina, siguiendo:
#   Hintzen, Aarts & Rijnsdorp (2019) ICES J. Mar. Sci. 76(4): 1072-1082
#   "Persistence in the fine-scale distribution and spatial aggregation of fishing"
#   (método basado en Ellis et al. 2014, CJFAS 71: 733-746)
#
# Adaptación HUMUH (Golfo de Cádiz, cajas verdes + notas de venta):
#   - Conteo  = nº de pings de pesca (o nº de mareas) por celda fina y año.
#               Con ping cada ~190 s, el ping equivale a la "pasada" de Hintzen.
#   - Celda fina   (gridcell)       = res_fina[arte]  (por defecto 0.0025°, ~250 m)
#   - Unidad grande = C-square de resolución configurable (res_unidad:
#                     0.1° por defecto; opciones 0.05°, 0.1° o 0.5°).
#                     Sustituye al rectángulo ICES del paper.
#
# Pasos (numeración de ecuaciones y figuras = paper):
#   1. Lectura, ANÁLISIS EXPLORATORIO (E1-E9) y conteos por celda fina × año
#      × arte (EN_PUERTO == 0 + filtro velocidad; pings o días por arte)
#   2. Diagnóstico de la distribución: sobredispersión, ley de Taylor, AIC
#      (Poisson / Geométrica / NB) y NB1 vs NB2 -> elección de la familia
#   3. Universo de celdas por C-square (mar y profundidad accesible al arte)
#   4. Ec. 1: Binomial Negativa truncada en cero, beta aleatorio por año
#      dentro de cada C-square -> TMB.
#      NB2: Var = mu + beta*mu^2 (beta = 1/theta, independiente de mu)
#      NB1: Var = mu(1+beta) (paper). familia = "auto" elige por AIC.
#   5. Ec. 2: ceros en exceso (P_false) ~ esfuerzo -> inaccesible / desfavorable
#   6. Variabilidad interanual por celda (theta NB, offset log esfuerzo total)
#      + test de permutación (100 reparticiones aleatorias)
#   7. Covariables: heterogeneidad batimétrica (BPI-5 km) y hábitat MSFD
#   8-10. Tablas, figuras emulando el paper (Fig. 1-8) y figuras ad hoc
#   11. Persistencia del caladero (beta-binomial de años pescados por celda)
#       para artes con metodo = "betabin" (p.ej. HMD)
#
# Entrada : outputs/CV_NV_2018_2024_ZE.fst   (salida de R/03_join_CV_NV.R)
# Salidas : outputs/hintzen_*.csv|fst   y   figs/hintzen_*.png
# Ejecutar desde la raíz del proyecto (HUMUH.Rproj)
# =============================================================================

library(data.table)
library(fst)
library(ggplot2)
library(sf)
library(TMB)

# =============================================================================
# 0. CONFIGURACIÓN  (editar aquí)
# =============================================================================

# --- Artes a analizar (códigos ya agrupados, ver agrupar_artes) -------------
# Dragas: HMD (draga hidráulica) y DRB (rastra). Disponible también: LHP
artes_sel <- c("FPO", "GTR", "GN", "LL", "HMD", "DRB")

# --- Agrupación de artes: código original = grupo ---------------------------
# Enmalles -> GN ; Trasmallos -> GTR ; Palangres -> LL. Cada ping conserva el
# filtro de velocidad de su arte ORIGINAL antes de agruparse (LE_MET4_orig).
# GTN (enmalle-trasmallo combinado) va con los trasmallos.
agrupar_artes <- c(GNS = "GN", GTN = "GTR", LLS = "LL")

anios <- 2018:2024

# --- Resolución de la celda fina (grados) por arte ---------------------------
# Debe dividir exactamente a res_unidad (res_unidad / res_fina = entero).
# Criterio: ~ distancia mediana entre pings de pesca del arte
#   HMD ~178 m, FPO ~135 m, DRB ~22 m, GTR ~81 m, GNS ~106 m
res_fina <- c(HMD = 0.0025, FPO = 0.0025, DRB = 0.001,
              GTR = 0.0025, GTN = 0.0025, GNS = 0.0025, GN = 0.0025,
              LL  = 0.0025, LLS = 0.0025, LHP = 0.0025)

# --- Unidad de agregación: C-square (en lugar del rectángulo ICES) ----------
res_unidad <- 0.1    # opciones: 0.05 (~5 km) | 0.1 (~10 km) | 0.5 (~50 km)

# --- Velocidades de pesca por arte (mismos rangos que R/04_effort_grid.R) ----
vel_pesca <- list(
  DRB = c(0.1, 1.0), HMD = c(1.0, 3.0), FPO = c(0.1, 3.0),
  GN  = c(0.1, 4.5), GNS = c(0.1, 2.0), GTN = c(0.1, 2.0),
  GTR = c(0.1, 2.0), LL  = c(0.1, 2.5), LLS = c(0.1, 2.5),
  LHP = c(0.1, 1.5)
)

# --- Profundidad máxima operativa por arte (m) para definir el universo ------
prof_max <- c(HMD = 30, FPO = 200, DRB = 30, GTR = 200, GTN = 200,
              GNS = 200, GN = 200, LL = 200, LLS = 200, LHP = 200)

# --- Unidad de conteo POR ARTE ----------------------------------------------
# "pings"  = nº de pings de pesca en la celda (≈ tiempo pescando en la celda)
# "mareas" = nº de días de pesca (barco × día) distintos en la celda
#            (≈ veces que se pescó la celda; recomendado si un mismo día deja
#             muchos pings en la misma celda, como la draga sobre un banco)
# Artes no listados usan "pings". Ver tabla eda_pings_por_marea_celda.
conteo <- c(FPO = "pings", GTR = "pings", GN = "pings", LL = "pings",
            HMD = "mareas", DRB = "mareas")

# --- Análisis exploratorio previo (sección 1b) -------------------------------
eda <- TRUE

# --- Método por arte -----------------------------------------------------------
# "hintzen" = NB truncada por C-square (secciones 2-10): agregación, zonas
#             inaccesibles/desfavorables y variabilidad interanual
# "betabin" = persistencia del caladero (sección 11): beta-binomial del nº de
#             años (de los disponibles) en que se pescó cada celda
#             Recomendado para caladeros fijos (p.ej. draga sobre bancos)
# "ambos"   = Hintzen (secciones 2-10) + beta-binomial (sección 11)
# Artes no listados usan "hintzen".
metodo <- c(FPO = "hintzen", GTR = "hintzen", GN = "hintzen", LL = "hintzen",
            HMD = "ambos", DRB = "ambos")
min_conteo_anio <- 1     # betabin: celda "pescada" en un año si conteo >= este valor
min_celdas_bb   <- 30    # betabin: mín. celdas de huella para ajustar un C-square

# --- Universo de celdas posibles dentro de cada C-square ---------------------
universo <- "mar"      # "mar"    = todas las celdas en mar y con prof <= prof_max
                       #            (+ todas las celdas pescadas alguna vez)
                       # "huella" = solo celdas pescadas algún año 2018-2024
                       #            (con "huella" la fracción inaccesible = 0)

min_celdas_pos  <- 30    # nº mínimo de celdas con conteo > 0 para ajustar un C-square-año
beta_max_valido <- 1e3   # beta > este valor = ajuste fallido (no converge) -> se descarta
beta_min_valido <- 1e-3  # beta < este valor = ajuste en el límite (sin sobredispersión
                         # estimable; suele ser C-square-año con pocos datos) -> se descarta
inacc_regla     <- "segundo_minimo"  # fracción inaccesible: "minimo" (paper) o
                                     # "segundo_minimo" (robusto a un año mal ajustado)
min_anios_theta <- 3     # nº mínimo de años con esfuerzo para test de persistencia
n_perm          <- 100   # permutaciones para el test de variabilidad interanual
umbral_hotspot  <- 0.5   # hotspot = celdas que acumulan el 50 % del esfuerzo del año

# --- Parametrización de la Binomial Negativa ---------------------------------
familia <- "auto"        # "auto" = la elige el diagnóstico de la sección 2 (AIC NB1 vs NB2)
                         # "NB2"  = Var = mu + beta*mu^2 (beta = 1/theta)
                         # "NB1"  = Var = mu*(1+beta)  (Hintzen et al. 2019)

# --- Diagnóstico de distribución (sección 2) ----------------------------------
diagnostico     <- TRUE  # FALSE para saltarlo (requiere familia "NB1" o "NB2")
min_celdas_diag <- 30    # mín. celdas > 0 por C-square-año para el diagnóstico
n_max_diag      <- 60    # máx. C-square-año por arte (los de más celdas) para AIC

motor <- "TMB"           # "TMB" (beta aleatorio por año, como el paper)
                         # "glmmTMB" (beta fijo por año; sin compilar C++)

usar_habitat <- TRUE     # proporción de sustrato MSFD por C-square (Fig. 7)
dir_habitat  <- "data/shp/capa_habitat_BHTY"

archivo_cvnv <- "outputs/CV_NV_2018_2024_ZE.fst"
archivo_bati <- "outputs/marmap_coord_-7.5;36;-5;37.3_res_1.csv"  # respaldo sin internet
pref_out <- "outputs/hintzen_"
pref_fig <- "figs/hintzen_"

set.seed(2019)

artes_label <- c(
  DRB = "DRB (Rastra)", HMD = "HMD (Draga hidráulica)", FPO = "FPO (Nasas/Trampas)",
  GN  = "GN (Enmalle s/e)", GNS = "GNS (Enmalle fondo)", GTN = "GTN (Enmalle-trasmallo)",
  GTR = "GTR (Trasmallo)", LL = "LL (Palangre s/e)", LLS = "LLS (Palangre calado)",
  LHP = "LHP (Línea de mano)"
)
artes_label_orig <- artes_label
for (g in unique(agrupar_artes)) {
  orig <- c(g, names(agrupar_artes)[agrupar_artes == g])
  nombre <- c(GN = "Enmalles", LL = "Palangres", GTR = "Trasmallos")[g]
  if (is.na(nombre)) nombre <- g
  artes_label[[g]] <- sprintf("%s (%s: %s)", g, nombre, paste(orig, collapse = "+"))
}

# Comprobaciones
stopifnot(all(artes_sel %in% names(res_fina)),
          all(artes_sel %in% names(vel_pesca)),
          all(artes_sel %in% names(prof_max)))
k_arte <- round(res_unidad / res_fina[artes_sel])
if (any(abs(res_unidad / res_fina[artes_sel] - k_arte) > 1e-6))
  stop("res_fina debe dividir exactamente a res_unidad")
if (!res_unidad %in% c(0.05, 0.1, 0.5))
  stop("res_unidad debe ser un C-square válido: 0.05, 0.1 o 0.5")

dir.create("figs", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)

# --- Tema gráfico (igual que R/05_maps.R) ------------------------------------
tema_base <- if (requireNamespace("ggthemes", quietly = TRUE)) ggthemes::theme_few() else theme_bw()
tema <- tema_base + theme(
  legend.position  = "bottom",
  plot.title       = element_text(face = "bold", size = 12),
  plot.strip.text       = element_text(size = 9, face = "bold"),
  legend.title     = element_text(size = 9, face = "bold")
)
tema_mapa <- tema + theme(panel.background = element_rect(fill = "#eaf4fb", colour = NA),
                          axis.text = element_text(size = 7))
# Colores por arte (Okabe-Ito, aptos para daltonismo). Los 4 artes del análisis
# (FPO, GTR, GN, LL) tienen colores bien distintos; el resto queda por si se añaden.
col_artes <- c(FPO = "#E69F00", GTR = "#CC79A7", GN = "#0072B2", LL = "#009E73",
               HMD = "#56B4E9", DRB = "#D55E00", LHP = "#999999",
               GNS = "#0072B2", GTN = "#0072B2", LLS = "#009E73")

guardar <- function(p, nombre, w = 10, h = 7) {
  f <- paste0(pref_fig, nombre, ".png")
  ggsave(f, p, width = w, height = h, dpi = 200)
  message("  figura: ", f)
}
escribir <- function(x, nombre) {
  f <- paste0(pref_out, nombre, ".csv")
  fwrite(x, f)
  message("  tabla : ", f)
}

# --- Unidad de conteo por arte -----------------------------------------------
conteo_de  <- function(a) if (a %in% names(conteo)) conteo[[a]] else "pings"
lab_conteo <- function(a) if (conteo_de(a) == "mareas") "días de pesca" else "pings"
metodo_de  <- function(a) if (a %in% names(metodo)) metodo[[a]] else "hintzen"
stopifnot(all(metodo %in% c("hintzen", "betabin", "ambos")))
stopifnot(all(conteo %in% c("pings", "mareas")),
          inacc_regla %in% c("minimo", "segundo_minimo"))

# --- Funciones de la NB según la parametrización elegida --------------------
stopifnot(familia %in% c("auto", "NB1", "NB2"))
if (familia == "auto" && !diagnostico) stop("familia = 'auto' requiere diagnostico = TRUE")
nb_size <- function(mu, beta) if (familia == "NB2") 1 / beta else mu / beta
p0_nb   <- function(mu, beta) dnbinom(0, size = nb_size(mu, beta), mu = mu)
dnb_tr  <- function(x, mu, beta)        # pmf truncada en cero
  dnbinom(x, size = nb_size(mu, beta), mu = mu) / (1 - p0_nb(mu, beta))
pnb_tr_cola <- function(q, mu, beta)    # P(Y >= q | Y > 0)
  pnbinom(q - 1, size = nb_size(mu, beta), mu = mu, lower.tail = FALSE) / (1 - p0_nb(mu, beta))

# =============================================================================
# 1. LECTURA, ANÁLISIS EXPLORATORIO Y CONTEOS POR CELDA FINA × AÑO × ARTE
# =============================================================================
message("1. Leyendo CV_NV...")

# --- 1a. Lectura ---------------------------------------------------------------
cols <- c("VE_REF", "FT_REF", "ANYO", "SI_TIMESTAMP", "SI_LATI", "SI_LONG",
          "SI_SPCA", "SI_TDIFF", "SI_DISTANCECA", "EN_PUERTO", "LE_MET4")
cv <- fst::read_fst(archivo_cvnv, columns = cols, as.data.table = TRUE)
# artes originales que entran (los seleccionados + los que se agrupan en ellos)
artes_orig <- unique(c(artes_sel, names(agrupar_artes)[agrupar_artes %in% artes_sel]))
stopifnot(all(artes_orig %in% names(vel_pesca)))
cv <- cv[LE_MET4 %in% artes_orig & ANYO %in% anios]

# filtro de velocidad según el arte ORIGINAL
vel_dt <- rbindlist(lapply(artes_orig, function(a)
  data.table(LE_MET4 = a, vmin = vel_pesca[[a]][1], vmax = vel_pesca[[a]][2])))
cv <- merge(cv, vel_dt, by = "LE_MET4")
cv[, pesca := EN_PUERTO == 0 & SI_SPCA >= vmin & SI_SPCA <= vmax]

# agrupación (enmalles -> GN, palangres -> LL)
cv[, LE_MET4_orig := LE_MET4]
cv[LE_MET4 %in% names(agrupar_artes), LE_MET4 := agrupar_artes[LE_MET4]]
message("   pings por arte original -> grupo:")
print(cv[, .(pings = .N, pings_pesca = sum(pesca)), by = .(LE_MET4_orig, grupo = LE_MET4)][order(grupo)])
cv[, rf := res_fina[LE_MET4]]
cv[, `:=`(ix = as.integer(floor(SI_LONG / rf + 1e-9)),
          iy = as.integer(floor(SI_LATI / rf + 1e-9)))]

# =============================================================================
# 1b. ANÁLISIS EXPLORATORIO DE LOS DATOS (antes de cualquier modelo)
# =============================================================================
# Preguntas:
#   E1 ¿Cuántos datos hay por arte y año? (pings, barcos, días, % puerto, % pesca)
#   E2 ¿El filtro de velocidad captura el modo de pesca? (histograma de velocidad)
#   E3 ¿La frecuencia de ping es regular? (intervalo entre pings)
#   E4 ¿Cómo se reparte la actividad en el año? (pings de pesca por mes)
#   E5 ¿Cuánto avanza el barco entre pings? (justifica el tamaño de celda fina)
#   E6 ¿Cuánto dura un día de pesca? (pings y horas por marea)
#   E7 ¿Cuántos pings deja un mismo día en una misma celda?
#       ~1 -> el ping equivale a una "pasada"; >> 1 -> mide tiempo -> usar "mareas"
#   E8 Relación pings vs días de pesca por celda-año
#   E9 Celdas con más pings: ¿pesca real o barcos parados/fondeados?
# Salidas: outputs/hintzen_eda_*.csv  y  figs/hintzen_eda_*.png
# =============================================================================
if (eda) {
  message("1b. Análisis exploratorio...")
  cv[, mes := as.integer(format(SI_TIMESTAMP, "%m"))]

  # --- E1. Volumen de datos por arte y año ----------------------------------
  e1 <- cv[, .(pings_total   = .N,
               pct_en_puerto = 100 * mean(EN_PUERTO == 1),
               pct_pesca     = 100 * mean(pesca),
               pings_pesca   = sum(pesca),
               barcos        = uniqueN(VE_REF),
               dias_pesca    = uniqueN(FT_REF[pesca]),
               horas_pesca   = sum(pmin(SI_TDIFF, 600)[pesca], na.rm = TRUE) / 3600),
           by = .(LE_MET4, ANYO)][order(LE_MET4, ANYO)]
  escribir(e1, "eda_resumen_arte_anio")
  print(e1)
  e1[, `:=`(pings_pesca = as.numeric(pings_pesca), barcos = as.numeric(barcos),
            dias_pesca = as.numeric(dias_pesca))]
  e1l <- melt(e1, id.vars = c("LE_MET4", "ANYO"),
              measure.vars = c("pings_pesca", "barcos", "dias_pesca", "horas_pesca"))
  e1l[, variable := factor(variable, levels = c("pings_pesca", "barcos", "dias_pesca", "horas_pesca"),
                           labels = c("Pings de pesca", "Barcos", "Días de pesca", "Horas de pesca"))]
  pE1 <- ggplot(e1l, aes(factor(ANYO), value, fill = LE_MET4)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = col_artes, labels = artes_label, name = NULL) +
    facet_wrap(~ variable, scales = "free_y") + tema +
    labs(title = "E1. Volumen de datos por arte y año",
         x = "Año", y = NULL)
  guardar(pE1, "eda_E1_volumen_anual", 11, 7)

  # --- E2. Velocidad (fuera de puerto) y rango de pesca ---------------------
  e2 <- cv[EN_PUERTO == 0 & SI_SPCA <= 10, .(LE_MET4 = LE_MET4_orig, SI_SPCA)]
  pE2 <- ggplot(e2, aes(SI_SPCA)) +
    geom_rect(data = vel_dt, aes(xmin = vmin, xmax = vmax, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "orange", alpha = 0.25) +
    geom_histogram(binwidth = 0.1, boundary = 0, fill = "grey25") +
    facet_wrap(~ LE_MET4, scales = "free_y", labeller = as_labeller(artes_label_orig)) + tema +
    labs(title = "E2. Velocidad de los pings fuera de puerto (por arte original)",
         x = "Velocidad (nudos)", y = "Nº de pings")
  guardar(pE2, "eda_E2_velocidad", 11, 5)
  e2t <- cv[pesca == TRUE, .(pct_menor_0.5kn = 100 * mean(SI_SPCA < 0.5),
                             vel_p10 = quantile(SI_SPCA, 0.1), vel_mediana = median(SI_SPCA),
                             vel_p90 = quantile(SI_SPCA, 0.9)), by = .(LE_MET4, LE_MET4_orig)]
  escribir(e2t, "eda_velocidad_pesca")
  rm(e2)

  # --- E3. Intervalo entre pings --------------------------------------------
  pE3 <- ggplot(cv[pesca == TRUE & SI_TDIFF <= 900], aes(SI_TDIFF)) +
    geom_histogram(binwidth = 10, boundary = 0, fill = "grey25") +
    facet_wrap(~ LE_MET4, scales = "free_y", labeller = as_labeller(artes_label)) +
    scale_y_sqrt() + tema +
    labs(title = "E3. Intervalo entre pings consecutivos (pings de pesca)",
         x = "Segundos", y = "Nº de pings (raíz)")
  guardar(pE3, "eda_E3_intervalo_pings", 11, 5)

  # --- E4. Estacionalidad ----------------------------------------------------
  e4 <- cv[pesca == TRUE, .N, by = .(LE_MET4, ANYO, mes)]
  pE4 <- ggplot(e4, aes(factor(mes), factor(ANYO), fill = N)) +
    geom_tile(colour = "white") +
    scale_fill_viridis_c(trans = "log10", name = "Pings de pesca") +
    facet_wrap(~ LE_MET4, labeller = as_labeller(artes_label)) + tema +
    labs(title = "E4. Actividad pesquera por mes y año", x = "Mes", y = "Año")
  guardar(pE4, "eda_E4_estacionalidad", 11, 5)

  # --- E5. Distancia entre pings de pesca vs tamaño de celda ----------------
  e5 <- cv[pesca == TRUE & SI_DISTANCECA > 0 & SI_DISTANCECA < 2000]
  celda_m <- data.table(LE_MET4 = artes_sel,
                        lado_m = res_fina[artes_sel] * 111320 * cos(36.8 * pi / 180))
  pE5 <- ggplot(e5, aes(SI_DISTANCECA)) +
    geom_histogram(bins = 80, fill = "grey25") +
    geom_vline(data = celda_m, aes(xintercept = lado_m), colour = "red", linetype = "dashed") +
    scale_x_log10() +
    facet_wrap(~ LE_MET4, scales = "free_y", labeller = as_labeller(artes_label)) + tema +
    labs(title = "E5. Distancia recorrida entre pings de pesca consecutivos",
         x = "Metros (log)", y = "Nº de pings")
  guardar(pE5, "eda_E5_distancia_pings", 11, 5)
  e5t <- e5[, .(dist_mediana_m = median(SI_DISTANCECA), dist_p90_m = quantile(SI_DISTANCECA, 0.9)),
            by = LE_MET4]
  e5t <- merge(e5t, celda_m, by = "LE_MET4")
  escribir(e5t, "eda_distancia_pings")
  rm(e5)

  # --- E6. Pings y horas por día de pesca (marea) ---------------------------
  e6 <- cv[pesca == TRUE, .(pings = .N, horas = sum(pmin(SI_TDIFF, 600), na.rm = TRUE) / 3600,
                            celdas = uniqueN(paste(ix, iy))), by = .(LE_MET4, FT_REF)]
  pE6 <- ggplot(e6, aes(horas)) +
    geom_histogram(binwidth = 0.25, boundary = 0, fill = "grey25") +
    coord_cartesian(xlim = c(0, 14)) +
    facet_wrap(~ LE_MET4, scales = "free_y", labeller = as_labeller(artes_label)) + tema +
    labs(title = "E6. Horas de pesca por día de pesca (barco × día)",
         x = "Horas de pesca en el día", y = "Nº de días de pesca")
  guardar(pE6, "eda_E6_horas_por_dia", 11, 5)

  # --- E7. Pings que deja un mismo día en una misma celda (clave) -----------
  e7 <- cv[pesca == TRUE, .N, by = .(LE_MET4, FT_REF, ix, iy)]
  e7t <- e7[, .(n = .N, mediana = as.numeric(median(N)), media = mean(N),
                pct_1_ping = 100 * mean(N == 1), pct_mas_10 = 100 * mean(N > 10),
                p90 = as.numeric(quantile(N, 0.9)), max = as.numeric(max(N))), by = LE_MET4]
  e7t[, conteo_usado := vapply(LE_MET4, conteo_de, character(1))]
  e7t[, sugerencia := fifelse(mediana <= 2, "pings ~ pasadas: 'pings' es adecuado",
                              "pings miden tiempo en la celda: usar 'mareas'")]
  escribir(e7t, "eda_pings_por_marea_celda")
  print(e7t)
  pE7 <- ggplot(e7, aes(pmin(N, 100))) +
    geom_histogram(binwidth = 1, boundary = 0.5, fill = "grey25") +
    scale_y_sqrt() +
    facet_wrap(~ LE_MET4, scales = "free_y", labeller = as_labeller(artes_label)) + tema +
    labs(title = "E7. Pings que deja un mismo día de pesca en una misma celda fina",
         x = "Pings por (día de pesca × celda)", y = "Frecuencia (raíz)")
  guardar(pE7, "eda_E7_pings_por_dia_celda", 11, 5)
  rm(e7)

  # --- E8. Pings vs días de pesca por celda-año ------------------------------
  e8 <- cv[pesca == TRUE, .(pings = .N, dias = uniqueN(FT_REF)), by = .(LE_MET4, ANYO, ix, iy)]
  pE8 <- ggplot(e8, aes(dias, pings)) +
    geom_bin2d(bins = 50) +
    geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
    scale_x_log10() + scale_y_log10() +
    scale_fill_viridis_c(trans = "log10", name = "Nº celdas-año") +
    facet_wrap(~ LE_MET4, labeller = as_labeller(artes_label)) + tema +
    labs(title = "E8. Pings vs días de pesca por celda fina × año",
         x = "Días de pesca en la celda (log)", y = "Pings en la celda (log)")
  guardar(pE8, "eda_E8_pings_vs_dias", 11, 5)
  rm(e8)

  # --- E9. Celdas con más pings: ¿pesca o barcos parados? --------------------
  e9 <- cv[pesca == TRUE, .(pings = .N, barcos = uniqueN(VE_REF), dias = uniqueN(FT_REF),
                            vel_mediana = as.numeric(median(SI_SPCA)),
                            lon = mean(SI_LONG), lat = mean(SI_LATI)),
           by = .(LE_MET4, ix, iy)][order(LE_MET4, -pings)][, head(.SD, 15), by = LE_MET4]
  e9[, pings_por_dia := pings / dias]
  escribir(e9, "eda_top_celdas")
  print(e9)
  rm(e6); gc()
}

# --- 1c. Conteos por celda fina × año × arte ----------------------------------
message("1c. Conteos por celda fina...")
cv <- cv[pesca == TRUE]
message(sprintf("   pings de pesca (fuera de puerto, rango de velocidad): %s",
                format(nrow(cv), big.mark = ",")))

cnt <- cv[, .(pings = .N, mareas = uniqueN(FT_REF)), by = .(LE_MET4, ANYO, ix, iy)]
rm(cv); gc()
cnt[, y := if (conteo_de(LE_MET4[1]) == "mareas") mareas else pings, by = LE_MET4]
message("   unidad de conteo: ",
        paste(sprintf("%s = %s", artes_sel, vapply(artes_sel, conteo_de, character(1))), collapse = ", "))

cnt[, rf := res_fina[LE_MET4]]
cnt[, k  := as.integer(round(res_unidad / rf))]
cnt[, `:=`(ux = ix %/% k, uy = iy %/% k)]
cnt[, `:=`(lon = (ix + 0.5) * rf, lat = (iy + 0.5) * rf)]

# --- Código C-square a la resolución res_unidad (función de R/04_effort_grid.R)
csquare.cod <- NULL
try({
  ex <- parse("R/04_effort_grid.R")
  i_fun <- which(vapply(ex, function(e) is.call(e) && length(e) >= 3 &&
                          identical(e[[2]], as.name("csquare.cod")), logical(1)))
  if (length(i_fun)) eval(ex[[i_fun[1]]])
}, silent = TRUE)

unidades <- unique(cnt[, .(LE_MET4, ux, uy, k, rf)])
unidades[, `:=`(lon_u = (ux + 0.5) * res_unidad, lat_u = (uy + 0.5) * res_unidad)]
unidades[, csq := if (is.function(csquare.cod))
  csquare.cod(lon_u, lat_u, degrees = res_unidad) else sprintf("u%d_%d", ux, uy)]
cnt <- merge(cnt, unidades[, .(LE_MET4, ux, uy, csq, lon_u, lat_u)],
             by = c("LE_MET4", "ux", "uy"))

fst::write_fst(cnt, paste0(pref_out, "conteos_celda.fst"), compress = 75)
message(sprintf("   celdas-año con esfuerzo: %s | C-squares: %s",
                format(nrow(cnt), big.mark = "."), nrow(unidades)))

# --- 1d. Reparto de artes por método ------------------------------------------
artes_todos <- artes_sel
met_sel  <- vapply(artes_sel, metodo_de, character(1))
artes_bb <- artes_sel[met_sel %in% c("betabin", "ambos")]   # sección 11
artes_sel <- artes_sel[met_sel %in% c("hintzen", "ambos")]  # secciones 2-10
cnt_bb   <- cnt[LE_MET4 %in% artes_bb]
cnt      <- cnt[LE_MET4 %in% artes_sel]
unidades <- unidades[LE_MET4 %in% artes_sel]
message("   método Hintzen (NB): ", if (length(artes_sel)) paste(artes_sel, collapse = ", ") else "ninguno")
message("   método beta-binomial (persistencia): ", if (length(artes_bb)) paste(artes_bb, collapse = ", ") else "ninguno")

if (length(artes_sel) > 0) {   # ===== inicio bloque Hintzen (secciones 2-10) =====

# =============================================================================
# 2. DIAGNÓSTICO DE LA DISTRIBUCIÓN DE LOS CONTEOS (elección de la familia)
# =============================================================================
# Preguntas:
#   a) ¿Hay sobredispersión?  Var/Media >> 1  -> descarta Poisson
#   b) ¿Cómo escala la varianza? Ley de Taylor Var = a*Media^b
#        b ~ 1 -> NB1 (paper) ; b ~ 2 -> NB2
#   c) AIC por C-square-año (truncadas): Poisson, Geométrica, NB
#   d) NB1 vs NB2: modelo conjunto por arte con media propia por C-square-año
#      y dispersión COMÚN. Dentro de un solo C-square-año NB1 y NB2 son la
#      misma distribución; solo se distinguen al compartir la dispersión
#      entre unidades con distinta media.
# Salidas: outputs/hintzen_dist_*.csv  y  figs/hintzen_dist_*.png
# =============================================================================
if (diagnostico) {
  message("2. Diagnóstico de la distribución de los conteos...")
  if (!requireNamespace("glmmTMB", quietly = TRUE)) stop("Instalar glmmTMB para el diagnóstico")

  cnt[, csq_anio := paste(csq, ANYO, sep = "_")]

  # --- a) resumen general ----------------------------------------------------
  d_res <- cnt[, .(n_celdas = .N, media = mean(y), mediana = as.numeric(median(y)),
                   var_media = var(y) / mean(y), p_uno = mean(y == 1),
                   p90 = as.numeric(quantile(y, 0.9)), max = as.numeric(max(y))), by = LE_MET4]
  escribir(d_res, "dist_resumen")
  print(d_res)

  p_h <- ggplot(cnt, aes(y)) +
    geom_histogram(bins = 60, fill = "grey30") +
    scale_x_log10() + scale_y_log10() +
    facet_wrap(~ LE_MET4, scales = "free", labeller = as_labeller(artes_label)) + tema +
    labs(title = "Distribución de la intensidad por celda fina (celda × año)",
         x = "Intensidad por celda (log)", y = "Nº de celdas (log)")
  guardar(p_h, "dist_histograma", 10, 5)

  # --- b) ley de Taylor ------------------------------------------------------
  mv <- cnt[, .(media = mean(y), varianza = var(y), n = .N),
            by = .(LE_MET4, ANYO, csq, csq_anio)][n >= min_celdas_diag]
  n_var0 <- mv[!(varianza > 0), .N, by = LE_MET4]
  mv <- mv[varianza > 0 & is.finite(varianza)]
  taylor <- mv[, {
    m <- lm(log(varianza) ~ log(media))
    ci <- confint(m)[2, ]
    list(n_csq_anio = .N, b_taylor = coef(m)[2], b_IC95_inf = ci[1], b_IC95_sup = ci[2],
         R2 = summary(m)$r.squared, mediana_var_media = median(varianza / media))
  }, by = LE_MET4]
  taylor <- merge(taylor, setnames(n_var0, "N", "n_excluidos_var0"), by = "LE_MET4", all.x = TRUE)
  escribir(taylor, "dist_taylor")
  print(taylor)

  p_t <- ggplot(mv, aes(media, varianza)) +
    geom_point(aes(colour = factor(ANYO)), alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black") +
    scale_x_log10() + scale_y_log10() +
    scale_colour_viridis_d(name = "Año") +
    facet_wrap(~ LE_MET4, scales = "free", labeller = as_labeller(artes_label)) + tema +
    labs(title = "Ley de Taylor: media vs varianza por C-square × año (log-log)",
         x = "Intensidad media por celda", y = "Varianza")
  guardar(p_t, "dist_taylor", 10, 5)

  # --- c) AIC por C-square-año ----------------------------------------------
  sel_d <- mv[order(-n)][, head(.SD, n_max_diag), by = LE_MET4]
  aic1 <- function(y, fam) {
    m <- tryCatch(suppressWarnings(glmmTMB::glmmTMB(y ~ 1, family = fam, data = data.frame(y = y))),
                  error = function(e) NULL)
    if (is.null(m)) return(c(NA_real_, NA_real_))
    c(AIC(m), tryCatch(sigma(m), error = function(e) NA_real_))
  }
  aic_geom <- function(y) {
    nll <- function(lm) { p <- 1 / (1 + exp(lm)); -sum(dgeom(y, p, log = TRUE) - log(1 - p)) }
    2 * optimize(nll, c(-5, 12))$objective + 2
  }
  d_aic <- rbindlist(lapply(seq_len(nrow(sel_d)), function(i) {
    s <- sel_d[i]
    y <- cnt[LE_MET4 == s$LE_MET4 & csq_anio == s$csq_anio, y]
    po <- aic1(y, glmmTMB::truncated_poisson())
    nb <- aic1(y, glmmTMB::truncated_nbinom2())
    data.table(LE_MET4 = s$LE_MET4, csq = s$csq, ANYO = s$ANYO, n = length(y),
               Poisson = po[1], Geometrica = aic_geom(y), NB = nb[1], theta_NB = nb[2])
  }))
  d_aicl <- melt(d_aic, id.vars = c("LE_MET4", "csq", "ANYO", "n"),
                 measure.vars = c("Poisson", "Geometrica", "NB"),
                 variable.name = "modelo", value.name = "AIC")
  d_aicl[, dAIC := AIC - min(AIC, na.rm = TRUE), by = .(LE_MET4, csq, ANYO)]
  d_gana <- d_aicl[, .(pct_mejor = 100 * mean(dAIC == 0, na.rm = TRUE),
                       pct_aceptable = 100 * mean(dAIC < 2, na.rm = TRUE),
                       dAIC_mediana = median(dAIC, na.rm = TRUE)), by = .(LE_MET4, modelo)]
  escribir(d_aic, "dist_aic_csq_anio")
  escribir(d_gana, "dist_aic_resumen")
  print(d_gana)

  p_a <- ggplot(d_aicl[is.finite(dAIC)], aes(modelo, pmin(dAIC, 500), fill = modelo)) +
    geom_boxplot(outlier.size = 0.5) + scale_y_sqrt() +
    facet_wrap(~ LE_MET4, labeller = as_labeller(artes_label)) + tema +
    theme(legend.position = "none") +
    labs(title = "ΔAIC respecto al mejor modelo, por C-square × año (truncadas en cero)",
         x = NULL, y = "ΔAIC")
  guardar(p_a, "dist_dAIC", 9, 5)

  # --- d) NB1 vs NB2 con dispersión común (modelo conjunto por arte) --------
  d_nb <- rbindlist(lapply(artes_sel, function(a) {
    ids <- sel_d[LE_MET4 == a, csq_anio]
    if (length(ids) < 2) return(NULL)
    dd <- cnt[LE_MET4 == a & csq_anio %in% ids, .(y, f = factor(csq_anio))]
    m1 <- tryCatch(suppressWarnings(glmmTMB::glmmTMB(y ~ 0 + f, family = glmmTMB::truncated_nbinom1(), data = dd)),
                   error = function(e) NULL)
    m2 <- tryCatch(suppressWarnings(glmmTMB::glmmTMB(y ~ 0 + f, family = glmmTMB::truncated_nbinom2(), data = dd)),
                   error = function(e) NULL)
    data.table(LE_MET4 = a, n_csq_anio = length(ids), n_celdas = nrow(dd),
               AIC_NB1 = if (!is.null(m1)) AIC(m1) else NA_real_,
               AIC_NB2 = if (!is.null(m2)) AIC(m2) else NA_real_)
  }))
  d_nb[, dAIC_NB1_menos_NB2 := AIC_NB1 - AIC_NB2]
  d_nb[, familia_AIC := fifelse(is.na(dAIC_NB1_menos_NB2), NA_character_,
                                fifelse(dAIC_NB1_menos_NB2 > 0, "NB2", "NB1"))]
  d_nb <- merge(d_nb, taylor[, .(LE_MET4, b_taylor, b_IC95_inf, b_IC95_sup)], by = "LE_MET4")
  escribir(d_nb, "dist_nb1_vs_nb2")
  print(d_nb)

  # --- e) observado vs esperado (6 C-square-año con más celdas por arte) ---
  for (a in artes_sel) {
    ej <- sel_d[LE_MET4 == a][1:min(6, .N)]
    if (!nrow(ej)) next
    f_oe <- rbindlist(lapply(seq_len(nrow(ej)), function(i) {
      y <- cnt[LE_MET4 == a & csq_anio == ej$csq_anio[i], y]
      mp <- glmmTMB::glmmTMB(y ~ 1, family = glmmTMB::truncated_poisson(), data = data.frame(y = y))
      mn <- glmmTMB::glmmTMB(y ~ 1, family = glmmTMB::truncated_nbinom2(), data = data.frame(y = y))
      lp <- exp(glmmTMB::fixef(mp)$cond); mu <- exp(glmmTMB::fixef(mn)$cond); th <- sigma(mn)
      ym <- max(y)                               # histograma completo, sin acumular cola
      x <- 1:ym
      fp <- dpois(x, lp) / (1 - dpois(0, lp))
      fn <- dnbinom(x, size = th, mu = mu) / (1 - dnbinom(0, size = th, mu = mu))
      data.table(panel = sprintf("%s · %d", ej$csq[i], ej$ANYO[i]), x = x,
                 obs = tabulate(y, nbins = ym) / length(y),
                 Poisson = fp, NB = fn)
    }))
    f_oel <- melt(f_oe, id.vars = c("panel", "x", "obs"), variable.name = "modelo", value.name = "esp")
    p_oe <- ggplot(f_oel, aes(x)) +
      geom_col(data = f_oe, aes(y = obs), fill = "grey70") +
      geom_line(aes(y = esp, colour = modelo)) +
      scale_colour_manual(values = c(Poisson = "#d7301f", NB = "#2171b5"), name = NULL) +
      scale_y_sqrt() + facet_wrap(~ panel, scales = "free") + tema +
      labs(title = paste("Observado (barras) vs esperado (truncadas) —", artes_label[[a]]),
           x = sprintf("Nº de %s por celda (> 0)", lab_conteo(a)), y = "Proporción")
    guardar(p_oe, paste0("dist_obs_vs_esperado_", a), 11, 7)
  }

  # --- f) elección de la familia --------------------------------------------
  if (familia == "auto") {
    votos <- d_nb[!is.na(familia_AIC), familia_AIC]
    familia <- if (length(votos) && all(votos == "NB1")) "NB1" else "NB2"
    message(sprintf("   familia elegida por AIC (NB1 vs NB2 con dispersión común): %s", familia))
    if (length(unique(votos)) > 1)
      message("   AVISO: los artes no coinciden en la familia; se usa NB2 para todos (comparables)")
  }
  cnt[, csq_anio := NULL]
}
lab_beta <- if (familia == "NB2") "β = 1/θ (NB2)" else "β (NB1)"
message(sprintf("   Familia usada en la Ec. 1: %s", familia))

# =============================================================================
# 3. BATIMETRÍA, COSTA Y UNIVERSO DE CELDAS
# =============================================================================
message("3. Batimetría, costa y universo de celdas...")

bb <- c(xmin = min(unidades$lon_u) - 0.15, xmax = max(unidades$lon_u) + 0.15,
        ymin = min(unidades$lat_u) - 0.15, ymax = max(unidades$lat_u) + 0.15)

# --- Batimetría (NOAA vía marmap; respaldo: csv del proyecto) ----------------
bati <- NULL
if (requireNamespace("marmap", quietly = TRUE)) {
  bati <- tryCatch({
    b <- marmap::getNOAA.bathy(lon1 = bb[["xmin"]], lon2 = bb[["xmax"]],
                               lat1 = bb[["ymin"]], lat2 = bb[["ymax"]], resolution = 1)
    as.data.table(marmap::fortify.bathy(b))[, .(x, y, z)]
  }, error = function(e) NULL)
}
if (is.null(bati) && file.exists(archivo_bati)) {
  bati <- fread(archivo_bati)
  setnames(bati, c("x", "y", "z"))
  message("   batimetría desde csv local (cobertura parcial)")
}
if (is.null(bati)) message("   AVISO: sin batimetría; universo solo por costa")

if (!is.null(bati)) {
  bx <- sort(unique(round(bati$x, 6))); by_ <- sort(unique(round(bati$y, 6)))
  bres <- median(diff(bx))
  bati[, `:=`(gx = as.integer(round((x - bx[1]) / bres)),
              gy = as.integer(round((y - by_[1]) / bres)))]
  bati <- bati[, .(x = mean(x), y = mean(y), z = mean(z)), by = .(gx, gy)]
  # BPI-5 km: profundidad menos la media en un radio ~5 km (ventana cuadrada)
  rad <- max(1L, as.integer(round(5 / (bres * 111.32 * cos(mean(bati$y) * pi / 180)))))
  M <- matrix(NA_real_, max(bati$gx) + 1, max(bati$gy) + 1)
  M[cbind(bati$gx + 1, bati$gy + 1)] <- bati$z
  S <- matrix(0, nrow(M), ncol(M)); Nn <- S
  for (dx in -rad:rad) for (dy in -rad:rad) {
    ii <- (1:nrow(M)) + dx; jj <- (1:ncol(M)) + dy
    okx <- ii >= 1 & ii <= nrow(M); oky <- jj >= 1 & jj <= ncol(M)
    sub <- M[ii[okx], jj[oky], drop = FALSE]
    v <- !is.na(sub)
    S[okx, oky] <- S[okx, oky] + ifelse(v, sub, 0)
    Nn[okx, oky] <- Nn[okx, oky] + v
  }
  bati[, bpi := z - (S / pmax(Nn, 1))[cbind(gx + 1, gy + 1)]]
  prof_de <- function(lon, lat) {
    q <- data.table(gx = as.integer(round((lon - bx[1]) / bres)),
                    gy = as.integer(round((lat - by_[1]) / bres)))
    bati[q, on = .(gx, gy), z]
  }
}

# --- Costa (rnaturalearth, como en R/05_maps.R) -----------------------------
costa <- NULL
if (requireNamespace("rnaturalearth", quietly = TRUE)) {
  costa <- tryCatch({
    esc <- if (requireNamespace("rnaturalearthhires", quietly = TRUE)) "large" else "medium"
    sf::sf_use_s2(FALSE)
    cc <- rnaturalearth::ne_countries(scale = esc, returnclass = "sf")
    cc <- suppressWarnings(sf::st_crop(sf::st_make_valid(cc),
                                       xmin = -9, xmax = -4, ymin = 34, ymax = 39))
    sf::sf_use_s2(TRUE)
    cc
  }, error = function(e) NULL)
}

# --- Universo: todas las celdas finas de cada C-square con esfuerzo ----------
U <- unidades[, {
  kk <- k[1]
  CJ(ix = ux * kk + 0:(kk - 1L), iy = uy * kk + 0:(kk - 1L))
}, by = .(LE_MET4, ux, uy, csq, lon_u, lat_u, k, rf)]
U[, `:=`(lon = (ix + 0.5) * rf, lat = (iy + 0.5) * rf)]

huella <- unique(cnt[, .(LE_MET4, ix, iy)])[, huella := TRUE]
U <- merge(U, huella, by = c("LE_MET4", "ix", "iy"), all.x = TRUE)
U[is.na(huella), huella := FALSE]

U[, prof := if (!is.null(bati)) -prof_de(lon, lat) else NA_real_]
if (!is.null(costa)) {
  sf::sf_use_s2(FALSE)
  costa_u <- sf::st_union(costa)
  pts <- sf::st_as_sf(U[, .(lon, lat)], coords = c("lon", "lat"), crs = 4326)
  U[, tierra := lengths(sf::st_intersects(pts, costa_u)) > 0]
  sf::sf_use_s2(TRUE)
  rm(pts)
} else {
  U[, tierra := !is.na(prof) & prof <= 0]
}
U[, pmax_arte := prof_max[LE_MET4]]
U[, mar := !tierra & (is.na(prof) | prof <= pmax_arte)]
U[, en_universo := if (universo == "huella") huella else (huella | mar)]

N_univ <- U[en_universo == TRUE, .(N_univ = .N), by = .(LE_MET4, ux, uy)]
escribir(U[en_universo == TRUE, .(LE_MET4, csq, ix, iy, lon, lat, prof, huella)], "universo_celdas")

# =============================================================================
# 4. EC. 1 — NB TRUNCADA EN CERO (NB2 o NB1) CON BETA ALEATORIO POR AÑO (TMB)
# =============================================================================
message("4. Ajuste NB truncada (Ec. 1) por C-square...")

# Tabla C-square × año
CA <- cnt[, .(n_pos = .N, E = sum(y)), by = .(LE_MET4, ux, uy, csq, lon_u, lat_u, ANYO)]
CA <- merge(CA, N_univ, by = c("LE_MET4", "ux", "uy"))
CA[, esfuerzo_medio := E / N_univ]   # mu de la Ec. 2: intensidad media en el universo
CA[, ajustable := n_pos >= min_celdas_pos]

# --- Plantilla TMB -----------------------------------------------------------
tmb_cpp <- '
#include <TMB.hpp>
// Hintzen et al. 2019, Ec. 1 — NB truncada en 0
// nb2 = 1: Var = mu + beta mu^2   (size = 1/beta)    [NB2]
// nb2 = 0: Var = mu (1 + beta)    (size = mu/beta)   [NB1, paper]
// log(beta_t) = mu_logbeta + eps_t ;  eps_t ~ N(0, sd)  (aleatorio por año)
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_IVECTOR(t);
  DATA_INTEGER(re);
  DATA_INTEGER(nb2);
  PARAMETER_VECTOR(log_mu);
  PARAMETER(mu_logbeta);
  PARAMETER(log_sd_logbeta);
  PARAMETER_VECTOR(eps);
  int nt = log_mu.size();
  Type nll = 0;
  vector<Type> log_beta(nt);
  for (int j = 0; j < nt; j++) log_beta(j) = mu_logbeta + eps(j);
  if (re == 1) for (int j = 0; j < nt; j++)
    nll -= dnorm(eps(j), Type(0), exp(log_sd_logbeta), true);
  for (int i = 0; i < y.size(); i++) {
    Type mu   = exp(log_mu(t(i)));
    Type beta = exp(log_beta(t(i)));
    Type var, logp0;
    if (nb2 == 1) {
      var   = mu + beta * mu * mu;
      logp0 = -(Type(1) / beta) * log(Type(1) + beta * mu);
    } else {
      var   = mu * (Type(1) + beta);
      logp0 = -(mu / beta) * log(Type(1) + beta);
    }
    Type ll = dnbinom2(y(i), mu, var, true);
    ll -= logspace_sub(Type(0), logp0);   // truncamiento: / (1 - P(Y=0))
    nll -= ll;
  }
  vector<Type> beta = exp(log_beta);
  vector<Type> mu = exp(log_mu);
  REPORT(beta); REPORT(mu);
  ADREPORT(log_beta);
  ADREPORT(log_mu);
  return nll;
}
'

if (motor == "TMB") {
  dir.create("code/tmb", showWarnings = FALSE, recursive = TRUE)
  # Nombre de DLL con versión: si en la sesión de R quedó cargada una versión
  # anterior de la plantilla, R reutilizaría la vieja aunque se recompile.
  dll_tmb  <- "hintzen_tnb_v2"
  cpp_file <- file.path("code/tmb", paste0(dll_tmb, ".cpp"))
  writeLines(tmb_cpp, cpp_file)
  for (d in c("hintzen_tnb1", dll_tmb))
    if (d %in% names(getLoadedDLLs())) try(dyn.unload(TMB::dynlib(file.path("code/tmb", d))), silent = TRUE)
  ok_tmb <- tryCatch({
    TMB::compile(cpp_file)
    dyn.load(TMB::dynlib(file.path("code/tmb", dll_tmb)))
    TRUE
  }, error = function(e) { message("   No compila TMB: ", conditionMessage(e)); FALSE })
  if (!ok_tmb) {
    motor <- "glmmTMB"
    message("   -> se usa motor glmmTMB (beta fijo por año)")
  }
}

ajustar_tmb <- function(y, ti, nt, re) {
  m0 <- as.numeric(tapply(y, factor(ti, levels = 1:nt), mean))
  v0 <- as.numeric(tapply(y, factor(ti, levels = 1:nt), var))
  b_ini <- if (familia == "NB2") (v0 - m0) / m0^2 else v0 / m0 - 1
  b0 <- log(max(median(b_ini, na.rm = TRUE), 0.05))
  dat <- list(y = as.numeric(y), t = as.integer(ti - 1L), re = as.integer(re),
              nb2 = as.integer(familia == "NB2"))
  if (re) {
    par <- list(log_mu = log(pmax(m0 - 0.5, 0.1)), mu_logbeta = b0,
                log_sd_logbeta = log(0.3), eps = rep(0, nt))
    obj <- MakeADFun(dat, par, random = "eps", DLL = dll_tmb, silent = TRUE)
  } else {
    par <- list(log_mu = log(pmax(m0 - 0.5, 0.1)), mu_logbeta = 0,
                log_sd_logbeta = 0, eps = rep(b0, nt))
    obj <- MakeADFun(dat, par, map = list(mu_logbeta = factor(NA), log_sd_logbeta = factor(NA)),
                     DLL = dll_tmb, silent = TRUE)
  }
  opt <- suppressWarnings(nlminb(obj$par, obj$fn, obj$gr,
                                 control = list(iter.max = 1000, eval.max = 2000)))
  if (opt$convergence != 0 || !is.finite(opt$objective)) stop("no converge")
  sdr <- suppressWarnings(sdreport(obj))
  s <- summary(sdr, "report")
  lb <- s[rownames(s) == "log_beta", , drop = FALSE]
  lm_ <- s[rownames(s) == "log_mu", , drop = FALSE]
  sd_re <- if (re) exp(opt$par[["log_sd_logbeta"]]) else NA_real_
  data.table(t = 1:nt, mu = exp(lm_[, 1]), log_beta = lb[, 1], se_log_beta = lb[, 2],
             sd_re = sd_re, re = re, nll = opt$objective)
}

ajustar_glmmtmb <- function(y, ti, nt) {
  d <- data.frame(y = y, fy = factor(ti, levels = 1:nt))
  f  <- if (nt > 1) y ~ 0 + fy else y ~ 1
  fd <- if (nt > 1) ~ 0 + fy else ~ 1
  fam <- if (familia == "NB2") glmmTMB::truncated_nbinom2() else glmmTMB::truncated_nbinom1()
  m <- glmmTMB::glmmTMB(f, dispformula = fd, family = fam, data = d)
  cf <- summary(m)$coefficients
  cm <- cf$cond; cd <- cf$disp
  # NB1: disp = log(beta) ; NB2: disp = log(theta) -> log(beta) = -log(theta)
  lb <- if (familia == "NB2") -cd[, 1] else cd[, 1]
  data.table(t = 1:nt, mu = exp(cm[, 1]), log_beta = lb, se_log_beta = cd[, 2],
             sd_re = NA_real_, re = 0L, nll = -as.numeric(logLik(m)))
}

ajustes <- list()
grupos <- unique(CA[ajustable == TRUE, .(LE_MET4, ux, uy)])
pb <- txtProgressBar(min = 0, max = max(1, nrow(grupos)), style = 3)
for (g in seq_len(nrow(grupos))) {
  gi <- grupos[g]
  yrs <- CA[gi, on = .(LE_MET4, ux, uy)][ajustable == TRUE, sort(ANYO)]
  d <- cnt[gi, on = .(LE_MET4, ux, uy)][ANYO %in% yrs]
  ti <- match(d$ANYO, yrs); nt <- length(yrs)
  res <- NULL
  if (motor == "TMB") {
    res <- tryCatch(ajustar_tmb(d$y, ti, nt, re = nt >= 3), error = function(e) NULL)
    if (is.null(res) && nt >= 3)
      res <- tryCatch(ajustar_tmb(d$y, ti, nt, re = FALSE), error = function(e) NULL)
  } else {
    res <- tryCatch(ajustar_glmmtmb(d$y, ti, nt), error = function(e) NULL)
  }
  if (!is.null(res)) {
    res[, `:=`(LE_MET4 = gi$LE_MET4, ux = gi$ux, uy = gi$uy, ANYO = yrs[t])]
    ajustes[[length(ajustes) + 1]] <- res
  }
  setTxtProgressBar(pb, g)
}
close(pb)
AJ <- rbindlist(ajustes)
# Control: en el C-square-año con más celdas, beta de TMB debe coincidir con glmmTMB
if (nrow(AJ) && requireNamespace("glmmTMB", quietly = TRUE)) {
  ck <- CA[ajustable == TRUE][which.max(n_pos)]
  yck <- cnt[LE_MET4 == ck$LE_MET4 & ux == ck$ux & uy == ck$uy & ANYO == ck$ANYO, y]
  fck <- if (familia == "NB2") glmmTMB::truncated_nbinom2() else glmmTMB::truncated_nbinom1()
  mck <- glmmTMB::glmmTMB(y ~ 1, family = fck, data = data.frame(y = yck))
  b_glmm <- if (familia == "NB2") 1 / sigma(mck) else sigma(mck)
  b_tmb  <- exp(AJ[LE_MET4 == ck$LE_MET4 & ux == ck$ux & uy == ck$uy & ANYO == ck$ANYO, log_beta])
  message(sprintf("   Control %s %s %d: beta TMB = %.3f | glmmTMB (sin efecto aleatorio) = %.3f",
                  ck$LE_MET4, ck$csq, ck$ANYO, b_tmb, b_glmm))
  if (abs(log(b_tmb / b_glmm)) > log(2))
    warning("beta de TMB difiere > 2x de glmmTMB: reiniciar R (Session > Restart R) y volver a correr")
}
message(sprintf("   C-square-año ajustados: %s de %s (motor %s)",
                nrow(AJ), nrow(CA), motor))

# Ajustes fallidos o en el límite: beta no finito, enorme (no converge) o ~0
# (el optimizador llega al borde: no hay información para estimar la dispersión)
AJ[, valido := is.finite(log_beta) & is.finite(mu) &
                exp(log_beta) < beta_max_valido & exp(log_beta) > beta_min_valido]
if (any(!AJ$valido)) {
  message(sprintf("   ajustes descartados (beta fuera de [%g, %g] o no finito): %d de %d",
                  beta_min_valido, beta_max_valido, sum(!AJ$valido), nrow(AJ)))
  escribir(AJ[valido == FALSE], "ajustes_fallidos")
}
AJ[valido == FALSE, `:=`(log_beta = NA_real_, se_log_beta = NA_real_, mu = NA_real_)]

CA <- merge(CA, AJ[, .(LE_MET4, ux, uy, ANYO, mu, log_beta, se_log_beta, sd_re, re)],
            by = c("LE_MET4", "ux", "uy", "ANYO"), all.x = TRUE)
CA[, `:=`(beta    = exp(log_beta),
          beta_lo = exp(log_beta - 1.96 * se_log_beta),
          beta_hi = exp(log_beta + 1.96 * se_log_beta))]

# =============================================================================
# 5. EC. 2 — CEROS EN EXCESO, ZONAS INACCESIBLES Y DESFAVORABLES
# =============================================================================
message("5. Ceros en exceso (Ec. 2)...")

# P(Y=0):  NB2 (1 + beta*mu)^(-1/beta)  |  NB1 (1 + beta)^(-mu/beta)
CA[, P0 := p0_nb(mu, beta)]
# fracción de celdas que pertenecen al proceso NB:  n_pos = pi * N * (1 - P0)
CA[, pi_nb := n_pos / (N_univ * (1 - P0))]
CA[, P_false := pmin(1, pmax(0, 1 - pi_nb))]
CA[, ceros_obs := (N_univ - n_pos) / N_univ]

# Ec. 2 por C-square: logit(P_false) = a + b * esfuerzo_medio
ec2 <- CA[!is.na(P_false), {
  if (.N >= 3) {
    z <- round(P_false * N_univ)
    m <- tryCatch(suppressWarnings(glm(cbind(z, N_univ - z) ~ esfuerzo_medio, family = binomial)),
                  error = function(e) NULL)
    if (is.null(m)) list(ec2_a = NA_real_, ec2_b = NA_real_, ec2_p_b = NA_real_)
    else list(ec2_a = coef(m)[1], ec2_b = coef(m)[2],
              ec2_p_b = summary(m)$coefficients[2, 4])
  } else list(ec2_a = NA_real_, ec2_b = NA_real_, ec2_p_b = NA_real_)
}, by = .(LE_MET4, ux, uy)]

# Inaccesible = mínimo de P_false en la serie (paper) o segundo mínimo (robusto
# a un único año mal ajustado); desfavorable = media - inaccesible
clasif <- CA[!is.na(P_false), {
  pf <- sort(P_false)
  pin <- if (inacc_regla == "segundo_minimo" && length(pf) >= 3) pf[2] else pf[1]
  list(n_anios_ajuste = .N, P_inacc = pin, P_desfav = max(0, mean(P_false) - pin))
}, by = .(LE_MET4, ux, uy)]
clasif[, P_fav := 1 - P_inacc - P_desfav]
clasif <- merge(clasif, ec2, by = c("LE_MET4", "ux", "uy"), all.x = TRUE)

# =============================================================================
# 6. VARIABILIDAD INTERANUAL POR CELDA (THETA) + TEST DE PERMUTACIÓN
# =============================================================================
message("6. Persistencia interanual (theta) y permutaciones...")

# theta NB por celda con offset = esfuerzo total del C-square en el año.
# E[y_it] = s_i * E_t ;  Var = m + m^2/theta.  Estimador de momentos
# (rápido, permite permutar; theta alto = poca variación entre años).
theta_celdas <- function(Y, E) {
  s <- rowSums(Y) / sum(E)
  Mu <- outer(s, E)
  num <- rowSums(Mu^2)
  den <- rowSums((Y - Mu)^2 - Mu)
  th <- ifelse(den <= 0, 100, num / den)
  pmin(pmax(th, 0.01), 100)
}

persist <- list()
theta_cel <- list()
grp_t <- CA[, .(n_anios = .N), by = .(LE_MET4, ux, uy)][n_anios >= min_anios_theta]
for (g in seq_len(nrow(grp_t))) {
  gi <- grp_t[g]
  d  <- cnt[gi, on = .(LE_MET4, ux, uy)]
  cel <- U[gi, on = .(LE_MET4, ux, uy)][en_universo == TRUE, .(ix, iy)]
  yrs <- sort(unique(d$ANYO))
  Y <- matrix(0, nrow(cel), length(yrs))
  pos <- cel[d, on = .(ix, iy), which = TRUE]
  Y[cbind(pos, match(d$ANYO, yrs))] <- d$y
  E <- colSums(Y)
  activas <- rowSums(Y) > 0
  if (sum(activas) < 5) next
  Ya <- Y[activas, , drop = FALSE]
  th_obs <- theta_celdas(Ya, E)
  stat_obs <- mean(log(th_obs))
  # Permutación: cada año el esfuerzo se reparte al azar entre las celdas de
  # la huella del C-square (celdas pescadas algún año) -> elección aleatoria
  # de la localización, sin mover esfuerzo a celdas nunca pescadas
  stat_perm <- replicate(n_perm, {
    Yp <- apply(Ya, 2, sample)
    a  <- rowSums(Yp) > 0
    mean(log(theta_celdas(Yp[a, , drop = FALSE], E)))
  })
  persist[[length(persist) + 1]] <- data.table(
    LE_MET4 = gi$LE_MET4, ux = gi$ux, uy = gi$uy,
    n_celdas_activas = sum(activas), n_anios = length(yrs),
    theta_media = exp(stat_obs),
    theta_perm_media = exp(mean(stat_perm)),
    p_valor = mean(stat_perm >= stat_obs)   # p bajo = más estable que el azar
  )
  theta_cel[[length(theta_cel) + 1]] <- data.table(
    LE_MET4 = gi$LE_MET4, ux = gi$ux, uy = gi$uy,
    ix = cel$ix[activas], iy = cel$iy[activas], theta = th_obs)
}
PER <- rbindlist(persist)
TH  <- rbindlist(theta_cel)

# =============================================================================
# 7. COVARIABLES POR C-SQUARE: BATIMETRÍA (BPI) Y HÁBITAT
# =============================================================================
message("7. Covariables batimétricas y de hábitat...")

cov_u <- unique(unidades[, .(LE_MET4, ux, uy, csq, lon_u, lat_u)])
if (!is.null(bati)) {
  bati[, `:=`(ux = as.integer(floor(x / res_unidad + 1e-9)),
              uy = as.integer(floor(y / res_unidad + 1e-9)))]
  bstat <- bati[z < 0, .(prof_media = mean(-z), prof_sd = sd(z),
                         var_bpi = var(bpi), n_bati = .N), by = .(ux, uy)]
  cov_u <- merge(cov_u, bstat, by = c("ux", "uy"), all.x = TRUE)
}
# profundidad desde el universo de celdas (respaldo)
cov_u <- merge(cov_u, U[en_universo == TRUE & !is.na(prof) & prof > 0,
                        .(prof_sd_celdas = sd(prof)), by = .(LE_MET4, ux, uy)],
               by = c("LE_MET4", "ux", "uy"), all.x = TRUE)

hab_u <- NULL
if (usar_habitat && dir.exists(dir_habitat)) {
  hab_u <- tryCatch({
    shp <- list.files(dir_habitat, pattern = "\\.shp$", full.names = TRUE)
    sf::sf_use_s2(FALSE)
    hab <- do.call(rbind, lapply(shp, function(f) {
      h <- sf::read_sf(f)
      h <- sf::st_transform(h[, "MSFD_BBHT"], 4326)
      h
    }))
    hab <- sf::st_make_valid(hab)
    pts <- sf::st_as_sf(U[en_universo == TRUE, .(LE_MET4, ux, uy, lon, lat)],
                        coords = c("lon", "lat"), crs = 4326, remove = FALSE)
    j <- sf::st_join(pts, hab, join = sf::st_within, left = TRUE)
    sf::sf_use_s2(TRUE)
    j <- as.data.table(sf::st_drop_geometry(j))
    j <- j[!duplicated(j[, .(LE_MET4, ux, uy, lon, lat)])]
    j[, sustrato := fcase(
      grepl("sand", MSFD_BBHT, ignore.case = TRUE), "Arena",
      grepl("mud", MSFD_BBHT, ignore.case = TRUE), "Fango",
      grepl("coarse|mixed", MSFD_BBHT, ignore.case = TRUE), "Grueso/mixto",
      grepl("rock|reef|biogenic", MSFD_BBHT, ignore.case = TRUE), "Roca/arrecife",
      default = "Otro/sin dato")]
    dcast(j[, .N, by = .(LE_MET4, ux, uy, sustrato)][, p := N / sum(N), by = .(LE_MET4, ux, uy)],
          LE_MET4 + ux + uy ~ sustrato, value.var = "p", fill = 0)
  }, error = function(e) { message("   hábitat no disponible: ", conditionMessage(e)); NULL })
}

# =============================================================================
# 8. TABLAS DE RESULTADOS
# =============================================================================
message("8. Guardando tablas...")

tend <- CA[!is.na(beta), {
  if (.N >= 3) {
    m <- lm(beta ~ ANYO)
    list(pend_beta = coef(m)[2], p_pend_beta = suppressWarnings(summary(m)$coefficients[2, 4]))
  } else list(pend_beta = NA_real_, p_pend_beta = NA_real_)
}, by = .(LE_MET4, ux, uy)]

RES <- CA[, .(beta_media = mean(beta, na.rm = TRUE),
              beta_mediana = median(beta, na.rm = TRUE),
              esfuerzo_total = sum(E),
              n_anios_esfuerzo = .N), by = .(LE_MET4, ux, uy, csq, lon_u, lat_u)]
RES <- merge(RES, N_univ, by = c("LE_MET4", "ux", "uy"))
RES <- merge(RES, tend,   by = c("LE_MET4", "ux", "uy"), all.x = TRUE)
RES <- merge(RES, clasif, by = c("LE_MET4", "ux", "uy"), all.x = TRUE)
RES <- merge(RES, PER,    by = c("LE_MET4", "ux", "uy"), all.x = TRUE)
RES <- merge(RES, cov_u[, !c("csq", "lon_u", "lat_u")], by = c("LE_MET4", "ux", "uy"), all.x = TRUE)
if (!is.null(hab_u)) RES <- merge(RES, hab_u, by = c("LE_MET4", "ux", "uy"), all.x = TRUE)
RES[, beta_media := fifelse(is.nan(beta_media), NA_real_, beta_media)]
RES[, estable := !is.na(p_valor) & p_valor < 0.05]

escribir(CA[order(LE_MET4, csq, ANYO)], "beta_csq_anio")
escribir(RES[order(LE_MET4, -esfuerzo_total)], "csq_resumen")
fst::write_fst(TH, paste0(pref_out, "theta_celdas.fst"))

# Resumen por arte (cifras equivalentes a las del paper)
RA <- RES[, .(
  n_csq              = .N,
  n_csq_ajustados    = sum(!is.na(beta_media)),
  pct_esf_ajustado   = 100 * sum(esfuerzo_total[!is.na(beta_media)]) / sum(esfuerzo_total),
  beta_mediana       = median(beta_media, na.rm = TRUE),
  beta_min           = suppressWarnings(min(beta_media, na.rm = TRUE)),
  beta_max           = suppressWarnings(max(beta_media, na.rm = TRUE)),
  n_tendencia_sig    = sum(p_pend_beta < 0.05, na.rm = TRUE),
  P_inacc_media      = weighted.mean(P_inacc, N_univ, na.rm = TRUE),
  P_desfav_media     = weighted.mean(P_desfav, N_univ, na.rm = TRUE),
  n_csq_test_theta   = sum(!is.na(p_valor)),
  n_csq_estables     = sum(estable),
  pct_esf_estable    = 100 * sum(esfuerzo_total[estable]) /
                       sum(esfuerzo_total[!is.na(p_valor)])
), by = LE_MET4]
escribir(RA, "resumen_arte")
print(RA)

# =============================================================================
# 9. FIGURAS EMULANDO EL PAPER
# =============================================================================
message("9. Figuras...")

lims <- function(a, marg = 0.06) {
  u <- unidades[LE_MET4 == a]
  list(x = range(u$lon_u) + c(-1, 1) * (res_unidad / 2 + marg),
       y = range(u$lat_u) + c(-1, 1) * (res_unidad / 2 + marg))
}
mapa <- function(a, ...) {
  l <- lims(a)
  list(
    if (!is.null(costa)) geom_sf(data = costa, fill = "grey85", colour = "grey45",
                                 linewidth = 0.2, inherit.aes = FALSE),
    coord_sf(xlim = l$x, ylim = l$y, expand = FALSE),
    labs(x = NULL, y = NULL),
    tema_mapa
  )
}
# mini-serie temporal dentro de cada C-square (como Fig. 2 y 5a del paper)
mini_serie <- function(dt, var, escala_max) {
  r <- range(anios); w <- res_unidad * 0.8
  dt[, `:=`(xs = lon_u - w / 2 + (ANYO - r[1]) / diff(r) * w,
            ys = lat_u - w / 2 + pmin(get(var) / escala_max, 1) * w)]
  dt
}

for (a in artes_sel) {
  lab <- artes_label[[a]]
  ca_a <- CA[LE_MET4 == a]; res_a <- RES[LE_MET4 == a]

  # --- Fig. 1: conteos por celda fina (suma 2018-2024) -----------------------
  c1 <- cnt[LE_MET4 == a, .(y = sum(y)), by = .(lon, lat)]
  grid_u <- unidades[LE_MET4 == a]
  p1 <- ggplot() +
    geom_tile(data = c1, aes(lon, lat, fill = y), width = res_fina[[a]], height = res_fina[[a]]) +
    geom_tile(data = grid_u, aes(lon_u, lat_u), width = res_unidad, height = res_unidad,
              fill = NA, colour = "grey30", linewidth = 0.15) +
    scale_fill_viridis_c(option = "inferno", trans = "log10", direction = -1,
                         name = sprintf("Nº de %s (2018-2024)", lab_conteo(a))) +
    mapa(a) +
    labs(title = paste("Distribución fina del esfuerzo —", lab))
  guardar(p1, paste0("fig1_conteos_mapa_", a), 10, 7)

  # Zoom a los 4 C-squares con más esfuerzo (insets de la Fig. 1)
  top4 <- res_a[order(-esfuerzo_total)][1:min(4, .N), .(ux, uy, csq)]
  z1 <- cnt[LE_MET4 == a][top4, on = .(ux, uy)][, .(y = sum(y)), by = .(csq, lon, lat)]
  p1z <- ggplot(z1, aes(lon, lat, fill = y)) +
    geom_tile(width = res_fina[[a]], height = res_fina[[a]]) +
    scale_fill_viridis_c(option = "inferno", trans = "log10", direction = -1,
                         name = sprintf("Nº de %s", lab_conteo(a))) +
    facet_wrap(~ csq, scales = "free") + tema +
    theme(axis.text = element_text(size = 6)) +
    labs(title = paste("Zoom a escala fina — 4 C-squares con más esfuerzo —", lab), x = NULL, y = NULL)
  guardar(p1z, paste0("fig1_zoom_csq_", a), 10, 9)

  # --- Fig. 2: beta por C-square-año (mapa + mini-series) --------------------
  if (any(!is.na(ca_a$beta))) {
    bmax <- quantile(ca_a$beta, 0.95, na.rm = TRUE)
    ms <- mini_serie(ca_a[!is.na(beta)], "beta", bmax)
    p2 <- ggplot() +
      geom_tile(data = res_a[!is.na(beta_media)], aes(lon_u, lat_u, fill = beta_media),
                width = res_unidad, height = res_unidad, colour = "white") +
      geom_tile(data = res_a[p_pend_beta < 0.05], aes(lon_u, lat_u),
                width = res_unidad, height = res_unidad, fill = NA, colour = "black", linewidth = 0.6) +
      geom_path(data = ms, aes(xs, ys, group = csq), linewidth = 0.25) +
      geom_point(data = ms, aes(xs, ys), size = 0.3) +
      scale_fill_distiller(palette = "YlOrRd", direction = 1, trans = "log10",
                           name = paste("Agregación media,", lab_beta, "(log)")) +
      mapa(a) +
      labs(title = paste("Agregación del esfuerzo (β) por C-square —", lab))
    guardar(p2, paste0("fig2_beta_mapa_", a), 10, 7)
  }
}

# Fig. 2b: heatmap beta C-square × año
if (any(!is.na(CA$beta))) {
  hb <- CA[!is.na(beta)]
  hb[, csq_f := factor(csq, levels = unique(hb[order(lat_u, lon_u), csq]))]
  p2b <- ggplot(hb, aes(factor(ANYO), csq_f, fill = beta)) +
    geom_tile(colour = "white") +
    scale_fill_distiller(palette = "YlOrRd", direction = 1, trans = "log10", name = paste(lab_beta, "(log)")) +
    facet_wrap(~ LE_MET4, scales = "free_y", labeller = as_labeller(artes_label)) +
    tema + theme(axis.text.y = element_text(size = 5)) +
    labs(title = "Agregación (β) por C-square y año", x = "Año", y = "C-square (ordenados S→N)")
  guardar(p2b, "fig2b_beta_heatmap", 12, 10)
}

# --- Fig. 3a: agregado vs aleatorio (C-square-año con más esfuerzo) ----------
# --- Fig. 3b: observado vs ajustado (6 C-square-año con más esfuerzo) -------
for (a in artes_sel) {
  sel <- CA[LE_MET4 == a & !is.na(beta)][order(-E)]
  if (!nrow(sel)) next
  s1 <- sel[1]
  d1 <- cnt[LE_MET4 == a & ux == s1$ux & uy == s1$uy & ANYO == s1$ANYO]
  ymax <- max(d1$y)                          # histograma completo, sin acumular cola
  obs <- c(s1$N_univ - nrow(d1), tabulate(d1$y, nbins = ymax))
  lam <- s1$E / s1$N_univ
  poi <- dpois(0:ymax, lam)
  f3a <- rbind(data.table(y = 0:ymax, prop = obs / sum(obs), tipo = "Observado (agregado)"),
               data.table(y = 0:ymax, prop = poi, tipo = "Aleatorio (Poisson, mismo esfuerzo)"))
  p3a <- ggplot(f3a, aes(y, prop, fill = tipo)) +
    geom_col(position = "dodge", width = 0.9) +
    scale_fill_manual(values = c("grey60", "black"), name = NULL) +
    scale_y_sqrt() +
    tema + labs(title = paste("Frecuencia de intensidad: agregado vs aleatorio —", artes_label[[a]]),
                x = sprintf("Nº de %s por celda", lab_conteo(a)), y = "Proporción de celdas")
  guardar(p3a, paste0("fig3a_agregado_vs_aleatorio_", a), 9, 6)

  top6 <- sel[1:min(6, .N)]
  f3b <- rbindlist(lapply(seq_len(nrow(top6)), function(i) {
    s <- top6[i]
    d <- cnt[LE_MET4 == a & ux == s$ux & uy == s$uy & ANYO == s$ANYO]
    ym <- max(d$y)                             # histograma completo, sin acumular cola
    o <- tabulate(d$y, nbins = ym) / nrow(d)
    f <- dnb_tr(1:ym, s$mu, s$beta)
    data.table(panel = sprintf("%s · %d · β=%.2f", s$csq, s$ANYO, s$beta),
               y = 1:ym, Observado = o, Ajustado = f)
  }))
  p3b <- ggplot(f3b, aes(y)) +
    geom_col(aes(y = Observado), fill = "black", width = 1) +
    geom_line(aes(y = Ajustado), colour = "red") +
    scale_y_sqrt() +
    facet_wrap(~ panel, scales = "free") + tema +
    labs(title = paste0(familia, " truncada: observado (barras) vs ajustado (rojo) — ", artes_label[[a]]),
         x = sprintf("Nº de %s por celda (> 0)", lab_conteo(a)), y = "Proporción")
  guardar(p3b, paste0("fig3b_observado_vs_ajustado_", a), 11, 7)
}

# --- Fig. 4: beta vs heterogeneidad batimétrica ------------------------------
if ("var_bpi" %in% names(RES)) {
  f4 <- melt(RES[!is.na(beta_media)],
             id.vars = c("LE_MET4", "csq", "beta_media"),
             measure.vars = intersect(c("var_bpi", "prof_sd", "prof_media"), names(RES)),
             variable.name = "covariable", value.name = "valor")
  f4[, covariable := factor(covariable, levels = c("var_bpi", "prof_sd", "prof_media"),
                            labels = c("Varianza BPI-5 km", "DE profundidad (m)", "Profundidad media (m)"))]
  p4 <- ggplot(f4[!is.na(valor)], aes(valor, beta_media, colour = LE_MET4)) +
    geom_point(size = 1.6) + scale_y_log10() +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span = 1) +
    scale_colour_manual(values = col_artes, labels = artes_label, name = NULL) +
    facet_grid(LE_MET4 ~ covariable, scales = "free") + tema +
    labs(title = "Agregación (β) vs heterogeneidad batimétrica por C-square",
         x = NULL, y = "Agregación (β)")
  guardar(p4, "fig4_beta_vs_batimetria", 11, 7)
}

# --- Fig. 5: P_false (inaccesible + desfavorable) -----------------------------
for (a in artes_sel) {
  ca_a <- CA[LE_MET4 == a & !is.na(P_false)]; res_a <- RES[LE_MET4 == a & !is.na(P_inacc)]
  if (!nrow(ca_a)) next
  ms <- mini_serie(copy(ca_a), "P_false", 1)
  pf <- ca_a[, .(P_false_media = mean(P_false)), by = .(csq, lon_u, lat_u)]
  p5a <- ggplot() +
    geom_tile(data = pf, aes(lon_u, lat_u, fill = P_false_media),
              width = res_unidad, height = res_unidad, colour = "white") +
    geom_path(data = ms, aes(xs, ys, group = csq), linewidth = 0.25) +
    geom_point(data = ms, aes(xs, ys), size = 0.3) +
    scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0, 1),
                         name = "Proporción inaccesible | desfavorable") +
    mapa(a) + labs(title = paste("Celdas no pescadas no explicadas por la NB (P_false) —", artes_label[[a]]))
  guardar(p5a, paste0("fig5a_pfalse_mapa_", a), 10, 7)

  f5b <- melt(res_a, id.vars = c("csq", "lon_u", "lat_u"),
              measure.vars = c("P_fav", "P_desfav", "P_inacc"),
              variable.name = "tipo", value.name = "prop")
  f5b[, tipo := factor(tipo, levels = c("P_fav", "P_desfav", "P_inacc"),
                       labels = c("Favorable", "Desfavorable", "Inaccesible"))]
  p5b <- ggplot() +
    geom_tile(data = f5b, aes(lon_u, lat_u, fill = prop), width = res_unidad, height = res_unidad,
              colour = "white") +
    scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, 1), name = "Proporción") +
    facet_wrap(~ tipo, ncol = 1) + mapa(a) +
    labs(title = paste("Proporción favorable / desfavorable / inaccesible —", artes_label[[a]]))
  guardar(p5b, paste0("fig5b_categorias_mapa_", a), 9, 12)

  res_a[, csq_f := factor(csq, levels = res_a[order(-P_inacc - P_desfav), csq])]
  f5c <- melt(res_a, id.vars = "csq_f", measure.vars = c("P_fav", "P_desfav", "P_inacc"),
              variable.name = "tipo", value.name = "prop")
  f5c[, tipo := factor(tipo, levels = c("P_fav", "P_desfav", "P_inacc"),
                       labels = c("Favorable", "Desfavorable", "Inaccesible"))]
  p5c <- ggplot(f5c, aes(csq_f, prop, fill = tipo)) + geom_col(width = 1) +
    scale_fill_manual(values = c(Favorable = "white", Desfavorable = "orange", Inaccesible = "red3"),
                      name = NULL) +
    tema + theme(axis.text.x = element_text(angle = 90, size = 5, vjust = 0.5),
                 panel.background = element_rect(fill = "grey95")) +
    labs(title = paste("Composición del universo de celdas por C-square —", artes_label[[a]]),
         x = "C-square", y = "Proporción de celdas")
  guardar(p5c, paste0("fig5c_categorias_barras_", a), 12, 6)
}

# --- Fig. 6: % de ceros vs esfuerzo ------------------------------------------
b_ref <- median(CA$beta, na.rm = TRUE)
mu_x  <- seq(0, 5, by = 0.05)
p0_x  <- p0_nb(pmax(mu_x, 1e-6), b_ref)
f6a <- rbind(
  data.table(mu = mu_x, ceros = 100 * p0_x, tipo = "Totalmente accesible"),
  data.table(mu = mu_x, ceros = 100 * (0.4 + 0.6 * p0_x), tipo = "Parcialmente inaccesible"),
  data.table(mu = mu_x, ceros = 100 * (1 - (1 - p0_x) * (1 - 0.6 * exp(-0.4 * mu_x))),
             tipo = "Parcialmente desfavorable"))
p6a <- ggplot(f6a, aes(mu, ceros, linetype = tipo)) + geom_line() + tema +
  labs(title = "Esquema: % de celdas con cero vs esfuerzo (Fig. 6a)",
       x = "Esfuerzo (conteo medio por celda)", y = "% de celdas con cero", linetype = NULL)
guardar(p6a, "fig6a_esquema_ceros", 8, 6)

f6b <- CA[!is.na(P_false)]
f6b <- merge(f6b, clasif[, .(LE_MET4, ux, uy, P_inacc, P_desfav, n_anios_ajuste, ec2_a, ec2_b)],
             by = c("LE_MET4", "ux", "uy"))
f6b <- f6b[n_anios_ajuste >= 4]
if (nrow(f6b)) {
  ej <- f6b[, .(P_inacc = P_inacc[1], P_desfav = P_desfav[1]), by = .(LE_MET4, csq)]
  ej <- ej[, .SD[c(which.min(P_inacc + P_desfav), which.max(P_inacc), which.max(P_desfav))],
           by = LE_MET4]
  ej[, tipo := rep(c("Accesible", "Inaccesible", "Desfavorable"), length.out = .N), by = LE_MET4]
  ej <- unique(ej, by = c("LE_MET4", "csq"))
  f6 <- merge(f6b, ej[, .(LE_MET4, csq, tipo)], by = c("LE_MET4", "csq"))
  f6[, panel := paste(LE_MET4, "·", tipo, "·", csq)]
  curv <- f6[!is.na(ec2_a), {
    xs <- seq(min(esfuerzo_medio), max(esfuerzo_medio), length.out = 50)
    list(esfuerzo_medio = xs, P_false = plogis(ec2_a[1] + ec2_b[1] * xs))
  }, by = panel]
  p6b <- ggplot(f6, aes(esfuerzo_medio, P_false)) +
    geom_line(data = curv, linetype = "dashed") +
    geom_point() +
    geom_hline(data = unique(f6[, .(panel, P_inacc)]), aes(yintercept = P_inacc),
               colour = "red", linetype = "dotted") +
    facet_wrap(~ panel, scales = "free_x") + ylim(0, 1) + tema +
    labs(title = "Ceros en exceso (P_false) vs esfuerzo — Ec. 2",
         x = "Esfuerzo (conteo medio por celda del universo)", y = expression(P[false]))
  guardar(p6b, "fig6b_pfalse_vs_esfuerzo", 11, 7)
}

# --- Fig. 7: inaccesible / desfavorable vs sustrato ---------------------------
if (!is.null(hab_u)) {
  sust <- setdiff(names(hab_u), c("LE_MET4", "ux", "uy"))
  sust <- intersect(sust, c("Arena", "Fango", "Grueso/mixto", "Roca/arrecife"))
  if (length(sust)) {
    f7 <- melt(RES[!is.na(P_inacc)], id.vars = c("LE_MET4", "csq", "P_inacc", "P_desfav"),
               measure.vars = sust, variable.name = "sustrato", value.name = "p_sustrato")
    f7 <- melt(f7, measure.vars = c("P_inacc", "P_desfav"), variable.name = "tipo", value.name = "prop")
    f7[, tipo := factor(tipo, labels = c("Inaccesible", "Desfavorable"))]
    p7 <- ggplot(f7, aes(p_sustrato, prop, colour = LE_MET4)) +
      geom_point(size = 1.2) +
      geom_smooth(method = "glm", formula = y ~ x, method.args = list(family = quasibinomial),
                  se = TRUE) +
      scale_colour_manual(values = col_artes, labels = artes_label, name = NULL) +
      facet_grid(sustrato ~ tipo) + tema +
      labs(title = "Proporción inaccesible / desfavorable vs sustrato (MSFD)",
           x = "Proporción de sustrato en el C-square", y = "Proporción de celdas")
    guardar(p7, "fig7_habitat", 9, 10)
  }
}

# --- Fig. 8: variabilidad interanual (theta) ----------------------------------
if (nrow(PER)) {
  for (a in artes_sel) {
    r8 <- RES[LE_MET4 == a & !is.na(p_valor)]
    if (!nrow(r8)) next
    r8[, p_cat := cut(p_valor, c(-Inf, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, Inf),
                      labels = c("< 0,05", "0,05–0,1", "0,1–0,15", "0,15–0,2", "0,2–0,25", "0,25–0,3", "> 0,3"))]
    p8a <- ggplot() +
      geom_tile(data = r8, aes(lon_u, lat_u, fill = p_cat), width = res_unidad, height = res_unidad,
                colour = "white") +
      geom_text(data = r8, aes(lon_u, lat_u, label = sprintf("%.2f", p_valor)), size = 1.8) +
      scale_fill_brewer(palette = "Blues", drop = FALSE, name = "p-valor") +
      mapa(a) +
      labs(title = paste("¿Variabilidad interanual menor que una elección aleatoria? —", artes_label[[a]]))
    guardar(p8a, paste0("fig8a_pvalor_theta_mapa_", a), 10, 7)
  }
  p8b <- ggplot(RES[!is.na(theta_media) & !is.na(beta_media)],
                aes(theta_media, beta_media, colour = LE_MET4)) +
    geom_point(aes(shape = estable), size = 2) + scale_y_log10() +
    geom_smooth(method = "loess", formula = y ~ x, span = 1, se = TRUE) +
    scale_x_log10() +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "Estable (p<0,05)") +
    scale_colour_manual(values = col_artes, labels = artes_label, name = NULL) +
    tema + labs(title = "Agregación (β) vs estabilidad interanual (θ)",
                x = "Variabilidad interanual θ (log)", y = "Agregación (β)")
  guardar(p8b, "fig8b_beta_vs_theta", 9, 7)
}

# =============================================================================
# 10. FIGURAS AD HOC
# =============================================================================
message("10. Figuras ad hoc...")

# --- A1. Persistencia de hotspots a escala fina ------------------------------
hs <- copy(cnt)[order(LE_MET4, ANYO, -y)]
hs[, cum := cumsum(y) / sum(y), by = .(LE_MET4, ANYO)]
hs[, hotspot := (cum - y / sum(y)) < umbral_hotspot, by = .(LE_MET4, ANYO)]
hot <- hs[hotspot == TRUE]
pers_cel <- hot[, .(n_anios_hot = uniqueN(ANYO)), by = .(LE_MET4, ix, iy, lon, lat, csq)]
fst::write_fst(pers_cel, paste0(pref_out, "hotspot_persistencia_celda.fst"))
escribir(pers_cel[, .N, by = .(LE_MET4, n_anios_hot)][order(LE_MET4, n_anios_hot)],
         "hotspot_persistencia_resumen")

for (a in artes_sel) {
  pc <- pers_cel[LE_MET4 == a]
  if (!nrow(pc)) next
  pA1 <- ggplot() +
    geom_tile(data = pc, aes(lon, lat, fill = factor(n_anios_hot)),
              width = res_fina[[a]], height = res_fina[[a]]) +
    scale_fill_viridis_d(option = "plasma", direction = -1, name = "Nº de años como hotspot") +
    mapa(a) +
    labs(title = paste("Persistencia de hotspots a escala fina —", artes_label[[a]]))
  guardar(pA1, paste0("adhoc_hotspot_persistencia_", a), 10, 7)
}

# --- A2. Similitud de Jaccard de hotspots entre años -------------------------
jac <- rbindlist(lapply(artes_sel, function(a) {
  h <- hot[LE_MET4 == a, .(id = paste(ix, iy)), by = ANYO]
  yy <- sort(unique(h$ANYO))
  CJ(a1 = yy, a2 = yy)[, .(LE_MET4 = a,
    jaccard = { s1 <- h[ANYO == a1, id]; s2 <- h[ANYO == a2, id]
                length(intersect(s1, s2)) / length(union(s1, s2)) }), by = .(a1, a2)]
}))
escribir(jac, "hotspot_jaccard")
pA2 <- ggplot(jac, aes(factor(a1), factor(a2), fill = jaccard)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", jaccard)), size = 2.5) +
  scale_fill_distiller(palette = "GnBu", direction = 1, limits = c(0, 1), name = "Jaccard") +
  facet_wrap(~ LE_MET4, labeller = as_labeller(artes_label)) + coord_equal() + tema +
  labs(title = "Solapamiento de hotspots entre años (índice de Jaccard)", x = NULL, y = NULL)
guardar(pA2, "adhoc_hotspot_jaccard", 11, 6)

# --- A3. Curvas de concentración (Lorenz) y Gini -----------------------------
lor <- cnt[order(LE_MET4, ANYO, -y)][, .(
  frac_celdas = seq_len(.N) / .N, frac_esf = cumsum(y) / sum(y)), by = .(LE_MET4, ANYO)]
gini <- cnt[, {
  v <- as.numeric(sort(y)); n <- length(v)      # numérico: evita desbordamiento de enteros
  list(n_celdas = n,
       gini = (2 * sum(seq_len(n) * v) / (n * sum(v))) - (n + 1) / n,
       pct_celdas_50 = 100 * mean(cumsum(sort(y, decreasing = TRUE)) / sum(y) <= 0.5),
       pct_celdas_90 = 100 * mean(cumsum(sort(y, decreasing = TRUE)) / sum(y) <= 0.9))
}, by = .(LE_MET4, ANYO)][order(LE_MET4, ANYO)]
escribir(gini, "concentracion_gini")
pA3 <- ggplot(lor, aes(frac_celdas, frac_esf, colour = factor(ANYO))) +
  geom_line(linewidth = 0.6) + geom_abline(linetype = "dotted") +
  scale_colour_viridis_d(name = "Año") +
  facet_wrap(~ LE_MET4, labeller = as_labeller(artes_label)) + tema +
  labs(title = "Concentración del esfuerzo a escala fina (curvas de Lorenz)",
       x = "Fracción acumulada de celdas pescadas (de mayor a menor)",
       y = "Fracción acumulada del esfuerzo")
guardar(pA3, "adhoc_lorenz", 11, 6)

# --- A4. beta por año y beta vs esfuerzo medio -------------------------------
if (any(!is.na(CA$beta))) {
  pA4 <- ggplot(CA[!is.na(beta)], aes(factor(ANYO), beta, fill = LE_MET4)) +
    geom_boxplot(outlier.size = 0.6) + scale_y_log10() +
    scale_fill_manual(values = col_artes, labels = artes_label, name = NULL) +
    tema + labs(title = "Distribución anual de la agregación (β) entre C-squares",
                x = "Año", y = paste(lab_beta, "(log)"))
  guardar(pA4, "adhoc_beta_anual", 10, 6)

  pA5 <- ggplot(CA[!is.na(beta)], aes(esfuerzo_medio, beta, colour = LE_MET4)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", formula = y ~ x, span = 1) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = col_artes, labels = artes_label, name = NULL) +
    tema + labs(title = "¿Depende la agregación del nivel de esfuerzo?",
                x = "Esfuerzo medio por celda (log)", y = paste(lab_beta, "(log)"))
  guardar(pA5, "adhoc_beta_vs_esfuerzo", 9, 6)
}

}   # ===== fin bloque Hintzen (secciones 2-10) =====

# =============================================================================
# 11. PERSISTENCIA DEL CALADERO — BETA-BINOMIAL DE AÑOS PESCADOS POR CELDA
# =============================================================================
# Para cada celda fina de la huella (pescada al menos un año):
#   k = nº de años en que se pescó (conteo >= min_conteo_anio), de n años con
#       actividad del arte. Como solo hay celdas con k >= 1 -> truncada en 0.
# Modelos:
#   Binomial(n, p)          -> cada año la flota elige celdas al azar
#                              (SIN persistencia: k se reparte "en campana")
#   Beta-binomial(n, mu, rho) -> p varía entre celdas (celdas "siempre" y "casi nunca")
#   rho (correlación intra-celda, 0-1) = ÍNDICE DE PERSISTENCIA:
#       rho ~ 0 -> sin persistencia (como binomial) ; rho -> 1 -> caladero fijo
#   mu = probabilidad media de que una celda de la huella se pesque un año dado
# Se ajusta por arte (global) y por C-square (mapa de persistencia).
# Salidas: outputs/hintzen_betabin_*.csv|fst  y  figs/hintzen_betabin_*.png
# =============================================================================
if (length(artes_bb) > 0) {
  message("11. Beta-binomial de persistencia: ", paste(artes_bb, collapse = ", "))

  # --- Funciones ----------------------------------------------------------------
  bb_ab <- function(mu, rho) { s <- (1 - rho) / rho; c(a = mu * s, b = (1 - mu) * s) }
  dbb_tr <- function(k, n, mu, rho) {               # beta-binomial truncada en 0
    ab <- bb_ab(mu, rho); a <- ab[["a"]]; b <- ab[["b"]]
    exp(lchoose(n, k) + lbeta(k + a, n - k + b) - lbeta(a, b)) /
      (1 - exp(lbeta(a, n + b) - lbeta(a, b)))
  }
  dbin_tr <- function(k, n, p) dbinom(k, n, p) / (1 - (1 - p)^n)
  nll_bb <- function(th, k, n) {
    mu <- plogis(th[1]); rho <- plogis(th[2]); ab <- bb_ab(mu, rho)
    a <- ab[["a"]]; b <- ab[["b"]]
    ll  <- lchoose(n, k) + lbeta(k + a, n - k + b) - lbeta(a, b)
    lp0 <- lbeta(a, n + b) - lbeta(a, b)
    -sum(ll - log(-expm1(lp0)))
  }
  nll_bin <- function(t, k, n) {
    p <- plogis(t); -sum(dbinom(k, n, p, log = TRUE) - log(1 - (1 - p)^n))
  }
  ajustar_bb <- function(k, n) {
    m0 <- min(max(mean(k / n), 0.05), 0.95)
    o1 <- optim(c(qlogis(m0), qlogis(0.3)), nll_bb, k = k, n = n,
                method = "Nelder-Mead", control = list(maxit = 5000, reltol = 1e-10))
    o0 <- optimize(nll_bin, c(-12, 12), k = k, n = n)
    data.table(mu = plogis(o1$par[1]), rho = plogis(o1$par[2]),
               p_binomial = plogis(o0$minimum),
               AIC_betabin = 2 * o1$value + 4, AIC_binomial = 2 * o0$objective + 2,
               convergencia = o1$convergence)
  }

  # --- Datos: años pescados por celda ----------------------------------------
  n_anios_arte <- cnt_bb[, .(n = uniqueN(ANYO)), by = LE_MET4]
  K <- cnt_bb[y >= min_conteo_anio,
              .(k = uniqueN(ANYO), esfuerzo = sum(y)),
              by = .(LE_MET4, ix, iy, lon, lat, ux, uy, csq, lon_u, lat_u)]
  K <- merge(K, n_anios_arte, by = "LE_MET4")
  fst::write_fst(K, paste0(pref_out, "betabin_celdas.fst"))

  # --- Ajuste por arte ---------------------------------------------------------
  BB <- K[, {
    f <- ajustar_bb(k, n)
    cbind(data.table(n_anios = n[1], n_celdas = .N,
                     pct_celdas_todos_anios = 100 * mean(k == n[1]),
                     pct_celdas_un_anio = 100 * mean(k == 1),
                     pct_esf_celdas_todos_anios = 100 * sum(esfuerzo[k == n[1]]) / sum(esfuerzo)),
          f)
  }, by = LE_MET4]
  BB[, dAIC_binomial_menos_betabin := AIC_binomial - AIC_betabin]
  escribir(BB, "betabin_resumen_arte")
  print(BB)

  # --- Ajuste por C-square -----------------------------------------------------
  BBc <- K[, if (.N >= min_celdas_bb) {
    f <- ajustar_bb(k, n)
    cbind(data.table(n_celdas = .N, k_medio = mean(k),
                     pct_celdas_todos_anios = 100 * mean(k == n[1]),
                     esfuerzo = sum(esfuerzo)), f)
  }, by = .(LE_MET4, csq, ux, uy, lon_u, lat_u)]
  escribir(BBc, "betabin_csq")

  # --- Mapas: costa y límites ---------------------------------------------------
  if (!exists("costa") || is.null(costa)) {
    costa <- NULL
    if (requireNamespace("rnaturalearth", quietly = TRUE)) {
      costa <- tryCatch({
        esc <- if (requireNamespace("rnaturalearthhires", quietly = TRUE)) "large" else "medium"
        sf::sf_use_s2(FALSE)
        cc <- rnaturalearth::ne_countries(scale = esc, returnclass = "sf")
        cc <- suppressWarnings(sf::st_crop(sf::st_make_valid(cc),
                                           xmin = -9, xmax = -4, ymin = 34, ymax = 39))
        sf::sf_use_s2(TRUE)
        cc
      }, error = function(e) NULL)
    }
  }
  mapa_bb <- function(a, marg = 0.06) {
    kk <- K[LE_MET4 == a]
    list(
      if (!is.null(costa)) geom_sf(data = costa, fill = "grey85", colour = "grey45",
                                   linewidth = 0.2, inherit.aes = FALSE),
      coord_sf(xlim = range(kk$lon) + c(-1, 1) * marg,
               ylim = range(kk$lat) + c(-1, 1) * marg, expand = FALSE),
      labs(x = NULL, y = NULL), tema_mapa)
  }

  for (a in artes_bb) {
    Ka <- K[LE_MET4 == a]; ba <- BB[LE_MET4 == a]; n_a <- ba$n_anios
    if (!nrow(Ka)) next

    # B1. Distribución de k: observado vs binomial (sin persistencia) vs beta-binomial
    obs <- Ka[, .N, by = k][, prop := N / sum(N)]
    esp <- rbind(
      data.table(k = 1:n_a, prop = dbin_tr(1:n_a, n_a, ba$p_binomial),
                 modelo = "Binomial (sin persistencia)"),
      data.table(k = 1:n_a, prop = dbb_tr(1:n_a, n_a, ba$mu, ba$rho),
                 modelo = "Beta-binomial"))
    pB1 <- ggplot() +
      geom_col(data = obs, aes(factor(k), prop), fill = "grey70") +
      geom_line(data = esp, aes(k, prop, colour = modelo), linewidth = 0.9) +
      geom_point(data = esp, aes(k, prop, colour = modelo)) +
      scale_colour_manual(values = c("Binomial (sin persistencia)" = "#d7301f",
                                     "Beta-binomial" = "#2171b5"), name = NULL) +
      tema +
      labs(title = paste("Años pescados por celda —", artes_label[[a]]),
           x = sprintf("Nº de años con pesca (de %d)", n_a), y = "Proporción de celdas")
    guardar(pB1, paste0("betabin_B1_anios_pescados_", a), 9, 6)

    # B2. Mapa: nº de años pescados por celda fina
    pB2 <- ggplot() +
      geom_tile(data = Ka, aes(lon, lat, fill = factor(k)),
                width = res_fina[[a]], height = res_fina[[a]]) +
      scale_fill_viridis_d(option = "plasma", direction = -1,
                           name = sprintf("Años con pesca (de %d)", n_a)) +
      mapa_bb(a) +
      labs(title = paste("Persistencia del caladero a escala fina —", artes_label[[a]]))
    guardar(pB2, paste0("betabin_B2_mapa_anios_", a), 10, 7)

    # B3. Mapa: persistencia (rho) por C-square
    ca <- BBc[LE_MET4 == a]
    if (nrow(ca)) {
      pB3 <- ggplot() +
        geom_tile(data = ca, aes(lon_u, lat_u, fill = rho),
                  width = res_unidad, height = res_unidad, colour = "white") +
        geom_text(data = ca, aes(lon_u, lat_u, label = sprintf("%.2f", rho)), size = 2) +
        scale_fill_distiller(palette = "PuBuGn", direction = 1, limits = c(0, 1),
                             name = "ρ persistencia") +
        mapa_bb(a) +
        labs(title = paste("Índice de persistencia (ρ beta-binomial) por C-square —", artes_label[[a]]))
      guardar(pB3, paste0("betabin_B3_mapa_rho_", a), 10, 7)
    }

    # B4. ¿Qué parte del esfuerzo viene de celdas persistentes?
    ek <- Ka[, .(esfuerzo = sum(esfuerzo), celdas = .N), by = k][order(k)]
    ek[, `:=`(pct_esfuerzo = 100 * esfuerzo / sum(esfuerzo), pct_celdas = 100 * celdas / sum(celdas))]
    ekl <- melt(ek, id.vars = "k", measure.vars = c("pct_celdas", "pct_esfuerzo"))
    ekl[, variable := factor(variable, labels = c("% de celdas", "% del esfuerzo"))]
    pB4 <- ggplot(ekl, aes(factor(k), value, fill = variable)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("grey60", "#2171b5"), name = NULL) +
      tema +
      labs(title = paste("Celdas y esfuerzo según años de pesca —", artes_label[[a]]),
           x = sprintf("Nº de años con pesca (de %d)", n_a), y = "%")
    guardar(pB4, paste0("betabin_B4_esfuerzo_por_anios_", a), 9, 6)
    escribir(cbind(LE_MET4 = a, ek), paste0("betabin_esfuerzo_por_anios_", a))
  }
}

message("Script hintzen2019_agregacion_persistencia.R completado.")

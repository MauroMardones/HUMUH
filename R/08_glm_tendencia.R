# =============================================================================
# 08_glm_tendencia.R
# Tendencia temporal del esfuerzo por celda — Nested GLM
# glm(horas_pesca ~ ANYO, family = gaussian) por celda × arte
# Año anidado dentro de celda × arte (conforme a Pennino et al. y WGSFD)
# Extrae slope, p-valor y pseudo-R² → mapa de tendencias
# =============================================================================
#
# OBJETO DE ENTRADA
# -----------------
# esfuerzo  [data.table]  → generado en 04_effort_grid.R
#
# Si retomas el análisis desde cero (nueva sesión R), carga así:
#   library(data.table)
#   library(fst)
#   esfuerzo <- fst::read_fst("outputs/esfuerzo_2018_2024.fst",
#                              as.data.table = TRUE)
#
# LÓGICA
# ------
# Para cada celda × arte que tuvo esfuerzo en al menos N_MIN años:
#   - Se completa el panel con 0 en años sin pesca (esfuerzo real = 0)
#   - Se ajusta GLM Gaussiano: glm(horas_pesca ~ ANYO)
#   - Se extrae: slope (tendencia), p-valor, pseudo-R² (1 - dev/null.dev)
# Slope > 0 → esfuerzo creciente; slope < 0 → esfuerzo decreciente
# =============================================================================

library(data.table)
library(ggplot2)
library(ggthemes)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(viridis)
library(marmap)

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

N_MIN  <- 3      # mínimo de años con datos para ajustar el modelo
# Esto es relevante para evitar ajustar modelos con muy pocos puntos (p.ej. 1 o 2 años con pesca)
ALPHA  <- 0.05   # umbral de significancia

# =============================================================================
# 1. COMPLETAR PANEL CON CEROS
# Para cada celda × arte que aparece en al menos 1 año, rellenar los años
# sin observación con horas_pesca = 0
# =============================================================================
message("Completando panel con ceros...")

esf <- esfuerzo[LE_MET4 %in% artes_obj & horas_pesca > 0,
                .(lon_cel, lat_cel, ANYO, LE_MET4, horas_pesca)]

anos    <- sort(unique(esf$ANYO))
n_anos  <- length(anos)

# Grid completo: todas las combinaciones celda × arte × año
# Solo para combinaciones celda × arte que existen (evita generar celdas vacías)
celdas_arte <- unique(esf[, .(lon_cel, lat_cel, LE_MET4)])
panel_full  <- celdas_arte[, CJ(ANYO = anos), by = .(lon_cel, lat_cel, LE_MET4)]
panel_full  <- merge(panel_full, esf,
                     by  = c("lon_cel", "lat_cel", "LE_MET4", "ANYO"),
                     all.x = TRUE)
panel_full[is.na(horas_pesca), horas_pesca := 0]

message(sprintf("Filas panel completo: %s", format(nrow(panel_full), big.mark = ".")))

# =============================================================================
# 2. FILTRAR CELDAS CON SUFICIENTES AÑOS CON DATOS (> 0)
# =============================================================================
n_anos_pos <- panel_full[horas_pesca > 0, .(n_pos = .N), by = .(lon_cel, lat_cel, LE_MET4)]
celdas_ok  <- n_anos_pos[n_pos >= N_MIN, .(lon_cel, lat_cel, LE_MET4)]
panel_ok   <- panel_full[celdas_ok, on = .(lon_cel, lat_cel, LE_MET4)]

message(sprintf("Celdas × arte con >= %d años activos: %s",
  N_MIN, format(nrow(celdas_ok), big.mark = ".")))

# =============================================================================
# 3. AJUSTAR lm POR CELDA × ARTE
# =============================================================================
message("Ajustando modelos de tendencia (nested GLM: año anidado en celda × arte)...")

# GLM Gaussiano con link identidad por celda × arte
# Equivale a lm() en estimación pero se enmarca como GLM nested (año dentro de celda)
# conforme a Pennino et al. y literatura WGSFD/ICES de tendencias espaciales de esfuerzo
# Pseudo-R² = 1 - (devianza residual / devianza nula)

tendencia <- panel_ok[, {
  fit    <- glm(horas_pesca ~ ANYO, family = gaussian(link = "identity"))
  cf     <- summary(fit)$coefficients
  dev    <- summary(fit)$deviance
  null_d <- summary(fit)$null.deviance
  r2     <- if (null_d > 0) 1 - dev / null_d else NA_real_
  list(
    slope   = cf["ANYO", "Estimate"],
    pval    = cf["ANYO", "Pr(>|t|)"],
    r2      = r2,
    n_anos  = .N
  )
}, by = .(lon_cel, lat_cel, LE_MET4)]

# Añadir código C-square (join desde esfuerzo; csq es 1:1 con lon_cel/lat_cel)
csq_lookup <- unique(esfuerzo[, .(lon_cel, lat_cel, csq)])
tendencia[csq_lookup, csq := i.csq, on = .(lon_cel, lat_cel)]

# Clasificar tendencia
tendencia[, sig      := pval < ALPHA]
tendencia[, direccion := fcase(
  slope > 0 & sig,  "Aumento significativo",
  slope < 0 & sig,  "Descenso significativo",
  slope > 0 & !sig, "Aumento no significativo",
  slope < 0 & !sig, "Descenso no significativo"
)]
tendencia[, direccion := factor(direccion, levels = c(
  "Aumento significativo",
  "Aumento no significativo",
  "Descenso no significativo",
  "Descenso significativo"
))]
tendencia[, arte_label := artes_label[LE_MET4]]

message(sprintf("Celdas con modelo ajustado: %s", format(nrow(tendencia), big.mark = ".")))
tendencia[, .N, by = .(LE_MET4, direccion)][order(LE_MET4, direccion)]

fst::write_fst(tendencia, "outputs/tendencia_glm_2018_2024.fst", compress = 75)
message("Guardado: outputs/tendencia_glm_2018_2024.fst")

# =============================================================================
# 4. MAPAS DE TENDENCIA POR ARTE
# =============================================================================

if (!exists("bathy_df")) {
  message("Descargando batimetría NOAA...")
  bathy <- marmap::getNOAA.bathy(lon1 = -8, lon2 = -5,
                                  lat1 = 35, lat2 = 38,
                                  resolution = 1)
  bathy_df <- marmap::fortify.bathy(bathy)
  bathy_df <- bathy_df[bathy_df$z <= 0, ]
}
if (!exists("costa_fill")) {
  sf::sf_use_s2(FALSE)
  costa_fill <- ne_countries(scale = "large", returnclass = "sf")
  costa_fill <- st_crop(costa_fill, xmin = -9, xmax = -4, ymin = 33, ymax = 39)
  costa_line <- ne_coastline(scale = "large", returnclass = "sf")
  costa_line <- st_crop(costa_line, xmin = -9, xmax = -4, ymin = 33, ymax = 39)
  sf::sf_use_s2(TRUE)
}

xlim <- c(-9, -5)
ylim <- c(35, 37.5)

colores_dir <- c(
  "Aumento significativo"     = "#d62728",
  "Aumento no significativo"  = "#f7b6b2",
  "Sin tendencia detectable"  = "#d9d9d9",
  "Descenso no significativo" = "#aec7e8",
  "Descenso significativo"    = "#1f77b4"
)

# Todas las celdas con esfuerzo (para capa de fondo)
celdas_esf <- unique(esf[, .(lon_cel, lat_cel, LE_MET4)])

tema_mapa <- theme_few() +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    legend.text      = element_text(size = 8),
    legend.title     = element_text(size = 9, face = "bold"),
    panel.background = element_rect(fill = "#eaf4fb", colour = NA)
  )

# --- Mapa por arte (dirección de tendencia + celdas sin tendencia detectable) -
for (art in artes_obj) {
  dat   <- tendencia[LE_MET4 == art]
  fondo <- celdas_esf[LE_MET4 == art]
  if (nrow(fondo) == 0) next

  p <- ggplot() +
    geom_contour(
      data = bathy_df, aes(x = x, y = y, z = z),
      breaks = c(-25, -50, -200),
      colour = "#005b96", linewidth = 0.2, alpha = 0.4
    ) +
    geom_tile(data = fondo,
              aes(x = lon_cel, y = lat_cel),
              fill = "#d9d9d9") +
    geom_tile(data = dat,
              aes(x = lon_cel, y = lat_cel, fill = direccion)) +
    scale_fill_manual(
      values = colores_dir,
      name   = "Tendencia",
      drop   = FALSE
    ) +
    geom_sf(data = costa_fill, fill = "white", colour = NA) +
    geom_sf(data = costa_line, colour = "grey40", linewidth = 0.1) +
    coord_sf(xlim = xlim, ylim = ylim) +
    xlab(expression(paste(Longitud^o, ~"O"))) +
    ylab(expression(paste(Latitud^o,  ~"N"))) +
    labs(title = paste0("Tendencia del esfuerzo — ", artes_label[art])) +
    tema_mapa

  fname <- paste0("figs/tendencia_", art, ".png")
  ggsave(fname, p, width = 8, height = 7, dpi = 200)
  message(sprintf("Guardado: %s", fname))
}

# --- Mapa del slope continuo (todas las celdas, alpha por significancia) -----
for (art in artes_obj) {
  dat <- tendencia[LE_MET4 == art]
  if (nrow(dat) == 0) next

  lim_slope <- max(abs(dat$slope)) * 1.05

  p <- ggplot() +
    geom_contour(
      data = bathy_df, aes(x = x, y = y, z = z),
      breaks = c(-25, -50, -200),
      colour = "#005b96", linewidth = 0.2, alpha = 0.4
    ) +
    geom_tile(data = dat[sig == FALSE],
              aes(x = lon_cel, y = lat_cel, fill = slope),
              alpha = 0.8) +
    geom_tile(data = dat[sig == TRUE],
              aes(x = lon_cel, y = lat_cel, fill = slope),
              alpha = 1.0) +
    scale_fill_gradient2(
      low      = "#1f77b4",
      mid      = "white",
      high     = "#d62728",
      midpoint = 0,
      limits   = c(-lim_slope, lim_slope),
      name     = "Slope\n(h/año)"
    ) +
    geom_sf(data = costa_fill, fill = "white", colour = NA) +
    geom_sf(data = costa_line, colour = "grey40", linewidth = 0.1) +
    coord_sf(xlim = xlim, ylim = ylim) +
    xlab(expression(paste(Longitud^o, ~"O"))) +
    ylab(expression(paste(Latitud^o,  ~"N"))) +
    labs(title = paste0("Slope del esfuerzo — ", artes_label[art])) +
    tema_mapa

  fname <- paste0("figs/slope_", art, ".png")
  ggsave(fname, p, width = 8, height = 7, dpi = 200)
  message(sprintf("Guardado: slope_%s.png", art))
}

# =============================================================================
# 5. RESUMEN: % CELDAS POR CATEGORÍA × ARTE
# =============================================================================

resumen_dir <- tendencia[, .(n = .N), by = .(LE_MET4, direccion)]
resumen_dir[, total := sum(n), by = LE_MET4]
resumen_dir[, pct   := round(n / total * 100, 1)]
resumen_dir[, arte_label := artes_label[LE_MET4]]

p_resumen <- ggplot(resumen_dir,
                    aes(x = arte_label, y = pct, fill = direccion)) +
  geom_col(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(pct >= 5, paste0(pct, "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 3, colour = "grey20") +
  scale_fill_manual(values = colores_dir, name = "Tendencia") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  coord_flip() +
  xlab(NULL) +
  ylab("% de celdas") +
  labs(title = "Distribución de tendencias por arte") +
  theme_few() +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 8))

ggsave("figs/tendencia_resumen_artes.png", p_resumen,
       width = 10, height = 6, dpi = 200)
message("Guardado: tendencia_resumen_artes.png")

message("Script 08 completado.")

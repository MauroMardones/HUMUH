# =============================================================================
# 06_core_fishing_ground.R
# Core fishing ground: celdas que concentran el 90% del esfuerzo acumulado
# por arte (2018-2024)
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
# =============================================================================

library(data.table)
library(ggplot2)
library(ggthemes)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(viridis)
library(marmap)

# --- Parámetros --------------------------------------------------------------
artes_obj <- c("DRB", "HMD", "FPO", "GN", "GNS", "GTN", "GTR", "LL", "LLS", "LHP")
umbral    <- c(0.75, 0.90)   # core 75% y 90%

# --- Acumulado 2018-2024 por celda × arte ------------------------------------
# csq es 1:1 con (lon_cel, lat_cel) → incluirlo en by no cambia la agregación
acumulado <- esfuerzo[LE_MET4 %in% artes_obj,
  .(horas_pesca = sum(horas_pesca, na.rm = TRUE)),
  by = .(lon_cel, lat_cel, csq, LE_MET4)]
acumulado <- acumulado[horas_pesca > 0]

# --- Función core fishing ground ---------------------------------------------
# Para cada arte:
#   1. Ordena celdas de mayor a menor esfuerzo
#   2. Calcula proporción acumulada
#   3. Marca las celdas hasta alcanzar cada umbral

core_fg <- function(dt, u) {
  dt <- dt[order(-horas_pesca)]
  dt[, cum_prop := cumsum(horas_pesca) / sum(horas_pesca)]
  # "core" = celdas que JUNTAS acumulan hasta el umbral u
  # La celda que supera el umbral también se incluye (primera fuera del core)
  dt[, paste0("core_", u * 100) := as.integer(data.table::shift(cum_prop, 1, fill = 0) < u)]
  dt
}

# Aplica para cada arte y umbral
cfg_list <- lapply(artes_obj, function(art) {
  dt <- copy(acumulado[LE_MET4 == art])
  for (u in umbral) dt <- core_fg(dt, u)
  dt
})
cfg <- rbindlist(cfg_list)

message(sprintf("Celdas totales en core_90 (todos los artes): %s",
  format(cfg[core_90 == 1, .N], big.mark = ".")))
cfg[, .(
  n_celdas_total = .N,
  n_core75  = sum(core_75),
  n_core90  = sum(core_90),
  pct_celdas_core90 = round(sum(core_90) / .N * 100, 1)
), by = LE_MET4]

# --- Guardar resultado -------------------------------------------------------
fst::write_fst(cfg, "outputs/core_fishing_ground_2018_2024.fst", compress = 75)
message("Guardado: outputs/core_fishing_ground_2018_2024.fst")

# =============================================================================
# MAPAS CORE FISHING GROUND
# =============================================================================

# Reutiliza batimetría y costa si están en memoria (desde 05_maps.R),
# de lo contrario las descarga/carga aquí.
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

breaks_bathy <- c(-10, -25, -50, -100, -200, -500)
xlim <- c(-9, -5)
ylim <- c(35, 37.5)

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

tema_mapa <- theme_few() +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    legend.text      = element_text(size = 9),
    legend.title     = element_text(size = 10, face = "bold"),
    panel.background = element_rect(fill = "#eaf4fb", colour = NA)
  )

# Paleta: fondo = todo el esfuerzo (gris), core75 = naranja, core90 = rojo
for (art in artes_obj) {
  dat <- cfg[LE_MET4 == art]
  if (nrow(dat) == 0) next

  # Asignar categoría para el color
  dat[, categoria := fcase(
    core_75 == 1, "Core 75%",
    core_90 == 1, "Core 90%",
    default      = "Resto"
  )]
  dat[, categoria := factor(categoria,
    levels = c("Core 75%", "Core 90%", "Resto"))]

  p <- ggplot() +
    geom_contour(
      data = bathy_df, aes(x = x, y = y, z = z),
      breaks = breaks_bathy, colour = "#005b96",
      linewidth = 0.2, alpha = 0.4
    ) +
    geom_tile(data = dat,
              aes(x = lon_cel, y = lat_cel, fill = categoria)) +
    scale_fill_manual(
      values = c(
        "Core 75%" = "#d62728",
        "Core 90%" = "#ff7f0e",
        "Resto"    = "#aec7e8"
      ),
      name = "Zona de pesca",
      drop = FALSE
    ) +
    geom_sf(data = costa_fill, fill = "white", colour = NA) +
    geom_sf(data = costa_line, colour = "grey40", linewidth = 0.1) +
    coord_sf(xlim = xlim, ylim = ylim) +
    xlab(expression(paste(Longitud^o, ~"O"))) +
    ylab(expression(paste(Latitud^o, ~"N"))) +
    labs(title = paste0("Core fishing ground — ", artes_label[art])) +
    tema_mapa

  fname <- paste0("figs/core_fg_", art, ".png")
  ggsave(fname, p, width = 8, height = 7, dpi = 200)
  message(sprintf("Guardado: %s", fname))
}

# =============================================================================
# ANÁLISIS ANUAL DEL CORE FISHING GROUND
# Para cada año × arte: cuántas celdas concentran el 90% del esfuerzo
# =============================================================================

anual_arte <- esfuerzo[LE_MET4 %in% artes_obj & horas_pesca > 0,
  .(horas_pesca = sum(horas_pesca, na.rm = TRUE)),
  by = .(lon_cel, lat_cel, csq, ANYO, LE_MET4)]

# Aplica el mismo algoritmo de core por año × arte
cfg_anual_list <- lapply(artes_obj, function(art) {
  lapply(sort(unique(anual_arte$ANYO)), function(yr) {
    dt <- copy(anual_arte[LE_MET4 == art & ANYO == yr])
    if (nrow(dt) == 0) return(NULL)
    for (u in umbral) dt <- core_fg(dt, u)
    dt
  })
})
cfg_anual <- rbindlist(unlist(cfg_anual_list, recursive = FALSE))

# Resumen: n celdas y % celdas en core 90% por año × arte
resumen_anual <- cfg_anual[, .(
  n_celdas_usadas = .N,
  n_core75        = sum(core_75),
  n_core90        = sum(core_90),
  pct_core75      = round(sum(core_75) / .N * 100, 1),
  pct_core90      = round(sum(core_90) / .N * 100, 1)
), by = .(ANYO, LE_MET4)]

resumen_anual[, arte_label := artes_label[LE_MET4]]

fst::write_fst(cfg_anual,      "outputs/core_fg_anual_2018_2024.fst",   compress = 75)
fst::write_fst(resumen_anual,  "outputs/resumen_core_fg_2018_2024.fst", compress = 75)
message("Guardados: core_fg_anual y resumen_core_fg")


# =============================================================================
# PLOTS DE BARRAS: CORE FISHING GROUND POR AÑO Y ARTE
# =============================================================================

# --- 1. N° de celdas en core 90% por año — facetado por arte ----------------
p_nceldas <- ggplot(resumen_anual, aes(x = factor(ANYO), y = n_core90, fill = arte_label)) +
  geom_col(show.legend = FALSE, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = n_core90), vjust = -0.4, size = 2.8, colour = "grey30") +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  facet_wrap(~ arte_label, scales = "free_y", ncol = 3) +
  xlab("Año") +
  ylab("N° de celdas (core 90%)") +
  labs(
      subtitle = "N° de celdas de 0.05° que concentran el 90% del esfuerzo anual"
  ) +
  theme_few() +
  theme(
    strip.text  = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

ggsave("figs/core_fg_nceldas_anual.png", p_nceldas, width = 13, height = 9, dpi = 200)
message("Guardado: core_fg_nceldas_anual.png")


# --- 2. % de celdas usadas que son core 90% — concentración relativa --------
p_pct <- ggplot(resumen_anual, aes(x = factor(ANYO), y = pct_core90, fill = arte_label)) +
  geom_col(show.legend = FALSE, colour = "white", linewidth = 0.3) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_text(aes(label = paste0(pct_core90, "%")),
            vjust = -0.4, size = 2.8, colour = "grey30") +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  facet_wrap(~ arte_label, ncol = 3) +
  xlab("Año") +
  ylab("% de celdas activas en core 90%") +
  labs(
    title    = "Concentración relativa del esfuerzo por año y arte",
    subtitle = "% de celdas activas que constituyen el core 90% | línea = 50%"
  ) +
  theme_few() +
  theme(
    strip.text  = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

ggsave("figs/core_fg_pct_anual.png", p_pct, width = 13, height = 9, dpi = 200)
message("Guardado: core_fg_pct_anual.png")


# --- 3. Panel doble: n celdas + % celdas, por arte (línea temporal) ----------
p_linea <- ggplot(resumen_anual, aes(x = ANYO)) +
  geom_col(aes(y = n_core90, fill = arte_label),
           alpha = 0.6, show.legend = FALSE) +
  geom_line(aes(y = pct_core90 * max(resumen_anual$n_core90) / 100,
                group = 1),
            colour = "#d62728", linewidth = 0.8) +
  geom_point(aes(y = pct_core90 * max(resumen_anual$n_core90) / 100),
             colour = "#d62728", size = 1.8) +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  scale_y_continuous(
    name     = "N° celdas core 90%",
    sec.axis = sec_axis(
      ~ . * 100 / max(resumen_anual$n_core90),
      name   = "% celdas activas (core 90%)",
      labels = function(x) paste0(round(x), "%")
    )
  ) +
  scale_x_continuous(breaks = sort(unique(resumen_anual$ANYO))) +
  facet_wrap(~ arte_label, scales = "free_y", ncol = 3) +
  xlab("Año") +
  labs(
    title    = "Tendencia del core fishing ground por arte",
    subtitle = "Barras = N° celdas | Línea roja = % de celdas activas en core 90%"
  ) +
  theme_few() +
  theme(
    strip.text        = element_text(size = 9, face = "bold"),
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 8),
    axis.title.y.right = element_text(colour = "#d62728"),
    axis.text.y.right  = element_text(colour = "#d62728")
  )

ggsave("figs/core_fg_tendencia_anual.png", p_linea, width = 13, height = 9, dpi = 200)
message("Guardado: core_fg_tendencia_anual.png")

message("Script 06 completado.")

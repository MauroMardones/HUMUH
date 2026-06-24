# =============================================================================
# 05_maps.R
# Mapas de esfuerzo pesquero (horas de pesca) por arte y año
# Celdas 0.05° x 0.05° (~5x5 km) con geom_tile
# Estética: batimetría NOAA + theme_few
# Golfo de Cádiz - Demarcación Sudatlántica
# =============================================================================

library(data.table)
library(ggplot2)
library(ggthemes)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(marmap)
library(metR)
library(sf)
library(scales)

# esfuerzo debe estar en memoria desde 04_effort_grid.R
# o bien: esfuerzo <- fst::read_fst("outputs/esfuerzo_2018_2024.fst", as.data.table = TRUE)

# --- Batimetría NOAA ---------------------------------------------------------
message("Descargando batimetría NOAA...")
bathy <- marmap::getNOAA.bathy(
  lon1 = -8, lon2 = -5,
  lat1 = 35, lat2 = 38,
  resolution = 1
)
bathy_df <- marmap::fortify.bathy(bathy)
bathy_df <- bathy_df[bathy_df$z <= 0, ]  # solo zona marina

breaks_bathy  <- c(-10, -25, -50, -100, -200, -500)
breaks_etiq   <- c(-25, -50, -200)

# Breaks de horas de pesca para la leyenda (escala log1p)
breaks_hrs <- c(1, 5, 10, 50, 100, 500, 1000, 5000)
labels_hrs <- c("1", "5", "10", "50", "100", "500", "1 000", "5 000")

# --- Costa -------------------------------------------------------------------
library(rnaturalearthhires)
sf::sf_use_s2(FALSE)
costa_fill <- ne_countries(scale = "large", returnclass = "sf")
costa_fill <- st_crop(costa_fill, xmin = -9, xmax = -4, ymin = 33, ymax = 39)
costa_line <- ne_coastline(scale = "large", returnclass = "sf")
costa_line <- st_crop(costa_line, xmin = -9, xmax = -4, ymin = 33, ymax = 39)
sf::sf_use_s2(TRUE)

# --- Extensión del mapa ------------------------------------------------------
xlim <- c(-9, -5)
ylim <- c(35.5, 37.5)

# --- Tema base ---------------------------------------------------------------
tema_mapa <- theme_few() +
  theme(
    legend.position  = "bottom",
    legend.key       = element_rect(fill = "white"),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    legend.text      = element_text(size = 9),
    legend.title     = element_text(size = 10, face = "bold"),
    strip.text       = element_text(size = 9, face = "bold"),
    panel.background = element_rect(fill = "#eaf4fb", colour = NA)
  )

# --- Artes y etiquetas -------------------------------------------------------
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

# --- Esfuerzo medio por celda (media entre años) -----------------------------
medio <- esfuerzo[LE_MET4 != "TOTAL",
  .(horas_pesca = mean(horas_pesca, na.rm = TRUE)),
  by = .(lon_cel, lat_cel, LE_MET4)]
medio <- medio[horas_pesca > 0]

# --- Datos anuales -----------------------------------------------------------
anual <- esfuerzo[LE_MET4 != "TOTAL" & horas_pesca > 0]


# =============================================================================
# FUNCIÓN BASE DEL MAPA (evita repetir capas)
# =============================================================================
base_map <- function(dat_tile, titulo) {
  lim_max <- max(dat_tile$horas_pesca, na.rm = TRUE)
  bk <- breaks_hrs[breaks_hrs <= lim_max * 1.5]
  lb <- labels_hrs[breaks_hrs <= lim_max * 1.5]

  ggplot() +
    geom_contour(
      data = bathy_df,
      aes(x = x, y = y, z = z),
      breaks    = breaks_bathy,
      colour    = "#005b96",
      linewidth = 0.2,
      alpha     = 0.5
    ) +
    metR::geom_text_contour(
      data  = bathy_df,
      aes(x = x, y = y, z = z),
      breaks       = breaks_etiq,
      colour       = "#005b96",
      size         = 2.2,
      skip         = 0,
      label.placer = metR::label_placer_fraction(frac = 0.5),
      rotate       = FALSE,
      formatter    = function(x) paste0(abs(x), " m")
    ) +
    geom_tile(data = dat_tile,
              aes(x = lon_cel, y = lat_cel, fill = horas_pesca)) +
    geom_sf(data = costa_fill, fill = "white", colour = NA) +
    geom_sf(data = costa_line, colour = "grey40", linewidth = 0.1) +
    scale_fill_viridis(
      option = "D",
      trans  = "log1p",
      breaks = bk,
      labels = lb,
      name   = "Horas pesca",
      guide  = guide_colorbar(barwidth = 12, barheight = 0.6,
                              title.position = "top", title.hjust = 0.5)
    ) +
    coord_sf(xlim = xlim, ylim = ylim) +
    xlab(expression(paste(Longitud^o, ~"O"))) +
    ylab(expression(paste(Latitud^o,  ~"N"))) +
    labs(title = titulo) +
    tema_mapa
}

# =============================================================================
# MAPAS DE ESFUERZO MEDIO 2018-2024 (media entre años por celda)
# =============================================================================
for (art in names(artes_label)) {
  dat <- medio[LE_MET4 == art]
  if (nrow(dat) == 0) next
  p <- base_map(dat, paste0("Esfuerzo medio 2018-2024 — ", artes_label[art]))
  fname <- paste0("figs/esfuerzo_medio_", art, ".png")
  ggsave(fname, p, width = 10, height = 6, dpi = 200)
  message(sprintf("Guardado: %s", fname))
}

# =============================================================================
# MAPAS DE ESFUERZO MEDIO ANUAL (facet_wrap por año)
# =============================================================================
base_map_anual <- function(dat_tile, titulo) {
  lim_max <- max(dat_tile$horas_pesca, na.rm = TRUE)
  bk <- breaks_hrs[breaks_hrs <= lim_max * 1.5]
  lb <- labels_hrs[breaks_hrs <= lim_max * 1.5]

  ggplot() +
    geom_contour(
      data = bathy_df,
      aes(x = x, y = y, z = z),
      breaks    = breaks_bathy,
      colour    = "#005b96",
      linewidth = 0.15,
      alpha     = 0.4
    ) +
    geom_tile(data = dat_tile,
              aes(x = lon_cel, y = lat_cel, fill = horas_pesca)) +
    geom_sf(data = costa_fill, fill = "white", colour = NA) +
    geom_sf(data = costa_line, colour = "grey40", linewidth = 0.1) +
    scale_fill_viridis(
      option = "inferno",
      trans  = "log1p",
      breaks = bk,
      labels = lb,
      name   = "Horas pesca",
      guide  = guide_colorbar(barwidth = 10, barheight = 0.5,
                              title.position = "top", title.hjust = 0.5)
    ) +
    coord_sf(xlim = xlim, ylim = ylim) +
    facet_wrap(~ ANYO, ncol = 4) +
    xlab(expression(paste(Longitud^o, ~"O"))) +
    ylab(expression(paste(Latitud^o,  ~"N"))) +
    labs(title = titulo) +
    tema_mapa
}

for (art in names(artes_label)) {
  dat <- anual[LE_MET4 == art]
  if (nrow(dat) == 0) next
  p <- base_map_anual(dat, paste0("Esfuerzo medio anual — ", artes_label[art]))
  fname <- paste0("figs/esfuerzo_anual_", art, ".png")
  ggsave(fname, p, width = 16, height = 7, dpi = 200)
  message(sprintf("Guardado: %s", fname))
}

# =============================================================================
# HISTOGRAMAS DE DISTRIBUCIÓN DE ESFUERZO POR ARTE
# Métrica: horas totales por celda (suma entre años)
# =============================================================================

total_hist <- esfuerzo[LE_MET4 != "TOTAL",
  .(horas_pesca = sum(horas_pesca, na.rm = TRUE)),
  by = .(lon_cel, lat_cel, LE_MET4)]
total_hist <- total_hist[horas_pesca > 0]

# --- Histograma conjunto (facetado) ------------------------------------------
total_artes <- total_hist[LE_MET4 %in% names(artes_label)]
total_artes[, arte_label := artes_label[LE_MET4]]

p_hist_all <- ggplot(total_artes, aes(x = horas_pesca)) +
  geom_histogram(
    aes(fill = LE_MET4),
    bins   = 40,
    colour = "white",
    linewidth = 0.2,
    show.legend = FALSE
  ) +
  scale_fill_viridis(option = "turbo", discrete = TRUE) +
  facet_wrap(~ arte_label, scales = "free_x", ncol = 3) +
  xlab("Horas de pesca totales por celda") +
  ylab("N° de celdas") +
  labs(title = "Distribución del esfuerzo total 2018-2024 por arte") +
  theme_few() +
  theme(
    strip.text   = element_text(size = 9, face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 8)
  )

ggsave("figs/hist_esfuerzo_all.png", p_hist_all, width = 13, height = 9, dpi = 200)
message("Guardado: hist_esfuerzo_all.png")

# --- Histogramas individuales -------------------------------------------------
for (art in names(artes_label)) {
  dat_h <- total_hist[LE_MET4 == art]
  if (nrow(dat_h) == 0) next

  p_h <- ggplot(dat_h, aes(x = horas_pesca)) +
    geom_histogram(
      fill   = "#E63946",
      colour = "white",
      bins   = 40,
      linewidth = 0.2
    ) +
    geom_vline(
      xintercept = median(dat_h$horas_pesca),
      linetype = "dashed", colour = "grey30", linewidth = 0.7
    ) +
    annotate("text",
      x     = median(dat_h$horas_pesca) * 1.15,
      y     = Inf,
      label = paste0("Mediana: ", round(median(dat_h$horas_pesca), 1), " h"),
      vjust = 1.5, hjust = 0, size = 3.2, colour = "grey30"
    ) +
    xlab("Horas de pesca totales por celda") +
    ylab("N° de celdas") +
    labs(title = paste0("Distribución del esfuerzo total — ", artes_label[art])) +
    theme_few() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  fname <- paste0("figs/hist_esfuerzo_", art, ".png")
  ggsave(fname, p_h, width = 7, height = 5, dpi = 200)
  message(sprintf("Guardado: hist_esfuerzo_%s.png", art))
}

message("Script 05 completado. Mapas e histogramas guardados en figs/")

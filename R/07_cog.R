# =============================================================================
# 07_cog.R
# Centro de gravedad (COG) del esfuerzo pesquero por arte y año
# COG_lat  = Σ(horas × lat)  / Σ(horas)
# COG_lon  = Σ(horas × lon)  / Σ(horas)
# COG_prof = Σ(horas × prof) / Σ(horas)   ← desplazamiento batimétrico
# =============================================================================
#
# OBJETO DE ENTRADA
# -----------------
# esfuerzo  [data.table]  → generado en 04_effort_grid.R
#
# Si retomas el análisis desde cero (nueva sesión R), carga así:
#   library(data.table)
#   library(fst)
  # esfuerzo <- fst::read_fst("outputs/esfuerzo_2018_2024.fst",
  #                            as.data.table = TRUE)
#
# =============================================================================


library(data.table)
library(ggplot2)
library(ggthemes)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(marmap)

artes_obj <- c("DRB", "HMD", "FPO", "GN", "GNS", "GTN", "GTR", "LL", "LLS", "LHP", "LHM")

artes_label <- c(
  DRB = "DRB (Rastra)",
  HMD = "HMD (Draga hidráulica)",
  FPO = "FPO (Nasas)",
  GN  = "GN (Enmalle s/e)",
  GNS = "GNS (Enmalle fondo)",
  GTN = "GTN (Enmalle-trasmallo)",
  GTR = "GTR (Trasmallo)",
  LL  = "LL (Palangre s/e)",
  LLS = "LLS (Palangre calado)",
  LHP = "LHP (Línea de mano)",
  LHM = "LHM (Anzuelos/Línea de mano)"
)

# =============================================================================
# 1. PROFUNDIDAD EN CADA CELDA
# Descarga batimetría NOAA (o reutiliza si está en memoria) y extrae
# la profundidad en el centroide de cada celda del grid de esfuerzo
# =============================================================================
message("Obteniendo profundidades por celda...")

if (!exists("bathy")) {
  bathy <- marmap::getNOAA.bathy(lon1 = -8, lon2 = -5,
                                  lat1 = 35, lat2 = 38,
                                  resolution = 1)
}

# Centroides únicos del grid de esfuerzo
celdas <- unique(esfuerzo[LE_MET4 %in% artes_obj,
                           .(lon_cel, lat_cel)])

# --- Profundidad: join por redondeo con bathy_df ----------------------------
# bathy_df viene de fortify.bathy(bathy): columnas x, y, z
# Resolución NOAA = 1 arc-minute = 1/60° ≈ 0.01667°
# Redondeamos centroides de celda a la misma resolución y hacemos merge

if (!exists("bathy_df")) {
  bathy_df <- marmap::fortify.bathy(bathy)
}

bathy_dt <- as.data.table(bathy_df)

# Agregar bathy a la misma grilla del esfuerzo (0.05°)
# Evita problemas de snap entre distintas resoluciones
res_grid <- 0.05
bathy_dt[, lon_cel := floor(x / res_grid) * res_grid + res_grid / 2]
bathy_dt[, lat_cel := floor(y / res_grid) * res_grid + res_grid / 2]
bathy_key <- bathy_dt[z < 0, .(prof_m = mean(abs(z))), by = .(lon_cel, lat_cel)]

# Join directo por lon_cel / lat_cel
celdas <- merge(celdas, bathy_key, by = c("lon_cel", "lat_cel"), all.x = TRUE)

message(sprintf("Celdas con profundidad asignada: %s  |  NA: %s",
  celdas[!is.na(prof_m), .N],
  celdas[is.na(prof_m),  .N]))

message(sprintf("Celdas con profundidad asignada: %s  |  NA: %s",
  celdas[!is.na(prof_m), .N],
  celdas[is.na(prof_m),  .N]))

# Unir profundidad al esfuerzo
esfuerzo_cog <- merge(
  esfuerzo[LE_MET4 %in% artes_obj & horas_pesca > 0],
  celdas,
  by = c("lon_cel", "lat_cel"),
  all.x = TRUE
)

# =============================================================================
# 2. CÁLCULO DEL COG POR AÑO × ARTE
# =============================================================================
message("Calculando centro de gravedad...")

# Filtro marino: excluye tierra, puertos y estuarios
# prof_m >= 20 m excluye zonas portuarias, estuarios (Odiel/Tinto, Guadalquivir)
# y zonas intermareales donde la batimetría NOAA (1 arc-min) no resuelve bien
esfuerzo_mar <- esfuerzo_cog[!is.na(prof_m) & prof_m >= 10]

cog_pos <- esfuerzo_mar[, .(
  COG_lat     = sum(horas_pesca * lat_cel, na.rm = TRUE) / sum(horas_pesca, na.rm = TRUE),
  COG_lon     = sum(horas_pesca * lon_cel, na.rm = TRUE) / sum(horas_pesca, na.rm = TRUE),
  horas_total = sum(horas_pesca, na.rm = TRUE),
  n_celdas    = .N
), by = .(ANYO, LE_MET4)]

cog_prof <- esfuerzo_mar[, .(
  COG_prof = sum(horas_pesca * prof_m, na.rm = TRUE) / sum(horas_pesca, na.rm = TRUE)
), by = .(ANYO, LE_MET4)]

cog <- merge(cog_pos, cog_prof, by = c("ANYO", "LE_MET4"), all.x = TRUE)

cog[, arte_label := artes_label[LE_MET4]]

fst::write_fst(cog, "outputs/cog_2018_2024.fst", compress = 75)
message("Guardado: outputs/cog_2018_2024.fst")

print(cog[order(LE_MET4, ANYO)])

# =============================================================================
# 3. PLOTS
# =============================================================================

tema_cog <- theme_few() +
  theme(
    strip.text   = element_text(size = 9,  face = "bold"),
    axis.text.x  = element_text(size = 8),
    plot.title   = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    legend.position = "none"
  )

# --- 3a. COG latitud por año -------------------------------------------------
p_lat <- ggplot(cog, aes(x = ANYO, y = COG_lat, colour = arte_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_colour_viridis_d(option = "turbo") +
  scale_x_continuous(breaks = sort(unique(cog$ANYO))) +
  facet_wrap(~ arte_label, ncol = 3, scales = "free_y") +
  xlab("Año") +
  ylab("Latitud COG (°N)") +
  labs(
    title    = "Centro de gravedad — Latitud",
    ) +
  tema_cog

ggsave("figs/cog_latitud.png", p_lat, width = 13, height = 9, dpi = 200)
message("Guardado: cog_latitud.png")

# --- 3b. COG longitud por año ------------------------------------------------
p_lon <- ggplot(cog, aes(x = ANYO, y = COG_lon, colour = arte_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_colour_viridis_d(option = "turbo") +
  scale_x_continuous(breaks = sort(unique(cog$ANYO))) +
  facet_wrap(~ arte_label, ncol = 3, scales = "free_y") +
  xlab("Año") +
  ylab("Longitud COG (°O)") +
  labs(
    title    = "Centro de gravedad — Longitud",
    ) +
  tema_cog

ggsave("figs/cog_longitud.png", p_lon, width = 13, height = 9, dpi = 200)
message("Guardado: cog_longitud.png")

# --- 3c. COG profundidad por año ---------------------------------------------
p_prof <- ggplot(cog, aes(x = ANYO, y = COG_prof, colour = arte_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_colour_viridis_d(option = "turbo") +
  scale_x_continuous(breaks = sort(unique(cog$ANYO))) +
  scale_y_reverse() +    # profundidad mayor hacia abajo
  facet_wrap(~ arte_label, ncol = 3, scales = "free_y") +
  xlab("Año") +
  ylab("Profundidad COG (m)") +
  labs(
    title    = "Centro de gravedad — Profundidad",
    ) +
  tema_cog

ggsave("figs/cog_profundidad.png", p_prof, width = 13, height = 9, dpi = 200)
message("Guardado: cog_profundidad.png")

# --- 3d. Trayectoria espacial del COG (lat × lon) ----------------------------
# scale = "large" para capturar bien la Bahía de Algeciras y el estrecho
if (!exists("costa_fill")) {
  sf::sf_use_s2(FALSE)
  costa_fill <- ne_countries(scale = "large", returnclass = "sf")
  costa_fill <- st_crop(costa_fill, xmin = -9, xmax = -4, ymin = 33, ymax = 39)
  costa_line <- ne_coastline(scale = "large", returnclass = "sf")
  costa_line <- st_crop(costa_line, xmin = -9, xmax = -4, ymin = 33, ymax = 39)
  sf::sf_use_s2(TRUE)
}
if (!exists("bathy_df")) {
  bathy_df <- marmap::fortify.bathy(bathy)
  bathy_df  <- bathy_df[bathy_df$z <= 0, ]
}

p_traj <- ggplot() +
  geom_contour(
    data = bathy_df, aes(x = x, y = y, z = z),
    breaks = c(-25, -50, -200),
    colour = "#005b96", linewidth = 0.2, alpha = 0.4
  ) +
  geom_sf(data = costa_fill, fill = "white", colour = NA) +
  geom_sf(data = costa_line, colour = "grey40", linewidth = 0.1) +
  geom_point(data = cog,
             aes(x = COG_lon, y = COG_lat,
                 colour = arte_label, size = horas_total),
             alpha = 0.8) +
  geom_text(data = cog,
            aes(x = COG_lon, y = COG_lat, label = ANYO),
            size = 2.5, vjust = -0.8, colour = "grey20") +
  scale_colour_viridis_d(option = "turbo", name = "Arte") +
  scale_size_continuous(name = "Horas totales", range = c(1, 5)) +
  coord_sf(xlim = c(-9, -5), ylim = c(35.5, 37.5)) +
  facet_wrap(~ arte_label, ncol = 5) +
  xlab(expression(paste(Longitud^o, ~"O"))) +
  ylab(expression(paste(Latitud^o,  ~"N"))) +
  labs(title = "Centro de gravedad del esfuerzo 2018-2024") +
  theme_few() +
  theme(
    strip.text       = element_text(size = 9, face = "bold"),
    legend.position  = "bottom",
    panel.background = element_rect(fill = "#eaf4fb", colour = NA)
  )

ggsave("figs/cog_trayectoria.png", p_traj, width = 13, height = 12, dpi = 200)
message("Guardado: cog_trayectoria.png")

message("Script 07 completado.")


# =============================================================================
# 10_habitat_effort.R
# Intersección esfuerzo × hábitat MSFD (BBHT)
# - Une los 5 shapefiles de hábitat del Golfo de Cádiz
# - Spatial join centroide de celda → hábitat
# - Cuenta celdas con esfuerzo por arte × año × hábitat
# - Gráficos de barras: esfuerzo acumulado y anual por hábitat
# =============================================================================

library(data.table)
library(sf)
library(ggplot2)
library(ggthemes)
library(viridis)
library(fst)

# =============================================================================
# 1. CARGAR Y UNIR SHAPEFILES DE HÁBITAT
# =============================================================================
message("Cargando hábitats MSFD...")

shp_dir   <- "data/shp/capa_habitat_BHTY"
shp_files <- list.files(shp_dir, pattern = "\\.shp$", full.names = TRUE)

habitat_sf <- do.call(rbind, lapply(shp_files, function(f) {
  d <- sf::read_sf(f)
  d[, c("MSFD_BBHT", "geometry")]
}))
habitat_sf <- sf::st_make_valid(habitat_sf)

message(sprintf("Hábitats únicos: %s", length(unique(habitat_sf$MSFD_BBHT))))

# =============================================================================
# 2. CARGAR ESFUERZO Y CONSTRUIR CENTROIDES SF
# =============================================================================
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

esfuerzo <- fst::read_fst("outputs/esfuerzo_2018_2024.fst", as.data.table = TRUE)

esf <- esfuerzo[LE_MET4 %in% artes_obj & horas_pesca > 0]

# Centroides únicos → sf
celdas_sf <- unique(esf[, .(lon_cel, lat_cel)])
celdas_sf <- sf::st_as_sf(celdas_sf,
                           coords = c("lon_cel", "lat_cel"),
                           crs    = 4326)

# =============================================================================
# 3. SPATIAL JOIN: centroide → hábitat
# =============================================================================
message("Spatial join centroides × hábitat...")

sf::sf_use_s2(FALSE)
join <- sf::st_join(celdas_sf, habitat_sf["MSFD_BBHT"], join = sf::st_within)
sf::sf_use_s2(TRUE)

# Recuperar coordenadas y convertir a data.table
coords  <- as.data.frame(sf::st_coordinates(join))
join_dt <- data.table(
  lon_cel   = coords$X,
  lat_cel   = coords$Y,
  MSFD_BBHT = join$MSFD_BBHT
)

# Celdas sin hábitat asignado → "Sin asignar"
join_dt[is.na(MSFD_BBHT), MSFD_BBHT := "Sin asignar"]

message(sprintf("Celdas con hábitat asignado: %s  |  Sin asignar: %s",
  join_dt[MSFD_BBHT != "Sin asignar", .N],
  join_dt[MSFD_BBHT == "Sin asignar", .N]))

# =============================================================================
# 4. UNIR HÁBITAT AL ESFUERZO
# =============================================================================
esf <- merge(esf, join_dt, by = c("lon_cel", "lat_cel"), all.x = TRUE)
esf[is.na(MSFD_BBHT), MSFD_BBHT := "Sin asignar"]

# =============================================================================
# 5. AGREGACIÓN
# =============================================================================

# --- 5a. Acumulado 2018-2024: horas y n_celdas por arte × hábitat ------------
acum_hab <- esf[, .(
  horas_pesca = sum(horas_pesca, na.rm = TRUE),
  n_celdas    = uniqueN(paste(lon_cel, lat_cel))
), by = .(LE_MET4, MSFD_BBHT)]
acum_hab[, arte_label := artes_label[LE_MET4]]

# --- 5b. Anual: horas y n_celdas por arte × año × hábitat -------------------
anual_hab <- esf[, .(
  horas_pesca = sum(horas_pesca, na.rm = TRUE),
  n_celdas    = uniqueN(paste(lon_cel, lat_cel))
), by = .(LE_MET4, ANYO, MSFD_BBHT)]
anual_hab[, arte_label := artes_label[LE_MET4]]

# Guardar
fst::write_fst(acum_hab,  "outputs/habitat_esfuerzo_acum.fst",  compress = 75)
fst::write_fst(anual_hab, "outputs/habitat_esfuerzo_anual.fst", compress = 75)
message("Guardados: habitat_esfuerzo_acum.fst y habitat_esfuerzo_anual.fst")

# =============================================================================
# 6. GRÁFICOS
# =============================================================================

# Orden habitats por horas acumuladas totales (todos los artes)
hab_order <- acum_hab[, .(total = sum(horas_pesca)), by = MSFD_BBHT][order(-total), MSFD_BBHT]
acum_hab[,  MSFD_BBHT := factor(MSFD_BBHT, levels = rev(hab_order))]
anual_hab[, MSFD_BBHT := factor(MSFD_BBHT, levels = rev(hab_order))]

# Porcentajes acumulados por arte
acum_hab[, pct_horas   := horas_pesca / sum(horas_pesca) * 100, by = LE_MET4]
acum_hab[, pct_celdas  := n_celdas    / sum(n_celdas)    * 100, by = LE_MET4]

# Porcentajes anuales por arte × año
anual_hab[, pct_horas  := horas_pesca / sum(horas_pesca) * 100, by = .(LE_MET4, ANYO)]
anual_hab[, pct_celdas := n_celdas    / sum(n_celdas)    * 100, by = .(LE_MET4, ANYO)]

tema <- theme_few() +
  theme(
    strip.text  = element_text(size = 8, face = "bold"),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8),
    plot.title  = element_text(face = "bold", size = 11)
  )

# --- 6a. Barras acumuladas: horas por hábitat, facetado por arte -------------
p_acum <- ggplot(acum_hab, aes(x = MSFD_BBHT, y = horas_pesca, fill = arte_label)) +
  geom_col(show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  coord_flip() +
  facet_wrap(~ arte_label, scales = "free_x", ncol = 3) +
  xlab("Hábitat MSFD") +
  ylab("Horas de pesca acumuladas 2018-2024") +
  labs(title = "Esfuerzo acumulado por hábitat MSFD y arte") +
  tema

ggsave("figs/habitat_esfuerzo_acum.png", p_acum, width = 14, height = 10, dpi = 200)
message("Guardado: figs/habitat_esfuerzo_acum.png")

# --- 6b. Barras acumuladas: n_celdas por hábitat, facetado por arte ----------
p_acum_n <- ggplot(acum_hab, aes(x = MSFD_BBHT, y = n_celdas, fill = arte_label)) +
  geom_col(show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  coord_flip() +
  facet_wrap(~ arte_label, scales = "free_x", ncol = 3) +
  xlab("Hábitat MSFD") +
  ylab("N° de celdas") +
  labs(title = "N° de celdas por hábitat MSFD y arte (acumulado)") +
  tema

ggsave("figs/habitat_nceldas_acum.png", p_acum_n, width = 14, height = 10, dpi = 200)
message("Guardado: figs/habitat_nceldas_acum.png")

# --- 6c. Barras anuales: horas por hábitat × año, facetado por arte ----------
p_anual <- ggplot(anual_hab, aes(x = MSFD_BBHT, y = horas_pesca, fill = factor(ANYO))) +
  geom_col(position = "dodge") +
  scale_fill_viridis(discrete = TRUE, option = "mako", name = "Año") +
  coord_flip() +
  facet_wrap(~ arte_label, scales = "free_x", ncol = 3) +
  xlab("Hábitat MSFD") +
  ylab("Horas de pesca") +
  labs(title = "Esfuerzo anual por hábitat MSFD y arte") +
  tema +
  theme(legend.position = "bottom")

ggsave("figs/habitat_esfuerzo_anual.png", p_anual, width = 16, height = 12, dpi = 200)
message("Guardado: figs/habitat_esfuerzo_anual.png")

# --- 6d. Panel resumen: todos los artes apilados por hábitat (acumulado) -----
acum_total <- acum_hab[, .(horas_pesca = sum(horas_pesca)), by = .(MSFD_BBHT, arte_label)]

p_total <- ggplot(acum_total, aes(x = MSFD_BBHT, y = horas_pesca, fill = arte_label)) +
  geom_col() +
  scale_fill_viridis(discrete = TRUE, option = "turbo", name = "Arte") +
  coord_flip() +
  xlab("Hábitat MSFD") +
  ylab("Horas de pesca acumuladas 2018-2024") +
  labs(title = "Esfuerzo total por hábitat MSFD (todos los artes)") +
  tema +
  theme(legend.position = "right")

ggsave("figs/habitat_esfuerzo_total.png", p_total, width = 12, height = 7, dpi = 200)
message("Guardado: figs/habitat_esfuerzo_total.png")

acum_total2 <- acum_hab[, .(MSFD_BBHT, arte_label, horas_pesca, n_celdas, pct_horas, pct_celdas)]
acum_total2[, pct_horas_hab  := horas_pesca / sum(horas_pesca) * 100, by = MSFD_BBHT]
acum_total2[, pct_celdas_hab := n_celdas    / sum(n_celdas)    * 100, by = MSFD_BBHT]

p_total_pct_horas <- ggplot(acum_total2, aes(x = MSFD_BBHT, y = pct_horas_hab, fill = arte_label)) +
  geom_col() +
  scale_fill_viridis(discrete = TRUE, option = "turbo", name = "Arte") +
  coord_flip() +
  xlab("Hábitat MSFD") +
  ylab("% de horas de pesca") +
  labs(title    = "% del esfuerzo por hábitat MSFD (todos los artes)",
       subtitle = "100% = total horas en cada hábitat") +
  tema +
  theme(legend.position = "right")

ggsave("figs/habitat_pct_horas_total.png", p_total_pct_horas, width = 12, height = 7, dpi = 200)
message("Guardado: figs/habitat_pct_horas_total.png")

p_total_pct_celdas <- ggplot(acum_total2, aes(x = MSFD_BBHT, y = pct_celdas_hab, fill = arte_label)) +
  geom_col() +
  scale_fill_viridis(discrete = TRUE, option = "turbo", name = "Arte") +
  coord_flip() +
  xlab("Hábitat MSFD") +
  ylab("% de celdas") +
  labs(title    = "% de celdas por hábitat MSFD (todos los artes)",
       subtitle = "100% = total celdas en cada hábitat") +
  tema +
  theme(legend.position = "right")

ggsave("figs/habitat_pct_celdas_total.png", p_total_pct_celdas, width = 12, height = 7, dpi = 200)
message("Guardado: figs/habitat_pct_celdas_total.png")

# % sobre el gran total (todas las horas / todas las celdas = 100%)
gran_total_horas  <- sum(acum_total2$horas_pesca)
gran_total_celdas <- sum(acum_total2$n_celdas)
acum_total2[, pct_horas_gt  := horas_pesca / gran_total_horas  * 100]
acum_total2[, pct_celdas_gt := n_celdas    / gran_total_celdas * 100]

p_gt_horas <- ggplot(acum_total2, aes(x = MSFD_BBHT, y = pct_horas_gt, fill = arte_label)) +
  geom_col() +
  scale_fill_viridis(discrete = TRUE, option = "turbo", name = "Arte") +
  coord_flip() +
  xlab("Hábitat MSFD") +
  ylab("% de horas de pesca") +
  labs(title    = "% del esfuerzo total por hábitat MSFD (todos los artes)",
       subtitle = "100% = total horas de toda la flota 2018-2024") +
  tema +
  theme(legend.position = "right")

ggsave("figs/habitat_pct_horas_grandtotal.png", p_gt_horas, width = 12, height = 7, dpi = 200)
message("Guardado: figs/habitat_pct_horas_grandtotal.png")

p_gt_celdas <- ggplot(acum_total2, aes(x = MSFD_BBHT, y = pct_celdas_gt, fill = arte_label)) +
  geom_col() +
  scale_fill_viridis(discrete = TRUE, option = "turbo", name = "Arte") +
  coord_flip() +
  xlab("Hábitat MSFD") +
  ylab("% de celdas") +
  labs(title    = "% de celdas totales por hábitat MSFD (todos los artes)",
       subtitle = "100% = total celdas de toda la flota 2018-2024") +
  tema +
  theme(legend.position = "right")

ggsave("figs/habitat_pct_celdas_grandtotal.png", p_gt_celdas, width = 12, height = 7, dpi = 200)
message("Guardado: figs/habitat_pct_celdas_grandtotal.png")

# --- 6d2. Panel inverso: flotas en eje X, hábitat en leyenda -----------------
p_total_inv <- ggplot(acum_total, aes(x = arte_label, y = horas_pesca, fill = MSFD_BBHT)) +
  geom_col() +
  scale_fill_viridis(discrete = TRUE, option = "turbo", name = "Hábitat MSFD") +
  coord_flip() +
  xlab("Arte") +
  ylab("Horas de pesca acumuladas 2018-2024") +
  labs(title = "Esfuerzo total por arte (desglose por hábitat MSFD)") +
  tema +
  theme(legend.position  = "right",
        legend.text      = element_text(size = 7),
        legend.title     = element_text(size = 8, face = "bold"))

ggsave("figs/habitat_esfuerzo_total_inv.png", p_total_inv, width = 14, height = 7, dpi = 200)
message("Guardado: figs/habitat_esfuerzo_total_inv.png")

# =============================================================================
# 7. TABLA RESUMEN: % horas y % celdas por arte × hábitat (acumulado)
# =============================================================================
tabla_pct <- acum_hab[, .(
  arte_label,
  MSFD_BBHT,
  horas_pesca,
  n_celdas,
  pct_horas  = round(pct_horas,  2),
  pct_celdas = round(pct_celdas, 2)
)]
setorder(tabla_pct, arte_label, -pct_horas)

fwrite(tabla_pct, "outputs/CSV/habitat_pct_resumen.csv")
message("Guardado: outputs/CSV/habitat_pct_resumen.csv")

# --- 6e. % horas acumuladas por hábitat, facetado por arte ------------------
p_pct_horas <- ggplot(acum_hab, aes(x = MSFD_BBHT, y = pct_horas, fill = arte_label)) +
  geom_col(show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  coord_flip() +
  facet_wrap(~ arte_label, scales = "free_x", ncol = 3) +
  xlab("Hábitat MSFD") +
  ylab("% de horas de pesca") +
  labs(title    = "% del esfuerzo acumulado por hábitat MSFD y arte",
       subtitle = "100% = total horas por arte 2018-2024") +
  tema

ggsave("figs/habitat_pct_horas_acum.png", p_pct_horas, width = 14, height = 10, dpi = 200)
message("Guardado: figs/habitat_pct_horas_acum.png")

# --- 6f. % celdas acumuladas por hábitat, facetado por arte -----------------
p_pct_celdas <- ggplot(acum_hab, aes(x = MSFD_BBHT, y = pct_celdas, fill = arte_label)) +
  geom_col(show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  coord_flip() +
  facet_wrap(~ arte_label, scales = "free_x", ncol = 3) +
  xlab("Hábitat MSFD") +
  ylab("% de celdas") +
  labs(title    = "% de celdas por hábitat MSFD y arte",
       subtitle = "100% = total celdas por arte 2018-2024") +
  tema

ggsave("figs/habitat_pct_celdas_acum.png", p_pct_celdas, width = 14, height = 10, dpi = 200)
message("Guardado: figs/habitat_pct_celdas_acum.png")

# --- 6g. % horas anuales por hábitat × año, facetado por arte ---------------
p_pct_horas_anual <- ggplot(anual_hab, aes(x = MSFD_BBHT, y = pct_horas, fill = factor(ANYO))) +
  geom_col(position = "dodge") +
  scale_fill_viridis(discrete = TRUE, option = "mako", name = "Año") +
  coord_flip() +
  facet_wrap(~ arte_label, scales = "free_x", ncol = 3) +
  xlab("Hábitat MSFD") +
  ylab("% de horas de pesca") +
  labs(title    = "% del esfuerzo anual por hábitat MSFD y arte",
       subtitle = "100% = total horas por arte y año") +
  tema +
  theme(legend.position = "bottom")

ggsave("figs/habitat_pct_horas_anual.png", p_pct_horas_anual, width = 16, height = 12, dpi = 200)
message("Guardado: figs/habitat_pct_horas_anual.png")

# --- 6h. % celdas anuales por hábitat × año, facetado por arte --------------
p_pct_celdas_anual <- ggplot(anual_hab, aes(x = MSFD_BBHT, y = pct_celdas, fill = factor(ANYO))) +
  geom_col(position = "dodge") +
  scale_fill_viridis(discrete = TRUE, option = "mako", name = "Año") +
  coord_flip() +
  facet_wrap(~ arte_label, scales = "free_x", ncol = 3) +
  xlab("Hábitat MSFD") +
  ylab("% de celdas") +
  labs(title    = "% de celdas anuales por hábitat MSFD y arte",
       subtitle = "100% = total celdas por arte y año") +
  tema +
  theme(legend.position = "bottom")

ggsave("figs/habitat_pct_celdas_anual.png", p_pct_celdas_anual, width = 16, height = 12, dpi = 200)
message("Guardado: figs/habitat_pct_celdas_anual.png")

# --- 6i. Mapa: hábitats MSFD + grilla de esfuerzo (todos los artes) ----------
library(rnaturalearth)
library(rnaturalearthdata)

if (!exists("costa_fill")) {
  sf::sf_use_s2(FALSE)
  costa_fill <- ne_countries(scale = "large", returnclass = "sf")
  costa_fill <- st_crop(costa_fill, xmin = -9, xmax = -4, ymin = 33, ymax = 39)
  costa_line <- ne_coastline(scale = "large", returnclass = "sf")
  costa_line <- st_crop(costa_line, xmin = -9, xmax = -4, ymin = 33, ymax = 39)
  sf::sf_use_s2(TRUE)
}

# Centroides de celdas con esfuerzo (todos los artes, acumulado)
celdas_esf <- esf[, .(horas_pesca = sum(horas_pesca, na.rm = TRUE)),
                  by = .(lon_cel, lat_cel, MSFD_BBHT)]

p_mapa <- ggplot() +
  geom_sf(data = habitat_sf, aes(fill = MSFD_BBHT),
          colour = NA, alpha = 0.6) +
  scale_fill_viridis(discrete = TRUE, option = "turbo",
                     name = "Hábitat MSFD",
                     guide = guide_legend(ncol = 1, keyheight = 0.8)) +
  geom_tile(data = celdas_esf,
            aes(x = lon_cel, y = lat_cel),
            fill = NA, colour = "black", linewidth = 0.15,
            width = 0.05, height = 0.05) +
  geom_sf(data = costa_fill, fill = "white", colour = NA) +
  geom_sf(data = costa_line, colour = "grey40", linewidth = 0.2) +
  coord_sf(xlim = c(-9, -5), ylim = c(35.5, 37.5)) +
  xlab(expression(paste(Longitud^o, ~"O"))) +
  ylab(expression(paste(Latitud^o,  ~"N"))) +
  labs(title   = "Hábitats MSFD y grilla de esfuerzo — Golfo de Cádiz 2018-2024") +
  theme_few() +
  theme(
    legend.position  = "right",
    legend.text      = element_text(size = 7),
    legend.title     = element_text(size = 8, face = "bold"),
    panel.background = element_rect(fill = "#eaf4fb", colour = NA),
    plot.title       = element_text(face = "bold", size = 11)
  )

ggsave("figs/mapa_habitats_grilla.png", p_mapa, width = 14, height = 8, dpi = 200)
message("Guardado: figs/mapa_habitats_grilla.png")

message("Script 10 completado.")

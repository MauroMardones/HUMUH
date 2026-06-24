# =============================================================================
# 09_mapa_habitats.R
# Mapa de hábitats bentónicos MSFD — Golfo de Cádiz
# Fuente: Habitats_region_IV.shp (MSFD BBHT classification)
# =============================================================================
#
# OBJETO DE ENTRADA
# -----------------
# data/shp/habitats/Habitats_region_IV.shp   (EPSG:3857)
#
# SALIDA
# ------
# figs/mapa_habitats_MSFD.png
# =============================================================================

library(sf)
library(ggplot2)
library(ggthemes)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)

# =============================================================================
# 1. CARGAR Y RECORTAR SHAPEFILE DE HÁBITATS
# =============================================================================
message("Cargando shapefile de hábitats...")

sf::sf_use_s2(FALSE)

hab <- st_read("data/shp/habitats/Habitats_region_IV.shp", quiet = TRUE)
hab <- st_zm(hab, drop = TRUE, what = "ZM")   # elimina coordenadas Z/M
hab <- st_transform(hab, crs = 4326)
hab <- st_buffer(hab, dist = 0)               # repara geometrías inválidas vía buffer nulo

xlim <- c(-9, -5)
ylim <- c(35, 37.5)

hab_clip <- st_crop(hab,
                    xmin = xlim[1], xmax = xlim[2],
                    ymin = ylim[1], ymax = ylim[2])

sf::sf_use_s2(TRUE)

message(sprintf("Polígonos tras recorte: %s", format(nrow(hab_clip), big.mark = ".")))

# =============================================================================
# 2. AGRUPAR CATEGORÍAS MSFD_BBHT EN GRUPOS BENTÓNICOS AMPLIOS
# =============================================================================
hab_clip <- hab_clip %>%
  mutate(grupo = case_when(
    grepl("Abyssal",         MSFD_BBHT, ignore.case = TRUE) ~ "Abisal",
    grepl("Lower bathyal",   MSFD_BBHT, ignore.case = TRUE) ~ "Batial inferior",
    grepl("Upper bathyal",   MSFD_BBHT, ignore.case = TRUE) ~ "Batial superior",
    grepl("Offshore",        MSFD_BBHT, ignore.case = TRUE) ~ "Circalitoral offshore",
    grepl("Circalittoral",   MSFD_BBHT, ignore.case = TRUE) ~ "Circalitoral",
    grepl("Infralittoral",   MSFD_BBHT, ignore.case = TRUE) ~ "Infralitoral",
    MSFD_BBHT == "Na"                                        ~ "Sin clasificar",
    TRUE                                                     ~ "Otro"
  )) %>%
  mutate(grupo = factor(grupo, levels = c(
    "Infralitoral",
    "Circalitoral",
    "Circalitoral offshore",
    "Batial superior",
    "Batial inferior",
    "Abisal",
    "Sin clasificar"
  )))

message("Distribución por grupo:")
print(table(hab_clip$grupo))

# =============================================================================
# 3. COSTA
# =============================================================================
if (!exists("costa_fill")) {
  sf::sf_use_s2(FALSE)
  costa_fill <- ne_countries(scale = "large", returnclass = "sf")
  costa_fill <- st_crop(costa_fill, xmin = -10, xmax = -4, ymin = 34, ymax = 39)
  costa_line <- ne_coastline(scale = "large", returnclass = "sf")
  costa_line <- st_crop(costa_line, xmin = -10, xmax = -4, ymin = 34, ymax = 39)
  sf::sf_use_s2(TRUE)
}

# =============================================================================
# 4. MAPA
# =============================================================================
colores_hab <- c(
  "Infralitoral"          = "#f4d166",
  "Circalitoral"          = "#f4a131",
  "Circalitoral offshore" = "#c25e1a",
  "Batial superior"       = "#6baed6",
  "Batial inferior"       = "#2171b5",
  "Abisal"                = "#08306b",
  "Sin clasificar"        = "#cccccc"
)

p_hab <- ggplot() +
  geom_sf(data = hab_clip,
          aes(fill = grupo),
          colour = NA, alpha = 0.92) +
  scale_fill_manual(
    values = colores_hab,
    name   = "Hábitat MSFD BBHT",
    drop   = FALSE
  ) +
  geom_sf(data = costa_fill, fill = "white",   colour = NA) +
  geom_sf(data = costa_line, colour = "grey40", linewidth = 0.1) +
  coord_sf(xlim = xlim, ylim = ylim) +
  xlab(expression(paste(Longitud^o, ~"O"))) +
  ylab(expression(paste(Latitud^o,  ~"N"))) +
  labs(title = "Hábitats bentónicos MSFD — Golfo de Cádiz") +
  theme_few() +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8),
    plot.title       = element_text(face = "bold", size = 12),
    panel.background = element_rect(fill = "#eaf4fb", colour = NA)
  )

ggsave("figs/mapa_habitats_MSFD.png", p_hab, width = 10, height = 6, dpi = 200)
message("Guardado: figs/mapa_habitats_MSFD.png")

# =============================================================================
# mapa_puntos_pesca_ultimo_anio.R
# Mapa de puntos de pesca (pings en estado de pesca) por pesquería (LE_MET4)
# para el último año disponible.
# Entrada: objeto de salida de R/03_join_CV_NV.R -> outputs/CV_NV_2018_2024_ZE.fst
# Estado de pesca: mismas ventanas de velocidad que R/04_effort_grid.R
# Correr desde la raíz del proyecto (HUMUH.Rproj)
# =============================================================================

library(data.table)
library(fst)
library(ggplot2)
library(ggthemes)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(marmap)
library(patchwork)

# --- 1. Leer solo las columnas necesarias -------------------------------------
cols <- c("VE_REF", "ANYO", "SI_TIMESTAMP", "SI_LATI", "SI_LONG", "SI_SPCA", "LE_MET4")
CV_NV <- read_fst("outputs/CV_NV_2018_2024_ZE.fst",
                  columns = cols, as.data.table = TRUE)
CV_NV[, SI_SPCA := as.numeric(SI_SPCA)]

# --- 2. Último año con datos --------------------------------------------------
ult_anio <- max(CV_NV$ANYO, na.rm = TRUE)
message("Último año: ", ult_anio)
pts <- CV_NV[ANYO == ult_anio]
rm(CV_NV); gc()

# --- 3. Estado de pesca: ventanas de velocidad por arte (igual que 04) -------
vel <- data.table(
  LE_MET4 = c("DRB","HMD","FPO","GN","GNS","GTN","GTR","LL","LLS","LHP"),
  vmin    = c(0.1, 1.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
  vmax    = c(1.0, 3.0, 3.0, 4.5, 2.0, 2.0, 2.0, 2.5, 2.5, 1.5)
)
pts <- pts[vel, on = "LE_MET4", nomatch = 0L][SI_SPCA >= vmin & SI_SPCA <= vmax]
message(sprintf("Pings de pesca %d: %s", ult_anio, format(nrow(pts), big.mark = ".")))
print(pts[, .(pings = .N, barcos = uniqueN(VE_REF)), by = LE_MET4][order(-pings)])

# --- 3b. Filtrar un mes (ejemplo) --------------------------------------------
# NULL = mes con más pings de pesca en el último año; o fija uno, p.ej. mes <- 5
mes <- NULL
if (!inherits(pts$SI_TIMESTAMP, c("POSIXct", "Date")))
  pts[, SI_TIMESTAMP := as.POSIXct(SI_TIMESTAMP, tz = "UTC")]
pts[, MES := month(SI_TIMESTAMP)]
if (is.null(mes)) mes <- pts[, .N, by = MES][which.max(N), MES]
pts <- pts[MES == mes]
mes_lab <- format(as.Date(sprintf("%d-%02d-01", ult_anio, mes)), "%B %Y")
message(sprintf("Mes: %s · pings de pesca: %s", mes_lab, format(nrow(pts), big.mark = ".")))
print(pts[, .(pings = .N, barcos = uniqueN(VE_REF)), by = LE_MET4][order(-pings)])

# --- 4. Submuestreo opcional (puntos grandes: conviene pocos) -----------------
max_pts_arte <- 3000   # pon Inf para usar todos
set.seed(123)
pts <- pts[, .SD[sample(.N, min(.N, max_pts_arte))], by = LE_MET4]

# --- 5. Etiquetas y orden de paneles -----------------------------------------
artes_label <- c(
  DRB = "DRB (Dredge)",                 HMD = "HMD (Hydraulic dredge)",
  FPO = "FPO (Pots & traps)",           GN  = "GN (Gillnet, unspecified)",
  GNS = "GNS (Set gillnet)",            GTN = "GTN (Combined gillnet-trammel)",
  GTR = "GTR (Trammel net)",            LL  = "LL (Longline, unspecified)",
  LLS = "LLS (Set longline)",           LHP = "LHP (Hand line)"
)
pts[, arte := factor(artes_label[LE_MET4], levels = artes_label)]
pts <- pts[!is.na(arte)]

# --- 6. Costa y batimetría ---------------------------------------------------
sf::sf_use_s2(FALSE)
costa <- ne_countries(scale = "large", returnclass = "sf")
costa_full <- costa
costa <- st_crop(costa, xmin = -9, xmax = -4, ymin = 33, ymax = 39)
sf::sf_use_s2(TRUE)

bathy <- getNOAA.bathy(lon1 = -8, lon2 = -5, lat1 = 35, lat2 = 38, resolution = 1)
bathy_df <- fortify.bathy(bathy)
bathy_df <- bathy_df[bathy_df$z <= 0, ]

xlim <- c(-7.9, -5.2)   # algo más amplio que antes
ylim <- c(35.7, 37.45)

# --- 7. Mapa único con puntos superpuestos -----------------------------------
pts <- pts[sample(.N)]   # orden aleatorio: ningún arte tapa sistemáticamente a otro

p <- ggplot() +
  geom_contour(data = bathy_df, aes(x = x, y = y, z = z),
               breaks = c(-25, -50, -100, -200),
               colour = "#005b96", linewidth = 0.2, alpha = 0.5) +
  geom_point(data = pts, aes(x = SI_LONG, y = SI_LATI, fill = arte),
             shape = 21, colour = "grey40", stroke = 0.25,
             size = 1.8, alpha = 0.85) +
  geom_sf(data = costa, fill = "grey90", colour = "grey40", linewidth = 0.1) +
  coord_sf(xlim = xlim, ylim = ylim) +
  scale_fill_viridis_d(option = "turbo", name = "Fishery (LE_MET4)") +
  guides(fill = guide_legend(override.aes = list(size = 6, alpha = 1), ncol = 1)) +
  labs(x = NULL, y = NULL) +
  theme_few() +
  theme(legend.position  = "right",
        legend.title     = element_text(face = "bold", size = 15),
        legend.text      = element_text(size = 13),
        legend.key.size  = unit(0.9, "cm"),
        axis.text        = element_text(size = 11),
        panel.background = element_rect(fill = "#eaf4fb", colour = NA))

# --- 7b. Mapa de ubicación (inset) -------------------------------------------
sf::sf_use_s2(FALSE)
iberia <- st_crop(costa_full, xmin = -12, xmax = 5, ymin = 33, ymax = 45)
sf::sf_use_s2(TRUE)
caja <- st_as_sfc(st_bbox(c(xmin = xlim[1], xmax = xlim[2],
                            ymin = ylim[1], ymax = ylim[2]), crs = 4326))
p_loc <- ggplot() +
  geom_sf(data = iberia, fill = "grey80", colour = "grey50", linewidth = 0.1) +
  geom_sf(data = caja, fill = NA, colour = "#c0392b", linewidth = 0.9) +
  coord_sf(xlim = c(-10.5, 4), ylim = c(34.5, 44.5), expand = FALSE) +
  theme_void() +
  theme(panel.background = element_rect(fill = "#eaf4fb", colour = "grey40", linewidth = 0.5))

# esquina inferior izquierda (mar abierto)
p <- p + inset_element(p_loc, left = 0.01, bottom = 0.01, right = 0.27, top = 0.33,
                       align_to = "panel")

print(p)

ggsave(sprintf("figs/puntos_pesca_%d_%02d.png", ult_anio, mes), p,
       width = 11, height = 7, dpi = 300, bg = "white")

# Copia con nombre fijo para la slide "Why the Gulf of Cádiz?" (docs/WGSFD_Turku_2026.Rmd)
ggsave("figs/puntos_pesca_mes.png", p,
       width = 11, height = 7, dpi = 300, bg = "white")

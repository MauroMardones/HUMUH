# =============================================================================
# fig_csquare_nested.R
# Nested C-square figure for the Gulf of Cádiz (WGSFD Turku 2026 slide)
# 10° -> 5° -> 1° -> 0.5° -> 0.1° -> 0.05°, over our own cumulative effort
# Output: figs/csquare_nested_GoC.png
# Run from the project root (HUMUH.Rproj)
# =============================================================================

library(sf)
library(ggplot2)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthhires)

sf_use_s2(FALSE)

navy <- "#0b2e4f"; teal <- "#0f7c8c"; sand <- "#e9c46a"; red <- "#c0392b"

# --- 1. C-square code (Rees, 2003), same logic as csquare.cod() in 04 --------
csq_code <- function(lon, lat, degrees) {
  g    <- 4 - ((2 * floor(1 + lon / 200) - 1) * (2 * floor(1 + lat / 200) + 1))
  code <- as.character(g * 1000 + floor(abs(lat) / 10) * 100 + floor(abs(lon) / 10))
  latR <- abs(lat) - 10 * floor(abs(lat) / 10)
  lonR <- abs(lon) - 10 * floor(abs(lon) / 10)
  steps <- c(5, 1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001)
  for (s in steps) {
    if (degrees > s) break
    q <- 2 * floor(latR * 0.2) + floor(lonR * 0.2) + 1
    if (s %in% c(5, 0.5, 0.05, 0.005)) {          # quadrant level
      code <- paste0(code, ":", q)
    } else {                                      # decimal level
      ld <- floor(latR); od <- floor(lonR)
      code <- sub(":[1-4]$", "", code)            # replace the quadrant digit
      code <- paste0(code, ":", q * 100 + ld * 10 + od)
      latR <- round((latR - ld) * 10, 7); lonR <- round((lonR - od) * 10, 7)
    }
  }
  code
}

# bounding box of the square of size `deg` that contains (lon, lat)
sq_box <- function(lon, lat, deg) {
  x0 <- floor(lon / deg) * deg; y0 <- floor(lat / deg) * deg
  st_as_sfc(st_bbox(c(xmin = x0, ymin = y0, xmax = x0 + deg, ymax = y0 + deg),
                    crs = 4326))
}

# regular grid lines of size `deg` inside a box
grid_lines <- function(xlim, ylim, deg) {
  xs <- seq(floor(xlim[1] / deg) * deg, ceiling(xlim[2] / deg) * deg, by = deg)
  ys <- seq(floor(ylim[1] / deg) * deg, ceiling(ylim[2] / deg) * deg, by = deg)
  xlim <- unname(xlim); ylim <- unname(ylim)
  rbind(data.frame(x = xs, xend = xs, y = ylim[1], yend = ylim[2]),
        data.frame(x = xlim[1], xend = xlim[2], y = ys, yend = ys))
}

# --- 2. Data -----------------------------------------------------------------
land <- ne_countries(scale = "large", returnclass = "sf")

f_eff  <- "outputs/shp/esfuerzo_acumulado_GoC_2018_2024.gpkg"
layers <- st_layers(f_eff)$name
eff <- do.call(rbind, lapply(layers, function(l)
  st_read(f_eff, layer = l, quiet = TRUE)[, c("lon_cel", "lat_cel", "csq", "horas_pesca")]))
eff <- aggregate(horas_pesca ~ lon_cel + lat_cel + csq, data = st_drop_geometry(eff), FUN = sum)

# sanity check: our code reproduces the codes stored by 04_effort_grid.R
stopifnot(all(csq_code(eff$lon_cel, eff$lat_cel, 0.05) == eff$csq))

eff_sf <- st_as_sf(do.call(rbind, lapply(seq_len(nrow(eff)), function(i)
  data.frame(eff[i, ], geometry = sq_box(eff$lon_cel[i], eff$lat_cel[i], 0.05)))))

# --- 3. Pick the example: the busiest 1° square, then the busiest cells inside
eff$sq1  <- csq_code(eff$lon_cel, eff$lat_cel, 1)
best1    <- names(which.max(tapply(eff$horas_pesca, eff$sq1, sum)))
in1      <- eff[eff$sq1 == best1, ]
in1$sq05 <- csq_code(in1$lon_cel, in1$lat_cel, 0.5)
# the 0.5° square with most fished cells (the busiest one is mostly land)
n05      <- table(in1$sq05)
best05   <- names(n05)[which.max(n05 + tapply(in1$horas_pesca, in1$sq05, sum)[names(n05)] / 1e9)]
in05     <- in1[in1$sq05 == best05, ]
# a well-fished cell near the centre of the 0.5° square, so it reads clearly
c05      <- c(mean(floor(in05$lon_cel / .5) * .5 + .25), mean(floor(in05$lat_cel / .5) * .5 + .25))
cand     <- in05[in05$horas_pesca >= median(in05$horas_pesca), ]
cell     <- cand[which.min((cand$lon_cel - c05[1])^2 + (cand$lat_cel - c05[2])^2), ]
lon0 <- cell$lon_cel; lat0 <- cell$lat_cel

levels_tab <- data.frame(
  deg  = c(10, 5, 1, 0.5, 0.1, 0.05),
  lab  = c("10°", "5°", "1°", "0.5°", "0.1°", "0.05°"))
levels_tab$code <- sapply(levels_tab$deg, function(d) csq_code(lon0, lat0, d))
levels_tab$size <- sprintf("%s × %s km",
  formatC(levels_tab$deg * 111.32 * cos(lat0 * pi / 180), format = "fg", digits = 2, big.mark = ","),
  formatC(levels_tab$deg * 111.32, format = "fg", digits = 2, big.mark = ","))
print(levels_tab)

box <- lapply(levels_tab$deg, function(d) sq_box(lon0, lat0, d)); names(box) <- levels_tab$lab

# --- 4. Panels -----------------------------------------------------------------
th <- theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "#eef7f8", colour = NA),
        panel.border = element_rect(fill = NA, colour = "grey60"),
        axis.title = element_blank(), axis.text = element_text(size = 7, colour = "grey40"),
        plot.title = element_text(colour = navy, face = "bold", size = 12),
        plot.subtitle = element_text(colour = "grey30", size = 9.5))

fill_sc <- scale_fill_gradient(low = "#cfe7ea", high = navy, trans = "log10",
                               name = "Hours fished\n2018–2024", labels = scales::comma)

# A: 10° and 5°
xA <- c(-12, 1); yA <- c(29, 44)
pA <- ggplot() +
  geom_sf(data = land, fill = "grey80", colour = "grey55", linewidth = .2) +
  geom_sf(data = box[["10°"]], fill = NA, colour = navy, linewidth = .8, linetype = "dashed") +
  geom_sf(data = box[["5°"]],  fill = alpha(teal, .15), colour = teal, linewidth = .8) +
  geom_sf(data = box[["1°"]],  fill = NA, colour = red, linewidth = .9) +
  annotate("label", x = -9.8, y = 39.6, label = levels_tab$code[1], hjust = 0, vjust = 1,
           size = 3.3, colour = navy, fill = "white", label.size = 0) +
  annotate("label", x = -9.8, y = 35.3, label = levels_tab$code[2], hjust = 0, vjust = 0,
           size = 3.3, colour = teal, fill = "white", label.size = 0) +
  coord_sf(xlim = xA, ylim = yA, expand = FALSE) +
  labs(title = "10° → 5°", subtitle = "Global and intermediate squares") + th

# B: 1° and 0.5° over effort
b1 <- st_bbox(box[["1°"]])
xB <- c(-7.8, -5.0); yB <- c(35.8, 37.4)
pB <- ggplot() +
  geom_sf(data = eff_sf, aes(fill = horas_pesca), colour = NA) +
  geom_segment(data = grid_lines(xB, yB, 1), aes(x, y, xend = xend, yend = yend),
               colour = "grey45", linewidth = .3, linetype = "dotted") +
  geom_sf(data = land, fill = "grey80", colour = "grey55", linewidth = .2) +
  geom_sf(data = box[["1°"]],   fill = NA, colour = red,  linewidth = .9) +
  geom_sf(data = box[["0.5°"]], fill = NA, colour = sand, linewidth = 1.1) +
  annotate("label", x = b1["xmin"] + .03, y = b1["ymin"] + .03, label = levels_tab$code[3],
           hjust = 0, vjust = 0, size = 3.3, colour = red, fill = "white", label.size = 0) +
  fill_sc +
  coord_sf(xlim = xB, ylim = yB, expand = FALSE) +
  labs(title = "1° → 0.5°", subtitle = "Gulf of Cádiz, our cumulative effort") + th +
  theme(legend.position = "none")

# C: 0.1° and 0.05° inside the 0.5° square
b05 <- st_bbox(box[["0.5°"]])
xC <- unname(c(b05["xmin"], b05["xmax"])); yC <- unname(c(b05["ymin"], b05["ymax"]))
pC <- ggplot() +
  geom_sf(data = eff_sf, aes(fill = horas_pesca), colour = "white", linewidth = .15) +
  geom_segment(data = grid_lines(xC, yC, 0.1), aes(x, y, xend = xend, yend = yend),
               colour = "grey35", linewidth = .35) +
  geom_sf(data = land, fill = "grey80", colour = "grey55", linewidth = .2) +
  geom_sf(data = box[["0.1°"]],  fill = NA, colour = teal, linewidth = 1) +
  geom_sf(data = box[["0.05°"]], fill = NA, colour = red,  linewidth = 1.1) +
  annotate("label", x = xC[1] + .015, y = yC[2] - .015, label = levels_tab$code[6],
           hjust = 0, vjust = 1, size = 3.3, colour = red, fill = "white", label.size = 0) +
  fill_sc +
  coord_sf(xlim = xC, ylim = yC, expand = FALSE) +
  labs(title = "0.1° → 0.05°", subtitle = "Analysis cell (~5 km)") + th +
  theme(legend.position = "right", legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))

p <- pA + pB + pC + plot_layout(widths = c(1, 1.45, 1.15))

dir.create("figs", showWarnings = FALSE)
ggsave("figs/csquare_nested_GoC.png", p, width = 13, height = 4.6, dpi = 200, bg = "white")
write.csv(levels_tab[, c("lab", "code", "size")], "outputs/csquare_levels_example.csv",
          row.names = FALSE)
message("Saved figs/csquare_nested_GoC.png")

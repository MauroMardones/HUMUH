# =============================================================================
# fig_dt_capping.R
# Histograma del intervalo entre pings (Δt) en estado de pesca, con umbrales
# propuestos para acotar el esfuerzo ("Δt capping").
# Entrada : outputs/CV_NV_2018_2024_ZE.fst (salida de R/03_join_CV_NV.R)
# Salida  : figs/dt_capping_hist.png  (slide "Open questions")
#           outputs/dt_capping_resumen.csv
# Correr desde la raíz del proyecto (HUMUH.Rproj)
# =============================================================================

library(data.table)
library(fst)
library(ggplot2)

umbrales <- c(10, 30, 60)   # minutos: propuestas a discutir en el grupo

# --- 1. Datos: solo las columnas necesarias ----------------------------------
x <- read_fst("outputs/CV_NV_2018_2024_ZE.fst",
              columns = c("SI_TDIFF", "SI_SPCA", "LE_MET4"),
              as.data.table = TRUE)
x[, SI_SPCA := as.numeric(SI_SPCA)]

# --- 2. Pings en estado de pesca (mismas ventanas que R/04_effort_grid.R) ----
vel <- data.table(
  LE_MET4 = c("DRB","HMD","FPO","GN","GNS","GTN","GTR","LL","LLS","LHP"),
  vmin    = c(0.1, 1.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
  vmax    = c(1.0, 3.0, 3.0, 4.5, 2.0, 2.0, 2.0, 2.5, 2.5, 1.5))
x <- x[vel, on = "LE_MET4", nomatch = 0L][SI_SPCA >= vmin & SI_SPCA <= vmax]
dt <- x[!is.na(SI_TDIFF) & SI_TDIFF > 0, SI_TDIFF / 60]   # minutos
rm(x); gc()

# --- 3. Resumen por umbral ----------------------------------------------------
tot_h <- sum(dt) / 60
res <- rbindlist(lapply(umbrales, function(u) data.table(
  umbral_min  = u,
  pct_pings   = 100 * mean(dt > u),                       # intervalos por encima
  pct_horas   = 100 * sum(dt[dt > u]) / sum(dt),          # horas que vienen de ellos
  pct_recorte = 100 * sum(pmax(dt - u, 0)) / sum(dt)      # horas que se pierden al truncar
)))
med <- median(dt)
message(sprintf("Pings de pesca: %s · horas: %s · mediana Δt: %.1f min",
                format(length(dt), big.mark = ","), format(round(tot_h), big.mark = ","), med))
print(res)
fwrite(res, "outputs/dt_capping_resumen.csv")

# --- 4. Histograma en escala log (bins precalculados: rápido con millones) ---
brk <- 10^seq(log10(0.5), log10(24 * 60), length.out = 60)
h <- data.table(dt = dt)[, bin := cut(dt, brk, include.lowest = TRUE)][, .N, by = bin]
h[, `:=`(xmin = brk[as.integer(bin)], xmax = brk[as.integer(bin) + 1])]
h <- h[!is.na(bin)]

ymax <- max(h$N)
lab <- res[, sprintf("%d min\n%.1f%% of intervals\n%.0f%% of hours",
                     umbral_min, pct_pings, pct_horas)]
cols <- c("#b8860b", "#c0392b", "#0b2e4f")

p <- ggplot(h) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 1, ymax = N),
            fill = "#0f7c8c", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = med, colour = "grey40", linetype = "dotted", linewidth = 0.6) +
  annotate("text", x = med, y = ymax, label = sprintf("median\n%.1f min", med),
           hjust = 1.05, vjust = 1, size = 3.6, colour = "grey30") +
  geom_vline(xintercept = umbrales, colour = cols, linetype = "dashed", linewidth = 0.9) +
  annotate("label", x = umbrales * 1.08, y = ymax * c(0.6, 0.04, 0.6), label = lab,
           hjust = 0, size = 3.4, colour = cols, fill = "white", label.size = 0.3) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 30, 60, 180, 600, 1440),
                labels = c("1", "2", "5", "10", "30", "60", "180", "600", "1440")) +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Δt between consecutive fishing pings (minutes, log scale)",
       y = "Number of intervals (log scale)") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", colour = "#0b2e4f"),
        plot.subtitle = element_text(colour = "grey30", size = 11))

ggsave("figs/dt_capping_hist.png", p, width = 11, height = 4.8, dpi = 200, bg = "white")
message("Guardado: figs/dt_capping_hist.png")

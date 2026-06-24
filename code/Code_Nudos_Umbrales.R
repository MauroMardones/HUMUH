library(fst)
library(data.table)
library(ggplot2)

CV_NV <- fst::read_fst("outputs/CV_NV_2018_2024_ZE.fst", as.data.table = TRUE)

umbrales <- data.table(
  LE_MET4 = c("DRB", "HMD", "FPO", "GN", "GNS", "GTN", "GTR", "LL", "LLS", "LHP"),
  lo      = c(0.1,   1.0,   0.1,   0.1,  0.1,   0.1,   0.1,   0.1,  0.1,   0.1),
  hi      = c(1.0,   3.0,   3.0,   2.0,  2.0,   2.0,   2.0,   2.5,  2.5,   1.5)
)

dat <- CV_NV[LE_MET4 %in% umbrales$LE_MET4 & !is.na(SI_SPCA)]

# percentil 99 por arte para truncar eje X
p99 <- dat[, .(p99 = quantile(SI_SPCA, 0.99, na.rm = TRUE)), by = LE_MET4]
dat <- dat[p99, on = "LE_MET4"][SI_SPCA <= p99]

ggplot(dat, aes(x = SI_SPCA)) +
  geom_histogram(bins = 80, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(data = umbrales, aes(xintercept = lo), color = "red",     linetype = "dashed", linewidth = 0.7) +
  geom_vline(data = umbrales, aes(xintercept = hi), color = "darkred", linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ LE_MET4, scales = "free", ncol = 3) +
  labs(x = "SI_SPCA (nudos)", y = "N pings") +
  theme_bw()

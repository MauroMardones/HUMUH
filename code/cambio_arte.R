# =============================================================================
# cambio_arte.R
# ¿Cuántas embarcaciones cambian de arte y cuántas veces?
# Por año (2018-2024) y para el período total.
#
# UNIDAD DE ANÁLISIS
#   El arte se asigna por marea diaria (FT_REF = barco + día), usando el arte
#   dominante en kilos de la nota de venta (ver R/03_join_CV_NV.R).
#   Para cada barco se ordena su secuencia de mareas por fecha y se mira si el
#   arte (LE_MET4) cambia entre una marea y la siguiente.
#
# MÉTRICAS (por barco, dentro de cada año y en el total 2018-2024)
#   n_viajes      : n° de mareas diarias con arte asignado
#   n_artes       : n° de artes distintos usados
#   n_cambios     : n° de veces que el arte cambia entre mareas consecutivas
#   cambia        : TRUE si n_artes >= 2                       (criterio estricto)
#   n_artes_rob   : n° de artes usados en >= min_viajes mareas
#   n_cambios_rob : cambios contados solo con esos artes
#   cambia_rob    : TRUE si n_artes_rob >= 2                   (criterio robusto)
#   arte_ppal     : arte con más mareas; pct_ppal = % de mareas en ese arte
#
#   El criterio robusto evita contar como "cambio" una marea suelta con otro
#   arte, que suele ser un artefacto del arte dominante en la nota de venta
#   (p. ej. un día en que se vendió más pulpo que chirla).
#
# DENOMINADORES
#   n_trackeadas : barcos con cajas verdes ese año (CV formateadas, antes del cruce)
#   n_con_arte   : barcos con al menos una marea con LE_MET4 (CV × NV)
#   La proporción principal es sobre n_con_arte (sólo en ellos se puede ver
#   un cambio); también se reporta sobre n_trackeadas.
#
# ENTRADAS : outputs/CV_NV_2018_2024_ZE.fst  (03_join_CV_NV.R)
#            outputs/CV_2017_2024_ZE_fmt.fst (02_format_CV.R; solo si denom_cv = TRUE)
# SALIDAS  : outputs/CSV/cambio_arte_resumen.csv      (año + Total)
#            outputs/CSV/cambio_arte_barco.csv        (barco × período)
#            outputs/CSV/cambio_arte_transiciones.csv (de → a, Total)
#            outputs/CSV/cambio_arte_interanual.csv   (arte principal año a año)
#            figs/cambio_arte_*.png
# Correr desde la raíz del proyecto (HUMUH.Rproj)
# =============================================================================

library(data.table)
library(fst)
library(ggplot2)
library(ggthemes)

# --- Parámetros ---------------------------------------------------------------
f_cvnv     <- "outputs/CV_NV_2018_2024_ZE.fst"
f_cv       <- "outputs/CV_2017_2024_ZE_fmt.fst"
anios      <- 2018:2024
min_viajes <- 3        # mareas mínimas con un arte para contarlo (criterio robusto)
denom_cv   <- TRUE     # TRUE = lee CV formateadas para el n° de barcos trackeados
chunk      <- 2e7      # filas por bloque al leer los .fst grandes

dir.create("outputs/CSV", showWarnings = FALSE, recursive = TRUE)
dir.create("figs",        showWarnings = FALSE)

# Lectura por bloques de un .fst grande, quedándose con filas únicas
leer_unico <- function(f, cols) {
  n <- fst::metadata_fst(f)$nrOfRows
  ini <- seq(1, n, by = chunk)
  out <- lapply(ini, function(i) {
    unique(read_fst(f, columns = cols, from = i, to = min(i + chunk - 1, n),
                    as.data.table = TRUE))
  })
  unique(rbindlist(out))
}

# =============================================================================
# 1. MAREAS CON ARTE (una fila por FT_REF)
# =============================================================================
message("Leyendo mareas con arte (CV × NV)...")
viajes <- leer_unico(f_cvnv, c("VE_REF", "FT_REF", "LE_MET4"))
viajes <- viajes[!is.na(LE_MET4)]

# Fecha desde FT_REF (VE_REF_YYYY_MM_DD); el año es el de la marea
viajes[, fecha := as.Date(sub(".*_(\\d{4}_\\d{2}_\\d{2})$", "\\1", FT_REF), "%Y_%m_%d")]
viajes[, ANYO := as.integer(format(fecha, "%Y"))]
viajes <- viajes[ANYO %in% anios]

# Control: cada marea debería tener un solo arte (03 colapsa NV por FT_REF)
dup <- viajes[, .N, by = FT_REF][N > 1]
if (nrow(dup)) warning(nrow(dup), " FT_REF con más de un arte; se usa el primero")
viajes <- unique(viajes, by = "FT_REF")
setorder(viajes, VE_REF, fecha)

message(sprintf("  %s mareas, %s barcos",
                format(nrow(viajes), big.mark = ".", decimal.mark = ","), uniqueN(viajes$VE_REF)))

# =============================================================================
# 2. BARCOS TRACKEADOS (denominador amplio)
# =============================================================================
if (denom_cv) {
  message("Leyendo barcos trackeados (CV formateadas)...")
  cv_barco <- leer_unico(f_cv, c("VE_REF", "ANYO"))[ANYO %in% anios]
} else {
  cv_barco <- unique(viajes[, .(VE_REF, ANYO)])   # sin CV: = barcos con arte
}
n_track <- rbind(cv_barco[, .(n_trackeadas = uniqueN(VE_REF)), by = .(periodo = as.character(ANYO))],
                 data.table(periodo = "Total", n_trackeadas = uniqueN(cv_barco$VE_REF)))

# =============================================================================
# 3. MÉTRICAS POR BARCO
# =============================================================================
# n° de cambios en una secuencia ordenada de artes
n_cambios_fun <- function(x) if (length(x) < 2) 0L else sum(x[-1] != x[-length(x)])

metricas_barco <- function(dt) {
  dt[, {
    tab      <- table(LE_MET4)
    artes_ok <- names(tab)[tab >= min_viajes]
    x_rob    <- LE_MET4[LE_MET4 %in% artes_ok]
    .(n_viajes      = .N,
      n_artes       = length(tab),
      n_cambios     = n_cambios_fun(LE_MET4),
      n_artes_rob   = length(artes_ok),
      n_cambios_rob = n_cambios_fun(x_rob),
      arte_ppal     = names(tab)[which.max(tab)],
      pct_ppal      = round(100 * max(tab) / .N, 1),
      artes         = paste(sort(names(tab)), collapse = "+"))
  }, by = VE_REF][, `:=`(cambia = n_artes >= 2, cambia_rob = n_artes_rob >= 2)]
}

# Por año: la secuencia se corta dentro de cada año
barco_anual <- rbindlist(lapply(anios, function(a) {
  if (!nrow(viajes[ANYO == a])) return(NULL)
  metricas_barco(viajes[ANYO == a])[, periodo := as.character(a)]
}))
# Total: secuencia completa 2018-2024 (incluye cambios entre años)
barco_total <- metricas_barco(viajes)[, periodo := "Total"]

barco <- rbind(barco_anual, barco_total)
setcolorder(barco, c("periodo", "VE_REF"))
fwrite(barco, "outputs/CSV/cambio_arte_barco.csv")

# =============================================================================
# 4. RESUMEN POR AÑO Y TOTAL
# =============================================================================
resumen <- barco[, .(
  n_con_arte        = .N,
  n_cambian         = sum(cambia),
  n_cambian_rob     = sum(cambia_rob),
  total_cambios     = sum(n_cambios),
  total_cambios_rob = sum(n_cambios_rob),
  # entre los barcos que cambian: cuántas veces lo hacen
  cambios_media     = round(mean(n_cambios[cambia]), 1),
  cambios_mediana   = as.numeric(median(n_cambios[cambia])),
  cambios_max       = if (any(cambia)) max(n_cambios) else 0L,
  artes_media       = round(mean(n_artes), 2)
), by = periodo]

resumen <- merge(n_track, resumen, by = "periodo", all.y = TRUE)

resumen[, `:=`(
  pct_cambian       = round(100 * n_cambian     / n_con_arte, 1),
  pct_cambian_rob   = round(100 * n_cambian_rob / n_con_arte, 1),
  pct_cambian_track = round(100 * n_cambian     / n_trackeadas, 1)
)]
resumen[, periodo := factor(periodo, levels = c(as.character(anios), "Total"))]
setorder(resumen, periodo)
setcolorder(resumen, c("periodo", "n_trackeadas", "n_con_arte",
                       "n_cambian", "pct_cambian", "n_cambian_rob", "pct_cambian_rob",
                       "pct_cambian_track"))
fwrite(resumen, "outputs/CSV/cambio_arte_resumen.csv")
print(resumen)

# =============================================================================
# 5. TRANSICIONES DE → A (período total)
# =============================================================================
trans <- viajes[, .(de = LE_MET4, a = shift(LE_MET4, type = "lead")), by = VE_REF][
  !is.na(a) & de != a]
trans <- trans[, .(n_cambios = .N, n_barcos = uniqueN(VE_REF)), by = .(de, a)][
  order(-n_cambios)]
fwrite(trans, "outputs/CSV/cambio_arte_transiciones.csv")

# =============================================================================
# 6. CAMBIO DEL ARTE PRINCIPAL ENTRE AÑOS
#    ¿El barco cambia de "flota" de un año a otro (no solo estacionalmente)?
# =============================================================================
ppal <- barco_anual[, .(VE_REF, ANYO = as.integer(periodo), arte_ppal)]
setorder(ppal, VE_REF, ANYO)
ppal[, `:=`(arte_prev = shift(arte_ppal), anio_prev = shift(ANYO)), by = VE_REF]
interanual <- ppal[!is.na(arte_prev), .(
  n_barcos         = .N,                         # barcos presentes en año y año anterior
  n_cambia_ppal    = sum(arte_ppal != arte_prev),
  pct_cambia_ppal  = round(100 * mean(arte_ppal != arte_prev), 1)
), by = ANYO][order(ANYO)]
interanual <- rbind(interanual[, .(periodo = as.character(ANYO), n_barcos, n_cambia_ppal, pct_cambia_ppal)],
  ppal[!is.na(arte_prev), .(periodo = "Total", n_barcos = uniqueN(VE_REF),
       n_cambia_ppal = uniqueN(VE_REF[arte_ppal != arte_prev]),
       pct_cambia_ppal = round(100 * uniqueN(VE_REF[arte_ppal != arte_prev]) / uniqueN(VE_REF), 1))])
fwrite(interanual, "outputs/CSV/cambio_arte_interanual.csv")
print(interanual)

# =============================================================================
# 7. FIGURAS
# =============================================================================
tema <- theme_few() +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        legend.position = "bottom",
        strip.text    = element_text(face = "bold"))
col_crit <- c("Estricto (≥ 2 artes)" = "#2a6f97",
              "Robusto (≥ 2 artes con ≥ min mareas)" = "#a9d6e5")
names(col_crit)[2] <- sprintf("Robusto (≥ 2 artes con ≥ %d mareas)", min_viajes)
niv_per <- c(as.character(anios), "Total")

# --- 7.1 Proporción de barcos que cambian de arte, por año y total ------------
d1 <- melt(resumen[, .(periodo, n_con_arte, n_cambian, n_cambian_rob)],
           id.vars = c("periodo", "n_con_arte"), value.name = "n")
d1[, criterio := factor(ifelse(variable == "n_cambian", names(col_crit)[1], names(col_crit)[2]),
                        levels = names(col_crit))]
d1[, prop := n / n_con_arte]
d1[, tipo := ifelse(periodo == "Total", "Total 2018-2024", "Anual")]

p1 <- ggplot(d1, aes(periodo, prop, fill = criterio)) +
  geom_col(position = position_dodge(0.8), width = 0.75, colour = "grey30", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%d/%d", n, n_con_arte)),
            position = position_dodge(0.8), vjust = -0.4, size = 2.6) +
  geom_vline(xintercept = length(anios) + 0.5, linetype = "dashed", colour = "grey50") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = col_crit, name = NULL) +
  labs(x = NULL, y = "% de embarcaciones con arte asignado",
       title = "Embarcaciones que cambian de arte")+
tema
ggsave("figs/cambio_arte_prop.png", p1, width = 9, height = 5.5, dpi = 200)

# --- 7.2 Cuántas veces cambian (solo barcos que cambian) ----------------------
d2 <- barco[cambia == TRUE]
d2[, periodo := factor(periodo, levels = niv_per)]
p2 <- ggplot(d2, aes(periodo, n_cambios)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.35, size = 1, colour = "#2a6f97") +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.55) +
  geom_vline(xintercept = length(anios) + 0.5, linetype = "dashed", colour = "grey50") +
  scale_y_continuous(trans = "log1p", breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500)) +
  labs(x = NULL, y = "N° de cambios de arte por barco (escala log)",
       title = "Frecuencia de cambio de arte",
       subtitle = "Solo barcos con ≥ 2 artes; cambio = arte distinto entre mareas consecutivas") +
  tema
ggsave("figs/cambio_arte_ncambios.png", p2, width = 9, height = 5.5, dpi = 200)

# --- 7.3 N° de artes por barco (composición) ----------------------------------
d3 <- barco[, .(n = .N), by = .(periodo, n_artes = factor(pmin(n_artes, 4)))]
levels(d3$n_artes)[levels(d3$n_artes) == "4"] <- "4+"
d3[, periodo := factor(periodo, levels = niv_per)]
p3 <- ggplot(d3, aes(periodo, n, fill = n_artes)) +
  geom_col(position = "fill", colour = "grey30", linewidth = 0.2) +
  geom_vline(xintercept = length(anios) + 0.5, linetype = "dashed", colour = "grey50") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_brewer(palette = "Blues", name = "N° artes") +
  labs(x = NULL, y = "% de embarcaciones",
       title = "N° de artes distintos usados por embarcación") +
  tema
ggsave("figs/cambio_arte_nartes.png", p3, width = 8, height = 5, dpi = 200)

# --- 7.4 Matriz de transiciones (Total) ---------------------------------------
artes_niv <- sort(unique(viajes$LE_MET4))
d4 <- CJ(de = artes_niv, a = artes_niv)
d4 <- merge(d4, trans, by = c("de", "a"), all.x = TRUE)
p4 <- ggplot(d4, aes(a, de, fill = n_cambios)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(is.na(n_cambios), "", sprintf("%d\n(%d)", n_cambios, n_barcos))),
            size = 2.4) +
  scale_fill_distiller(palette = "Blues", direction = 1, trans = "log10",
                       na.value = "grey95", name = "N° cambios") +
  labs(x = "Arte siguiente (a)", y = "Arte anterior (de)",
       title = "Transiciones de arte 2018-2024",
       subtitle = "N° de cambios entre mareas consecutivas (n° de barcos entre paréntesis)") +
  coord_equal() + tema + theme(legend.position = "right")
ggsave("figs/cambio_arte_transiciones.png", p4, width = 8, height = 7, dpi = 200)

# --- 7.5 % de barcos que cambian, según su arte principal (Total) -------------
d5 <- barco[periodo == "Total", .(n = .N, n_cambian = sum(cambia),
                                  n_cambian_rob = sum(cambia_rob)), by = arte_ppal]
d5 <- melt(d5, id.vars = c("arte_ppal", "n"), value.name = "k")
d5[, criterio := factor(ifelse(variable == "n_cambian", names(col_crit)[1], names(col_crit)[2]),
                        levels = names(col_crit))]
p5 <- ggplot(d5, aes(reorder(arte_ppal, -k / n), k / n, fill = criterio)) +
  geom_col(position = position_dodge(0.8), width = 0.75, colour = "grey30", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%d/%d", k, n)), position = position_dodge(0.8),
            vjust = -0.4, size = 2.6) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = col_crit, name = NULL) +
  labs(x = "Arte principal del barco (2018-2024)", y = "% de embarcaciones que cambian",
       title = "Cambio de arte según arte principal") +
  tema
ggsave("figs/cambio_arte_por_arte.png", p5, width = 9, height = 5.5, dpi = 200)

# --- 7.6 Cambio del arte principal entre años ---------------------------------
d6 <- interanual[periodo != "Total"]
p6 <- ggplot(d6, aes(periodo, pct_cambia_ppal / 100)) +
  geom_col(fill = "#2a6f97", colour = "grey30", linewidth = 0.2, width = 0.7) +
  geom_text(aes(label = sprintf("%d/%d", n_cambia_ppal, n_barcos)), vjust = -0.4, size = 2.8) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "% de embarcaciones",
       title = "Cambio del arte principal respecto al año anterior",
       subtitle = sprintf("Barcos presentes en dos años consecutivos. Total período: %d de %d barcos (%.1f %%) cambiaron alguna vez",
                        interanual[periodo == "Total", n_cambia_ppal],
                        interanual[periodo == "Total", n_barcos],
                        interanual[periodo == "Total", pct_cambia_ppal])) +
  tema
ggsave("figs/cambio_arte_interanual.png", p6, width = 8, height = 5, dpi = 200)

message("Listo: outputs/CSV/cambio_arte_*.csv y figs/cambio_arte_*.png")

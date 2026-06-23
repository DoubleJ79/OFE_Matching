## =============================================================================
## make_deck_figures.R
##   ONE reproducible script for every figure used in deck/build.js
##   (CEM_for_Agronomists.pptx). Replaces the scattered inline scripts that
##   originally produced fig2,3,5-15 — those were never saved, so the deck was
##   not reproducible. This rebuilds all of them from Home.gpkg with ONE coherent
##   confounder story and prints a VALIDATION block of every number the slides cite.
##
## Confounder story (kept consistent across the deck):
##   * primary estimate / balance / binning : RSP + ApDepth        (figs 5,6,7,8,14,15)
##   * "add LS to show dropping"             : RSP + ApDepth + LS   (figs 2,3  — slides 6 & 8)
##   * model robustness sweep                : {RSP, RSP+ApDepth, RSP+ApDepth+LS, LS}
##                                                                  (figs 10,11)
##
## Run:  "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" make_deck_figures.R
## Figures are written straight into this folder, where build.js reads them.
## =============================================================================

suppressPackageStartupMessages({
  library(sf); library(terra); library(lwgeom); library(nngeo)
  library(MatchIt); library(sandwich); library(cobalt); library(ggplot2)
})

OUT <- "C:/Users/crop/OneDrive - University of Guelph/R code/Caleb/OFE_Stats/CEMvsPSM"
gp  <- file.path(OUT, "Home", "Home.gpkg")
rd  <- function(L){ r <- tryCatch(terra::rast(gp, subds = L), error = function(e) NULL)
  if (is.null(r) || terra::ncell(r) == 0) r <- terra::rast(sprintf("GPKG:%s:%s", gp, L)); r }
GREEN <- "#1b9e77"; ORANGE <- "#d95f02"; PURPLE <- "#7570b3"

## =============================================================================
## SHARED PREAMBLE — build the 48 trial plots once, extract every raster once
## =============================================================================
ab <- st_zm(st_read(gp, "ABLine", quiet = TRUE))
if (any(st_geometry_type(ab) == "MULTILINESTRING")) ab <- st_cast(ab, "LINESTRING")
ML  <- as.numeric(st_length(ab) / ceiling(st_length(ab) / 15))
Seg <- nngeo::st_segments(st_segmentize(ab, dfMaxLength = ML))
PP  <- st_collection_extract(lwgeom::st_split(st_sf(geometry = st_buffer(Seg, 24, endCapStyle = "FLAT")), st_geometry(ab)), "POLYGON")
off <- st_geometry(st_cast(st_boundary(st_buffer(ab, 12)), "MULTILINESTRING"))
PP  <- st_collection_extract(lwgeom::st_split(PP, off), "POLYGON")

co <- st_coordinates(st_geometry(ab)); s <- co[1, ]; e <- co[nrow(co), ]; lv <- c(e[1] - s[1], e[2] - s[2])
ctr <- st_coordinates(suppressWarnings(st_centroid(PP)))
PP$along <- ((ctr[, 1] - s[1]) * lv[1] + (ctr[, 2] - s[2]) * lv[2]) / sqrt(sum(lv^2))
PP$cx <- ctr[, 1]; PP$cy <- ctr[, 2]

yt <- st_transform(st_read(gp, "Yield_treat", quiet = TRUE), st_crs(PP))
ix <- st_intersects(PP, yt)
PP$Treat <- vapply(ix, function(j) if (length(j) == 0) NA_integer_ else as.integer(mean(yt$Rate[j] == "high") >= 0.5), integer(1))
PP$Yield <- vapply(ix, function(j) if (length(j) == 0) NA_real_  else mean(yt$Yield[j], na.rm = TRUE), numeric(1))

for (nm in c("RSP", "ApDepth", "LS", "Slope", "OM", "Clay")) {
  v <- tolower(gsub("ApDepth", "apdepth", gsub("Disp-CEC", "cec", nm)))
  r <- rd(nm); PP[[v]] <- terra::extract(r[[1]], terra::project(terra::vect(PP), terra::crs(r)), fun = mean, na.rm = TRUE)[, 2]
}

soil_vars <- c("rsp", "apdepth", "ls", "slope", "om", "clay")
keep <- is.finite(PP$Treat) & is.finite(PP$Yield) & Reduce(`&`, lapply(soil_vars, function(v) is.finite(PP[[v]])))
PP <- PP[keep, ]
d  <- as.data.frame(PP)[, c("Treat", "Yield", soil_vars)]
N  <- nrow(d); tr <- which(d$Treat == 1); ct <- which(d$Treat == 0)
naive <- mean(d$Yield[d$Treat == 1]) - mean(d$Yield[d$Treat == 0])
cq <- function(v, k) unique(quantile(v, seq(0, 1, length.out = k + 1), na.rm = TRUE))

cat(sprintf("\n================ DATA: %d plots (high %d / low %d).  Naive N response = %.1f bu/ac ================\n",
            N, length(tr), length(ct), naive))
if (N != 48) cat(sprintf("  !! WARNING: expected 48 plots, got %d — slide text that says '48' / '24 vs 24' may be off.\n", N))

## ---- shared estimator helpers ----------------------------------------------
pps    <- function(v) as.numeric(predict(glm(reformulate(v, "Treat"), d, family = binomial("probit")), type = "response"))
ate_ci <- function(md){ md$Tn <- as.integer(as.character(md$Treat)); f <- lm(Yield ~ Tn, md, weights = md$weights)
  cl <- if (!is.null(md$subclass)) md$subclass else seq_len(nrow(md)); V <- sandwich::vcovCL(f, cluster = cl)
  e <- unname(coef(f)["Tn"]); s <- sqrt(V["Tn", "Tn"]); c(ATE = e, SE = s, lo = e - 1.96*s, hi = e + 1.96*s) }
SETS <- list("RSP_only" = "rsp", "RSP+ApDepth" = c("rsp","apdepth"),
             "RSP+ApDepth+LS" = c("rsp","apdepth","ls"), "LS_only" = "ls"); lev <- names(SETS)

## =============================================================================
## FIG 14 — confounder screen  (slide 3).  PRINTS the SMD / partial-cor table so
##          the slide-3 wording ("which vars are confounders") is verifiable.
## =============================================================================
scr <- do.call(rbind, lapply(soil_vars, function(v){ x <- d[[v]]
  smd <- abs((mean(x[d$Treat==1]) - mean(x[d$Treat==0])) / sqrt((var(x[d$Treat==1]) + var(x[d$Treat==0]))/2))
  pc  <- abs(cor(residuals(lm(d$Yield ~ d$Treat)), residuals(lm(x ~ d$Treat))))
  data.frame(var = v, smd = smd, pc = pc) }))
scr$type <- ifelse(scr$smd>=0.1 & scr$pc>=0.1, "confounder",
            ifelse(scr$pc>=0.1, "yield-related, balanced",
            ifelse(scr$smd>=0.1, "imbalanced, not yield-related", "neither")))
cat("\n---- (fig14) Confounder screen: confounder = |SMD|>=0.1 AND |partial cor|>=0.1 ----\n")
print(transform(scr, smd = round(smd,3), pc = round(pc,3)), row.names = FALSE)

fs <- ggplot(scr, aes(smd, pc, colour = type)) +
  annotate("rect", xmin = 0.1, xmax = Inf, ymin = 0.1, ymax = Inf, fill = GREEN, alpha = 0.08) +
  geom_hline(yintercept = 0.1, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0.1, linetype = "dashed", colour = "grey60") +
  geom_point(size = 4) + ggrepel::geom_text_repel(aes(label = toupper(var)), size = 4.2, seg.color = "grey70") +
  scale_colour_manual(values = c(confounder = ORANGE, `yield-related, balanced` = PURPLE,
                                 `imbalanced, not yield-related` = "#999999", neither = "#bbbbbb"), name = NULL) +
  labs(x = "imbalance between strips  (|SMD|)", y = "relation to yield  (|partial cor|)",
       title = "What actually confounds this trial: RSP and ApDepth",
       subtitle = "A confounder must be BOTH imbalanced and yield-related (top-right).") +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(OUT, "fig14_screen.png"), fs, width = 9, height = 5.6, dpi = 130)

## =============================================================================
## FIG 5 — 3x3 binning grid on RSP x ApDepth  (slide 5b)
## =============================================================================
rq <- quantile(d$rsp, c(0,1/3,2/3,1)); aq <- quantile(d$apdepth, c(0,1/3,2/3,1))
d$rb <- as.integer(cut(d$rsp, rq, include.lowest = TRUE)); d$ab <- as.integer(cut(d$apdepth, aq, include.lowest = TRUE))
g <- expand.grid(rb = 1:3, ab = 1:3)
g$nT <- mapply(function(r,a) sum(d$Treat==1 & d$rb==r & d$ab==a), g$rb, g$ab)
g$nC <- mapply(function(r,a) sum(d$Treat==0 & d$rb==r & d$ab==a), g$rb, g$ab)
g$status <- ifelse(g$nT>0 & g$nC>0, "matched", ifelse(g$nT+g$nC>0, "dropped", "empty"))
g$xmin <- rq[g$rb]; g$xmax <- rq[g$rb+1]; g$ymin <- aq[g$ab]; g$ymax <- aq[g$ab+1]
g$xc <- (g$xmin+g$xmax)/2; g$yc <- (g$ymin+g$ymax)/2; g$lab <- sprintf("%dT %dC", g$nT, g$nC)
mc5 <- paste(g$rb, g$ab)[g$status=="matched"]; kept5 <- sum(paste(d$rb, d$ab) %in% mc5)
cat(sprintf("\n---- (fig5) 3x3 grid: matched=%d dropped=%d empty=%d | %d of %d plots in matched cells ----\n",
            sum(g$status=="matched"), sum(g$status=="dropped"), sum(g$status=="empty"), kept5, N))
d$grp <- factor(ifelse(d$Treat==1, "Treated (high N)", "Control (low N)"), levels = c("Control (low N)","Treated (high N)"))
p5 <- ggplot() +
  geom_rect(data = g, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=status), colour = "grey55") +
  scale_fill_manual(values = c(matched="#CDEBD6", dropped="#F7DAC6", empty="#F2F2F2"),
                    breaks = c("dropped","empty","matched"), labels = c("dropped (one side)","empty","matched subclass"), name = NULL) +
  geom_text(data = g, aes(xc, yc, label = lab), size = 4.1, colour = "grey25") +
  geom_point(data = d, aes(rsp, apdepth, shape = grp, colour = grp), size = 3, stroke = 1.1, fill = NA) +
  scale_shape_manual(values = c("Control (low N)"=21, "Treated (high N)"=16), name = NULL) +
  scale_colour_manual(values = c("Control (low N)"="#6E62A6", "Treated (high N)"="#1B9E77"), name = NULL) +
  labs(x = "RSP  ->  3 bins (low | mid | high)", y = "ApDepth  ->  3 bins (shallow | mid | deep)",
       title = "Binning the 2 confounders into 3 bins each = a 3 x 3 grid of subclasses",
       subtitle = sprintf("Each cell = one subclass (one RSP bin x one ApDepth bin); kept only if it holds BOTH a treated and a control.\nGreen = matched, orange = one side only (dropped), grey = empty.   %d of %d plots land in matched cells.", kept5, N)) +
  theme_minimal(base_size = 13) + theme(legend.position = "top", panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 11))
ggsave(file.path(OUT, "fig5_binning.png"), p5, width = 9.7, height = 7, dpi = 120)

## =============================================================================
## FIG 15 — the 3x3 RSP x ApDepth subclasses mapped onto the plots  (slide 5b2)
## =============================================================================
rbl <- cut(PP$rsp, rq, include.lowest = TRUE, labels = c("low","mid","high"))
abl <- cut(PP$apdepth, aq, include.lowest = TRUE, labels = c("shallow","mid","deep"))
pal <- c("low·shallow"="#C6DBEF","low·mid"="#6BAED6","low·deep"="#08519C",
         "mid·shallow"="#C7E9C0","mid·mid"="#74C476","mid·deep"="#006D2C",
         "high·shallow"="#FDD0A2","high·mid"="#FD8D3C","high·deep"="#A63603")
PP15 <- PP; PP15$cell <- factor(paste0(rbl, "·", abl), levels = names(pal))
cen0 <- st_centroid(st_union(st_geometry(PP15)))
vv <- c(co[nrow(co),1]-co[1,1], co[nrow(co),2]-co[1,2]); ang <- atan2(vv[2], vv[1])
Rm <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), 2, 2); rg <- function(g) (g - cen0) * Rm + cen0
st_geometry(PP15) <- st_sfc(rg(st_geometry(PP15))); abr <- st_sfc(rg(st_geometry(ab)))
p15 <- ggplot(PP15) + geom_sf(aes(fill = cell), colour = "white", linewidth = 0.18) +
  geom_sf(data = st_sf(geometry = abr), colour = "grey15", linewidth = 0.6) +
  scale_fill_manual(values = pal, drop = TRUE, name = "Subclass  (RSP bin · ApDepth bin)") +
  guides(fill = guide_legend(nrow = 1, title.position = "top", title.hjust = 0.5)) +
  theme_void(base_size = 12) + theme(legend.position = "bottom", legend.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"), plot.margin = margin(2, 8, 2, 8))
ggsave(file.path(OUT, "fig15_binmap.png"), p15, width = 12, height = 2.2, dpi = 150)

## =============================================================================
## FIG 7 + FIG 8 — balance on RSP+ApDepth, BOTH engines (slide 5c)
##   fig7 love plot shows Before / After CEM (3 bins, the headline method) /
##   After full matching: CEM balances ApDepth and improves RSP, but its coarse
##   bins leave RSP ~0.13; full matching tightens BOTH under 0.1. fig8 eCDF shows
##   ApDepth overlap after CEM (ApDepth does balance under CEM).
## =============================================================================
m_cem  <- matchit(Treat ~ rsp + apdepth, data = d, method = "cem",  estimand = "ATE",
                  cutpoints = list(rsp = cq(d$rsp, 3), apdepth = cq(d$apdepth, 3)))
m_full <- matchit(Treat ~ rsp + apdepth, data = d, method = "full", estimand = "ATE",
                  distance = pps(c("rsp","apdepth")))
.smdw <- function(w, v){ t <- d$Treat; x <- d[[v]]; s <- sqrt((var(x[t==1]) + var(x[t==0]))/2)
  abs(weighted.mean(x[t==1], w[t==1]) - weighted.mean(x[t==0], w[t==0])) / s }
cat(sprintf("\n---- (fig7) |SMD| RSP/ApDepth — before %.3f/%.3f | after CEM(3 bins) %.3f/%.3f (kept %d) | after full %.3f/%.3f (kept %d) ----\n",
            .smdw(rep(1,N),"rsp"), .smdw(rep(1,N),"apdepth"),
            .smdw(m_cem$weights,"rsp"),  .smdw(m_cem$weights,"apdepth"),  sum(m_cem$weights>0),
            .smdw(m_full$weights,"rsp"), .smdw(m_full$weights,"apdepth"), sum(m_full$weights>0)))
lp <- love.plot(Treat ~ rsp + apdepth, data = d,
                weights = list("After CEM (3 bins)" = m_cem$weights, "After full matching" = m_full$weights),
                stats = "mean.diffs", abs = TRUE, thresholds = c(m = .1), s.d.denom = "pooled",
                var.names = c(rsp = "RSP (slope position)", apdepth = "ApDepth (depth to layer)"),
                colors = c("#9A9A9A", ORANGE, GREEN), shapes = c("circle","triangle","diamond"), size = 4,
                sample.names = c("Before matching","After CEM (3 bins)","After full matching"),
                title = "Matching recovers the balance randomization would give") +
  theme(legend.position = "top")
ggsave(file.path(OUT, "fig7_loveplot.png"), lp, width = 8.4, height = 4.0, dpi = 130)
bp <- bal.plot(m_cem, var.name = "apdepth", which = "both", type = "ecdf", colors = c(PURPLE, GREEN)) +
  labs(title = "ApDepth, treated vs control",
       subtitle = "Before: the high-N strip sits on shallower soil. After CEM matching: the curves overlap, as a randomized trial would.") +
  theme_minimal(base_size = 12)
ggsave(file.path(OUT, "fig8_ecdf.png"), bp, width = 9, height = 4.5, dpi = 130)

## =============================================================================
## FIG 11 + FIG 10 — four-model comparison & bin-sensitivity  (slides 7 & 7b)
## =============================================================================
cem_k <- function(v, k){ cps <- setNames(lapply(v, function(x) unique(quantile(d[[x]], seq(0,1,length.out = k+1)))), v)
  mm <- matchit(reformulate(v, "Treat"), d, method = "cem", estimand = "ATE", cutpoints = cps)
  c(ate_ci(match.data(mm)), ret = 100*sum(mm$weights > 0)/N) }
psm <- function(v) ate_ci(match.data(matchit(reformulate(v, "Treat"), d, method = "full", estimand = "ATE", distance = pps(v))))

ca <- do.call(rbind, lapply(lev, function(nm) do.call(rbind, lapply(3:7, function(k){
  a <- cem_k(SETS[[nm]], k); data.frame(model = nm, bins = k, ATE = a["ATE"], lo = a["lo"], hi = a["hi"], ret = a["ret"]) }))))
ca$model <- factor(ca$model, levels = lev)
cmp <- do.call(rbind, lapply(lev, function(nm){ cm <- ca[ca$model==nm, ]; pr <- psm(SETS[[nm]])
  rbind(data.frame(model = nm, engine = "CEM (mean of sweep)", ATE = mean(cm$ATE), lo = mean(cm$lo), hi = mean(cm$hi)),
        data.frame(model = nm, engine = "PSM (full match)",    ATE = pr["ATE"], lo = pr["lo"], hi = pr["hi"])) }))
cmp$model <- factor(cmp$model, levels = rev(lev))
cat("\n---- (fig11) ATE per model x engine ----\n"); print(transform(cmp, ATE = round(ATE,1), lo = round(lo,1), hi = round(hi,1)), row.names = FALSE)
cat(sprintf("---- (fig10) CEM bin-sweep ATE range across all models/bins: %.1f to %.1f bu/ac ----\n",
            min(ca$ATE), max(ca$ATE)))

p11 <- ggplot(cmp, aes(ATE, model, colour = engine)) +
  geom_vline(xintercept = naive, linetype = "dashed", colour = "grey55") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  scale_colour_manual(values = c("CEM (mean of sweep)" = GREEN, "PSM (full match)" = ORANGE), name = NULL) +
  labs(x = "N response (bu/ac, 95% CI)", y = NULL, title = "The four competing models: N response, CEM vs PSM",
       subtitle = sprintf("Dashed = naive (%.0f). All four land in the same window regardless of covariates or engine.", naive)) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(OUT, "fig11_models5.png"), p11, width = 9, height = 4.8, dpi = 130)

p10 <- ggplot(ca, aes(bins, ATE)) +
  geom_hline(yintercept = naive, linetype = "dotted", colour = "grey55") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.13, fill = GREEN) + geom_line(colour = GREEN) +
  geom_point(aes(size = ret), colour = GREEN) +
  scale_size_continuous(range = c(1,5), name = "% plots\nretained", limits = c(0,100), breaks = c(25,50,75,100)) +
  facet_wrap(~model, nrow = 1) + scale_x_continuous(breaks = 3:7) +
  labs(x = "CEM bins per variable", y = "N response (bu/ac, 95% CI)",
       title = "More bins & more confounders shrink the matched sample — the estimand drifts",
       subtitle = "Point size = % retained. Every model & engine holds ~42-46; individual estimates get noisier as bins rise and retention falls.") +
  theme_minimal(base_size = 12) + theme(axis.text.x = element_text(size = 8))
ggsave(file.path(OUT, "fig10_binsens.png"), p10, width = 11, height = 4.6, dpi = 130)

## =============================================================================
## FIG 2 + FIG 3 — transects & 3D drops, CEM on RSP+ApDepth+LS @ 2 bins  (slides 8 & 6)
##   FIG 2 legend fixed: dropped plots render as an × (shape+colour both mapped).
## =============================================================================
m23 <- matchit(Treat ~ rsp + apdepth + ls, data = d, method = "cem", estimand = "ATE",
               cutpoints = list(rsp = cq(d$rsp,2), apdepth = cq(d$apdepth,2), ls = cq(d$ls,2)))
sc23 <- m23$subclass; w23 <- m23$weights; PP$drop <- w23 == 0
cat(sprintf("\n---- (fig3) RSP+ApDepth+LS @ 2 bins: %d of %d plots dropped (no like-for-like) ----\n", sum(PP$drop), N))

ni  <- st_nearest_feature(st_centroid(PP[tr, ]), st_centroid(PP[ct, ])); raw <- PP$Yield[tr] - PP$Yield[ct][ni]
mat <- vapply(tr, function(i){ cc <- ct[which(sc23[ct]==sc23[i] & !is.na(sc23[ct]) & w23[ct] > 0)]
  if (!length(cc)) return(NA_real_); dd <- (PP$cx[cc]-PP$cx[i])^2 + (PP$cy[cc]-PP$cy[i])^2
  PP$Yield[i] - PP$Yield[cc[which.min(dd)]] }, numeric(1))
td <- rbind(data.frame(panel = "Raw: high-N plot vs adjacent low-N control", along = PP$along[tr], delta = raw),
            data.frame(panel = "Matched: nearest RSP+ApDepth+LS control (× = none in stratum)", along = PP$along[tr], delta = mat))
td$panel  <- factor(td$panel, levels = unique(td$panel))
td$status <- ifelse(is.na(td$delta), "dropped", ifelse(td$delta < 0, "negative", "positive"))
td$y_plot <- ifelse(is.na(td$delta), 0, td$delta)
cat(sprintf("---- (fig2) dropped in bottom panel = %d ; negative deltas = %d ----\n",
            sum(td$status=="dropped"), sum(td$status=="negative")))
p2 <- ggplot(td, aes(along, y_plot)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(data = subset(td, status != "dropped"), colour = "grey75") +
  geom_point(aes(colour = status, shape = status), size = 2.4) +
  facet_wrap(~panel, ncol = 1) +
  scale_colour_manual(values = c(positive = GREEN, negative = ORANGE, dropped = "firebrick"), name = NULL) +
  scale_shape_manual(values = c(positive = 16, negative = 16, dropped = 4), name = NULL) +
  labs(x = "distance along the field (m)", y = "treated - control delta (bu/ac)",
       title = "N response along the field: raw neighbour vs matched control",
       subtitle = "Top: each high-N plot vs the low-N plot across the line. Bottom: vs its nearest RSP+ApDepth+LS-matched control.") +
  theme_minimal(base_size = 12)
ggsave(file.path(OUT, "fig2_transects.png"), p2, width = 9.5, height = 7, dpi = 130)

## fig3 — extruded 3D DEM with kept/dropped plots
dem <- rd("LiDAR DEM")[[1]]
demc <- terra::crop(terra::project(dem, terra::crs(PP)), terra::ext(terra::vect(PP)) * 1.25)
demc <- terra::aggregate(demc, fact = max(1, floor(min(dim(demc)[1:2]) / 110)), fun = "mean", na.rm = TRUE)
M <- as.matrix(demc, wide = TRUE); Z <- t(M[nrow(M):1, ])
xs <- seq(terra::xmin(demc), terra::xmax(demc), length.out = nrow(Z)); ys <- seq(terra::ymin(demc), terra::ymax(demc), length.out = ncol(Z))
zmin <- min(Z, na.rm = TRUE); EX <- 7; ez <- function(z) (z - zmin) * EX + zmin; Ze <- ez(Z); dz <- 0.32 * diff(range(Ze, na.rm = TRUE))
zf <- (Ze[-1,-1] + Ze[-1,-ncol(Ze)] + Ze[-nrow(Ze),-1] + Ze[-nrow(Ze),-ncol(Ze)]) / 4
fcol <- terrain.colors(80)[cut(as.vector(zf), 80)]
png(file.path(OUT, "fig3_dem3d.png"), width = 1150, height = 860, res = 130)
pm <- persp(xs, ys, Ze, theta = 40, phi = 28, expand = 1, scale = FALSE, border = NA, col = fcol, shade = 0.45, ltheta = -55, box = FALSE,
            main = "The 48-plot trial draped on the LiDAR DEM  -  RSP+ApDepth+LS CEM: which plots drop")
cen <- st_coordinates(suppressWarnings(st_centroid(PP))); evz <- terra::extract(demc, cen[, 1:2, drop = FALSE]); zc <- ez(evz[[ncol(evz)]])
proj <- trans3d(cen[,1], cen[,2], zc, pm); ord <- order(proj$y, decreasing = TRUE)
rings <- lapply(st_geometry(PP), function(g){ cc <- st_coordinates(g); cc[, 1:2, drop = FALSE] })
for (i in ord) { r <- rings[[i]]; n <- nrow(r)
  b <- trans3d(r[,1], r[,2], rep(zc[i], n), pm); t <- trans3d(r[,1], r[,2], rep(zc[i]+dz, n), pm)
  fl <- if (PP$drop[i]) "#d6453b" else GREEN; wl <- if (PP$drop[i]) "#9e2d25" else "#147a59"
  for (k in 1:(n-1)) polygon(c(b$x[k], b$x[k+1], t$x[k+1], t$x[k]), c(b$y[k], b$y[k+1], t$y[k+1], t$y[k]), col = wl, border = NA)
  polygon(t$x, t$y, col = fl, border = "grey15", lwd = 0.4) }
legend("topright", legend = c("kept (comparable)", "dropped (no like-for-like)"), fill = c(GREEN, "#d6453b"), border = "grey15", bty = "n", cex = 0.95)
dev.off()

## =============================================================================
## FIG 9 — the field around the trial: yield points + plots  (slide 2b)
## =============================================================================
yl  <- st_transform(st_read(gp, "Yield", quiet = TRUE), st_crs(PP))
yf  <- intersect(c("Yield","yield","Yld"), names(yl))[1]
bb  <- st_bbox(st_buffer(st_union(PP), 70)); ylc <- suppressWarnings(st_crop(yl, bb))
ym <- ggplot() + geom_sf(data = ylc, aes(colour = .data[[yf]]), size = 0.6) +
  scale_colour_viridis_c(option = "viridis", name = "Yield") +
  geom_sf(data = st_geometry(PP), fill = NA, colour = "red", linewidth = 0.5) +
  labs(title = "The field around the trial: yield tracks the ground",
       subtitle = "Combine yield (points) with the 48 trial plots (red). The strips cross a lower-yielding eroded rise.") +
  theme_minimal(base_size = 12) + theme(axis.text = element_text(size = 7))
ggsave(file.path(OUT, "fig9_yieldmap.png"), ym, width = 8.6, height = 6.4, dpi = 130)

## =============================================================================
## FIG 12 + FIG 13 — cross-strip elevation gradient = ApDepth  (slide 7c)
## =============================================================================
trI <- which(PP$Treat == 1); ctI <- which(PP$Treat == 0)
niE <- st_nearest_feature(st_centroid(PP[trI, ]), st_centroid(PP[ctI, ]))
# elevation per plot from DEM (project/extract on plot polygons)
PP$elev <- terra::extract(terra::project(dem, terra::crs(PP)), terra::vect(PP), fun = mean, na.rm = TRUE)[, 2]
PP$sdiff <- NA_real_; PP$sdiff[trI] <- PP$elev[trI] - PP$elev[ctI][niE]
hi <- PP[trI, ]; hi$adiff <- PP$apdepth[trI] - PP$apdepth[ctI][niE]; hi$dyield <- PP$Yield[trI] - PP$Yield[ctI][niE]
demc2 <- terra::crop(terra::project(dem, terra::crs(PP)), terra::vect(st_buffer(st_union(PP), 45)))
cont  <- sf::st_as_sf(terra::as.contour(demc2, nlevels = 16))
m12 <- ggplot() + geom_sf(data = cont, colour = "grey72", linewidth = 0.3) +
  geom_sf(data = PP[ctI, ], fill = "grey88", colour = "grey55", linewidth = 0.3) +
  geom_sf(data = hi, aes(fill = sdiff), colour = "grey25", linewidth = 0.35) +
  scale_fill_gradient2(low = "#2c7bb6", mid = "grey92", high = "#d7191c", midpoint = 0, name = "high-N elev.\nminus control (m)") +
  coord_sf(expand = FALSE) +
  labs(title = "Where the strips aren't iso-elevation, high-N sits upslope (red)",
       subtitle = "LiDAR contours + trial plots; red high-N plots sit above their nearest control.") +
  theme_minimal(base_size = 12) + theme(axis.text = element_text(size = 6))
ggsave(file.path(OUT, "fig12_gradientmap.png"), m12, width = 7.6, height = 6.6, dpi = 130)

mu <- mean(hi$dyield[hi$sdiff > 0]); md <- mean(hi$dyield[hi$sdiff <= 0]); ccor <- cor(hi$sdiff, hi$dyield)
n_up <- sum(hi$sdiff > 0)
cat(sprintf("\n---- (fig13) %d of %d high-N plots upslope of control; cor(elev gap, delta)=%.2f; upslope mean=%.1f, iso/down mean=%.1f ----\n",
            n_up, length(trI), ccor, mu, md))
q13 <- ggplot(hi, aes(sdiff, dyield)) + geom_vline(xintercept = 0, colour = "grey70") +
  geom_smooth(method = "lm", se = FALSE, colour = "grey40", linewidth = 0.7) + geom_point(aes(colour = adiff), size = 3.3) +
  scale_colour_gradient2(low = "#d7191c", mid = "grey80", high = "#2c7bb6", midpoint = 0, name = "high-N ApDepth\nminus control") +
  annotate("text", x = max(hi$sdiff), y = max(hi$dyield), hjust = 1, vjust = 1, size = 4.1, colour = "grey20",
           label = sprintf("cor = %.2f\nupslope = %.1f\niso/down = %.1f", ccor, mu, md)) +
  labs(x = "cross-strip elevation gap: high-N minus its control (m)   [>0 = high-N upslope]", y = "raw delta yield (bu/ac)",
       title = "Smaller delta where high-N is upslope — it's on shallower soil than its control",
       subtitle = "Colour: red = high-N shallower than its control. Upslope pairs are red and low-delta.") +
  theme_minimal(base_size = 12)
ggsave(file.path(OUT, "fig13_gradient_delta.png"), q13, width = 8.8, height = 5.4, dpi = 130)

## =============================================================================
## FIG 6 — sensitivity to hidden bias (E-value & Rosenbaum critical Gamma) slide bk
## =============================================================================
fig6_done <- tryCatch({
  suppressPackageStartupMessages({ library(EValue); library(sensitivityfull); library(senstrat); library(gridExtra); library(grid) })
  sd_y <- sd(d$Yield)
  ev_pt <- function(est, se) round(unname(evalues.OLS(est = est, se = se, sd = sd_y, delta = 1, true = 0)["E-values","point"]), 1)
  bsenfm <- function(mdf){ rows <- list(); t1 <- logical(0)
    for (sc in unique(mdf$subclass)) { ss <- mdf[mdf$subclass==sc, ]; trr <- as.integer(as.character(ss$Treat))
      ytv <- ss$Yield[trr==1]; ycv <- ss$Yield[trr==0]; if (!length(ytv) || !length(ycv)) next
      if (length(ytv)==1){ rows[[length(rows)+1]] <- c(ytv, ycv); t1 <- c(t1, TRUE) }
      else if (length(ycv)==1){ rows[[length(rows)+1]] <- c(ycv, ytv); t1 <- c(t1, FALSE) } else next }
    if (length(rows) < 2) return(NULL); J <- max(lengths(rows))
    list(y = t(sapply(rows, function(r) c(r, rep(NA, J-length(r))))), treated1 = t1) }
  crit <- function(pf){ if (pf(1) > 0.05) return("1"); lo <- 1; hi <- 2
    while (pf(hi) <= 0.05){ lo <- hi; hi <- hi*2; if (hi > 1000) return(">1e3") }
    while (hi-lo > 0.02){ mm <- (lo+hi)/2; if (pf(mm) <= 0.05) lo <- mm else hi <- mm }; as.character(round((lo+hi)/2, 1)) }
  g_psm <- function(mdf){ b <- bsenfm(mdf); if (is.null(b)) return(NA); crit(function(g) senfm(b$y, b$treated1, gamma = g, alternative = "greater")$pval) }
  g_cem <- function(mdf){ z <- as.integer(as.character(mdf$Treat)); st <- as.integer(mdf$subclass)
    sco <- mscores(mdf$Yield, z, st); crit(function(g) as.numeric(senstrat(sco, z, st, gamma = g, alternative = "greater")$Result["P-value"])) }
  S6 <- list("RSP only"="rsp","RSP+ApDepth"=c("rsp","apdepth"),"RSP+ApDepth+LS"=c("rsp","apdepth","ls"),"LS only"="ls")
  res <- do.call(rbind, lapply(names(S6), function(nm){ v <- S6[[nm]]
    mp <- match.data(matchit(reformulate(v,"Treat"), d, method = "full", estimand = "ATE", distance = pps(v))); ap <- ate_ci(mp)
    cps <- setNames(lapply(v, function(x) unique(quantile(d[[x]], seq(0,1,length.out = 4)))), v)
    mcd <- match.data(matchit(reformulate(v,"Treat"), d, method = "cem", estimand = "ATE", cutpoints = cps)); acv <- ate_ci(mcd)
    rbind(data.frame(Model = nm, Engine = "PSM", ATE = sprintf("%.1f", ap["ATE"]), `95% CI` = sprintf("(%.0f, %.0f)", ap["lo"], ap["hi"]),
                     `E-value` = ev_pt(ap["ATE"], ap["SE"]), `Crit. Γ` = g_psm(mp), check.names = FALSE),
          data.frame(Model = "", Engine = "CEM", ATE = sprintf("%.1f", acv["ATE"]), `95% CI` = sprintf("(%.0f, %.0f)", acv["lo"], acv["hi"]),
                     `E-value` = ev_pt(acv["ATE"], acv["SE"]), `Crit. Γ` = g_cem(mcd), check.names = FALSE)) }))
  cat("\n---- (fig6) Sensitivity to hidden bias ----\n"); print(res, row.names = FALSE)
  gam <- suppressWarnings(as.numeric(res$`Crit. Γ`)); gam <- gam[is.finite(gam)]
  cat(sprintf("---- (fig6) smallest critical Gamma across all models/engines = %.1f ----\n", min(gam)))
  prim <- which(res$Model == "RSP+ApDepth"); prim <- c(prim, prim+1)
  fillm <- ifelse(res$Engine == "PSM", "#e8f3ec", "#fdeee4"); fillm[prim] <- ifelse(res$Engine[prim]=="PSM", "#cfe8d6", "#fbdcc4")
  ff <- matrix("plain", nrow(res), ncol(res)); ff[prim, ] <- "bold"
  tt <- ttheme_minimal(base_size = 12, core = list(bg_params = list(fill = fillm, col = "grey85"),
        fg_params = list(hjust = 0.5, x = 0.5, fontface = ff)), colhead = list(fg_params = list(fontface = "bold")))
  gt <- tableGrob(res, rows = NULL, theme = tt)
  png(file.path(OUT, "fig6_sensitivity.png"), width = 1320, height = 560, res = 132); grid.newpage(); grid.draw(gt)
  grid.text("RSP+ApDepth (bold) targets the actual confounders — the primary estimate.", y = 0.045,
            gp = gpar(fontsize = 10, col = "grey35", fontface = "italic")); dev.off()
  TRUE
}, error = function(e){ cat("\n!! (fig6) skipped:", conditionMessage(e), "\n"); FALSE })

cat("\n================ ALL FIGURES WRITTEN TO PROJECT FOLDER ================\n")
cat(sprintf("fig6 (sensitivity) generated: %s\n", fig6_done))

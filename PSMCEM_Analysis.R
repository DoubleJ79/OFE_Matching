## =============================================================================
## PSMCEM_Total  —  one script, two stages
## -----------------------------------------------------------------------------
## STAGE 1 (data prep + main analysis, formerly "PSMCEM_Process"):
##   * load AB line, the Yield_treat trial points, soil + covariate rasters from
##     ONE GeoPackage
##   * build 15 m plot polygons; assign Treat/Control and yield from Yield_treat
##     (Rate: high N = 210 lb/ac = Treated, low N = 30 lb/ac = Control), keeping
##     only the ~48 plots that are actually part of the rate trial
##   * average every covariate raster (RSP, LS, OM, Disp-CEC, Slope) into the plots
##   * fit a propensity score (PS_METHOD: non-spatial probit or SAR-probit) and
##     estimate the treatment effect with PSM full matching (doubly-robust ATE)
##
## STAGE 2 (design sensitivity, formerly "PSM_vs_CEM__Home"):
##   * compare the ATE and its STABILITY across four competing confounder SETS
##     ( RSP | RSP+ApDepth | RSP+ApDepth+LS | LS -- no Convergence Index ):
##       - PSM : full matching on a propensity score (PS_METHOD)   (reference)
##       - CEM : coarsened exact matching, swept over the number of quantile bins
##
## ESTIMAND CAVEAT (Stage 2): even with estimand = "ATE", CEM DROPS treated units
## in non-overlapping coarsened strata, so it estimates the ATE for the RETAINED
## subpopulation, which shifts with the bin count. PSM full matching drops nothing
## and keeps the full-sample ATE. A drifting ATE can be a shifting target
## population, not just noise — read it alongside n_dropped / n_drop_treated.
## =============================================================================


## =============================================================================
## 0.  PACKAGES
## =============================================================================
pkgs <- c(
  "MatchIt","spdep","spatialprobit","spatialreg","sf","nngeo","terra",
  "gstat","sp","ggplot2","dplyr","lwgeom","optmatch","marginaleffects","emmeans",
  "broom","cobalt","viridis","sandwich","gridExtra"
)
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
suppressPackageStartupMessages(invisible(lapply(pkgs, require, character.only = TRUE)))


## =============================================================================
## 1.  CONFIG  — the only block you should normally need to edit
## =============================================================================

# Project root that contains the per-site folders (each holding a <Site>.gpkg).
setwd("CEMvsPSM")

SITE <- "Home"                                 # "Home" or "Hunter"
SITE <- match.arg(SITE, c("Hunter", "Home"))

# Everything for a site now lives in a single GeoPackage: vectors as feature
# layers, rasters as 2d-gridded-coverage layers. `gpkg_layers()` returns the
# logical -> actual-layer-name mapping for the chosen site. Edit the right-hand
# names to match your GeoPackage; add raster layers as you finish importing them.
gpkg_layers <- function(SITE) {
  switch(
    SITE,
    "Home" = list(
      gpkg    = file.path("Home", "Home.gpkg"),
      abline  = "ABLine",
      yield   = "Yield",                       # yield point cloud (has a "Yield" field)
      soil    = "Soil Type",                   # vector polygon layer (handled in Stage 1)
      # covariate rasters: name on the LEFT becomes the plot column name.
      # Confounder = imbalanced across strips (|SMD| >= 0.1) AND yield-related
      # (|partial cor| >= 0.1). On this trial only RSP and ApDepth qualify -- see
      # the Stage-2 confounder screen. (Convergence Index deliberately NOT used.)
      rasters = c(
        rsp     = "RSP",                  # relative slope position - imbalanced + yield-related (a confounder)
        ls      = "LS",                   # length-slope factor     - balanced, not yield-related (inert here)
        dem     = "LiDAR DEM",            # elevation               - mildly imbalanced + yield-related
        apdepth = "ApDepth",              # depth to layer          - imbalanced + strongly yield-related (the PRIMARY confounder; |SMD| 0.14, partial cor +0.46)
        om      = "OM",                   # organic matter          - imbalanced but NOT yield-related here (not a confounder)
        cec     = "Disp-CEC",             # dispersible CEC         - not used in the current models
        slope   = "Slope",                # terrain slope           - balanced, not yield-related (inert here)
        clay    = "Clay"                  # % clay                  - yield-related but balanced (precision covariate; partial cor +0.49)
      )
    ),
    "Hunter" = list(
      gpkg    = file.path("Hunter", "Hunter.gpkg"),
      abline  = "NStrip_Centerline",           # treatment centre line
      yield   = "Yield",
      soil    = "SoilRaster",                  # categorical soil-zone raster
      rasters = c(
        rsp = "RSP",
        ls  = "LS",
        dem = "DEM"
      )
    )
  )
}
LYR <- gpkg_layers(SITE)

# Plot-polygon geometry.
PLOT_LENGTH  <- 15                             # along-line plot length (m)
BUFFER_WIDTH <- 24                             # full strip width (m); each side is 12 m

# Covariates to mean-CENTRE before matching.
CENTER_VARS <- c("ls", "rsp", "dem","apdepth","om","slope","clay")

# Stage-1 propensity-score covariates = the PSM baseline (causal anchor).
# Drive both the SAR-probit score and the doubly-robust outcome model.
# Each must be present as a raster in gpkg_layers()$rasters above.
PS_COVARS <- c("rsp", "apdepth", "clay")
PS_MODEL  <- reformulate(PS_COVARS, response = "Treat")

# Propensity-score engine, used by BOTH the Stage-1 main PSM and every Stage-2
# per-model PSM baseline:
#   "probit"    = non-spatial probit GLM  -> robust on small samples (default)
#   "sarprobit" = spatial SAR-probit      -> use on larger data sets; its MCMC
#                 needs an adequate N or it returns NaN betas (e.g. the 48-plot trial)
PS_METHOD <- "probit"
PS_METHOD <- match.arg(PS_METHOD, c("probit", "sarprobit"))

# Stage-2: the CEM models to compare. Each is matched on (multivariate CEM + PSM);
# the design-sensitivity (2) and adjacency-trap (3) are run per model so they can
# be compared head-to-head. Each variable must be an extracted column.
CONF_SETS <- list(
  "RSP_only"       = "rsp",
  "RSP+ApDepth"    = c("rsp", "apdepth"),
  "RSP+ApDepth+LS" = c("rsp", "apdepth", "ls"),
  "LS_only"        = "ls"
)
BIN_COUNTS <- 3:7                              # CEM bin-count sweep


## =============================================================================
## 2.  HELPERS
## =============================================================================

# --- Read one raster layer out of a GeoPackage --------------------------------
# terra addresses a gpkg raster table as a named subdataset; fall back to the
# GDAL "GPKG:file:table" connection string if the named lookup misses.
read_gpkg_raster <- function(gpkg, layer) {
  r <- tryCatch(terra::rast(gpkg, subds = layer), error = function(e) NULL)
  if (is.null(r) || terra::ncell(r) == 0)
    r <- tryCatch(terra::rast(sprintf("GPKG:%s:%s", gpkg, layer)), error = function(e) NULL)
  r                                                       # NULL if `layer` is not a raster
}

# --- Stop early with an actionable message if needed columns are missing ------
require_cols <- function(data, vars, what) {
  miss <- setdiff(vars, names(data))
  if (length(miss))
    stop(sprintf("%s needs column(s) [%s] which are not present.\n  Add the matching raster layer(s) to %s and list them in gpkg_layers()$rasters.",
                 what, paste(miss, collapse = ", "), LYR$gpkg), call. = FALSE)
  invisible(TRUE)
}

# --- Soil-zone code -> human-readable soil name (site specific) ---------------
soil_map_for <- function(SITE) {
  if (SITE == "Hunter")
    c(`1` = "London Loam degraded", `2` = "London Loam", `3` = "Guelph degraded",
      `4` = "Guelph", `5` = "Colwood", `6` = "Parkhill")
  else # "Home"
    c(`1` = "London Loam degraded", `2` = "London Loam", `3` = "Guelph degraded",
      `4` = "Parkhill")
}
recode_soil <- function(df, SITE) {
  lut <- soil_map_for(SITE)
  dplyr::mutate(df, SoilType = dplyr::recode(as.character(Soil_Zone), !!!lut,
                                             .default = as.character(Soil_Zone),
                                             .missing = NA_character_))
}


## =============================================================================
## STAGE 1 — DATA PREP
## =============================================================================

## ---- 3. Load inputs from the GeoPackage ------------------------------------
abLine <- sf::st_read(LYR$gpkg, layer = LYR$abline, quiet = TRUE)
yield  <- sf::st_read(LYR$gpkg, layer = LYR$yield,  quiet = TRUE)

# AB line must be a single 2-D LINESTRING for the segment/side logic below.
abLine <- sf::st_zm(abLine)
if (any(sf::st_geometry_type(abLine) == "MULTILINESTRING"))
  abLine <- sf::st_cast(abLine, "LINESTRING")

# Soil zone source. At Hunter this is a categorical raster (modal-extracted into
# each plot below); at Home it is a vector polygon layer, so a raster read fails
# and we fall back to a spatial join. Either way soil is informational only — it
# is not used by the matching — so a failure just leaves SoilType = NA.
soil_rast <- suppressWarnings(tryCatch(read_gpkg_raster(LYR$gpkg, LYR$soil), error = function(e) NULL))
soil_vec  <- if (is.null(soil_rast))
  tryCatch(sf::st_read(LYR$gpkg, layer = LYR$soil, quiet = TRUE), error = function(e) NULL) else NULL


## ---- 4. Build 15 m plot polygons from the AB line --------------------------
abLength <- sf::st_length(abLine)
NumSplits <- ceiling(abLength / PLOT_LENGTH)
MLength   <- as.numeric(abLength / NumSplits)          # actual (near-15 m) plot length

# Segmentize -> per-segment lines -> flat buffer = one long treated+control strip.
SegLine  <- nngeo::st_segments(sf::st_segmentize(abLine, dfMaxLength = MLength))
BuffLine <- sf::st_buffer(SegLine, dist = BUFFER_WIDTH, endCapStyle = "FLAT")
BuffLine <- sf::st_sf(geometry = BuffLine)
BuffLine$ID <- seq_len(nrow(BuffLine))

# Split the strip lengthwise by the AB line (treated vs control halves) ...
abLine_geom  <- sf::st_geometry(abLine)
PlotPolygons <- sf::st_collection_extract(lwgeom::st_split(BuffLine, abLine_geom), "POLYGON")

# ... then split again by a line offset half a buffer width to either side,
# giving the individual plot polygons.
offset_line <- sf::st_buffer(abLine, dist = BUFFER_WIDTH / 2) |>
  sf::st_boundary() |> sf::st_cast("MULTILINESTRING") |> sf::st_geometry()
PlotPolygons <- sf::st_collection_extract(
  lwgeom::st_split(PlotPolygons, offset_line), "POLYGON")

PlotPolygons$PID     <- seq_len(nrow(PlotPolygons))
PlotPolygons$MLength <- round(MLength, 1)

## ---- 5. Assign treatment + yield from the Yield_treat trial layer ----------
# The real rate trial is recorded in the 'Yield_treat' point layer: each point
# carries a Rate (high N = 210 lb/ac = Treated, low N = 30 lb/ac = Control) and a
# Yield. A plot's treatment = the majority Rate of the points inside it; its yield
# = the mean of those points. Plots with no trial points are not part of the
# ~48-plot rate trial and are dropped downstream (section 7 filter).
# (Treatment is therefore data-driven, NOT inferred from side of the AB line.)
yt <- sf::st_read(LYR$gpkg, "Yield_treat", quiet = TRUE)
yt <- sf::st_transform(yt, sf::st_crs(PlotPolygons))
ix <- sf::st_intersects(PlotPolygons, yt)
PlotPolygons$Treat       <- vapply(ix, function(j) if (length(j) == 0) NA_integer_ else as.integer(mean(yt$Rate[j] == "high") >= 0.5), integer(1))
PlotPolygons$Yield_treat <- vapply(ix, function(j) if (length(j) == 0) NA_real_  else mean(yt$Yield[j], na.rm = TRUE), numeric(1))
PlotPolygons$Treatment   <- ifelse(is.na(PlotPolygons$Treat), NA_character_,
                                   ifelse(PlotPolygons$Treat == 1L, "Treated", "Control"))
PlotPolygons$Distance <- cumsum(PlotPolygons$MLength)
PlotPolygons <- PlotPolygons[order(PlotPolygons$PID), ]

## ---- 6. Yield comes from the Yield_treat trial layer (no IDW needed) --------
# Section 5 already averaged the Yield_treat point yields into each plot; build a
# SpatVector carrying that yield so the covariate rasters can be extracted onto it.
plots_extracted <- terra::vect(PlotPolygons)
plots_extracted$Yield <- PlotPolygons$Yield_treat

## ---- 7. Average soil + every covariate raster into the plots ---------------
# modal soil class, then each gpkg covariate raster (rsp, ls, om, cec, slope).
# Reproject the plots to each raster's CRS before extracting (the covariate
# rasters may not share the vector CRS, and terra::extract does not align it).
plots_extracted$Soil_Zone <- if (!is.null(soil_rast))
  as.integer(terra::extract(soil_rast, terra::project(PlotPolygons, terra::crs(soil_rast)),
                            fun = modal, na.rm = TRUE, bind = FALSE)[, 2]) else NA_integer_

for (col in names(LYR$rasters)) {
  r    <- read_gpkg_raster(LYR$gpkg, LYR$rasters[[col]])
  pv   <- terra::project(plots_extracted, terra::crs(r))  # align CRS before extracting
  vals <- terra::extract(r, pv, fun = mean, na.rm = TRUE)[, 2]
  plots_extracted[[col]] <- vals                          # column name = names(rasters)
}

# Back to sf, attach soil name, add centroid coords + numeric Treat (1/0).
plots_sf <- PlotPolygons
plots_sf <- cbind(plots_sf, plots_extracted[, setdiff(names(plots_extracted), names(plots_sf))])
# Soil name per plot: recode raster codes (Hunter), else nearest soil polygon (Home), else NA.
if (!is.null(soil_rast)) {
  plots_sf <- recode_soil(plots_sf, SITE)
} else if (!is.null(soil_vec)) {
  fld  <- intersect(c("Soil Type", "Soil_Type", "SoilType"), names(soil_vec))[1]
  cent <- sf::st_transform(suppressWarnings(sf::st_centroid(plots_sf)), sf::st_crs(soil_vec))
  plots_sf$SoilType <- as.character(soil_vec[[fld]])[sf::st_nearest_feature(cent, soil_vec)]
} else {
  plots_sf$SoilType <- NA_character_
}
cc <- sf::st_coordinates(suppressWarnings(sf::st_centroid(plots_sf)))
plots_sf$centroid_x <- cc[, 1]
plots_sf$centroid_y <- cc[, 2]
plots_sf$Treat <- ifelse(plots_sf$Treatment == "Treated", 1L,
                  ifelse(plots_sf$Treatment == "Control", 0L, NA_integer_))
plots_sf$id <- NULL

# Keep only the plots that are part of the rate trial (have a Rate + a yield).
map_sf <- plots_sf[is.finite(plots_sf$Treat) & is.finite(plots_sf$Yield), ]
cat(sprintf("Rate-trial plots kept: %d (Treated/high-N %d, Control/low-N %d) of %d built.\n",
            nrow(map_sf), sum(map_sf$Treat == 1), sum(map_sf$Treat == 0), nrow(plots_sf)))
# Persist the processed plots back to their own GeoPackage (vector + attributes).
sf::st_write(map_sf, sprintf("%s_plots.gpkg", SITE), layer = "plots",
             delete_dsn = TRUE, quiet = TRUE)


## =============================================================================
## STAGE 1 — SHARED SPATIAL SETUP  (used by both matching stages)
## =============================================================================
# spdep / sarprobit operate on a Spatial* object and contiguity weights.
map_sp <- as(map_sf, "Spatial")
coords <- sp::coordinates(map_sp)

neighbors     <- spdep::poly2nb(map_sp, queen = FALSE)        # rook  (effect estimation)
wAdj          <- spdep::nb2listw(neighbors, style = "W")
neighbors_SAR <- spdep::poly2nb(map_sp, queen = TRUE)         # queen (SAR probit)
wAdj_SAR      <- spdep::nb2listw(neighbors_SAR, style = "W")
w_Adj         <- as(spatialreg::as_dgRMatrix_listw(wAdj_SAR), "CsparseMatrix")

# Mean-centre covariates (centre only, original units kept).
# as.numeric() flattens scale()'s 1-column matrix, which otherwise trips cut()/quantile().
for (v in intersect(CENTER_VARS, names(map_sp@data)))
  map_sp@data[[v]] <- as.numeric(scale(map_sp@data[[v]], center = TRUE, scale = FALSE))

write.csv(map_sp@data, sprintf("%s_map_sp_data.csv", SITE))

## --- Propensity-score engine: spatial SAR-probit OR non-spatial probit -------
# One place that builds the propensity score, dispatched by PS_METHOD. Returns the
# fitted score `p`, the model `mod`, and the `method` used. SAR-probit borrows the
# queen-contiguity weights w_Adj; the probit GLM ignores them.
fit_pscore <- function(formula, data) {
  if (identical(PS_METHOD, "sarprobit")) {
    mod <- spatialprobit::sarprobit(formula, w_Adj, data)
    list(p = as.numeric(1 / (1 + exp(-fitted(mod)))), mod = mod, method = "sarprobit")
  } else {
    mod <- glm(formula, data = data, family = binomial("probit"))
    list(p = as.numeric(predict(mod, type = "response")), mod = mod, method = "probit")
  }
}


## =============================================================================
## STAGE 1 — MAIN ANALYSIS: SAR-probit PSM + doubly-robust ATE
## =============================================================================
# Needs the PS_COVARS columns. If any covariate raster is missing this block is
# skipped with a note, so the Stage-2 confounder-set comparison can still run.
if (!all(all.vars(PS_MODEL) %in% names(map_sp@data))) {
  message(sprintf(
    "Stage-1 main PSM SKIPPED: missing covariate column(s) [%s].\n  Add the matching raster layer(s) to %s and list them in gpkg_layers()$rasters.",
    paste(setdiff(all.vars(PS_MODEL), names(map_sp@data)), collapse = ", "), LYR$gpkg))
} else {

  # Propensity score; engine selected by PS_METHOD (non-spatial probit / SAR-probit).
  model   <- Treat ~ rsp + apdepth + clay
  ps_fit  <- fit_pscore(model, map_sp@data)
  mod_SAR <- ps_fit$mod; p_score <- ps_fit$p
  message("  Stage-1 propensity engine: ", ps_fit$method)
  summary(mod_SAR); if (identical(PS_METHOD, "sarprobit")) impacts(mod_SAR)
  logLik(mod_SAR); AIC(mod_SAR); BIC(mod_SAR)

  # Common-support check: propensity-score overlap between treated and control.
  # Good overlap (curves cover the same range) = the PSM comparison is supported.
  plot(density(p_score[map_sp@data$Treat == 1]), col = "blue", lwd = 2,
       main = "Propensity overlap (blue = treated, red = control)", xlab = "propensity score")
  lines(density(p_score[map_sp@data$Treat == 0]), col = "red", lwd = 2)

  # PSM full matching on that score (keeps every unit -> full-sample ATE).
  # NB: MatchIt's argument is `distance` (a vector of propensity scores); `pscores`
  # is NOT a matchit argument and would be silently ignored, leaving matchit to
  # estimate its own default logistic GLM. Pass the PS_METHOD score as `distance`.
  matchit_test <- MatchIt::matchit(PS_MODEL, data = map_sp@data,
                                   method = "full", distance = p_score)
  summary(matchit_test)
  plot(matchit_test, type = "jitter", interactive = FALSE)
  plot(summary(matchit_test))
  plot(matchit_test, type = "qq")
  plot(matchit_test, type = "ecdf")
  cobalt::bal.tab(matchit_test, un = TRUE)
  cobalt::love.plot(matchit_test, thresholds = 0.1)

  matched_data <- MatchIt::match.data(matchit_test)
  matched_data$Treat <- factor(matched_data$Treat)
  for (v in intersect(CENTER_VARS, names(matched_data)))
    matched_data[[v]] <- as.numeric(matched_data[[v]])

  # Doubly-robust effect: outcome model interacts Treat with the covariates, then
  # average the Treat contrast (cluster-robust by matched subclass). The formula
  # is built from PS_COVARS so it tracks whatever covariate set you configured.
  dr_formula <- reformulate(sprintf("Treat * (%s)", paste(PS_COVARS, collapse = " + ")),
                            response = "Yield")
  fit1 <- lm(dr_formula, data = matched_data, weights = weights)
  treat_effect_avgcomp <- marginaleffects::avg_comparisons(
    fit1, variables = "Treat", vcov = ~subclass, wts = "weights")
  print(treat_effect_avgcomp)
  summary(fit1)

  # Matched data as points (for mapping / export).
  md_df <- as.data.frame(matched_data)
  write.csv(md_df, sprintf("%s_matched_data.csv", SITE), row.names = FALSE)
  matched_data_sp <- sp::SpatialPointsDataFrame(
    coords      = as.matrix(md_df[, c("centroid_x", "centroid_y")]),
    data        = md_df,
    proj4string = sp::CRS(sf::st_crs(yield)$wkt))

  # Nested outcome models for AIC/BIC comparison (which covariates earn their keep):
  # Treat-only, each covariate added singly, the full additive set, and the full
  # Treat x covariate interaction. Built from PS_COVARS, so it adapts automatically.
  cand <- c(list(character(0)), as.list(PS_COVARS))
  if (length(PS_COVARS) > 1) cand <- c(cand, list(PS_COVARS))
  aic_fits <- lapply(cand, function(vars)
    lm(reformulate(c("Treat", vars), response = "Yield"), data = matched_data, weights = weights))
  names(aic_fits) <- vapply(cand, function(v) paste(c("Treat", v), collapse = " + "), character(1))
  aic_fits[["Treat * covars (full)"]] <- fit1
  aic_table <- data.frame(
    model = names(aic_fits),
    AIC   = vapply(aic_fits, AIC, numeric(1)),
    BIC   = vapply(aic_fits, BIC, numeric(1)), row.names = NULL)
  print(aic_table[order(aic_table$AIC), ])
}


## =============================================================================
## STAGE 2 — COMPETING CONFOUNDER SETS
##   (1) confounder screen + head-to-head set comparison
##   (2) design sensitivity (PSM vs CEM, ATE vs bins) per set
##   (3) adjacency trap per set
## =============================================================================
CONF_SETS <- CONF_SETS[vapply(CONF_SETS, function(v) all(v %in% names(map_sp@data)), logical(1))]
if (!length(CONF_SETS)) stop("No CONF_SETS fully present; add the raster(s) to gpkg_layers()$rasters.")

## ---- Correlation matrix: outcome + all covariates --------------------------
cor_vars <- intersect(c("Yield", "apdepth", "cec", "om", "slope", "ls", "rsp", "dem", "clay"), names(map_sp@data))
cat("\n=== Correlation matrix (Yield + covariates) ===\n")
print(round(cor(map_sp@data[, cor_vars], use = "complete.obs"), 3))

## ---- Confounder screen: imbalance x yield relevance x effect on the estimate
# A variable biases the strip comparison only if it is BOTH imbalanced across the
# strips AND related to yield:
#   SMD_before  = strip imbalance (treated vs control, unmatched)
#   SMD_after   = imbalance after CEM on that variable (matching balances it)
#   partial_cor = relation to yield, ADJUSTED for Treatment (the N rate)
#   shift       = how much matching on it moves the ATE (the bias it carries)
# So: imbalanced but shift ~ 0 (e.g. RSP) -> not a confounder; high shift -> it is.
naive_ate <- mean(map_sp@data$Yield[map_sp@data$Treat == 1]) - mean(map_sp@data$Yield[map_sp@data$Treat == 0])
.smd <- function(x, tr, w) { kt <- tr == 1 & w > 0 & is.finite(x); kc <- tr == 0 & w > 0 & is.finite(x)
  s <- sqrt((var(x[tr == 1], na.rm = TRUE) + var(x[tr == 0], na.rm = TRUE)) / 2)
  if (is.na(s) || s == 0) return(NA_real_); (weighted.mean(x[kt], w[kt]) - weighted.mean(x[kc], w[kc])) / s }
screen_vars <- intersect(c("apdepth", "cec", "om", "slope", "ls", "rsp", "dem", "clay"), names(map_sp@data))
confounder_screen <- do.call(rbind, lapply(screen_vars, function(v) {
  x <- map_sp@data[[v]]; tr <- map_sp@data$Treat; ok <- is.finite(x)
  pr <- tryCatch(cor(residuals(lm(map_sp@data$Yield[ok] ~ tr[ok])), residuals(lm(x[ok] ~ tr[ok]))), error = function(e) NA_real_)
  brk <- unique(quantile(x[ok], seq(0, 1, length.out = 6), na.rm = TRUE))
  ate <- tryCatch({
    m  <- MatchIt::matchit(Treat ~ x, data = data.frame(Treat = tr, x = x, Yield = map_sp@data$Yield)[ok, ],
                           method = "cem", estimand = "ATE", cutpoints = list(x = brk))
    wf <- rep(0, length(x)); wf[ok] <- m$weights
    list(ate = unname(coef(lm(Yield ~ Treat, MatchIt::match.data(m), weights = MatchIt::match.data(m)$weights))["Treat"]),
         smd_after = .smd(x, tr, wf))
  }, error = function(e) list(ate = NA_real_, smd_after = NA_real_))
  data.frame(variable = v, SMD_before = round(.smd(x, tr, rep(1, length(x))), 3),
             SMD_after = round(ate$smd_after, 3), partial_cor_yield = round(pr, 3),
             ATE_if_matched = round(ate$ate, 2), shift_vs_naive = round(ate$ate - naive_ate, 2), row.names = NULL)
}))
confounder_screen <- confounder_screen[order(-abs(confounder_screen$shift_vs_naive)), ]
cat(sprintf("\n=== Confounder screen (naive ATE = %.2f bu/ac). Confounder = imbalanced AND yield-related ===\n", naive_ate))
print(confounder_screen); gridExtra::grid.table(confounder_screen)

# Auto-detect the ACTUAL confounders from the screen, so the balance yardstick
# below never has to be remembered or hand-curated: a variable counts only if it
# is imbalanced across the strips (|SMD_before| >= 0.1) AND related to yield
# (|partial_cor_yield| >= 0.1). Bias comes only from these; balancing anything
# else buys comparability/precision, not bias removal.
.imb <- abs(confounder_screen$SMD_before)        >= 0.1
.rel <- abs(confounder_screen$partial_cor_yield) >= 0.1
CONFOUNDERS <- intersect(confounder_screen$variable[which(.imb & .rel)], names(map_sp@data))
cat(sprintf("Auto-detected confounders (imbalanced & yield-related): %s\n",
            if (length(CONFOUNDERS)) paste(CONFOUNDERS, collapse = ", ") else "(none cleared both thresholds)"))

## ---- Matching diagnostics (shared by PSM and every CEM variant) ------------
# ATE + 95% CI: the Treat coefficient in the weighted lm IS the weighted
# difference in means; the SE is cluster-robust by matched set (sandwich::vcovCL).
ate_ci <- function(m) {
  m$Tnum <- ifelse(as.character(m$Treat) %in% c("1", "treated"), 1, 0)
  fit <- lm(Yield ~ Tnum, data = m, weights = m$weights)
  cl  <- if (!is.null(m$subclass)) m$subclass else seq_len(nrow(m))
  V   <- sandwich::vcovCL(fit, cluster = cl)
  est <- unname(coef(fit)["Tnum"]); se <- sqrt(V["Tnum", "Tnum"]); z <- qnorm(0.975)
  c(ATE = est, lo = est - z * se, hi = est + z * se)
}

# n_dropped (weight-0 units), n_drop_treated (shifts the estimand), ess_control
# (control effective sample size = how heavily controls are reused).
match_diag <- function(m) {
  tr <- as.integer(as.character(m$treat)); w <- m$weights
  wc <- w[tr == 0 & w > 0]; ess_c <- if (length(wc)) sum(wc)^2 / sum(wc^2) else 0
  c(n_dropped = sum(w == 0), n_drop_treated = sum(w == 0 & tr == 1),
    n_drop_control = sum(w == 0 & tr == 0), pct_retained = round(100 * mean(w > 0), 1),
    n_control_used = length(wc), ess_control = round(ess_c, 1))
}

# Multivariate L1 imbalance on a variable SET, on a FIXED reference coarsening
# (5 bins/var, computed once on the full sample -> independent of the matching bins).
l1_ref <- setNames(lapply(unique(unlist(CONF_SETS)), function(v)
  unique(quantile(map_sp@data[[v]], seq(0, 1, length.out = 6), na.rm = TRUE))), unique(unlist(CONF_SETS)))
l1_set <- function(md, vars) {
  cells <- do.call(interaction, c(lapply(vars, function(v) cut(md[[v]], l1_ref[[v]], include.lowest = TRUE)), list(drop = FALSE)))
  tr <- as.integer(as.character(md$Treat)); w <- md$weights
  wt <- sum(w[tr == 1]); wc <- sum(w[tr == 0]); if (wt == 0 || wc == 0) return(NA_real_)
  ft <- tapply(w[tr == 1], cells[tr == 1], sum); ft[is.na(ft)] <- 0; ft <- ft / wt
  fc <- tapply(w[tr == 0], cells[tr == 0], sum); fc[is.na(fc)] <- 0; fc <- fc / wc
  round(0.5 * sum(abs(ft - fc)), 3)
}

## ---- CEM on a variable SET, k bins per variable ----------------------------
cem_set <- function(vars, k, return_match = FALSE) {
  cps <- setNames(lapply(vars, function(v) unique(quantile(map_sp@data[[v]], seq(0, 1, length.out = k + 1), na.rm = TRUE))), vars)
  m <- tryCatch(MatchIt::matchit(reformulate(vars, "Treat"), data = map_sp@data, method = "cem",
                                 estimand = "ATE", cutpoints = cps), error = function(e) NULL)
  row <- tryCatch({
    md <- MatchIt::match.data(m); s <- ate_ci(md); dc <- match_diag(m)
    data.frame(bins = k, ATE = unname(s["ATE"]), lo = unname(s["lo"]), hi = unname(s["hi"]),
               pct_retained = unname(dc["pct_retained"]), L1 = l1_set(md, vars),
               n_dropped = unname(dc["n_dropped"]), ess_control = unname(dc["ess_control"]), row.names = NULL)
  }, error = function(e) data.frame(bins = k, ATE = NA_real_, lo = NA_real_, hi = NA_real_,
               pct_retained = 0, L1 = NA_real_, n_dropped = NA_integer_, ess_control = 0))
  if (return_match) list(row = row, m = m) else row
}

## ---- PSM full matching on a variable SET (propensity engine = PS_METHOD) ----
psm_set <- function(vars) {
  ps <- fit_pscore(reformulate(vars, "Treat"), map_sp@data)$p
  m  <- MatchIt::matchit(reformulate(vars, "Treat"), data = map_sp@data, method = "full", estimand = "ATE", distance = ps)
  list(ate = ate_ci(MatchIt::match.data(m)), diag = match_diag(m), L1 = l1_set(MatchIt::match.data(m), vars), m = m)
}

## ---- Run all competing confounder sets -------------------------------------
k_delta <- min(BIN_COUNTS)                               # coarsest bins for the per-set delta match (most retained)
sets <- setNames(lapply(names(CONF_SETS), function(nm) {
  vars <- CONF_SETS[[nm]]
  cem  <- do.call(rbind, lapply(BIN_COUNTS, function(k) cem_set(vars, k))); cem$set <- nm
  psm  <- psm_set(vars)
  mid  <- cem_set(vars, k_delta, return_match = TRUE)$m
  cat(sprintf("Set %-14s: PSM ATE=%.1f  CEM ATE=%.1f (sd %.2f)  retained=%.0f%%  L1=%.3f\n",
              nm, psm$ate["ATE"], mean(cem$ATE, na.rm = TRUE), sd(cem$ATE, na.rm = TRUE),
              mean(cem$pct_retained, na.rm = TRUE), psm$L1))
  list(name = nm, vars = vars, cem = cem, psm = psm, mid = mid)
}), names(CONF_SETS))

## ---- (1) HEAD-TO-HEAD comparison of the competing models -------------------
# ONE ROW PER MODEL x ENGINE, so PSM and CEM metrics are never conflated (each
# model gets a PSM row and a CEM row; the model name is printed once per pair):
#   engine   = PSM (full matching, keeps every unit) or CEM (coarsened strata, bin sweep)
#   ATE/shift= that engine's ATE and its move from the naive estimate
#   bal_all  = residual max|SMD| over ALL covariates (BAL_VARS) after THAT engine's
#              matching -> fair "comparable on everything" score
#   bal_conf = residual max|SMD| over the AUTO-DETECTED confounders -> bias-relevant balance
#   ATE_sd   = CEM-only spread across the bin sweep (NA for PSM)
#   retained = % units kept (PSM = 100; CEM drops non-overlap strata); ESS = effective control N
# Lower is better for both balances.
naive_set  <- mean(map_sp@data$Yield[map_sp@data$Treat == 1]) - mean(map_sp@data$Yield[map_sp@data$Treat == 0])
BAL_VARS   <- intersect(c("rsp", "om", "apdepth", "clay", "ls", "slope"), names(map_sp@data))   # union of all candidate covariates (neutral common yardstick)
resid_maxSMD <- function(w, vars) if (!length(vars)) NA_real_ else
  round(max(abs(vapply(vars, function(v) .smd(map_sp@data[[v]], map_sp@data$Treat, w), numeric(1))), na.rm = TRUE), 3)
engine_row <- function(S, engine) {
  if (engine == "PSM") { w <- S$psm$m$weights; a <- unname(S$psm$ate["ATE"])
    data.frame(set = S$name, engine = "PSM", ATE = round(a, 2), shift = round(a - naive_set, 2),
      bal_all = resid_maxSMD(w, BAL_VARS), bal_conf = resid_maxSMD(w, CONFOUNDERS),
      ATE_sd = NA_real_, retained = 100, ESS = round(unname(S$psm$diag["ess_control"]), 1), row.names = NULL)
  } else { w <- S$mid$weights; a <- mean(S$cem$ATE, na.rm = TRUE)
    data.frame(set = S$name, engine = "CEM", ATE = round(a, 2), shift = round(a - naive_set, 2),
      bal_all = resid_maxSMD(w, BAL_VARS), bal_conf = resid_maxSMD(w, CONFOUNDERS),
      ATE_sd = round(sd(S$cem$ATE, na.rm = TRUE), 2),
      retained = round(mean(S$cem$pct_retained, na.rm = TRUE), 1),
      ESS = round(unname(match_diag(S$mid)["ess_control"]), 1), row.names = NULL)
  }
}
ord <- order(vapply(sets, function(S) resid_maxSMD(S$psm$m$weights, CONFOUNDERS), numeric(1)))   # best PSM confounder balance first
compare_sets <- do.call(rbind, lapply(sets[ord], function(S) rbind(engine_row(S, "PSM"), engine_row(S, "CEM"))))
compare_sets$set[c(FALSE, TRUE)] <- ""   # show each model name once (its PSM row), blank the CEM row -> visual 2-row blocks
cat(sprintf("\n=== (1) Competing models x engine (naive=%.1f). PSM = full matching; CEM = mean over %s-bin sweep.
  bal_all  = balance on ALL covariates [%s] after that engine (fair across models).
  bal_conf = balance on auto-detected confounders [%s] -> the bias that matters.
  ATE_sd is CEM-only; retained is 100 for PSM (keeps all units). Ordered by PSM bal_conf (least bias first). ===\n",
  naive_set, paste(range(BIN_COUNTS), collapse = "-"),
  paste(BAL_VARS, collapse = ", "),
  if (length(CONFOUNDERS)) paste(CONFOUNDERS, collapse = ", ") else "none"))
print(compare_sets); gridExtra::grid.table(compare_sets)


## =============================================================================
## (2) DESIGN SENSITIVITY: ATE vs bins, per confounder set
## =============================================================================
lev     <- names(CONF_SETS)
cem_all <- do.call(rbind, lapply(sets, `[[`, "cem")); cem_all$set <- factor(cem_all$set, levels = lev)
psm_ref <- do.call(rbind, lapply(sets, function(S) data.frame(set = factor(S$name, levels = lev),
  ATE = unname(S$psm$ate["ATE"]), lo = unname(S$psm$ate["lo"]), hi = unname(S$psm$ate["hi"]), row.names = NULL)))
print(
  ggplot(cem_all, aes(bins, ATE)) +
    geom_rect(data = psm_ref, aes(xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi), inherit.aes = FALSE, alpha = 0.10) +
    geom_hline(data = psm_ref, aes(yintercept = ATE), linetype = "dashed", colour = "grey30") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, fill = "#1b9e77") +
    geom_line(colour = "#1b9e77") + geom_point(colour = "#1b9e77", size = 2) +
    facet_wrap(~set) + scale_x_continuous(breaks = BIN_COUNTS) +
    labs(x = "CEM bins per variable", y = "ATE with 95% CI",
         title = "(2) Design sensitivity: ATE vs bins, per confounder set",
         subtitle = "Green = CEM sweep; dashed line/grey band = that set's PSM full-matching ATE. Flatter = steadier.") +
    theme_minimal())
print(
  ggplot(cem_all, aes(bins, L1, colour = set)) +
    geom_line() + geom_point(size = 2) + scale_x_continuous(breaks = BIN_COUNTS) +
    labs(x = "CEM bins per variable", y = "multivariate L1 imbalance", colour = "set",
         title = "(2) Balance (L1) vs bins, per confounder set  (lower = better)") +
    theme_minimal())

## =============================================================================
## (3) ADJACENCY TRAP, per confounder set
## =============================================================================
# Spatial-neighbour pairing pairs each treated plot with its NEAREST control across
# the AB line, ignoring covariates. Matching restricts the comparison to the treated
# plot's own stratum and then pairs it with the GEOGRAPHICALLY CLOSEST control IN
# that stratum (not the stratum average) -- so the only thing that changes from the
# raw pairing is the like-for-like (same-subclass) constraint. CEM and PSM strata
# are both shown, per model.
ad_tr <- which(map_sf$Treat == 1); ad_ct <- which(map_sf$Treat == 0)
ni <- sf::st_nearest_feature(sf::st_centroid(map_sf[ad_tr, ]), sf::st_centroid(map_sf[ad_ct, ]))
adj_delta <- map_sf$Yield[ad_tr] - map_sf$Yield[ad_ct][ni]
# treated - the geographically CLOSEST retained control within its stratum of match `m`
strat_delta <- function(m) { sc <- m$subclass; w <- m$weights
  cx <- map_sf$centroid_x; cy <- map_sf$centroid_y
  vapply(ad_tr, function(i) {
    cc <- ad_ct[which(sc[ad_ct] == sc[i] & !is.na(sc[ad_ct]) & w[ad_ct] > 0)]
    if (!length(cc)) return(NA_real_)
    d <- (cx[cc] - cx[i])^2 + (cy[cc] - cy[i])^2          # squared distance, same stratum
    map_sf$Yield[i] - map_sf$Yield[cc[which.min(d)]]
  }, numeric(1)) }
dd <- rbind(
  data.frame(set = "(raw neighbour)", method = "Spatial", delta = adj_delta),
  do.call(rbind, lapply(sets, function(S) data.frame(set = S$name, method = "CEM", delta = strat_delta(S$mid)))),
  do.call(rbind, lapply(sets, function(S) data.frame(set = S$name, method = "PSM", delta = strat_delta(S$psm$m)))))
dd <- dd[is.finite(dd$delta), ]
dd$set    <- factor(dd$set, levels = c("(raw neighbour)", names(CONF_SETS)))
dd$method <- factor(dd$method, levels = c("Spatial", "CEM", "PSM"))
delta_summ <- dplyr::summarise(dplyr::group_by(dd, set, method),
  mean_delta = round(mean(delta), 1), n_negative = sum(delta < 0), n = dplyr::n(), .groups = "drop")
cat("\n=== (3) Across-strip delta by pairing (mean, # negative) ===\n"); print(as.data.frame(delta_summ))
print(
  ggplot(dd, aes(set, delta, fill = method)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_violin(position = position_dodge2(preserve = "single"), alpha = 0.4, colour = NA, scale = "width") +
    geom_boxplot(position = position_dodge2(preserve = "single", padding = 0.2), width = 0.6, outlier.size = 0.4) +
    scale_fill_manual(values = c(Spatial = "grey70", CEM = "#1b9e77", PSM = "#d95f02"), name = NULL) +
    labs(x = NULL, y = "treated - control delta (bu/ac)",
         title = "(3) Adjacency trap: raw vs CEM vs PSM pairing, per model",
         subtitle = "Each treated plot vs the geographically closest control IN its own stratum.\nRaw neighbour = nearest control across the strip, ignoring covariates.") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 20, hjust = 1)))


## =============================================================================
## (4) SENSITIVITY TO HIDDEN BIAS: E-value + Rosenbaum critical Gamma
##   Per model x engine: how strong an UNMEASURED confounder would have to be to
##   overturn the result. E-value (VanderWeele & Ding 2017, EValue::evalues.OLS) on
##   the risk-ratio scale; Rosenbaum critical Gamma (the bias odds ratio at which
##   the one-sided p first exceeds 0.05) via Huber M-statistics -- senfm for PSM
##   full matching, senstrat for CEM m:n strata. Both RISE with effect size, so read
##   them WITH the balance in (1), not instead of it. Also saved as <Site>_sensitivity.png.
## =============================================================================
sens_table <- tryCatch({
  for (p in c("EValue", "sensitivityfull", "senstrat")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages({ library(EValue); library(sensitivityfull); library(senstrat) })
  sd_y   <- sd(map_sp@data$Yield, na.rm = TRUE)
  ate_se <- function(md) { md$Tn <- as.integer(as.character(md$Treat)); f <- lm(Yield ~ Tn, md, weights = md$weights)
    cl <- if (!is.null(md$subclass)) md$subclass else seq_len(nrow(md)); V <- sandwich::vcovCL(f, cluster = cl)
    c(ATE = unname(coef(f)["Tn"]), SE = sqrt(V["Tn", "Tn"])) }
  ev_pt  <- function(est, se) round(unname(evalues.OLS(est = est, se = se, sd = sd_y, delta = 1, true = 0)["E-values", "point"]), 2)
  build_senfm <- function(mdf) { rows <- list(); t1 <- logical(0)
    for (sc in unique(mdf$subclass)) { s <- mdf[mdf$subclass == sc, ]; tr <- as.integer(as.character(s$Treat))
      yt <- s$Yield[tr == 1]; yc <- s$Yield[tr == 0]; if (!length(yt) || !length(yc)) next
      if (length(yt) == 1) { rows[[length(rows) + 1]] <- c(yt, yc); t1 <- c(t1, TRUE) }
      else if (length(yc) == 1) { rows[[length(rows) + 1]] <- c(yc, yt); t1 <- c(t1, FALSE) } else next }
    if (length(rows) < 2) return(NULL); J <- max(lengths(rows))
    list(y = t(sapply(rows, function(r) c(r, rep(NA, J - length(r))))), treated1 = t1) }
  crit <- function(pf) { if (pf(1) > 0.05) return("1"); lo <- 1; hi <- 2
    while (pf(hi) <= 0.05) { lo <- hi; hi <- hi * 2; if (hi > 1000) return(">1000") }
    while (hi - lo > 0.02) { m <- (lo + hi) / 2; if (pf(m) <= 0.05) lo <- m else hi <- m }; as.character(round((lo + hi) / 2, 2)) }
  g_psm <- function(mdf) { b <- build_senfm(mdf); if (is.null(b)) return(NA_character_)
    crit(function(g) senfm(b$y, b$treated1, gamma = g, alternative = "greater")$pval) }
  g_cem <- function(mdf) { z <- as.integer(as.character(mdf$Treat)); st <- as.integer(mdf$subclass); sc <- mscores(mdf$Yield, z, st)
    crit(function(g) as.numeric(senstrat(sc, z, st, gamma = g, alternative = "greater")$Result["P-value"])) }
  st <- do.call(rbind, lapply(sets[ord], function(S) {
    mdp <- MatchIt::match.data(S$psm$m); ap <- ate_se(mdp); mdc <- MatchIt::match.data(S$mid); ac <- ate_se(mdc)
    rbind(data.frame(Model = S$name, Engine = "PSM", ATE = round(unname(ap["ATE"]), 1),
                     `E-value` = ev_pt(ap["ATE"], ap["SE"]), `Crit. Gamma` = g_psm(mdp), row.names = NULL, check.names = FALSE),
          data.frame(Model = "",      Engine = "CEM", ATE = round(unname(ac["ATE"]), 1),
                     `E-value` = ev_pt(ac["ATE"], ac["SE"]), `Crit. Gamma` = g_cem(mdc), row.names = NULL, check.names = FALSE)) }))
  cat("\n=== (4) Sensitivity to hidden bias (E-value & Rosenbaum critical Gamma) ===\n")
  cat("E-value  = min confounder assoc. (RR scale) with BOTH treatment & yield needed to explain away the ATE.\n")
  cat("Crit.Gamma = hidden-bias odds ratio at which one-sided p>0.05 (senfm=PSM full match, senstrat=CEM strata).\n")
  cat("Higher = more robust. BOTH rise with effect size -> read with the balance table (1), not instead of it.\n")
  print(st)
  cols <- ifelse(st$Engine == "PSM", "#e8f3ec", "#fdeee4")
  tt <- gridExtra::ttheme_minimal(base_size = 12, core = list(bg_params = list(fill = cols, col = "grey85"),
    fg_params = list(hjust = 0.5, x = 0.5)), colhead = list(fg_params = list(fontface = "bold")))
  gt <- gridExtra::tableGrob(st, rows = NULL, theme = tt); grid::grid.newpage(); grid::grid.draw(gt)
  ggplot2::ggsave(sprintf("%s_sensitivity.png", SITE), gt, width = 8, height = 4.6, dpi = 130)
  st
}, error = function(e) { message("(4) sensitivity skipped: ", conditionMessage(e)); NULL })

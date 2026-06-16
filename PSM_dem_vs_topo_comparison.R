## =============================================================================
## DEM vs. Topo-Derivatives Covariate Comparison
## Goal: test whether matching on ec+dem or dem-only gives a significantly
##       different treatment effect than matching on ec+ls+rsp.
## Requires: map_sp, w_Adj, wAdj, yield objects from the main script.
## Models compared:
##   Model 1 (topo): ec + ls + rsp   (main script)
##   Model 2 (dem):  ec + dem
##   Model 3 (dem-only): dem
## =============================================================================
##
## PIPELINE OVERVIEW (read once; the per-model blocks below all repeat it):
##   1. Build spatial weights matrices from polygon adjacency (spdep).
##   2. Mean-center the covariates (ec, ls, rsp, dem).
##   3. For each candidate covariate set:
##        a. Fit a spatial autoregressive (SAR) probit of Treat on the
##           covariates  -> sarprobit() {spatialprobit}.  Its only job is to
##           produce propensity scores; its coefficients are NOT interpreted.
##        b. Convert fitted values to propensity scores via the logistic CDF.
##        c. Full matching on those scores -> matchit(method="full") {MatchIt}.
##        d. Estimate the treatment effect on the matched sample with a
##           doubly-robust g-computation -> avg_comparisons() {marginaleffects},
##           with cluster-robust SEs by matched subclass and matching weights.
##   4. Summary tables A-F rendered to the RStudio Plots panel.
##
## PACKAGES (load these before running): sf, sp, spdep, spatialprobit,
##   MatchIt, cobalt, marginaleffects, gridExtra, grid, EValue, sensitivityfull,
##   senstrat.
##   Upstream objects required from the main script: map_sf, yield, Treat/Yield
##   columns in the data.
##
## METHOD REFERENCES (verify against the originals before citing):
##   - Full matching:        Hansen 2004; MatchIt vignette (Ho et al.).
##   - Doubly-robust effect:  any DR / g-computation reference; here it is an
##                            outcome model + averaging via marginaleffects.
##   - E-value:               VanderWeele & Ding 2017, Ann Intern Med.
##   - Rosenbaum Gamma:       Rosenbaum 2007 (Biometrics 63:456-464);
##                            full-matching models -> sensitivityfull::senfm;
##                            CEM (m:n strata)      -> senstrat::senstrat.
## =============================================================================
# Convert the SpatialPolygonsDataFrame to a Spatial object (needed by spdep/sp).
map_sp <- as(map_sf, "Spatial")


## Extract the coordinates from the spatial object 'map_sp' and assigns them to the variable 'coords'
coords <- coordinates(map_sp)  
#Rook contiguity neighbors for effect estimation
neighbors<- poly2nb(map_sp, queen=FALSE) # find all adjacent neighbors based on rook contiguity
wAdj <- nb2listw(neighbors, style="W")  # create a row standardized weights matrix based on adjacency
#Queen contiguity neighbors for SAR probit calculation of propensity scores

neighbors_SAR <- poly2nb(map_sp, queen=TRUE)
wAdj_SAR <- nb2listw( neighbors_SAR, style="W")  # use queen contiguity for SAR due to lower Morans I for fit1 residuals
# These W matrices can be modified for use in the SAR Probit according to https://journal.r-project.org/archive/2013/RJ-2013-013/RJ-2013-013.pdf
w_Adj <- as(as_dgRMatrix_listw(wAdj_SAR), "CsparseMatrix")

# Mean-center each covariate (center = TRUE, scale = FALSE => subtract the mean,
# do NOT divide by SD). Centering only; the covariates keep their original units.
map_sp@data$ls <- scale(map_sp@data$ls, center = TRUE, scale = FALSE)
map_sp@data$rsp <- scale(map_sp@data$rsp, center = TRUE, scale = FALSE)
map_sp@data$ec <- scale(map_sp@data$ec, center = TRUE, scale = FALSE)
map_sp@data$dem <- scale(map_sp@data$dem, center = TRUE, scale = FALSE)

# Snapshot the prepared analysis data to disk for inspection.
write.csv(map_sp@data, "map_sp_data.csv")

# ---- MODEL 1 (focal "topo" set): propensity model = Treat ~ rsp + ls + ec ----
# SAR probit treatment-assignment model. Output used ONLY to derive propensity
# scores; rho captures spatial autocorrelation in assignment.
model <- Treat ~ rsp + ls + ec
#mod1 <- probit(formula = model, data = map_sp@data, method = "ML")  # non-spatial alt (unused)
mod_SAR <- sarprobit(model, w_Adj, map_sp@data)   # w_Adj = queen-contiguity weights
summary(mod_SAR)     # coefficients (not interpreted) + rho
impacts(mod_SAR)     # direct/indirect/total marginal effects (diagnostic only)
logLik(mod_SAR)      # \
AIC(mod_SAR)         #  > fit statistics, collected later in TABLE D
BIC(mod_SAR)         # /


#### PROPENSITY SCORE CALCULATION
# fitted() returns the linear predictor; the logistic CDF 1/(1+e^-x) maps it to
# a 0-1 propensity score. (Probit would use pnorm(); logistic is used here.)
predicted_values <- fitted(mod_SAR)
p_score <- 1 / (1 + exp(-predicted_values))


# MATCHIT: full matching on the propensity score. Every unit is retained and
# placed in a matched subclass (variable treated:control ratios); "full" gave
# the best balance among the methods tried.
matchit_test <- matchit(model, data = map_sp@data, method = "full", pscores = p_score)
#matchit_test <- matchit(model, data = map_sp@data, exact = ~SoilTyp, method = "nearest", distance = "logit", pscores = p_score) #by soil type
summary(matchit_test)

bal.tab(matchit_test, un = TRUE)  # balance table (un=TRUE shows pre-match SMDs too)

# Diagnostics for this model (also produced for each model below):
plot(matchit_test, type = "jitter", interactive = FALSE)  # PS overlap
plot(summary(matchit_test))                               # love plot of SMDs
plot(matchit_test, type = "qq")                           # covariate QQ by group
for (v in c("rsp", "ls", "ec")) {                         # per-covariate eCDF overlap
  print(bal.plot(matchit_test, var.name = v, which = "both",
                 type = "ecdf", mirror = FALSE) +
          ggplot2::ggtitle(paste0("eCDF: ", v, " — ec+ls+rsp matching")))
}
love.plot(matchit_test, thresholds = 0.1)  # SMDs vs the 0.1 balance threshold


# Extract the matched sample (adds 'weights' and 'subclass' columns).
matched_data <- match.data(matchit_test)
# Coerce types so the outcome model and contrasts behave (Treat as factor;
# covariates back to plain numeric after scale() returned a matrix column).
matched_data$Treat <- factor(matched_data$Treat)
matched_data$rsp <- as.numeric(matched_data$rsp)
matched_data$ls <- as.numeric(matched_data$ls)
matched_data$ec <- as.numeric(matched_data$ec)


##### EFFECT ESTIMATION (doubly robust)
# Outcome model includes Treat * covariates, fit with the matching weights. The
# treatment effect then comes from avg_comparisons() below. "Doubly robust":
# correct if EITHER the matching OR this outcome model is right.
#fit1 <- lm(Yield ~ Treat * SoilType, data = matched_data, weights = weights) # for by soil type
fit1 <- lm(Yield ~ Treat * (rsp+ls+ec), data = matched_data, weights = weights)
# avg_comparisons: average marginal effect of Treat (ATE = E[Y|1] - E[Y|0]),
#   vcov = ~subclass -> SEs clustered by matched subclass,
#   wts  = "weights" -> use the full-matching weights in the averaging.
treat_effect_avgcomp <- avg_comparisons(fit1, variables = "Treat",
                                        vcov = ~subclass, wts = "weights")
#newdata = subset(matched_data, Treat == 1),# restrict to treated for ATT instead of ATE
treat_effect_avgcomp
# Summary of the lm model
summary(fit1)

# Export the focal matched sample as a spatial points layer (centroids) for
# mapping/QGIS. Not used in the effect estimates below.
md_df <- as.data.frame(matched_data)          # drop the 'matchdata' class
coords_mat <- as.matrix(md_df[, c("centroid_x", "centroid_y")])
write.csv(md_df, file = "RAinbarrel-matched_data09Oct2025.csv", row.names = FALSE)
matched_data_sp <- sp::SpatialPointsDataFrame(
  coords      = coords_mat,
  data        = md_df,
  proj4string = sp::CRS(sf::st_crs(yield)$wkt)   # reuse the yield layer's CRS
)

## =============================================================================
## MODELS 2-8: each block below repeats the SAME pattern as Model 1 above, just
## with a different covariate set:
##     model_*  <- Treat ~ <covariates>
##     mod_SAR_* <- sarprobit(...)             # propensity model
##     p_score_* <- 1/(1+exp(-fitted(...)))    # logistic CDF -> propensity score
##     matchit_* <- matchit(..., method="full")# full matching
##     matched_* <- match.data(...)            # matched sample (+weights,+subclass)
##     fit_*     <- lm(Yield ~ Treat*cov, weights=weights)   # outcome model
##     te_*      <- avg_comparisons(fit_*, "Treat", vcov=~subclass, wts="weights")
##     est_* / se_*  <- point estimate and SE, collected for the tables.
## Only the covariate set changes; read Model 1 to understand them all.
## =============================================================================

# ---- 1. SAR Probit for ec + dem model ------------------------------------

model_dem <- Treat ~ ec + dem

mod_SAR_dem <- sarprobit(model_dem, w_Adj, map_sp@data)
cat("\n=== SAR Probit: ec + dem ===\n")
summary(mod_SAR_dem)
cat("LogLik:", logLik(mod_SAR_dem),
    "  AIC:", AIC(mod_SAR_dem),
    "  BIC:", BIC(mod_SAR_dem), "\n")

# Propensity scores
p_score_dem <- 1 / (1 + exp(-fitted(mod_SAR_dem)))

# ---- 2. MatchIt with ec + dem -------------------------------------------

matchit_dem <- matchit(model_dem, data = map_sp@data,
                       method = "full", pscores = p_score_dem)
cat("\n=== Balance: ec + dem matching ===\n")
print(summary(matchit_dem))
plot(matchit_dem, type = "jitter", interactive = FALSE)
plot(summary(matchit_dem))
plot(matchit_dem, type = "qq")
for (v in c("ec", "dem")) {
  print(bal.plot(matchit_dem, var.name = v, which = "both",
                 type = "ecdf", mirror = FALSE) +
          ggplot2::ggtitle(paste0("eCDF: ", v, " — ec+dem matching")))
}
bal.tab(matchit_dem, un = TRUE)
love.plot(matchit_dem, thresholds = 0.1,
          title = "Balance: ec + dem")

# ---- 3. Effect estimation: ec + dem -------------------------------------

matched_dem <- match.data(matchit_dem)
matched_dem$Treat <- factor(matched_dem$Treat)
matched_dem$ec    <- as.numeric(matched_dem$ec)
matched_dem$dem   <- as.numeric(matched_dem$dem)

fit_dem <- lm(Yield ~ Treat * (ec + dem), data = matched_dem,
              weights = weights)

te_dem <- avg_comparisons(
  fit_dem,
  variables = "Treat",
  vcov      = ~subclass,
  wts       = "weights"
)
cat("\n=== Treatment effect (ec + dem matching) ===\n")
print(te_dem)

# For reference, re-state the original ec+ls+rsp effect
# (treat_effect_avgcomp must exist from the main script)
cat("\n=== Treatment effect (ec + ls + rsp matching) ===\n")
print(treat_effect_avgcomp)

# ---- 4. SAR Probit for dem-only model ------------------------------------

model_dem_only <- Treat ~ dem

mod_SAR_dem_only <- sarprobit(model_dem_only, w_Adj, map_sp@data)
cat("\n=== SAR Probit: dem only ===\n")
summary(mod_SAR_dem_only)
cat("LogLik:", logLik(mod_SAR_dem_only),
    "  AIC:", AIC(mod_SAR_dem_only),
    "  BIC:", BIC(mod_SAR_dem_only), "\n")

# Propensity scores
p_score_dem_only <- 1 / (1 + exp(-fitted(mod_SAR_dem_only)))

# ---- 5. MatchIt with dem only --------------------------------------------

matchit_dem_only <- matchit(model_dem_only, data = map_sp@data,
                            method = "full", pscores = p_score_dem_only)
cat("\n=== Balance: dem-only matching ===\n")
print(summary(matchit_dem_only))
plot(matchit_dem_only, type = "jitter", interactive = FALSE)
plot(summary(matchit_dem_only))
plot(matchit_dem_only, type = "qq")
print(bal.plot(matchit_dem_only, var.name = "dem", which = "both",
               type = "ecdf", mirror = FALSE) +
        ggplot2::ggtitle("eCDF: dem — dem-only matching"))
bal.tab(matchit_dem_only, un = TRUE)
love.plot(matchit_dem_only, thresholds = 0.1,
          title = "Balance: dem only")

# ---- 6. Effect estimation: dem only -------------------------------------

matched_dem_only <- match.data(matchit_dem_only)
matched_dem_only$Treat <- factor(matched_dem_only$Treat)
matched_dem_only$dem   <- as.numeric(matched_dem_only$dem)

fit_dem_only <- lm(Yield ~ Treat * dem, data = matched_dem_only,
                   weights = weights)

te_dem_only <- avg_comparisons(
  fit_dem_only,
  variables = "Treat",
  vcov      = ~subclass,
  wts       = "weights"
)
cat("\n=== Treatment effect (dem-only matching) ===\n")
print(te_dem_only)

# ---- 7. Unmatched baseline -----------------------------------------------
# Run on full unmatched data to show what matching is actually adding

fit_unmatched <- lm(Yield ~ Treat, data = map_sp@data)
te_unmatched <- avg_comparisons(fit_unmatched, variables = "Treat")
cat("\n=== Treatment effect (no matching) ===\n")
print(te_unmatched)

est_unmatched <- te_unmatched$estimate
se_unmatched  <- te_unmatched$std.error

# ---- 8. SAR Probit for ec-only model -------------------------------------

model_ec_only <- Treat ~ ec

mod_SAR_ec_only <- sarprobit(model_ec_only, w_Adj, map_sp@data)
cat("\n=== SAR Probit: ec only ===\n")
summary(mod_SAR_ec_only)
cat("LogLik:", logLik(mod_SAR_ec_only),
    "  AIC:", AIC(mod_SAR_ec_only),
    "  BIC:", BIC(mod_SAR_ec_only), "\n")

p_score_ec_only <- 1 / (1 + exp(-fitted(mod_SAR_ec_only)))

# ---- 9. MatchIt with ec only ---------------------------------------------

matchit_ec_only <- matchit(model_ec_only, data = map_sp@data,
                           method = "full", pscores = p_score_ec_only)
cat("\n=== Balance: ec-only matching ===\n")
print(summary(matchit_ec_only))
plot(matchit_ec_only, type = "jitter", interactive = FALSE)
plot(summary(matchit_ec_only))
plot(matchit_ec_only, type = "qq")
print(bal.plot(matchit_ec_only, var.name = "ec", which = "both",
               type = "ecdf", mirror = FALSE) +
        ggplot2::ggtitle("eCDF: ec — ec-only matching"))
bal.tab(matchit_ec_only, un = TRUE)
love.plot(matchit_ec_only, thresholds = 0.1,
          title = "Balance: ec only")

# ---- 10. Effect estimation: ec only --------------------------------------

matched_ec_only <- match.data(matchit_ec_only)
matched_ec_only$Treat <- factor(matched_ec_only$Treat)
matched_ec_only$ec    <- as.numeric(matched_ec_only$ec)

fit_ec_only <- lm(Yield ~ Treat * ec, data = matched_ec_only,
                  weights = weights)

te_ec_only <- avg_comparisons(
  fit_ec_only,
  variables = "Treat",
  vcov      = ~subclass,
  wts       = "weights"
)
cat("\n=== Treatment effect (ec-only matching) ===\n")
print(te_ec_only)

est_ec_only <- te_ec_only$estimate
se_ec_only  <- te_ec_only$std.error

# ---- 11. SAR Probit + matching: rsp only ------------------------------------

model_rsp_only <- Treat ~ rsp

mod_SAR_rsp_only <- sarprobit(model_rsp_only, w_Adj, map_sp@data)
cat("\n=== SAR Probit: rsp only ===\n")
summary(mod_SAR_rsp_only)
cat("LogLik:", logLik(mod_SAR_rsp_only),
    "  AIC:", AIC(mod_SAR_rsp_only),
    "  BIC:", BIC(mod_SAR_rsp_only), "\n")

p_score_rsp_only <- 1 / (1 + exp(-fitted(mod_SAR_rsp_only)))

matchit_rsp_only <- matchit(model_rsp_only, data = map_sp@data,
                             method = "full", pscores = p_score_rsp_only)
cat("\n=== Balance: rsp-only matching ===\n")
print(summary(matchit_rsp_only))
plot(matchit_rsp_only, type = "jitter", interactive = FALSE)
plot(summary(matchit_rsp_only))
plot(matchit_rsp_only, type = "qq")
print(bal.plot(matchit_rsp_only, var.name = "rsp", which = "both",
               type = "ecdf", mirror = FALSE) +
        ggplot2::ggtitle("eCDF: rsp — rsp-only matching"))
bal.tab(matchit_rsp_only, un = TRUE)
love.plot(matchit_rsp_only, thresholds = 0.1, title = "Balance: rsp only")

matched_rsp_only <- match.data(matchit_rsp_only)
matched_rsp_only$Treat <- factor(matched_rsp_only$Treat)
matched_rsp_only$rsp   <- as.numeric(matched_rsp_only$rsp)

fit_rsp_only <- lm(Yield ~ Treat * rsp, data = matched_rsp_only, weights = weights)

te_rsp_only <- avg_comparisons(
  fit_rsp_only, variables = "Treat",
  vcov = ~subclass, wts = "weights"
)
cat("\n=== Treatment effect (rsp-only matching) ===\n")
print(te_rsp_only)

est_rsp_only <- te_rsp_only$estimate
se_rsp_only  <- te_rsp_only$std.error

# ---- 12. SAR Probit + matching: ls only -------------------------------------

model_ls_only <- Treat ~ ls

mod_SAR_ls_only <- sarprobit(model_ls_only, w_Adj, map_sp@data)
cat("\n=== SAR Probit: ls only ===\n")
summary(mod_SAR_ls_only)
cat("LogLik:", logLik(mod_SAR_ls_only),
    "  AIC:", AIC(mod_SAR_ls_only),
    "  BIC:", BIC(mod_SAR_ls_only), "\n")

p_score_ls_only <- 1 / (1 + exp(-fitted(mod_SAR_ls_only)))

matchit_ls_only <- matchit(model_ls_only, data = map_sp@data,
                            method = "full", pscores = p_score_ls_only)
cat("\n=== Balance: ls-only matching ===\n")
print(summary(matchit_ls_only))
plot(matchit_ls_only, type = "jitter", interactive = FALSE)
plot(summary(matchit_ls_only))
plot(matchit_ls_only, type = "qq")
print(bal.plot(matchit_ls_only, var.name = "ls", which = "both",
               type = "ecdf", mirror = FALSE) +
        ggplot2::ggtitle("eCDF: ls — ls-only matching"))
bal.tab(matchit_ls_only, un = TRUE)
love.plot(matchit_ls_only, thresholds = 0.1, title = "Balance: ls only")

matched_ls_only <- match.data(matchit_ls_only)
matched_ls_only$Treat <- factor(matched_ls_only$Treat)
matched_ls_only$ls    <- as.numeric(matched_ls_only$ls)

fit_ls_only <- lm(Yield ~ Treat * ls, data = matched_ls_only, weights = weights)

te_ls_only <- avg_comparisons(
  fit_ls_only, variables = "Treat",
  vcov = ~subclass, wts = "weights"
)
cat("\n=== Treatment effect (ls-only matching) ===\n")
print(te_ls_only)

est_ls_only <- te_ls_only$estimate
se_ls_only  <- te_ls_only$std.error

# ---- 13. CEM with DEM quartiles ------------------------------------------
# CEM bins dem into quartile strata and matches exactly within each bin.
# Balance is guaranteed by construction; no propensity score needed.

dem_quartile_cuts <- quantile(map_sp@data$dem, probs = c(0, 0.25, 0.5, 0.75, 1))

matchit_cem <- matchit(
  Treat ~ dem,
  data      = map_sp@data,
  method    = "cem",
  cutpoints = list(dem = dem_quartile_cuts)
)
cat("\n=== Balance: CEM dem (quartiles) ===\n")
print(summary(matchit_cem))
plot(summary(matchit_cem))  # jitter and qq omitted — no propensity score in CEM
print(bal.plot(matchit_cem, var.name = "dem", which = "both",
               type = "ecdf", mirror = FALSE) +
        ggplot2::ggtitle("eCDF: dem — CEM (quartiles)"))
bal.tab(matchit_cem, un = TRUE)
love.plot(matchit_cem, thresholds = 0.1,
          title = "Balance: CEM dem (quartiles)")

# ---- 12. Effect estimation: CEM dem --------------------------------------

matched_cem <- match.data(matchit_cem)
matched_cem$Treat <- factor(matched_cem$Treat)
matched_cem$dem   <- as.numeric(matched_cem$dem)

fit_cem <- lm(Yield ~ Treat * dem, data = matched_cem, weights = weights)

te_cem <- avg_comparisons(
  fit_cem,
  variables = "Treat",
  vcov      = ~subclass,
  wts       = "weights"
)
cat("\n=== Treatment effect (CEM dem quartiles) ===\n")
print(te_cem)

est_cem <- te_cem$estimate
se_cem  <- te_cem$std.error

# ---- 13–17. Summary tables (all displayed in RStudio Plots panel) --------

if (!requireNamespace("gridExtra",       quietly = TRUE)) install.packages("gridExtra")
if (!requireNamespace("EValue",          quietly = TRUE)) install.packages("EValue")
if (!requireNamespace("sensitivityfull", quietly = TRUE)) install.packages("sensitivityfull")
library(gridExtra)
library(grid)
library(EValue)
library(sensitivityfull)

# Helper: draw a data frame as a table in the Plots panel
# note: pass NULL to omit the footnote (use for tables with no significance stars).
# Default is the significance legend, appropriate only for tables with a Sig column.
plot_table <- function(df, title,
                       note = "*** p<0.001  ** p<0.01  * p<0.05  ns = not significant") {
  tt <- ttheme_default(
    core    = list(fg_params = list(cex = 0.8)),
    colhead = list(fg_params = list(cex = 0.85, fontface = "bold"))
  )
  tbl <- tableGrob(df, rows = NULL, theme = tt)
  title_grob <- textGrob(title, gp = gpar(fontsize = 11, fontface = "bold"))
  grid.newpage()
  if (is.null(note) || !nzchar(note)) {
    grid.draw(arrangeGrob(title_grob, tbl, ncol = 1, heights = c(0.10, 0.90)))
  } else {
    note_grob <- textGrob(note, gp = gpar(fontsize = 8, col = "grey40"))
    grid.draw(arrangeGrob(title_grob, tbl, note_grob,
                          ncol = 1, heights = c(0.08, 0.84, 0.08)))
  }
}

# --- Shared objects ---
# Two-sample z-test for the difference between two INDEPENDENT estimates.
# These ATEs come from separate matched datasets (no shared error term), so a
# z-test on (est1-est2)/sqrt(se1^2+se2^2) is used rather than ANOVA. Two-sided p.
z_test <- function(est1, se1, est2, se2) {
  z <- (est1 - est2) / sqrt(se1^2 + se2^2)
  p <- 2 * pnorm(-abs(z))
  c(z = z, p = p)
}
# Significance-star helper for the table p-value columns.
sig_star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                        ifelse(p < 0.05,  "*",   "ns")))

# Pull point estimate + SE from each model's avg_comparisons result.
# (est_unmatched/se_unmatched, est_ec_only/..., est_*_only, est_cem were set in
#  their blocks above; the three below were not, so set them here.)
est_topo     <- treat_effect_avgcomp$estimate; se_topo     <- treat_effect_avgcomp$std.error
est_dem      <- te_dem$estimate;               se_dem      <- te_dem$std.error
est_dem_only <- te_dem_only$estimate;          se_dem_only <- te_dem_only$std.error

# Stack all eight models into parallel vectors used by every table below.
# Order is fixed and shared: keep all_est, all_se, all_mod aligned.
all_est <- c(est_unmatched, est_topo, est_ec_only, est_dem,
             est_dem_only, est_rsp_only, est_ls_only, est_cem)
all_se  <- c(se_unmatched,  se_topo,  se_ec_only,  se_dem,
             se_dem_only,  se_rsp_only,  se_ls_only,  se_cem)
all_mod <- c("No matching", "PSM: ec+ls+rsp", "PSM: ec only", "PSM: ec+dem",
             "PSM: dem only", "PSM: rsp only", "PSM: ls only",
             "CEM: dem (quartiles)")

# --- TABLE A: ATE per model ---
ate_table <- data.frame(
  Model   = all_mod,
  ATE     = round(all_est, 3),
  SE      = round(all_se,  3),
  CI_low  = round(all_est - 1.96 * all_se, 3),
  CI_high = round(all_est + 1.96 * all_se, 3)
)
plot_table(ate_table, "TABLE A: Treatment Effect Estimates by Model", note = NULL)

# --- TABLE B: Each model vs no-matching baseline (Bonferroni corrected) ---
# z-test every matched model (all_est[-1], i.e. excluding "No matching") against
# the unmatched baseline. Tests whether matching changed the estimate.
zp_list <- lapply(seq_along(all_est[-1]), function(i)
  z_test(all_est[i + 1], all_se[i + 1], est_unmatched, se_unmatched))

raw_p   <- sapply(zp_list, `[`, "p")          # raw two-sided p-values
bonf_p  <- pmin(raw_p * length(raw_p), 1)     # Bonferroni: p * (#comparisons), capped at 1

baseline_table <- data.frame(
  Model        = all_mod[-1],
  Difference   = round(all_est[-1] - est_unmatched, 3),
  z            = round(sapply(zp_list, `[`, "z"), 3),
  p_raw        = signif(raw_p,  3),
  p_Bonferroni = signif(bonf_p, 3),
  Sig          = sig_star(bonf_p)
)
plot_table(baseline_table,
           "TABLE B: Each Matched Model vs No-Matching Baseline\n(Bonferroni corrected, 5 comparisons)")

# --- TABLE C: PSM models vs CEM dem ---
# Same z-test machinery, but the reference is the CEM (coarsened exact matching)
# estimate instead of the unmatched baseline: does any PSM model differ from CEM?
psm_models <- c("PSM: ec+ls+rsp", "PSM: ec only", "PSM: ec+dem",
                "PSM: dem only", "PSM: rsp only", "PSM: ls only")
psm_est    <- c(est_topo, est_ec_only, est_dem, est_dem_only, est_rsp_only, est_ls_only)
psm_se     <- c(se_topo,  se_ec_only,  se_dem,  se_dem_only,  se_rsp_only,  se_ls_only)

zp_cem <- lapply(seq_along(psm_est), function(i)
  z_test(psm_est[i], psm_se[i], est_cem, se_cem))

raw_p_cem  <- sapply(zp_cem, `[`, "p")
bonf_p_cem <- pmin(raw_p_cem * length(raw_p_cem), 1)   # Bonferroni over these comparisons

cem_table <- data.frame(
  Comparison   = paste(psm_models, "vs CEM"),
  Difference   = round(psm_est - est_cem, 3),
  z            = round(sapply(zp_cem, `[`, "z"), 3),
  p_raw        = signif(raw_p_cem,  3),
  p_Bonferroni = signif(bonf_p_cem, 3),
  Sig          = sig_star(bonf_p_cem)
)
plot_table(cem_table,
           "TABLE C: PSM Models vs CEM dem (quartiles)\n(Bonferroni corrected, 4 comparisons)")

# --- TABLE D: SAR Probit model fit ---
psm_fit_table <- data.frame(
  Model  = c("ec+ls+rsp", "ec only", "ec+dem", "dem only", "rsp only", "ls only"),
  AIC    = round(c(AIC(mod_SAR), AIC(mod_SAR_ec_only), AIC(mod_SAR_dem),
                   AIC(mod_SAR_dem_only), AIC(mod_SAR_rsp_only), AIC(mod_SAR_ls_only)), 2),
  BIC    = round(c(BIC(mod_SAR), BIC(mod_SAR_ec_only), BIC(mod_SAR_dem),
                   BIC(mod_SAR_dem_only), BIC(mod_SAR_rsp_only), BIC(mod_SAR_ls_only)), 2),
  LogLik = round(c(logLik(mod_SAR), logLik(mod_SAR_ec_only), logLik(mod_SAR_dem),
                   logLik(mod_SAR_dem_only), logLik(mod_SAR_rsp_only), logLik(mod_SAR_ls_only)), 2)
)
plot_table(psm_fit_table, "TABLE D: SAR Probit Propensity Score Model Fit", note = NULL)

# --- TABLE E: SMD balance comparison ---
bal_topo     <- bal.tab(matchit_test,     un = TRUE)
bal_dem      <- bal.tab(matchit_dem,      un = TRUE)
bal_dem_only <- bal.tab(matchit_dem_only, un = TRUE)
bal_ec_only  <- bal.tab(matchit_ec_only,  un = TRUE)
bal_rsp_only <- bal.tab(matchit_rsp_only, un = TRUE)
bal_ls_only  <- bal.tab(matchit_ls_only,  un = TRUE)
bal_cem      <- bal.tab(matchit_cem,      un = TRUE)

# Pull BOTH unmatched (before) and adjusted (after) SMDs per model.
# Each model's own balance table supplies the before-value for its covariates,
# so dem gets a real before-value (was NA when pulled only from the topo table).
pull_ba <- function(bal, covs) {
  idx <- match(covs, rownames(bal$Balance))
  data.frame(B = round(bal$Balance$Diff.Un[idx],  3),
             A = round(bal$Balance$Diff.Adj[idx], 3))
}
all_covs <- unique(c(rownames(bal_topo$Balance), rownames(bal_dem$Balance),
                     rownames(bal_dem_only$Balance), rownames(bal_ec_only$Balance),
                     rownames(bal_rsp_only$Balance), rownames(bal_ls_only$Balance),
                     rownames(bal_cem$Balance)))

ba_topo <- pull_ba(bal_topo,     all_covs)   # ec+ls+rsp  (focal)
ba_demO <- pull_ba(bal_dem_only, all_covs)   # dem only   (focal)
ba_ecO  <- pull_ba(bal_ec_only,  all_covs)   # ec only
ba_ecd  <- pull_ba(bal_dem,      all_covs)   # ec+dem
ba_rspO <- pull_ba(bal_rsp_only, all_covs)   # rsp only
ba_lsO  <- pull_ba(bal_ls_only,  all_covs)   # ls only
ba_cem  <- pull_ba(bal_cem,      all_covs)   # CEM dem quartiles

# Reusable renderer: greyscale label column, optional shaded model pair.
render_smd <- function(df, title, note, shade_cols = integer(0), cex = 0.7) {
  nr <- nrow(df); nc <- ncol(df)
  fill <- matrix("white", nr, nc); fill[, 1] <- "grey95"
  if (length(shade_cols)) fill[, shade_cols] <- "grey82"
  hdr <- rep("grey90", nc); if (length(shade_cols)) hdr[shade_cols] <- "grey70"
  tt <- ttheme_default(
    core    = list(fg_params = list(cex = cex),
                   bg_params = list(fill = fill, col = "grey60")),
    colhead = list(fg_params = list(cex = cex, fontface = "bold"),
                   bg_params = list(fill = hdr, col = "grey60")))
  grid.newpage()
  grid.draw(arrangeGrob(
    textGrob(title, gp = gpar(fontsize = 11, fontface = "bold")),
    tableGrob(df, rows = NULL, theme = tt),
    textGrob(note, gp = gpar(fontsize = 7.5, col = "grey40")),
    ncol = 1, heights = c(0.10, 0.74, 0.16)))
}

# TABLE E1 — focal comparison: full topo model vs DEM only
smd_focal <- data.frame(
  Covariate     = all_covs,
  `ec+ls+rsp.b` = ba_topo$B, `ec+ls+rsp.a` = ba_topo$A,
  `dem.b`       = ba_demO$B, `dem.a`       = ba_demO$A,
  check.names   = FALSE
)
render_smd(smd_focal,
  "TABLE E1: SMD before (.b) vs after (.a) — focal comparison",
  paste(".b = before matching (unmatched SMD)   .a = after matching   NA = covariate not in model",
        "ec+ls+rsp = full topo-derivative model;   dem = DEM only (shaded).", sep = "\n"),
  shade_cols = 4:5, cex = 0.75)

# TABLE E2 — supporting models
smd_other <- data.frame(
  Covariate  = all_covs,
  `ec.b`     = ba_ecO$B,  `ec.a`     = ba_ecO$A,
  `ec+dem.b` = ba_ecd$B,  `ec+dem.a` = ba_ecd$A,
  `rsp.b`    = ba_rspO$B, `rsp.a`    = ba_rspO$A,
  `ls.b`     = ba_lsO$B,  `ls.a`     = ba_lsO$A,
  `cem.b`    = ba_cem$B,  `cem.a`    = ba_cem$A,
  check.names = FALSE
)
render_smd(smd_other,
  "TABLE E2: SMD before (.b) vs after (.a) — supporting models",
  paste(".b = before matching (unmatched SMD)   .a = after matching   NA = covariate not in model",
        "ec = EC only;   ec+dem;   rsp = RSP only;   ls = LS only;   cem = CEM dem quartiles.", sep = "\n"),
  cex = 0.62)

# --- TABLE F: Sensitivity — E-values + Rosenbaum Gamma ---
tryCatch({

  sd_yield <- sd(map_sp@data$Yield, na.rm = TRUE)

  # E-VALUES via EValue::evalues.OLS (VanderWeele & Ding 2017), the reference
  # implementation for a continuous (OLS) outcome. It standardizes est/sd
  # internally, converts to an approximate RR, and returns a 2x3 matrix:
  #   row "E-values": column "point" = E-value for the estimate;
  #   the non-NA of "lower"/"upper" = E-value for the CI limit closer to null.
  # delta = 1 (binary exposure contrast, treated vs control); true = 0 (null).
  if (!requireNamespace("EValue", quietly = TRUE)) install.packages("EValue")
  library(EValue)
  ev_est <- numeric(length(all_est))
  ev_ci  <- numeric(length(all_est))
  for (i in seq_along(all_est)) {
    ev <- evalues.OLS(est = all_est[i], se = all_se[i], sd = sd_yield,
                      delta = 1, true = 0)
    ev_est[i] <- round(ev["E-values", "point"], 2)
    ci <- stats::na.omit(c(ev["E-values", "lower"], ev["E-values", "upper"]))
    ev_ci[i]  <- round(if (length(ci)) ci[1] else 1, 2)
  }

  # ROSENBAUM GAMMA via sensitivityfull::senfm() — the method designed for FULL
  # matching (Rosenbaum 2007, Biometrics 63:456-464; Huber's M-statistic).
  if (!requireNamespace("sensitivityfull", quietly = TRUE)) install.packages("sensitivityfull")
  library(sensitivityfull)

  # Convert a MatchIt full-match matched dataset to senfm's input.
  # In full matching every subclass is either 1 treated : k controls, OR
  # m treated : 1 control. senfm wants the SINGLETON in column 1 of each row,
  # the rest in the remaining columns (NA-padded), and a logical treated1 that
  # is TRUE when that singleton is the treated unit, FALSE when it is the control.
  build_senfm <- function(mdf, outcome = "Yield", treat = "Treat") {
    rows <- list(); treated1 <- logical(0)
    for (sc in unique(mdf$subclass)) {
      s  <- mdf[mdf$subclass == sc, ]
      tr <- as.integer(as.character(s[[treat]]))
      yt <- s[[outcome]][tr == 1]; yc <- s[[outcome]][tr == 0]
      if (!length(yt) || !length(yc)) next
      if (length(yt) == 1) {            # 1 treated : k controls -> treated first
        rows[[length(rows) + 1]] <- c(yt, yc); treated1 <- c(treated1, TRUE)
      } else if (length(yc) == 1) {     # m treated : 1 control  -> control first
        rows[[length(rows) + 1]] <- c(yc, yt); treated1 <- c(treated1, FALSE)
      } else next                       # m:n set (e.g., CEM stratum) -> not full matching
    }
    # If few/no sets are valid full-matching sets (e.g. CEM, whose strata are
    # m treated : n controls), senfm does not apply. Signal that to the caller.
    if (length(rows) < 2) return(NULL)
    J <- max(lengths(rows))
    y <- t(sapply(rows, function(r) c(r, rep(NA, J - length(r)))))
    list(y = y, treated1 = treated1)
  }

  # Critical Gamma = the hidden-bias odds ratio at which the (one-sided, positive
  # effect) sensitivity p-value first exceeds 0.05. senfm's p-value increases
  # monotonically in gamma, so bisection finds the crossing efficiently.
  gamma_senfm <- function(mdf, g_max = 1000, tol = 0.01) {
    d <- build_senfm(mdf)
    if (is.null(d) || nrow(d$y) < 2) return(NA_real_)   # not full matching (e.g. CEM)
    pval_at <- function(g) senfm(d$y, d$treated1, gamma = g, alternative = "greater")$pval
    if (pval_at(1) > 0.05) return(1)              # not significant even at Gamma = 1
    lo <- 1; hi <- 2
    while (pval_at(hi) <= 0.05) {                 # expand until the p-value crosses 0.05
      lo <- hi; hi <- hi * 2
      if (hi > g_max) return(paste0(">", g_max))
    }
    while (hi - lo > tol) {                        # bisect to locate the crossing
      mid <- (lo + hi) / 2
      if (pval_at(mid) <= 0.05) lo <- mid else hi <- mid
    }
    round((lo + hi) / 2, 2)
  }

  # CEM strata are m treated : n controls, which senfm does NOT handle. The
  # design-appropriate Rosenbaum analysis for stratified m:n data is senstrat()
  # (Rosenbaum/Gastwirth-Krieger), using the same Huber M-score family (mscores).
  if (!requireNamespace("senstrat", quietly = TRUE)) install.packages("senstrat")
  library(senstrat)

  # Critical Gamma for a stratified (CEM) design, same bisection logic as senfm.
  # mscores() turns the outcome into within-stratum Huber M-scores; senstrat()
  # returns the worst-case one-sided P-value bound in $Result["P-value"].
  gamma_senstrat <- function(mdf, outcome = "Yield", treat = "Treat",
                             g_max = 1000, tol = 0.01) {
    z  <- as.integer(as.character(mdf[[treat]]))
    st <- as.integer(mdf$subclass)
    sc <- mscores(mdf[[outcome]], z, st)
    pval_at <- function(g)
      as.numeric(senstrat(sc, z, st, gamma = g, alternative = "greater")$Result["P-value"])
    if (pval_at(1) > 0.05) return(1)              # not significant even at Gamma = 1
    lo <- 1; hi <- 2
    while (pval_at(hi) <= 0.05) {
      lo <- hi; hi <- hi * 2
      if (hi > g_max) return(paste0(">", g_max))
    }
    while (hi - lo > tol) {
      mid <- (lo + hi) / 2
      if (pval_at(mid) <= 0.05) lo <- mid else hi <- mid
    }
    round((lo + hi) / 2, 2)
  }

  # Sanity check: at Gamma = 1 both p-values should be small for the focal/CEM
  # models (the effect is strongly significant); print them for auditing.
  cat(sprintf("senfm    Gamma=1 p (ec+ls+rsp): %.4g\n",
              senfm(build_senfm(matched_data)$y,
                    build_senfm(matched_data)$treated1, gamma = 1,
                    alternative = "greater")$pval))
  cat(sprintf("senstrat Gamma=1 p (CEM dem):   %.4g\n",
              gamma_senstrat_p1 <- as.numeric(senstrat(
                mscores(matched_cem$Yield, as.integer(as.character(matched_cem$Treat)),
                        as.integer(matched_cem$subclass)),
                as.integer(as.character(matched_cem$Treat)),
                as.integer(matched_cem$subclass), gamma = 1)$Result["P-value"])))

  # 6 PSM models use senfm (full matching); CEM uses senstrat (stratified m:n).
  gamma_psm <- sapply(
    list(matched_data, matched_ec_only, matched_dem,
         matched_dem_only, matched_rsp_only, matched_ls_only),
    gamma_senfm)
  gamma_cem <- gamma_senstrat(matched_cem)
  gamma_vals <- c(gamma_psm, gamma_cem)   # order matches all_mod[-1]

  sensitivity_table <- data.frame(
    Model   = all_mod,
    ATE     = round(all_est, 2),
    SE      = round(all_se,  2),
    EV.est  = ev_est,
    EV.CI   = ev_ci,
    Gamma   = c("N/A", gamma_vals),
    stringsAsFactors = FALSE
  )

  plot_table(sensitivity_table,
             "TABLE F: Sensitivity — E-values and Rosenbaum Gamma",
             note = paste("EV = E-value (VanderWeele & Ding 2017, via EValue::evalues.OLS): RR-scale confounding to explain away the effect; higher = more robust.",
                          "Gamma = Rosenbaum bias odds ratio at which the one-sided p first exceeds 0.05; higher = more robust. PSM models use senfm (full matching);",
                          "CEM uses senstrat (stratified m:n) — both Rosenbaum M-statistic analyses. Gamma = 1 means not significant even with no hidden bias.", sep = "\n"))

}, error = function(e) {
  grid.newpage()
  grid.draw(textGrob(
    paste0("TABLE F ERROR — check console:\n\n", conditionMessage(e)),
    gp = gpar(fontsize = 10, col = "darkred")
  ))
})
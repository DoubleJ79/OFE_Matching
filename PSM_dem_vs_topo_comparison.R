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
# Convert the SpatialPolygonsDataFrame to a Spatial object
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

#center and scale the covariatesn
map_sp@data$ls <- scale(map_sp@data$ls, center = TRUE, scale = FALSE)
map_sp@data$rsp <- scale(map_sp@data$rsp, center = TRUE, scale = FALSE)
map_sp@data$ec <- scale(map_sp@data$ec, center = TRUE, scale = FALSE)
map_sp@data$dem <- scale(map_sp@data$dem, center = TRUE, scale = FALSE)

#write map_sp@data to .csv
write.csv(map_sp@data, "map_sp_data.csv")

#Probit and SAR/SEM Probit // Mean 

model <- Treat ~ rsp + ls + ec
#mod1 <- probit(formula = model, data = map_sp@data, method = "ML")
mod_SAR <- sarprobit(model, w_Adj, map_sp@data)
summary(mod_SAR) 
impacts(mod_SAR)
logLik(mod_SAR)
AIC(mod_SAR)
BIC(mod_SAR)


#### PROPENSITY SCORE CALCULATION

predicted_values <- fitted(mod_SAR)
p_score <- 1 / (1 + exp(-predicted_values))


# MATCHIT method 
matchit_test <- matchit(model, data = map_sp@data, method = "full", pscores = p_score) # full and quick methods seems to preform best and have best balance
#matchit_test <- matchit(model, data = map_sp@data, exact = ~SoilTyp, method = "nearest", distance = "logit", pscores = p_score) #by soil type
summary(matchit_test)

bal.tab(matchit_test, un = TRUE)  # balance table

plot(matchit_test, type = "jitter", interactive = FALSE)
plot(summary(matchit_test))
plot(matchit_test, type = "qq")
for (v in c("rsp", "ls", "ec")) {
  print(bal.plot(matchit_test, var.name = v, which = "both",
                 type = "ecdf", mirror = FALSE) +
          ggplot2::ggtitle(paste0("eCDF: ", v, " — ec+ls+rsp matching")))
}
love.plot(matchit_test, thresholds = 0.1)  # visual SMDs


matched_data <- match.data(matchit_test)
#relevel the Treatment factor to ensure it is ordered correctly for contrasts
matched_data$Treat <- factor(matched_data$Treat)
matched_data$rsp <- as.numeric(matched_data$rsp)
matched_data$ls <- as.numeric(matched_data$ls)
matched_data$ec <- as.numeric(matched_data$ec)


##### EFFECT ESTIMATION
#This is the doubly robust treatment effect estimate. Including covar in this model covers you for imbalance twice
#this method is identical to t test without covar
#fit1 <- lm(Yield ~ Treat * SoilType, data = matched_data, weights = weights) # for by soil type
fit1 <- lm(Yield ~ Treat * (rsp+ls+ec), data = matched_data, weights = weights) # for marginal effects
treat_effect_avgcomp <- avg_comparisons(fit1, variables = "Treat", # for treatment effect
                                        vcov = ~subclass, wts = "weights")
#newdata = subset(matched_data, Treat == 1),# remove this line if estimating ATE rather than ATT
treat_effect_avgcomp
# Summary of the lm model
summary(fit1)

# Create a SpatialPointsDataFrame using centroid_x and centroid_y as coordinates.
# Drop the 'matchdata' class
md_df <- as.data.frame(matched_data)

# Build coords matrix
coords_mat <- as.matrix(md_df[, c("centroid_x", "centroid_y")])

write.csv(md_df, file = "RAinbarrel-matched_data09Oct2025.csv", row.names = FALSE)

# Create the SpatialPointsDataFrame
matched_data_sp <- sp::SpatialPointsDataFrame(
  coords      = coords_mat,
  data        = md_df,
  proj4string = sp::CRS(sf::st_crs(yield)$wkt)
)

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
plot_table <- function(df, title) {
  tt <- ttheme_default(
    core    = list(fg_params = list(cex = 0.8)),
    colhead = list(fg_params = list(cex = 0.85, fontface = "bold"))
  )
  tbl <- tableGrob(df, rows = NULL, theme = tt)
  title_grob <- textGrob(title, gp = gpar(fontsize = 11, fontface = "bold"))
  note_grob  <- textGrob(
    "*** p<0.001  ** p<0.01  * p<0.05  ns = not significant",
    gp = gpar(fontsize = 8, col = "grey40")
  )
  grid.newpage()
  grid.draw(arrangeGrob(title_grob, tbl, note_grob,
                        ncol = 1, heights = c(0.08, 0.84, 0.08)))
}

# --- Shared objects ---
z_test <- function(est1, se1, est2, se2) {
  z <- (est1 - est2) / sqrt(se1^2 + se2^2)
  p <- 2 * pnorm(-abs(z))
  c(z = z, p = p)
}
sig_star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                        ifelse(p < 0.05,  "*",   "ns")))

est_topo     <- treat_effect_avgcomp$estimate; se_topo     <- treat_effect_avgcomp$std.error
est_dem      <- te_dem$estimate;               se_dem      <- te_dem$std.error
est_dem_only <- te_dem_only$estimate;          se_dem_only <- te_dem_only$std.error

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
plot_table(ate_table, "TABLE A: Treatment Effect Estimates by Model")

# --- TABLE B: Each model vs no-matching baseline (Bonferroni corrected) ---
zp_list <- lapply(seq_along(all_est[-1]), function(i)
  z_test(all_est[i + 1], all_se[i + 1], est_unmatched, se_unmatched))

raw_p   <- sapply(zp_list, `[`, "p")
bonf_p  <- pmin(raw_p * length(raw_p), 1)

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
psm_models <- c("PSM: ec+ls+rsp", "PSM: ec only", "PSM: ec+dem",
                "PSM: dem only", "PSM: rsp only", "PSM: ls only")
psm_est    <- c(est_topo, est_ec_only, est_dem, est_dem_only, est_rsp_only, est_ls_only)
psm_se     <- c(se_topo,  se_ec_only,  se_dem,  se_dem_only,  se_rsp_only,  se_ls_only)

zp_cem <- lapply(seq_along(psm_est), function(i)
  z_test(psm_est[i], psm_se[i], est_cem, se_cem))

raw_p_cem  <- sapply(zp_cem, `[`, "p")
bonf_p_cem <- pmin(raw_p_cem * length(raw_p_cem), 1)

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
plot_table(psm_fit_table, "TABLE D: SAR Probit Propensity Score Model Fit")

# --- TABLE E: SMD balance comparison ---
bal_topo     <- bal.tab(matchit_test,     un = TRUE)
bal_dem      <- bal.tab(matchit_dem,      un = TRUE)
bal_dem_only <- bal.tab(matchit_dem_only, un = TRUE)
bal_ec_only  <- bal.tab(matchit_ec_only,  un = TRUE)
bal_rsp_only <- bal.tab(matchit_rsp_only, un = TRUE)
bal_ls_only  <- bal.tab(matchit_ls_only,  un = TRUE)
bal_cem      <- bal.tab(matchit_cem,      un = TRUE)

pull_smd <- function(bal, covs) {
  idx <- match(covs, rownames(bal$Balance))
  round(bal$Balance$Diff.Adj[idx], 3)
}
all_covs <- unique(c(rownames(bal_topo$Balance), rownames(bal_dem$Balance),
                     rownames(bal_dem_only$Balance), rownames(bal_ec_only$Balance),
                     rownames(bal_rsp_only$Balance), rownames(bal_ls_only$Balance),
                     rownames(bal_cem$Balance)))

balance_compare <- data.frame(
  Covariate     = all_covs,
  SMD_before    = round(bal_topo$Balance$Diff.Un[
                    match(all_covs, rownames(bal_topo$Balance))], 3),
  PSM_ec_ls_rsp = pull_smd(bal_topo,     all_covs),
  PSM_ec_only   = pull_smd(bal_ec_only,  all_covs),
  PSM_ec_dem    = pull_smd(bal_dem,      all_covs),
  PSM_dem_only  = pull_smd(bal_dem_only, all_covs),
  PSM_rsp_only  = pull_smd(bal_rsp_only, all_covs),
  PSM_ls_only   = pull_smd(bal_ls_only,  all_covs),
  CEM_dem_Q4    = pull_smd(bal_cem,      all_covs)
)
plot_table(balance_compare,
           "TABLE E: Standardised Mean Differences\n(NA = covariate not in that model)")

# --- TABLE F: Sensitivity — E-values + Rosenbaum Gamma ---
tryCatch({

  sd_yield <- sd(map_sp@data$Yield, na.rm = TRUE)

  # E-values: VanderWeele & Ding 2017, continuous outcome, manual calculation
  ev_est <- numeric(length(all_est))
  ev_ci  <- numeric(length(all_est))
  for (i in seq_along(all_est)) {
    d1  <- all_est[i] / sd_yield
    d2  <- max(0, abs(all_est[i]) - 1.96 * all_se[i]) / sd_yield
    rr1 <- exp(0.91 * abs(d1))
    rr2 <- exp(0.91 * abs(d2))
    ev_est[i] <- round(rr1 + sqrt(rr1 * (rr1 - 1)), 2)
    ev_ci[i]  <- round(if (rr2 <= 1) 1 else rr2 + sqrt(rr2 * (rr2 - 1)), 2)
  }

  # Rosenbaum Gamma via rbounds::psens(x, y)
  if (!requireNamespace("rbounds", quietly = TRUE)) install.packages("rbounds")
  library(rbounds)

  get_xy <- function(mdf) {
    x <- y <- numeric(0)
    for (sc in unique(mdf$subclass)) {
      s   <- mdf[mdf$subclass == sc, ]
      tr  <- as.integer(as.character(s$Treat))
      t_m <- mean(s$Yield[tr == 1])
      c_m <- mean(s$Yield[tr == 0])
      if (is.finite(t_m) && is.finite(c_m)) { x <- c(x, t_m); y <- c(y, c_m) }
    }
    list(x = x, y = y)
  }

  # psens() returns a list; the bounds table is in res$bounds with columns
  # "Gamma", "Lower bound", "Upper bound". Critical Gamma = first value where
  # the Upper bound (conservative) p-value exceeds 0.05.
  gamma_for <- function(mdf, max_g = 100, step = 0.1) {
    xy <- get_xy(mdf)
    if (length(xy$x) < 2) return(NA_real_)
    res <- tryCatch(psens(xy$x, xy$y, Gamma = max_g, GammaInc = step),
                    error = function(e) NULL)
    if (is.null(res)) return(NA_real_)
    b   <- res$bounds
    ub  <- b[["Upper bound"]]
    g   <- b[["Gamma"]]
    idx <- which(ub > 0.05)
    if (length(idx) == 0) NA_real_ else round(g[idx[1]], 2)
  }

  gamma_vals <- sapply(
    list(matched_data, matched_ec_only, matched_dem,
         matched_dem_only, matched_rsp_only, matched_ls_only, matched_cem),
    gamma_for)

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
             "TABLE F: Sensitivity — E-values and Rosenbaum Gamma")

}, error = function(e) {
  grid.newpage()
  grid.draw(textGrob(
    paste0("TABLE F ERROR — check console:\n\n", conditionMessage(e)),
    gp = gpar(fontsize = 10, col = "darkred")
  ))
})
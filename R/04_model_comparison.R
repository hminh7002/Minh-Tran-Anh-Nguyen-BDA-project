# =============================================================================
# 04_model_comparison.R  -  predictive model comparison with cross-validation
#
# Two kinds of cross-validation, because the prediction task matters:
#   (a) PSIS-LOO, leave-one-EYE-out: "predict a new eye of a patient we already know".
#       Cheap (no refitting), but with only 2 eyes per patient the hierarchical models
#       can have many high Pareto-k values -> estimates then unreliable.
#   (b) Grouped 10-fold CV, leave-PATIENTS-out: "predict both eyes of a NEW patient".
#       Refits each model 10 times; this is the comparison we report.
#
# Outputs
#   results/loo_eye.csv          elpd_loo, p_loo, Pareto-k summary per model
#   results/compare_eye.csv      loo_compare table (leave-one-eye-out)
#   results/compare_patient.csv  loo_compare table (leave-patients-out, 10-fold)
#   figures/fig6_pareto_k.png    Pareto-k diagnostics per model
#   fits/kfold_*.rds             cached K-fold results (git-ignored)
#
# Needs the fits from 03_fit_models.R. Run:  source("R/04_model_comparison.R")
# First run ~5-10 min (30 refits), afterwards cached.
# =============================================================================

source("R/priors.R")
library(loo)
library(posterior)

options(brms.backend = "cmdstanr", mc.cores = 4)
m1 <- readRDS("fits/m1.rds"); m2 <- readRDS("fits/m2.rds"); m3 <- readRDS("fits/m3.rds")
models <- list(M1 = m1, M2 = m2, M3 = m3)

# ---- (a) PSIS-LOO, leave-one-eye-out --------------------------------------------
loos <- lapply(models, loo)
k_tab <- do.call(rbind, lapply(names(loos), function(nm) {
  l <- loos[[nm]]; k <- l$diagnostics$pareto_k
  data.frame(model = nm,
             elpd_loo = round(l$estimates["elpd_loo", "Estimate"], 1),
             se = round(l$estimates["elpd_loo", "SE"], 1),
             p_loo = round(l$estimates["p_loo", "Estimate"], 1),
             n_params = sum(grepl("^(b_|sd_|shape|r_id)", variables(models[[nm]]))),
             k_good = sum(k <= 0.7), k_bad = sum(k > 0.7 & k <= 1), k_very_bad = sum(k > 1))
}))
write.csv(k_tab, "results/loo_eye.csv", row.names = FALSE)
cat("\n==== PSIS-LOO, leave-one-eye-out ====\n"); print(k_tab, row.names = FALSE)

cmp_eye <- loo_compare(loos)
write.csv(as.data.frame(cmp_eye), "results/compare_eye.csv")
cat("\nComparison (leave-one-eye-out):\n"); print(cmp_eye, simplify = FALSE, digits = 1)

png("figures/fig6_pareto_k.png", width = 1500, height = 420, res = 110)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
for (nm in names(loos)) plot(loos[[nm]], main = paste(nm, "Pareto k (leave-one-eye-out)"),
                             label_points = FALSE)
dev.off()

# ---- (b) grouped 10-fold CV, leave-patients-out --------------------------------
# The same folds for all models, so the comparison is paired.
set.seed(2026)
folds <- loo::kfold_split_grouped(K = 10, x = d$id)
kf <- lapply(names(models), function(nm) {
  f <- file.path("fits", paste0("kfold_", tolower(nm), ".rds"))
  if (file.exists(f)) return(readRDS(f))
  r <- kfold(models[[nm]], folds = folds, chains = 4, refresh = 0, seed = 2026)
  saveRDS(r, f); r
})
names(kf) <- names(models)

kf_tab <- do.call(rbind, lapply(names(kf), function(nm) {
  e <- kf[[nm]]$estimates
  data.frame(model = nm, elpd_kfold = round(e["elpd_kfold", "Estimate"], 1),
             se = round(e["elpd_kfold", "SE"], 1))
}))
cmp_pat <- loo_compare(kf)
write.csv(as.data.frame(cmp_pat), "results/compare_patient.csv")
cat("\n==== Grouped 10-fold CV, leave-patients-out (predicting new patients) ====\n")
print(kf_tab, row.names = FALSE)
cat("\nComparison (leave-patients-out):\n"); print(cmp_pat, simplify = FALSE, digits = 1)
cat("\nRule of thumb: |elpd_diff| < 4 = negligible; otherwise compare with its SE.\n")

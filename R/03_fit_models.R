# =============================================================================
# 03_fit_models.R  -  fit M1-M3, convergence diagnostics, posterior predictive checks
#
# Models (Weibull survival with right-censoring, priors from R/priors.R)
#   M1 pooled:          trt + type + laser + risk_c              (eyes independent)
#   M2 hierarchical:    M1 + (1 | id)                            (patient effect / frailty)
#   M3 interactions:    trt * type + trt * laser + risk_c + (1 | id)
#
# Outputs
#   fits/m1.rds, fits/m2.rds, fits/m3.rds     cached fits (git-ignored; re-used by later scripts)
#   results/convergence.csv                   R-hat, ESS, divergences, tree-depth hits per model
#   results/posterior_summary.csv             time ratios and other parameters per model
#   figures/fig3_ppc_km.png                   KM overlay: observed vs replicated, per model
#   figures/fig4_ppc_km_groups.png            same, split by treatment x diabetes type
#   figures/fig5_ppc_pairs.png                check of between-eye dependence
#
# Run from the project folder:  source("R/03_fit_models.R")   (first run ~3-5 min, then cached)
# =============================================================================

source("R/priors.R")
library(posterior)
library(bayesplot)
library(dplyr)

options(brms.backend = "cmdstanr", mc.cores = 4)
for (dd in c("fits", "results", "figures")) dir.create(dd, showWarnings = FALSE)
color_scheme_set("blue")

# ---- 1. fit the three models ------------------------------------------------------
# file = ... caches each fit; it is only refitted when the model, data or priors change.
fit <- function(formula, prior, name) {
  brm(formula, data = d, family = weibull(), prior = prior,
      chains = 4, iter = 2000, warmup = 1000, seed = 2026, refresh = 0,
      file = file.path("fits", name), file_refit = "on_change")
}
m1 <- fit(f_m1, priors_main, "m1")
m2 <- fit(f_m2, c(priors_main, prior_sd), "m2")
m3 <- fit(f_m3, c(priors_main, prior_sd), "m3")
models <- list(M1 = m1, M2 = m2, M3 = m3)
labels <- c(M1 = "M1: pooled", M2 = "M2: hierarchical (patient effect)",
            M3 = "M3: hierarchical + interactions")

# ---- 2. convergence diagnostics ----------------------------------------------------
conv_row <- function(m, name) {
  s  <- summarise_draws(as_draws_df(m), "rhat", "ess_bulk", "ess_tail")
  s  <- s[!grepl("^(lp__|lprior)", s$variable), ]
  np <- nuts_params(m)
  data.frame(
    model           = name,
    parameters      = nrow(s),
    max_rhat        = round(max(s$rhat, na.rm = TRUE), 3),
    min_bulk_ess    = round(min(s$ess_bulk, na.rm = TRUE)),
    min_tail_ess    = round(min(s$ess_tail, na.rm = TRUE)),
    divergences     = sum(np$Value[np$Parameter == "divergent__"]),
    max_treedepth   = sum(np$Value[np$Parameter == "treedepth__"] >= 10)
  )
}
conv <- do.call(rbind, Map(conv_row, models, names(models)))
write.csv(conv, "results/convergence.csv", row.names = FALSE)
cat("\n==== Convergence diagnostics (all parameters incl. 197 patient effects) ====\n")
print(conv, row.names = FALSE)
cat("Targets: max R-hat < 1.01, ESS > 400, 0 divergences, 0 max-treedepth hits.\n")

# ---- 3. posterior summaries (population-level parameters) ------------------------
par_summary <- function(m, name) {
  s <- summarise_draws(as_draws_df(m), "median", ~quantile(.x, c(0.05, 0.95)))
  s <- s[grepl("^b_|^shape$|^sd_", s$variable), ]
  s$model <- name
  # coefficients as time ratios (exp scale) for interpretation
  is_b <- grepl("^b_", s$variable) & s$variable != "b_Intercept"
  s$time_ratio_median <- ifelse(is_b, exp(s$median), NA)
  s$time_ratio_q5     <- ifelse(is_b, exp(s$`5%`), NA)
  s$time_ratio_q95    <- ifelse(is_b, exp(s$`95%`), NA)
  s
}
ps <- do.call(rbind, Map(par_summary, models, names(models)))
ps[, sapply(ps, is.numeric)] <- round(ps[, sapply(ps, is.numeric)], 3)
write.csv(ps, "results/posterior_summary.csv", row.names = FALSE)
cat("\n==== Posterior medians and 90% intervals (time ratio = exp(coefficient)) ====\n")
print(as.data.frame(ps[, c("model", "variable", "median", "5%", "95%",
                           "time_ratio_median", "time_ratio_q5", "time_ratio_q95")]),
      row.names = FALSE)

# ---- 4. posterior predictive checks ----------------------------------------------
# Kaplan-Meier of the observed (censored) data vs KM curves of replicated event times.
set.seed(2026)
yrep <- lapply(models, posterior_predict, ndraws = 100)

ppc_km <- lapply(names(models), function(nm) {
  ppc_km_overlay(d$futime, yrep[[nm]][1:50, ], status_y = d$status) +
    coord_cartesian(xlim = c(0, 80)) +
    labs(title = labels[[nm]], x = "months", y = "share of eyes without vision loss") +
    theme_project() + theme(legend.position = c(0.85, 0.85), legend.title = element_blank())
})
g3 <- bayesplot_grid(plots = ppc_km, grid_args = list(nrow = 1))
ggsave("figures/fig3_ppc_km.png", g3, width = 15, height = 4.5, dpi = 150, bg = "white")

# by treatment x diabetes type: can each model reproduce the four groups?
grp <- interaction(ifelse(d$trt == 1, "treated", "untreated"), d$type, sep = ", ")
ppc_grp <- lapply(names(models), function(nm) {
  ppc_km_overlay_grouped(d$futime, yrep[[nm]][1:50, ], group = grp, status_y = d$status) +
    coord_cartesian(xlim = c(0, 80)) +
    labs(title = labels[[nm]], x = "months", y = "without vision loss") +
    theme_project() + theme(legend.position = "none")
})
g4 <- bayesplot_grid(plots = ppc_grp, grid_args = list(ncol = 1))
ggsave("figures/fig4_ppc_km_groups.png", g4, width = 13, height = 11, dpi = 150, bg = "white")

# ---- 5. targeted check: are the two eyes of a patient dependent? ----------------
# Statistic T = share of patients in whom BOTH eyes lost vision within 36 months,
# among patients whose status at 36 months is known for both eyes
# (each eye: vision loss by 36 months, or followed beyond 36 months).
t0 <- 36
known <- with(d, (status == 1 & futime <= t0) | futime > t0)
ok_id <- names(which(tapply(known, d$id, all)))
sel <- d$id %in% ok_id
both_lost <- function(time_vec) {
  lost <- time_vec[sel] <= t0
  mean(tapply(lost, d$id[sel], all))
}
T_obs <- both_lost(ifelse(d$status == 1, d$futime, Inf))
T_rep <- lapply(models, function(m) {
  yr <- posterior_predict(m, ndraws = 1000)
  apply(yr, 1, both_lost)
})
ppp <- sapply(T_rep, function(x) mean(x >= T_obs))
cat(sprintf("\n==== Between-eye dependence: share of patients (n = %d) with both eyes lost by %d months ====\n",
            length(ok_id), t0))
cat(sprintf("observed: %.3f\n", T_obs))
print(data.frame(model = names(T_rep),
                 replicated_median = sapply(T_rep, median),
                 q5 = sapply(T_rep, quantile, 0.05), q95 = sapply(T_rep, quantile, 0.95),
                 p_value = round(ppp, 3)), row.names = FALSE, digits = 3)

df5 <- data.frame(model = rep(labels[names(T_rep)], each = 1000), T = unlist(T_rep))
p5 <- ggplot(df5, aes(T)) +
  geom_histogram(bins = 30, fill = "#2a78d6", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = T_obs, colour = "#e34948", linewidth = 1) +
  facet_wrap(~ model, nrow = 1) +
  labs(title = sprintf("Both eyes lost vision within %d months: replicated (bars) vs observed (red line)", t0),
       x = "share of patients", y = "replicated data sets") +
  theme_project()
ggsave("figures/fig5_ppc_pairs.png", p5, width = 13, height = 3.8, dpi = 150, bg = "white")

cat("\nSaved results/convergence.csv, results/posterior_summary.csv,\n",
    "figures/fig3_ppc_km.png, figures/fig4_ppc_km_groups.png, figures/fig5_ppc_pairs.png\n")

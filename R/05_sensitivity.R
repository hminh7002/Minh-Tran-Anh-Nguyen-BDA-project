# =============================================================================
# 05_sensitivity.R  -  prior sensitivity analysis
#
# Part A  power-scaling sensitivity with priorsense (Kallioinen et al. 2023):
#         each prior / the likelihood is raised to a power close to 1 and the change
#         in the posterior is measured. Diagnosis per parameter:
#           prior-data conflict  = sensitive to both prior and likelihood
#           weak likelihood      = sensitive to prior, not to likelihood
# Part B  refit M3 with the WIDER priors from R/priors.R and compare key results.
#
# Outputs
#   results/priorsense_m2.csv, results/priorsense_m3.csv
#   results/sensitivity_refit.csv
#   figures/fig7_powerscale_m3.png
#   figures/fig8_sensitivity_refit.png
#
# Needs the fits from 03_fit_models.R. Run:  source("R/05_sensitivity.R")
# =============================================================================

source("R/priors.R")
library(posterior)
library(priorsense)
library(dplyr)

options(brms.backend = "cmdstanr", mc.cores = 4)
m2 <- readRDS("fits/m2.rds"); m3 <- readRDS("fits/m3.rds")

# ---- A. power-scaling sensitivity ------------------------------------------------
key_vars_m2 <- c("b_trt", "b_typeadult", "b_laserargon", "b_risk_c", "shape", "sd_id__Intercept")
key_vars_m3 <- c(key_vars_m2, "b_trt:typeadult", "b_trt:laserargon")

ps2 <- powerscale_sensitivity(m2, variable = key_vars_m2)
ps3 <- powerscale_sensitivity(m3, variable = key_vars_m3)
cat("\n==== priorsense: M2 ====\n"); print(ps2)
cat("\n==== priorsense: M3 ====\n"); print(ps3)
write.csv(as.data.frame(ps2), "results/priorsense_m2.csv", row.names = FALSE)
write.csv(as.data.frame(ps3), "results/priorsense_m3.csv", row.names = FALSE)

p7 <- powerscale_plot_dens(m3, variable = c("b_trt", "b_trt:typeadult", "sd_id__Intercept", "shape")) +
  theme_project()
ggsave("figures/fig7_powerscale_m3.png", p7, width = 12, height = 6, dpi = 150, bg = "white")

# ---- B. refit M3 with wider priors -------------------------------------------------
m3_wide <- brm(f_m3, data = d, family = weibull(),
               prior = c(priors_wide, prior_sd_wide),
               chains = 4, iter = 2000, warmup = 1000, seed = 2026, refresh = 0,
               file = "fits/m3_wide", file_refit = "on_change")

summ <- function(m, label) {
  dr <- as_draws_df(m)
  q <- function(x) quantile(x, c(0.05, 0.5, 0.95))
  rbind(
    data.frame(prior = label, quantity = "time ratio, juvenile onset",
               t(q(exp(dr$b_trt)))),
    data.frame(prior = label, quantity = "time ratio, adult onset",
               t(q(exp(dr$b_trt + dr$`b_trt:typeadult`)))),
    data.frame(prior = label, quantity = "extra effect in adults (ratio)",
               t(q(exp(dr$`b_trt:typeadult`)))),
    data.frame(prior = label, quantity = "argon vs xenon (ratio)",
               t(q(exp(dr$`b_trt:laserargon`)))),
    data.frame(prior = label, quantity = "sd of patient effects",
               t(q(dr$sd_id__Intercept))),
    data.frame(prior = label, quantity = "Weibull shape",
               t(q(dr$shape)))
  )
}
sens <- rbind(summ(m3, "main priors"), summ(m3_wide, "wide priors"))
names(sens)[3:5] <- c("q5", "median", "q95")
sens[, 3:5] <- round(sens[, 3:5], 2)
sens$quantity <- factor(sens$quantity, levels = unique(sens$quantity))  # keep the order above
write.csv(sens, "results/sensitivity_refit.csv", row.names = FALSE)
cat("\n==== M3 with main vs wide priors (median and 90% interval) ====\n")
print(sens, row.names = FALSE)

p8 <- ggplot(sens, aes(x = median, y = prior, colour = prior)) +
  geom_pointrange(aes(xmin = q5, xmax = q95)) +
  facet_wrap(~ quantity, scales = "free_x", ncol = 3) +
  scale_colour_manual(values = c("main priors" = "#2a78d6", "wide priors" = "#eb6834")) +
  labs(title = "Sensitivity of M3 results to the priors (median and 90% interval)",
       x = NULL, y = NULL) +
  theme_project() + theme(legend.position = "none")
ggsave("figures/fig8_sensitivity_refit.png", p8, width = 12, height = 5, dpi = 150, bg = "white")
cat("\nSaved results/*.csv and figures/fig7_powerscale_m3.png, figures/fig8_sensitivity_refit.png\n")

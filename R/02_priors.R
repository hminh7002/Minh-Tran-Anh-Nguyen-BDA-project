# =============================================================================
# 02_priors.R  -  priors and prior predictive checks
#
# What it does
#   1. defines the priors for models M1-M3 (one place, reused by later scripts)
#   2. samples from the PRIOR ONLY of the hierarchical model M2
#   3. checks what the priors imply in months and probabilities, before seeing data
#   4. saves figures/fig2_prior_predictive.png and prints a summary table
#
# Run from the project folder:  source("R/02_priors.R")
# =============================================================================

library(survival)
library(brms)
library(posterior)
library(bayesplot)
library(ggplot2)

options(brms.backend = "cmdstanr", mc.cores = 4)
dir.create("figures", showWarnings = FALSE)

# ---- data (same preprocessing as everywhere else) ----------------------------
data(retinopathy, package = "survival")
d <- within(retinopathy, {
  cens   <- 1 - status            # 1 = right-censored, 0 = vision loss observed
  risk_c <- risk - mean(risk)     # centred risk score (mean 0)
  type   <- factor(type)          # levels: juvenile (reference), adult
  laser  <- factor(laser)         # levels: xenon (reference), argon
})

# ---- 1. priors ---------------------------------------------------------------
# brms Weibull: log(mean time in months) = Intercept + b * predictors;
# exp(b) is a TIME RATIO (b > 0 means vision loss happens later).
# brms centres the predictors internally, so the Intercept prior is for an
# "average" eye.
priors_main <- c(
  prior(normal(4.1, 1),   class = Intercept),  # log(60): mean time ~5 years; 95%: ~8 months to ~40 years
  prior(normal(0, 0.75),  class = b),          # time ratios mostly within 1/4.4 to 4.4 (95%)
  prior(gamma(4, 4),      class = shape)       # Weibull shape: mean 1 (constant hazard), 90% within 0.34-1.94
)
prior_sd <- prior(exponential(2), class = sd)  # patient variation on log-time scale, mean 0.5

# wider alternative priors for the sensitivity analysis later (report section 9)
priors_wide <- c(
  prior(normal(4.1, 2),  class = Intercept),
  prior(normal(0, 2),    class = b),
  prior(gamma(1, 0.5),   class = shape)
)
prior_sd_wide <- prior(exponential(0.5), class = sd)

saveRDS(list(main = priors_main, sd = prior_sd, wide = priors_wide, sd_wide = prior_sd_wide),
        "R/priors.rds")

# ---- 2. sample from the prior only (model M2) ---------------------------------
m2_prior <- brm(
  futime | cens(cens) ~ trt + type + laser + risk_c + (1 | id),
  data = d, family = weibull(),
  prior = c(priors_main, prior_sd),
  sample_prior = "only",
  chains = 4, iter = 1000, seed = 2026, refresh = 0
)

# ---- 3. what do the priors imply? --------------------------------------------
dr <- as_draws_df(m2_prior)
k  <- dr$shape

# typical (random effect = 0) untreated juvenile, xenon, average risk eye
eta_untr <- dr$b_Intercept                      # b_Intercept is on the un-centred scale in draws
eta_trt  <- dr$b_Intercept + dr$b_trt
weib_S   <- function(t, eta, k) {               # survival probability S(t)
  scale <- exp(eta) / gamma(1 + 1 / k)
  exp(-(t / scale)^k)
}
median_time <- function(eta, k) exp(eta) / gamma(1 + 1 / k) * log(2)^(1 / k)

summ <- function(x) round(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), 2)
tab <- rbind(
  "median time to vision loss, untreated (months)" = summ(median_time(eta_untr, k)),
  "P(vision loss by 60 months), untreated"         = summ(1 - weib_S(60, eta_untr, k)),
  "P(vision loss by 60 months), treated"           = summ(1 - weib_S(60, eta_trt, k)),
  "time ratio treated / untreated"                 = summ(exp(dr$b_trt)),
  "Weibull shape"                                  = summ(k),
  "sd of patient effects (log scale)"              = summ(dr$sd_id__Intercept)
)
cat("\nPrior implications (5%, 25%, 50%, 75%, 95% quantiles):\n")
print(tab)

# ---- 4. prior predictive check figure -----------------------------------------
yrep <- posterior_predict(m2_prior, ndraws = 100)
yrep_cap <- pmin(yrep, 200)                     # cap extreme draws so the plot is readable

p1 <- ppc_km_overlay(d$futime, yrep_cap[1:50, ], status_y = d$status) +
  coord_cartesian(xlim = c(0, 80)) +
  labs(title = "Prior predictive: survival curves vs observed data",
       x = "months", y = "share of eyes without vision loss")

p2 <- ggplot(data.frame(m = median_time(eta_untr, k)), aes(m)) +
  geom_histogram(bins = 60, fill = "#2a78d6", alpha = 0.7) +
  scale_x_log10(breaks = c(1, 12, 60, 600),
                labels = c("1 month", "1 year", "5 years", "50 years")) +
  labs(title = "Prior: median time to vision loss (untreated eye)", x = NULL, y = "draws")

p3 <- ggplot(data.frame(tr = exp(dr$b_trt)), aes(tr)) +
  geom_histogram(bins = 60, fill = "#eb6834", alpha = 0.7) +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10)) +
  geom_vline(xintercept = 1, linetype = 2) +
  labs(title = "Prior: time ratio of laser treatment", x = "time ratio (log scale)", y = "draws")

fig <- bayesplot::bayesplot_grid(p1, p2, p3, grid_args = list(nrow = 1))
ggsave("figures/fig2_prior_predictive.png", fig, width = 15, height = 4.5, dpi = 150)
cat("\nSaved figures/fig2_prior_predictive.png\n")

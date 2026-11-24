# =============================================================================
# 02_priors.R  -  prior predictive checks
#
# What it does
#   1. loads data + priors from R/priors.R
#   2. samples from the PRIOR ONLY of the hierarchical model M2
#   3. checks what the priors imply in months and probabilities, before seeing data
#   4. saves figures/fig2_prior_predictive.png and prints a summary table
#
# Run from the project folder:  source("R/02_priors.R")
# =============================================================================

source("R/priors.R")
library(posterior)
library(bayesplot)

options(brms.backend = "cmdstanr", mc.cores = 4)
dir.create("figures", showWarnings = FALSE)

# ---- sample from the prior only (model M2) ------------------------------------
m2_prior <- brm(f_m2, data = d, family = weibull(),
                prior = c(priors_main, prior_sd),
                sample_prior = "only",
                chains = 4, iter = 1000, seed = 2026, refresh = 0)

# ---- what do the priors imply? --------------------------------------------------
dr <- as_draws_df(m2_prior)
k  <- dr$shape

# typical (random effect = 0) untreated juvenile, xenon, average-risk eye
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

# ---- prior predictive check figure ----------------------------------------------
yrep <- posterior_predict(m2_prior, ndraws = 100)
yrep_cap <- pmin(yrep, 200)                     # cap extreme draws so the plot is readable

p1 <- ppc_km_overlay(d$futime, yrep_cap[1:50, ], status_y = d$status) +
  coord_cartesian(xlim = c(0, 80)) +
  labs(title = "Prior predictive: survival curves vs observed data",
       x = "months", y = "share of eyes without vision loss") +
  theme_project() + theme(legend.position = c(0.85, 0.85), legend.title = element_blank())

# median time on a log axis; show 1 week - 100 years and report what falls outside
mt <- median_time(eta_untr, k)
lo <- 0.25; hi <- 1200                          # months: 1 week, 100 years
share_out <- mean(mt < lo | mt > hi)
p2 <- ggplot(data.frame(m = pmin(pmax(mt, lo), hi)), aes(m)) +
  geom_histogram(bins = 50, fill = "#2a78d6", colour = "white", linewidth = 0.2) +
  scale_x_log10(breaks = c(0.25, 1, 12, 120, 1200),
                labels = c("1 week", "1 month", "1 year", "10 years", "100 years"),
                limits = c(lo * 0.8, hi * 1.25)) +
  labs(title = "Prior: median time to vision loss (untreated eye)",
       subtitle = sprintf("%.0f%% of draws fall outside this range (piled at the edges)",
                          100 * share_out),
       x = "time (log scale)", y = "draws") +
  theme_project()

p3 <- ggplot(data.frame(tr = exp(dr$b_trt)), aes(tr)) +
  geom_histogram(bins = 50, fill = "#eb6834", colour = "white", linewidth = 0.2) +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10)) +
  geom_vline(xintercept = 1, linetype = 2) +
  labs(title = "Prior: time ratio of laser treatment", x = "time ratio (log scale)", y = "draws") +
  theme_project()

fig <- bayesplot::bayesplot_grid(p1, p2, p3, grid_args = list(nrow = 1))
ggsave("figures/fig2_prior_predictive.png", fig, width = 16, height = 4.8, dpi = 150, bg = "white")
cat("\nSaved figures/fig2_prior_predictive.png\n")

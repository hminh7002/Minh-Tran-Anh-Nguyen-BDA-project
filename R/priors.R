# =============================================================================
# priors.R  -  data preparation, priors and plot theme shared by all scripts
# (sourced by 02_priors.R, 03_fit_models.R, ...; not meant to be run alone)
# =============================================================================

library(survival)
library(brms)
library(ggplot2)

# ---- data --------------------------------------------------------------------
data(retinopathy, package = "survival")
d <- within(retinopathy, {
  cens   <- 1 - status            # brms convention: 1 = right-censored, 0 = vision loss observed
  risk_c <- risk - mean(risk)     # centred risk score (mean 0)
  type   <- factor(type)          # levels: juvenile (reference), adult
  laser  <- factor(laser)         # levels: xenon (reference), argon
})

# ---- priors ------------------------------------------------------------------
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

# wider alternative priors for the sensitivity analysis (report section 9)
priors_wide <- c(
  prior(normal(4.1, 2),  class = Intercept),
  prior(normal(0, 2),    class = b),
  prior(gamma(1, 0.5),   class = shape)
)
prior_sd_wide <- prior(exponential(0.5), class = sd)

# ---- model formulas ------------------------------------------------------------
f_m1 <- bf(futime | cens(cens) ~ trt + type + laser + risk_c)                     # pooled
f_m2 <- bf(futime | cens(cens) ~ trt + type + laser + risk_c + (1 | id))          # + patient effect
f_m3 <- bf(futime | cens(cens) ~ trt * type + trt * laser + risk_c + (1 | id))   # + interactions

# ---- one plot style for every figure: white background, black text -----------
theme_project <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(text = element_text(colour = "black"),
          axis.text = element_text(colour = "black"),
          plot.title = element_text(colour = "black", face = "bold", size = base_size),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          strip.background = element_rect(fill = "grey92", colour = "grey60"))
}

# =============================================================================
# 01_smoke_test.R  -  check that the whole toolchain works end to end
#
# Loads the retinopathy data, fits the pooled Weibull model (M1) with brms +
# cmdstanr, and prints diagnostics. First run compiles the model (1-2 minutes);
# sampling then takes a few seconds.
# =============================================================================

library(survival)
library(brms)
library(posterior)
library(bayesplot)

options(brms.backend = "cmdstanr", mc.cores = 4)

# ---- data --------------------------------------------------------------------
data(retinopathy, package = "survival")
d <- within(retinopathy, {
  cens   <- 1 - status            # brms convention: 1 = right-censored, 0 = event
  risk_c <- risk - mean(risk)     # centred risk score
  type   <- factor(type)
  laser  <- factor(laser)
})
cat("Rows:", nrow(d), " patients:", length(unique(d$id)), "\n")
print(table(treated = d$trt, vision_loss = d$status))

# ---- M1: pooled Weibull model with right-censoring ----------------------------
priors <- c(prior(normal(4.1, 1.5), class = Intercept),  # log(60 months) ~ 4.1
            prior(normal(0, 1),     class = b),
            prior(gamma(2, 1),      class = shape))

m1 <- brm(futime | cens(cens) ~ trt + type + laser + risk_c,
          data = d, family = weibull(), prior = priors,
          chains = 4, iter = 2000, seed = 2026, refresh = 0)

print(summary(m1))

# ---- quick diagnostics ---------------------------------------------------------
np <- nuts_params(m1)
cat("Divergent transitions:", sum(subset(np, Parameter == "divergent__")$Value), "\n")
s <- summarise_draws(as_draws_df(m1), "rhat", "ess_bulk", "ess_tail")
print(s[startsWith(s$variable, "b_") | s$variable == "shape", ])
cat("Max Rhat:", round(max(s$rhat, na.rm = TRUE), 3), "(should be < 1.01)\n")

# time ratio for laser treatment: how many times longer until vision loss
tr <- exp(as_draws_df(m1)$b_trt)
cat(sprintf("Time ratio treated/untreated: median %.2f, 90%% interval [%.2f, %.2f]\n",
            median(tr), quantile(tr, 0.05), quantile(tr, 0.95)))

cat("\nSmoke test finished: R, brms, cmdstanr and CmdStan all work.\n")

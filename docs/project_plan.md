# Project plan

**Project:** Laser treatment and vision loss in diabetic retinopathy: a Bayesian survival analysis
**Authors:** Minh Tran, Anh Nguyen
**Course:** Bayesian Data Analysis, Aalto University, 2026

> Check the official project requirements and deadlines on the course page before submitting; this plan
> follows the usual BDA project requirements (at least two models, priors justified, convergence
> diagnostics, posterior predictive checks, LOO comparison, prior sensitivity, discussion, reflection).

---

## 1. Research questions

1. **Treatment effect.** How much does laser treatment delay vision loss? Reported as
   - the ratio of typical time to vision loss (treated vs untreated),
   - the probability of vision loss within 36 and 60 months for treated and untreated eyes,
   - the absolute risk reduction and the number of eyes needed to treat to prevent one case.
2. **Who benefits most?** Does the benefit differ between juvenile- and adult-onset diabetes, and
   between xenon and argon lasers?
3. **Between-eye dependence.** How strongly are the two eyes of one patient correlated (size of the
   patient-level variation), and does accounting for it change the conclusions?

## 2. Data

`survival::retinopathy`: 197 patients, 394 eyes (one treated, one control per patient).

First-look facts (from our exploration):

- Vision loss in 101 of 197 untreated eyes and 54 of 197 treated eyes.
- Within patients: only the untreated eye lost vision in 63 patients, only the treated eye in 16,
  both in 38, neither in 80.
- 97 argon / 100 xenon patients; 114 juvenile / 83 adult onset.
- Follow-up 0.3–75 months; about 61% of eyes are censored.
- Kaplan–Meier curves suggest a larger benefit in adult-onset patients (to be tested with M3).

Preprocessing: convert `type`, `laser` to factors; centre `risk` (and `age` if used);
`cens = 1 - status` (brms convention: 1 = right-censored).

## 3. How this project differs from existing work

The data set is a classic benchmark. Earlier analyses mostly use it to demonstrate statistical
methods for paired survival data and report a single treatment coefficient:

- Huster, Brookmeyer & Self (1989): original analysis of paired survival with covariates
  (treatment effect, diabetes type).
- Cox models with clustered standard errors (the `survival` package example).
- Shared frailty models (e.g. the `frailtySurv` package paper; treatment log-hazard coefficient
  around −0.9).
- Copula and bivariate parametric models (e.g. conditional copulas for onset age, bivariate power
  generalized Weibull, cure-frailty copula models).
- Some Bayesian analyses of new distribution families (e.g. Bayesian non-parametric survival
  estimation; bivariate defective Gompertz with Bayesian estimates).

**Our contribution** (to state explicitly in the report):

1. A complete **Bayesian workflow** applied to a clinical question: prior predictive checks,
   convergence diagnostics, posterior predictive checks for censored data, leave-one-patient-out
   cross-validation, and prior sensitivity analysis.
2. **Clinically interpretable outputs** computed from posterior draws: absolute risks at 3 and 5
   years, absolute risk reduction, number needed to treat, with full uncertainty.
3. **Effect heterogeneity**: posterior estimates and probabilities for different treatment benefit by
   diabetes type and laser type (interactions), rather than only adjusting for these variables.
4. **Explicit comparison of pooled vs hierarchical** modelling of the eye pairs, showing what is
   gained by modelling the dependence.

## 4. Models

Observation model for all models: Weibull with right-censoring.
In brms, the Weibull family uses a log link on the **mean** time; positive coefficients mean
**longer** time to vision loss (protective), and `exp(b)` is a time ratio.

| | brms formula | what it adds |
|---|---|---|
| **M1 pooled** | `futime \| cens(cens) ~ trt + type + laser + risk_c` | baseline, eyes independent |
| **M2 hierarchical** | `futime \| cens(cens) ~ trt + type + laser + risk_c + (1 \| id)` | patient-level random effect (frailty): eyes of the same patient correlated |
| **M3 interactions** | `futime \| cens(cens) ~ trt * type + trt * laser + risk_c + (1 \| id)` | treatment effect varies by diabetes and laser type |

Optional extensions if time allows:

- M2 with an **exponential** observation model (shape fixed to 1) to test whether the Weibull shape is needed.
- A simplified **logistic** model of vision loss (yes/no) with `(1 | id)` as a sanity check.

### Priors (weakly informative, to be checked with prior predictive simulation)

| parameter | prior | reasoning |
|---|---|---|
| Intercept (log mean months, average eye) | `normal(4.1, 1)` | log(60): mean time ~5 years; 95% between ~8 months and ~40 years |
| coefficients `b` | `normal(0, 0.75)` | time ratios mostly within 1/4.4 to 4.4 (95%); no direction assumed |
| Weibull shape | `gamma(4, 4)` | mean 1 (constant hazard); 90% within 0.34–1.94; looser priors gave absurd median times (< 1 day) |
| `sd(id)` | `exponential(2)` | patient variation on log-time scale, mean 0.5 |

Prior predictive check (`R/02_priors.R`, `figures/fig2_prior_predictive.png`): implied median time to vision
loss for an untreated eye has 90% interval ~2–270 months; implied P(vision loss within 5 years) ranges
from 0.09 to 1.00; the treatment time ratio is centred at 1 (no effect assumed). The simulated survival
curves surround the observed Kaplan–Meier curve without being implausible.

Alternative (wider) priors for sensitivity, also defined in `R/02_priors.R`: intercept `normal(4.1, 2)`,
`b ~ normal(0, 2)`, shape `gamma(1, 0.5)`, `sd ~ exponential(0.5)`; plus `priorsense` power-scaling.

## 5. Workflow and code sketch

```r
library(brms); library(survival); library(bayesplot); library(loo); library(priorsense)
data(retinopathy)
d <- within(retinopathy, {
  cens   <- 1 - status
  risk_c <- risk - mean(risk)
  type   <- factor(type); laser <- factor(laser)
})

priors <- c(prior(normal(4.1, 1), class = Intercept),     # log(60) ≈ 4.1  (see R/02_priors.R)
            prior(normal(0, 0.75), class = b),
            prior(gamma(4, 4), class = shape))

# prior predictive check
m2_prior <- brm(futime | cens(cens) ~ trt + type + laser + risk_c + (1 | id),
                data = d, family = weibull(),
                prior = c(priors, prior(exponential(2), class = sd)),
                sample_prior = "only", backend = "cmdstanr")

# fits
m1 <- brm(futime | cens(cens) ~ trt + type + laser + risk_c, data = d,
          family = weibull(), prior = priors, backend = "cmdstanr",
          chains = 4, iter = 2000, seed = 1)
m2 <- update(m1, formula. = ~ . + (1 | id), newdata = d,
             prior = c(priors, prior(exponential(2), class = sd)))
m3 <- update(m2, formula. = futime | cens(cens) ~ trt * type + trt * laser + risk_c + (1 | id),
             newdata = d)

# diagnostics: Rhat, ESS, divergences
summary(m2)                                   # Rhat, Bulk_ESS, Tail_ESS (+ divergence warnings)
np <- nuts_params(m2)
sum(subset(np, Parameter == "divergent__")$Value)   # number of divergent transitions

# posterior predictive check for censored data (Kaplan-Meier overlay)
yrep <- posterior_predict(m2, ndraws = 100)
ppc_km_overlay(d$futime, yrep, status_y = d$status)

# cross-validation: LOO and leave-one-patient-out
loo_compare(loo(m1), loo(m2), loo(m3))
kf2 <- kfold(m2, folds = "grouped", group = "id", K = 10)

# prior sensitivity
powerscale_sensitivity(m2)
```

Note on M1: plain PSIS-LOO treats eyes as independent units; because the two eyes of a patient are
correlated, we also report **leave-one-patient-out** (grouped K-fold) for the final comparison.

### Clinical quantities from posterior draws

For each draw: Weibull shape `k`, mean `mu = exp(eta)`, scale `lambda = mu / gamma(1 + 1/k)`,
survival `S(t) = exp(-(t / lambda)^k)`.

- `P(vision loss by t) = 1 - S(t)` for treated and untreated eyes at t = 36 and 60 months,
  for a typical patient (random effect 0) and for a new patient (random effect drawn from N(0, sd)).
- Absolute risk reduction ARR = P_untreated − P_treated; NNT = 1 / ARR (report median and interval).
- Time ratio `exp(b_trt)`; in M3 the time ratio for each diabetes type and laser type, and
  `P(benefit larger in adult-onset | data)`.

## 6. Checks we must report

- Convergence: R-hat < 1.01, bulk/tail ESS > 400, no divergences (else non-centred / adapt_delta).
- Posterior predictive checks: Kaplan–Meier overlay overall and by `trt` × `type`.
- Model comparison: elpd differences with SE (LOO and grouped K-fold), Pareto-k diagnostics.
- Prior sensitivity: `priorsense` table plus the refit with alternative priors.

## 7. Timeline (adjust to the official deadline)

| week | tasks | owner |
|---|---|---|
| 1 | Set up repo and R environment; exploratory plots; prior predictive checks; fit M1 | both / Minh |
| 2 | Fit M2, M3; diagnostics; posterior predictive checks | Anh |
| 2 | LOO + grouped K-fold comparison; prior sensitivity | Minh |
| 3 | Clinical quantities and figures; write report sections | both |
| 3 | Presentation slides, final checks, submit | both |

Suggested split: Minh — data, M1, model comparison, sensitivity; Anh — hierarchical models,
diagnostics, posterior predictive checks. Both write the discussion and the reflection.

## 8. Risks and fallbacks

| risk | fallback |
|---|---|
| divergences in M2/M3 | `adapt_delta = 0.95`; tighter `sd` prior; check non-identifiability |
| Weibull fits poorly (PPC) | log-normal family or piecewise-exponential hazard |
| high Pareto-k in LOO | grouped K-fold CV only |
| time runs short | drop M3's laser interaction; keep M1 vs M2 + `trt:type` |

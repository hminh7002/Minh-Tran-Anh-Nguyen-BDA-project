# Laser treatment and vision loss in diabetic retinopathy: a Bayesian survival analysis

Course project for **Bayesian Data Analysis (Aalto University, 2026)**.
Authors: **Minh Tran** and **Anh Nguyen**.

## Question

Does laser photocoagulation delay vision loss in patients with high-risk diabetic retinopathy, by how
much in absolute terms, and does the benefit differ between juvenile- and adult-onset diabetes or between
xenon and argon lasers?

We answer this with Bayesian Weibull survival models that account for right-censoring and for the fact
that both eyes of a patient are correlated (a patient-level random effect, also called a frailty).

## Data

The `retinopathy` data set from the R package `survival` (Therneau), originally from the Diabetic
Retinopathy Study as analysed by Huster, Brookmeyer & Self (1989), *Biometrics* 45:145–156.

- 197 patients with high-risk diabetic retinopathy, 394 eyes (one row per eye).
- In each patient one eye was randomised to laser treatment; the other eye is an untreated control.
- Outcome: months until visual acuity fell below 5/200 at two consecutive visits ("vision loss"),
  minus the minimum possible time of 6.5 months; right-censored by death, dropout or end of study.

| column | meaning |
|---|---|
| `id` | patient id (2 rows per patient) |
| `laser` | laser type used on the treated eye: `xenon` / `argon` |
| `eye` | which eye was treated: `left` / `right` |
| `age` | age at diagnosis of diabetes (years) |
| `type` | `juvenile` (diagnosed before 20) / `adult` |
| `trt` | 1 = treated eye, 0 = control eye |
| `futime` | follow-up time in months (time to vision loss or censoring) |
| `status` | 1 = vision loss, 0 = censored |
| `risk` | risk score of the eye (6–12) |

Get the data:

```r
library(survival)
data(retinopathy)
```

```python
import pandas as pd
url = "https://raw.githubusercontent.com/vincentarelbundock/Rdatasets/master/csv/survival/retinopathy.csv"
d = pd.read_csv(url).drop(columns="rownames")
```

## Models (summary)

| model | description |
|---|---|
| M1 pooled | Weibull model, `trt + type + laser + risk`, eyes treated as independent |
| M2 hierarchical | M1 + patient random intercept `(1 | id)` |
| M3 hierarchical + interactions | M2 + `trt:type` and `trt:laser` |

All models are fitted with `brms` (Stan backend `cmdstanr`) using weakly informative priors, and
compared with leave-one-patient-out cross-validation. Details: [`docs/project_plan.md`](docs/project_plan.md).

## Repository structure

```
.
├── README.md
├── docs/
│   ├── project_plan.md      # questions, models, priors, workflow, timeline, roles
│   └── report_outline.md    # section-by-section outline of the final report
├── R/                       # analysis scripts (to be added)
├── figures/                 # generated figures (to be added)
└── report/                  # report source and PDF (to be added)
```

## Software

R (≥ 4.3) with `brms`, `cmdstanr` (and CmdStan), `survival`, `posterior`, `bayesplot`, `loo`,
`priorsense`, `ggplot2`, `dplyr`.

```r
install.packages(c("brms", "survival", "posterior", "bayesplot", "loo", "priorsense", "ggplot2", "dplyr"))
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

## Status

- [x] Data set chosen and explored
- [x] Project plan and report outline
- [ ] Models M1–M3 fitted, diagnostics
- [ ] Posterior predictive checks, LOO comparison, prior sensitivity
- [ ] Report and presentation

## References

- Huster, W. J., Brookmeyer, R., & Self, S. G. (1989). Modelling paired survival data with covariates.
  *Biometrics*, 45, 145–156.
- Therneau, T. M. `survival`: Survival Analysis. R package.
- Bürkner, P.-C. (2017). brms: An R package for Bayesian multilevel models using Stan.
  *Journal of Statistical Software*, 80(1).
- Gelman, A., et al. (2013). *Bayesian Data Analysis*, 3rd ed. CRC Press.

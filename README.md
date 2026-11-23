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
├── R/
│   ├── 00_setup.R           # one-time environment setup
│   └── 01_smoke_test.R      # checks the toolchain by fitting model M1
├── figures/                 # generated figures (to be added)
└── report/                  # report source and PDF (to be added)
```

## Setting up R

Software: R (≥ 4.3) with `brms`, `cmdstanr` (+ CmdStan), `survival`, `posterior`, `bayesplot`, `loo`,
`priorsense`, `ggplot2`, `dplyr`.

1. **Install R** from <https://cran.r-project.org> (Windows: "Download R for Windows" → "base").
2. **Windows only: install Rtools** from <https://cran.r-project.org/bin/windows/Rtools/>, choosing the
   version that matches your R version (e.g. Rtools46 for R 4.6.x). Stan needs its C++ compiler.
3. **Install an editor**: RStudio or Positron (<https://posit.co/download/>), or VS Code with the
   "R" extension (then also run `install.packages("languageserver")`).
4. **Open the project folder** in the editor and run, from the project root:
   ```r
   source("R/00_setup.R")      # installs packages, checks the toolchain, installs CmdStan (5-15 min)
   source("R/01_smoke_test.R") # fits the pooled Weibull model; should end with "Smoke test finished"
   ```

### Running the scripts in VS Code (Windows)

1. Install the **R** extension (publisher: REditorSupport) and, in R, `install.packages("languageserver")`.
2. If R is not found, open Settings (`Ctrl+,`), search `r.rterm.windows`, and set it to your R path, e.g.
   `C:\Program Files\R\R-4.6.2\bin\R.exe` (check the version folder in `C:\Program Files\R`).
3. *File → Open Folder…* → select the project folder (the one containing `README.md`).
4. `Ctrl+Shift+P` → **R: Create R terminal**. An R console opens in the bottom panel, already in the project
   folder (check with `getwd()`).
5. Type the `source(...)` lines above at the `>` prompt and press Enter.
   Alternatively open a script and press `Ctrl+Shift+S` (*R: Source file*).

If installing a package fails with *"lib ... is not writable"*, create your personal package library once
(in any R console), then reload VS Code (`Ctrl+Shift+P` → *Developer: Reload Window*):
```r
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)
.libPaths(Sys.getenv("R_LIBS_USER")); install.packages("languageserver")
```

Without the extension, use a normal terminal (`` Ctrl+` ``) in the project folder:
`& "C:\Program Files\R\R-4.6.2\bin\Rscript.exe" R/00_setup.R` (adjust the version).

Expected smoke-test result: 0 divergences, max R-hat < 1.01, and a time ratio (treated vs untreated)
around 2.6 with a 90% interval of roughly 1.9–3.7.

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

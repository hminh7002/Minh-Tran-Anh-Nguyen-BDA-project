# Report outline

**Title:** Laser treatment and vision loss in diabetic retinopathy: a Bayesian survival analysis
**Authors:** Minh Tran, Anh Nguyen

Target length: about 15–20 pages including figures; code in an appendix or linked repository.
Every section below lists **what to write**, **figures/tables**, and **who drafts it**.

---

## 1. Introduction (≈ 0.5–1 page) — Minh
- Diabetic retinopathy as a leading cause of vision loss; laser photocoagulation as treatment.
- Research questions (treatment effect in absolute terms, heterogeneity by diabetes/laser type,
  between-eye dependence).
- One-paragraph summary of the approach and the **main result** (fill in last).
- Short paragraph on related work and **what is new here**: earlier analyses (Huster et al. 1989;
  Cox/frailty models; copula and Bayesian distribution papers) mostly demonstrate methods and report one
  treatment coefficient; we give a full Bayesian workflow with clinically interpretable quantities and
  effect heterogeneity.

## 2. Data and analysis problem (≈ 1–1.5 pages) — Minh
- Source, study design (one eye randomised, other control), 197 patients / 394 eyes.
- Variables table; definition of the outcome and the 6.5-month shift; censoring (≈61%).
- **Figure 1:** Kaplan–Meier curves by treatment; by treatment × diabetes type.
- **Table 1:** events by treatment, by laser type, by diabetes type; paired outcome table.
- Why the data are hierarchical (two eyes per patient) and why censoring must be modelled.

## 3. Models (≈ 2 pages) — Anh
- Weibull survival model with right-censoring: likelihood contributions `f(t)` for events and `S(t)`
  for censored eyes; meaning of shape `k` (increasing/decreasing hazard).
- M1 pooled; M2 hierarchical (patient random intercept = frailty); M3 with `trt:type`, `trt:laser`.
- Mathematical notation for each model, plus the brms formulas.
- Interpretation of coefficients as log time ratios.

## 4. Priors (≈ 1 page) — Anh
- Table of priors with justification in units of months / time ratios.
- **Figure 2:** prior predictive check (simulated Kaplan–Meier curves or time distributions vs plausible range).

## 5. Implementation and computation (≈ 1 page) — Minh
- brms + cmdstanr; chains, iterations, warm-up, seed; `adapt_delta` if changed.
- Stan/brms code (short in text, full in appendix).

## 6. Convergence diagnostics (≈ 0.5–1 page) — Anh
- **Table 2:** R-hat, bulk-ESS, tail-ESS for key parameters of each model; number of divergences,
  tree-depth warnings.
- Interpretation in one paragraph; actions taken if problems appeared.

## 7. Posterior predictive checks (≈ 1–1.5 pages) — Anh
- **Figure 3:** Kaplan–Meier overlay of observed data vs replicated data for each model.
- **Figure 4:** same, split by treatment × diabetes type (does M3 capture the heterogeneity?).
- Interpretation: where each model fits / misfits.

## 8. Model comparison (≈ 1 page) — Minh
- **Table 3:** elpd_loo, p_loo, Pareto-k summary; leave-one-patient-out (grouped K-fold) elpd;
  elpd differences with SE.
- Interpretation: is the hierarchical model better? Are interactions supported?

## 9. Prior sensitivity analysis (≈ 0.5–1 page) — Minh
- **Table 4:** `priorsense` power-scaling diagnostics.
- Refit with alternative priors; how much key estimates change.

## 10. Results (≈ 2 pages) — both
- **Table 5:** posterior summaries (median, 90% intervals) of time ratios and `sd(id)`.
- **Figure 5:** posterior probability of vision loss by 36 and 60 months, treated vs untreated,
  by diabetes type (point + interval).
- Absolute risk reduction and number needed to treat with intervals.
- `P(benefit larger in adult-onset | data)`, laser-type comparison.
- Plain-language answers to the three research questions.

## 11. Discussion (≈ 1 page) — both
- Limitations: 50% subsample of the original study; 1970s treatment technology; shifted time scale;
  censoring assumed non-informative; Weibull assumption; only patient-level covariates.
- Possible improvements: flexible baseline hazard (splines), joint model of both eyes with correlation
  on the hazard scale, time-varying effects, more covariates.

## 12. Conclusion (≈ 0.3 page) — both
- Two or three sentences answering the research questions with uncertainty.

## 13. Self-reflection (≈ 0.3 page) — both
- What we learned (e.g. censored likelihoods, hierarchical survival models, grouped CV, sensitivity).

## References
- Huster, Brookmeyer & Self (1989); Therneau `survival`; Bürkner (2017) brms; Gelman et al. BDA3;
  Vehtari, Gelman & Gabry (2017) LOO; Kallioinen et al. (2023) priorsense; related analyses cited in §1.

## Appendix
- Full brms/Stan code (or link to repository), session info, additional diagnostics plots.

---

### Figure and table checklist

| # | content | section | status |
|---|---|---|---|
| Fig 1 | KM by treatment and by treatment × type | 2 | draft done (exploration) |
| Fig 2 | prior predictive check | 4 | |
| Fig 3 | PPC KM overlay per model | 7 | |
| Fig 4 | PPC by treatment × type | 7 | |
| Fig 5 | risk of vision loss at 36/60 months | 10 | |
| Tab 1 | descriptive counts | 2 | draft done |
| Tab 2 | convergence diagnostics | 6 | |
| Tab 3 | LOO / grouped K-fold comparison | 8 | |
| Tab 4 | prior sensitivity | 9 | |
| Tab 5 | posterior summaries | 10 | |

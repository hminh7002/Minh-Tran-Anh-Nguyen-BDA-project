# =============================================================================
# 00_setup.R  -  one-time setup of the R environment for the BDA project
#
# Run this ONCE on each computer (Windows, macOS or Linux), from the project
# folder, e.g. in RStudio: open the project folder, then  source("R/00_setup.R")
#
# Before running on Windows: install R (>= 4.3) and the matching Rtools
# (see README "Setting up R"). Rtools provides the C++ compiler that Stan needs.
# =============================================================================

cat("R version:", R.version.string, "\n")
# use a fixed CRAN mirror so installation never stops to ask (also works with Rscript)
if (is.null(getOption("repos")) || getOption("repos")[["CRAN"]] %in% c("@CRAN@", "")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
# install packages into a personal library (C:/Program Files/R is not writable without admin rights)
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))
cat("Installing packages into:", .libPaths()[1], "\n")
if (getRversion() < "4.3.0") stop("Please install R 4.3 or newer from https://cran.r-project.org")

# ---- 1. CRAN packages --------------------------------------------------------
cran_pkgs <- c(
  "brms",        # Bayesian regression models (writes Stan code for us)
  "survival",    # the retinopathy data + Kaplan-Meier
  "posterior",   # summaries and convergence diagnostics (Rhat, ESS)
  "bayesplot",   # MCMC and posterior predictive check plots
  "loo",         # leave-one-out cross-validation
  "priorsense",  # prior sensitivity (power-scaling)
  "ggplot2", "dplyr", "tidyr",
  "knitr", "rmarkdown"   # report
)
missing <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing) > 0) {
  cat("Installing from CRAN:", paste(missing, collapse = ", "), "\n")
  install.packages(missing)
} else {
  cat("All CRAN packages already installed.\n")
}

# ---- 2. cmdstanr (not on CRAN; from the Stan r-universe) ---------------------
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages("cmdstanr",
                   repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
}
library(cmdstanr)

# ---- 3. C++ toolchain + CmdStan ----------------------------------------------
# check_cmdstan_toolchain(fix = TRUE) also fixes common Rtools path problems on Windows.
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)

if (is.null(tryCatch(cmdstan_path(), error = function(e) NULL))) {
  cat("Installing CmdStan (downloads and compiles; takes 5-15 minutes)...\n")
  install_cmdstan(cores = max(1, parallel::detectCores() - 1))
}
cat("CmdStan version:", cmdstan_version(), "at", cmdstan_path(), "\n")

# ---- 4. Summary --------------------------------------------------------------
ok <- sapply(c(cran_pkgs, "cmdstanr"), requireNamespace, quietly = TRUE)
print(ok)
if (all(ok)) {
  cat("\nSetup complete. Next: source('R/01_smoke_test.R') to check that a model compiles and samples.\n")
} else {
  cat("\nSome packages failed to install:", paste(names(ok)[!ok], collapse = ", "), "\n")
}

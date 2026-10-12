# =============================================================================
# config.R
# Every tunable setting for the project lives here, so a reader can reproduce
# (or deliberately change) the whole study by editing one file.
# =============================================================================

# Reproducibility ------------------------------------------------------------
# The course requires the seed to be the last four digits of the student ID.
STUDENT_ID_LAST4 <- 6877L   # last four digits of student ID
SEED <- STUDENT_ID_LAST4

# Monte Carlo design ---------------------------------------------------------
# REPS  = number of simulated samples (how many sample means we compute).
# N_GRID = sample sizes n (how many observations go into ONE sample mean).
# These are different quantities: 10,000 reps of n = 30 means 10,000 sample
# means, each built from 30 observations. The sample size is still 30.
REPS   <- 10000L
N_GRID <- c(2, 5, 10, 20, 30, 50, 100, 250, 500, 1000, 2000)

# Sample sizes shown in the histogram / Q-Q panels (a subset of N_GRID).
N_SHOW <- c(2, 5, 10, 30, 100, 1000)

# "Approximately normal" criterion ---------------------------------------------
# A sampling distribution counts as approximately normal at a given n only if
# ALL four checks pass. Thresholds are chosen to be several Monte Carlo
# standard errors wide, so that a pass/fail is not caused by simulation noise
# (with 10,000 reps: SE(skewness) ~ 0.025, SE(excess kurtosis) ~ 0.05,
# SE(tail rate) ~ 0.002, relative SE(sd) ~ 0.7%).
CRITERIA <- list(
  skew     = 0.20,          # |skewness| <= 0.20
  exkurt   = 0.30,          # |excess kurtosis| <= 0.30
  se_ratio = c(0.97, 1.03), # simulated SD / correct theoretical SE
  tail     = c(0.04, 0.06)  # share of means beyond mu +/- 1.96 SE (normal: 0.05)
)

# Seeds used for the robustness check (does the conclusion depend on the seed?)
ROBUSTNESS_SEEDS <- SEED + c(0, 1000, 2000, 3000, 4000)

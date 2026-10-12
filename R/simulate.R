# =============================================================================
# simulate.R
# The Monte Carlo engine and the normality diagnostics.
# =============================================================================

# Moment-based shape statistics (population formulas, as in the class notes) --
skewness <- function(x) {
  d <- x - mean(x)
  mean(d^3) / mean(d^2)^(3 / 2)
}

excess_kurtosis <- function(x) {
  d <- x - mean(x)
  mean(d^4) / mean(d^2)^2 - 3
}

# Shapiro-Wilk accepts at most 5,000 values, so test a random subsample.
shapiro_subsample <- function(x, size = 5000) {
  x <- sample(x, min(size, length(x)))
  if (length(unique(x)) < 3) return(c(W = NA_real_, p = NA_real_))
  r <- shapiro.test(x)
  c(W = unname(r$statistic), p = r$p.value)
}

# Simulate the sampling distribution of the mean for every n -------------------
# Returns a named list: one vector of `reps` sample means per sample size.
simulate_means <- function(pop, n_grid = N_GRID, reps = REPS, seed = SEED) {
  set.seed(seed)
  out <- lapply(n_grid, function(n) pop$means(n, reps))
  names(out) <- paste0("n", n_grid)
  out
}

# Diagnose one vector of sample means -----------------------------------------
diagnose_means <- function(means, pop, n, criteria = CRITERIA) {
  se_true  <- pop$se(n)
  se_naive <- pop$se_naive(n)
  z <- (means - pop$mu) / se_true
  sw <- shapiro_subsample(means)

  sk <- skewness(means)
  ku <- excess_kurtosis(means)
  ratio <- sd(means) / se_true
  tail  <- mean(abs(z) > 1.96)

  checks <- c(
    skew     = abs(sk) <= criteria$skew,
    exkurt   = abs(ku) <= criteria$exkurt,
    se_ratio = isTRUE(ratio >= criteria$se_ratio[1] & ratio <= criteria$se_ratio[2]),
    tail     = isTRUE(tail >= criteria$tail[1] & tail <= criteria$tail[2])
  )
  checks[is.na(checks)] <- FALSE

  data.frame(
    population     = pop$label,
    n              = n,
    sim_mean       = mean(means),
    theo_mean      = pop$mu,
    sim_sd         = sd(means),
    theo_se        = se_true,
    naive_se       = se_naive,
    se_ratio       = ratio,
    naive_ratio    = sd(means) / se_naive,
    skewness       = sk,
    pred_skewness  = pop$skew_n(n),
    exkurt         = ku,
    pred_exkurt    = pop$exkurt_n(n),
    tail_rate      = tail,
    sw_W           = sw[["W"]],
    sw_p           = sw[["p"]],
    iqr            = IQR(means),
    pass_skew      = checks[["skew"]],
    pass_exkurt    = checks[["exkurt"]],
    pass_se        = checks[["se_ratio"]],
    pass_tail      = checks[["tail"]],
    approx_normal  = all(checks),
    row.names = NULL
  )
}

diagnose_all <- function(means_list, pop, n_grid = N_GRID) {
  do.call(rbind, Map(function(m, n) diagnose_means(m, pop, n), means_list, n_grid))
}

# Smallest n from which the criterion passes and KEEPS passing ----------------
# (a single lucky pass followed by failures does not count).
first_stable_n <- function(diag) {
  ok <- diag$approx_normal
  if (!ok[length(ok)]) return(NA_real_)
  last_fail <- max(c(0, which(!ok)))
  diag$n[last_fail + 1]
}

# Theory-based prediction of the n needed for the skewness / kurtosis checks --
# For i.i.d. data: skew(X-bar) = skew / sqrt(n), exkurt(X-bar) = exkurt / n.
# For the fatigue model the exact Poisson-binomial formulas are used instead.
# Returns NA when no theory is available (Cauchy, dependent process).
predict_n <- function(pop, criteria = CRITERIA, n_max = 5000) {
  ns <- 1:n_max
  sk <- vapply(ns, pop$skew_n, numeric(1))
  ku <- vapply(ns, pop$exkurt_n, numeric(1))
  if (all(is.na(sk))) return(NA_real_)
  ok <- abs(sk) <= criteria$skew & abs(ku) <= criteria$exkurt
  ok[is.na(ok)] <- FALSE
  if (!ok[n_max]) return(NA_real_)
  ns[max(c(0, which(!ok))) + 1]
}

# Run everything for one population ------------------------------------------
run_population <- function(pop, seed_offset = 0, n_grid = N_GRID, reps = REPS) {
  means <- simulate_means(pop, n_grid, reps, seed = SEED + seed_offset)
  diag  <- diagnose_all(means, pop, n_grid)
  list(pop = pop, means = means, diag = diag,
       stable_n = first_stable_n(diag), predicted_n = predict_n(pop))
}

# =============================================================================
# populations.R
# Each population is a list with the same fields, so one simulation engine
# can treat all of them alike:
#
#   id, label        short key and display name
#   kind             "continuous" or "discrete"
#   assumptions      which CLT conditions hold (independent / identical / finite var)
#   draw(m)          m individual observations (used to plot the population)
#   means(n, reps)   `reps` sample means, each from a sample of size n
#   mu               mean of the sample mean (the CLT centre)
#   sigma            population SD (NA if it does not exist)
#   se(n)            the CORRECT theoretical SD of the sample mean
#   se_naive(n)      the textbook sigma / sqrt(n) (differs from se() only when
#                    the i.i.d. assumption fails)
#   skew, exkurt     population skewness and excess kurtosis (NA if undefined)
#   skew_n(n), exkurt_n(n)  predicted skewness / kurtosis of the sample mean
# =============================================================================

# Sample means for i.i.d. data ------------------------------------------------
# Fill a (reps x n) matrix where each ROW is one sample, then take rowMeans.
# Large n is processed in blocks of rows to keep memory use modest.
iid_means <- function(draw, n, reps, max_cells = 2e6) {
  block <- max(1L, floor(max_cells / n))
  starts <- seq(1L, reps, by = block)
  unlist(lapply(starts, function(s) {
    r <- min(block, reps - s + 1L)
    rowMeans(matrix(draw(n * r), nrow = r, ncol = n))
  }))
}

# Helper: build an i.i.d. population from a random-number generator ----------
# For i.i.d. data with finite moments, the sample mean of n observations has
#   SE = sigma / sqrt(n),  skewness = skew / sqrt(n),  excess kurtosis = exkurt / n.
# These formulas give a PREDICTION of how quickly each population converges.
make_iid_population <- function(id, label, kind, draw, mu, sigma, skew, exkurt,
                                notes = "") {
  finite_var <- is.finite(sigma)
  list(
    id = id, label = label, kind = kind, notes = notes,
    assumptions = c(independent = TRUE, identical = TRUE, finite_variance = finite_var),
    draw  = draw,
    means = function(n, reps) iid_means(draw, n, reps),
    mu = mu, sigma = sigma, skew = skew, exkurt = exkurt,
    se       = function(n) sigma / sqrt(n),
    se_naive = function(n) sigma / sqrt(n),
    skew_n   = function(n) skew / sqrt(n),
    exkurt_n = function(n) exkurt / n
  )
}

# 1. Standard normal ---------------------------------------------------------
pop_normal <- make_iid_population(
  id = "normal", label = "Standard Normal", kind = "continuous",
  draw = function(m) rnorm(m, 0, 1),
  mu = 0, sigma = 1, skew = 0, exkurt = 0,
  notes = "Control case: the mean of normals is exactly normal for every n."
)

# 2. Fair six-sided die (discrete uniform on 1..6) -------------------------
die_values <- 1:6
pop_die <- make_iid_population(
  id = "die", label = "Six-Sided Die", kind = "discrete",
  draw = function(m) sample.int(6L, m, replace = TRUE),
  mu = 3.5,
  sigma = sqrt(35 / 12),                      # (k^2 - 1) / 12 with k = 6
  skew = 0,
  exkurt = -6 * (6^2 + 1) / (5 * (6^2 - 1)),  # discrete uniform: about -1.269
  notes = "Symmetric, bounded, flat (light tails)."
)

# 3. Exponential(rate = 1) -------------------------------------------------
pop_exp <- make_iid_population(
  id = "exponential", label = "Exponential (rate = 1)", kind = "continuous",
  draw = function(m) rexp(m, rate = 1),
  mu = 1, sigma = 1, skew = 2, exkurt = 6,
  notes = "Strong right skew; the classic slow-converging example."
)

# 4. Binomial(size = 5, p = 0.10) ------------------------------------------
b_size <- 5; b_p <- 0.10
b_var  <- b_size * b_p * (1 - b_p)
pop_binom <- make_iid_population(
  id = "binomial", label = "Binomial (5, 0.10)", kind = "discrete",
  draw = function(m) rbinom(m, size = b_size, prob = b_p),
  mu = b_size * b_p,
  sigma = sqrt(b_var),
  skew = (1 - 2 * b_p) / sqrt(b_var),
  exkurt = (1 - 6 * b_p * (1 - b_p)) / b_var,
  notes = "Highly discrete: 59% of draws are 0; strong right skew."
)

# 5. Cauchy(0, 1) -----------------------------------------------------------
# Mean and variance do not exist, so every moment-based formula is undefined.
pop_cauchy <- make_iid_population(
  id = "cauchy", label = "Cauchy (0, 1)", kind = "continuous",
  draw = function(m) rcauchy(m, location = 0, scale = 1),
  mu = 0,           # location / median, NOT a mean
  sigma = NA_real_, skew = NA_real_, exkurt = NA_real_,
  notes = "No finite mean or variance: the CLT's assumptions fail."
)

# =============================================================================
# Extra credit populations
# =============================================================================

# 6. Dependent machine-failure process --------------------------------------
#   lambda_t = lambda0 * exp(beta * (X_{t-1} - mu_ref)),   X_t | X_{t-1} ~ Exp(lambda_t)
# Simulation is vectorised ACROSS replications: `x` holds the current
# operating time for all `reps` independent machines, and we step every
# machine forward one failure at a time. A burn-in lets each machine forget
# its starting value so that samples come from the long-run (stationary)
# behaviour of the process rather than from its start-up period.
dep_params <- list(lambda0 = 1, beta = 1.25, mu_ref = 1, burn_in = 50)

dependent_step <- function(x, p = dep_params) {
  rexp(length(x), rate = p$lambda0 * exp(p$beta * (x - p$mu_ref)))
}

dependent_means <- function(n, reps, p = dep_params) {
  x <- rexp(reps, rate = p$lambda0)                # first operating time
  for (b in seq_len(p$burn_in)) x <- dependent_step(x, p)
  total <- numeric(reps)
  for (t in seq_len(n)) {
    x <- dependent_step(x, p)
    total <- total + x
  }
  total / n
}

dependent_chain <- function(len, p = dep_params) {
  x <- numeric(len + p$burn_in)
  x[1] <- rexp(1, p$lambda0)
  for (t in 2:length(x)) x[t] <- dependent_step(x[t - 1], p)
  x[-seq_len(p$burn_in)]
}

# The stationary moments of this process have no simple closed form, so we
# estimate them ONCE from a very long chain. Because the observations are
# correlated, Var(mean) is NOT sigma^2 / n. For a stationary series,
#   n * Var(X-bar) -> sigma^2 * (1 + 2 * sum_k rho_k)  (the long-run variance).
estimate_dependent_moments <- function(len = 1e6, max_lag = 30, seed = 1) {
  set.seed(seed)
  chain <- dependent_chain(len)
  rho <- acf(chain, lag.max = max_lag, plot = FALSE)$acf[-1]
  s   <- sd(chain)
  z   <- (chain - mean(chain)) / s
  list(
    chain_mean = mean(chain), chain_sd = s,
    skew = mean(z^3), exkurt = mean(z^4) - 3,
    rho = rho,
    lrv_factor = 1 + 2 * sum(rho),          # long-run variance / sigma^2
    chain_head = chain[1:200]
  )
}

# Exact variance of the mean of n consecutive values of a stationary series:
#   Var(X-bar) = (sigma^2 / n) * [1 + 2 * sum_{k=1}^{n-1} (1 - k/n) * rho_k]
# Autocorrelations beyond the estimated lags are treated as zero (they decay
# geometrically and are below 0.01 by lag ~20).
dependent_se <- function(n, mom) {
  k <- seq_len(min(n - 1, length(mom$rho)))
  factor <- 1 + 2 * sum((1 - k / n) * mom$rho[k])
  mom$chain_sd * sqrt(factor / n)
}

make_dependent_population <- function(mom) {
  list(
    id = "dependent", label = "Dependent Machine Failure", kind = "continuous",
    notes = "Exponential conditional on the previous failure; negative dependence.",
    assumptions = c(independent = FALSE, identical = TRUE, finite_variance = TRUE),
    draw  = function(m) dependent_chain(m),
    means = dependent_means,
    mu = mom$chain_mean, sigma = mom$chain_sd,
    skew = mom$skew, exkurt = mom$exkurt,
    se       = function(n) dependent_se(n, mom),    # accounts for autocorrelation
    se_naive = function(n) mom$chain_sd / sqrt(n),  # pretends data are i.i.d.
    skew_n   = function(n) NA_real_,
    exkurt_n = function(n) NA_real_,
    moments = mom
  )
}

# 7. Non-identical Bernoulli (basketball fatigue) ---------------------------
# A sample of size n is n shots whose make probability falls linearly from
# 0.90 (first shot) to 0.50 (last shot). Shots are independent but NOT
# identically distributed. Average p is exactly 0.70 for every n.
fatigue_p <- function(n, p_start = 0.90, p_end = 0.50) {
  seq(p_start, p_end, length.out = n)
}

fatigue_means <- function(n, reps) {
  p <- fatigue_p(n)
  # matrix fills column-wise, so column j holds shot j for every replication
  made <- matrix(rbinom(n * reps, size = 1, prob = rep(p, each = reps)),
                 nrow = reps, ncol = n)
  rowMeans(made)
}

# Exact moments of a Poisson-binomial proportion (sum of independent Bernoullis)
fatigue_moments <- function(n) {
  p <- fatigue_p(n); v <- p * (1 - p)
  list(
    mean = mean(p),
    se = sqrt(sum(v)) / n,
    skew = sum(v * (1 - 2 * p)) / sum(v)^1.5,
    exkurt = sum(v * (1 - 6 * v)) / sum(v)^2
  )
}

pop_fatigue <- list(
  id = "fatigue", label = "Non-identical Bernoulli (Fatigue)", kind = "discrete",
  notes = "Independent shots, but p falls from 0.90 to 0.50.",
  assumptions = c(independent = TRUE, identical = FALSE, finite_variance = TRUE),
  draw  = function(m) rbinom(m, 1, fatigue_p(m)),
  means = fatigue_means,
  mu = 0.70, sigma = sqrt(0.7 * 0.3),
  skew = NA_real_, exkurt = NA_real_,  # no single population: p changes per shot
  se       = function(n) fatigue_moments(n)$se,
  se_naive = function(n) sqrt(0.7 * 0.3 / n),       # pretends every p = 0.70
  skew_n   = function(n) fatigue_moments(n)$skew,
  exkurt_n = function(n) fatigue_moments(n)$exkurt
)

# 8. Custom age-at-death distribution (two-component mixture) ---------------
age <- 0:100
age_weight <- 0.045 * exp(-age / 5) +                 # early-death component
              0.030 * exp(-0.5 * ((age - 65) / 10)^2)  # adult-death component
age_prob <- age_weight / sum(age_weight)

age_mu     <- sum(age * age_prob)
age_sigma  <- sqrt(sum((age - age_mu)^2 * age_prob))
age_skew   <- sum(((age - age_mu) / age_sigma)^3 * age_prob)
age_exkurt <- sum(((age - age_mu) / age_sigma)^4 * age_prob) - 3

rdeath <- function(m) sample(age, size = m, replace = TRUE, prob = age_prob)

pop_age <- make_iid_population(
  id = "age", label = "Age at Death (Mixture)", kind = "discrete",
  draw = rdeath,
  mu = age_mu, sigma = age_sigma, skew = age_skew, exkurt = age_exkurt,
  notes = "Bimodal: infant-mortality spike plus a bell centred at 65."
)

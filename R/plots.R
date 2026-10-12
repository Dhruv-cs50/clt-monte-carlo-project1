# =============================================================================
# plots.R
# Base-R plotting helpers, so every population gets identical, comparable
# figures.
# =============================================================================

PAL <- list(
  pop    = "#E07B39",  # population (individual values)
  bars   = "#9CB8D6",  # sampling-distribution histograms
  normal = "#B22222",  # CLT-predicted normal curve
  pass   = "#2E7D5B",
  fail   = "#B22222",
  grid   = "gray88"
)

# Population: what do the individual observations look like? ---------------
plot_population <- function(pop, m = 10000, xlim = NULL, main = NULL, ...) {
  x <- pop$draw(m)
  main <- if (is.null(main)) paste0(pop$label, ": ", format(m, big.mark = ","),
                                    " individual values") else main
  if (pop$kind == "discrete") {
    tab <- table(factor(x, levels = sort(unique(x))))
    plot(as.numeric(names(tab)), as.numeric(tab) / m, type = "h", lwd = 6,
         lend = 1, col = PAL$pop, xlab = "Individual value",
         ylab = "Proportion", main = main, xlim = xlim, ...)
  } else {
    shown <- if (is.null(xlim)) x else x[x >= xlim[1] & x <= xlim[2]]
    if (is.null(xlim)) xlim <- range(shown)
    hist(shown, breaks = 60, probability = TRUE, col = PAL$pop, border = "white",
         xlab = "Individual value", main = main, xlim = xlim, ...)
  }
  abline(v = pop$mu, lty = 2, lwd = 2)
  invisible(x)
}

# Sampling distributions with the CLT normal curve, one panel per n -------------
plot_sampling_grid <- function(res, n_show = N_SHOW, common_x = TRUE,
                               xlim = NULL, ncol = 3) {
  pop <- res$pop
  nr <- ceiling(length(n_show) / ncol)
  old <- par(mfrow = c(nr, ncol), mar = c(4, 4, 3, 1))
  on.exit(par(old))
  if (is.null(xlim) && common_x) {
    xlim <- range(unlist(res$means[paste0("n", n_show)]), na.rm = TRUE)
  }
  for (n in n_show) {
    m <- res$means[[paste0("n", n)]]
    xl <- if (common_x || !is.null(xlim)) xlim else range(m)
    shown <- m[m >= xl[1] & m <= xl[2]]
    hist(shown, breaks = 50, probability = TRUE, col = PAL$bars,
         border = "white", xlim = xl, main = paste("n =", n),
         xlab = "Sample mean")
    se <- pop$se(n)
    if (is.finite(se)) {
      curve(dnorm(x, mean = pop$mu, sd = se), add = TRUE,
            col = PAL$normal, lwd = 2, n = 400)
    }
  }
}

# Normal Q-Q plots, one panel per n ------------------------------------------
plot_qq_grid <- function(res, n_show = N_SHOW, ncol = 3) {
  nr <- ceiling(length(n_show) / ncol)
  old <- par(mfrow = c(nr, ncol), mar = c(4, 4, 3, 1))
  on.exit(par(old))
  for (n in n_show) {
    m <- res$means[[paste0("n", n)]]
    qqnorm(m, main = paste("Q-Q plot: n =", n), pch = 16, cex = 0.3,
           col = adjustcolor("#2B5C8A", 0.5))
    qqline(m, col = PAL$normal, lwd = 2)
  }
}

# Skewness and excess kurtosis against n, with theory and thresholds --------
plot_shape_trend <- function(res, criteria = CRITERIA) {
  d <- res$diag; pop <- res$pop
  old <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  on.exit(par(old))
  ns_fine <- exp(seq(log(min(d$n)), log(max(d$n)), length.out = 200))
  for (what in c("skewness", "exkurt")) {
    y  <- d[[what]]
    th <- if (what == "skewness") criteria$skew else criteria$exkurt
    pred <- vapply(ns_fine, if (what == "skewness") pop$skew_n else pop$exkurt_n,
                   numeric(1))
    yl <- range(c(y, pred, th, -th), na.rm = TRUE)
    plot(d$n, y, log = "x", pch = 19, ylim = yl,
         col = ifelse(abs(y) <= th, PAL$pass, PAL$fail),
         xlab = "Sample size n (log scale)",
         ylab = if (what == "skewness") "Skewness of sample means" else "Excess kurtosis of sample means",
         main = if (what == "skewness") "Skewness vs n" else "Excess kurtosis vs n")
    rect(1e-3, -th, 1e6, th, col = adjustcolor(PAL$pass, 0.08), border = NA)
    abline(h = 0, col = "gray50")
    abline(h = c(-th, th), lty = 3, col = PAL$pass)
    if (!all(is.na(pred))) lines(ns_fine, pred, lwd = 2, col = "gray30", lty = 2)
    points(d$n, y, pch = 19, col = ifelse(abs(y) <= th, PAL$pass, PAL$fail))
    legend("topright", bty = "n", cex = 0.8,
           legend = c("Simulated", if (!all(is.na(pred))) "Theory", "Threshold band"),
           pch = c(19, if (!all(is.na(pred))) NA, 15),
           lty = c(NA, if (!all(is.na(pred))) 2, NA),
           col = c("gray20", if (!all(is.na(pred))) "gray30", adjustcolor(PAL$pass, 0.3)))
  }
}

# Simulated SD vs theoretical SE ------------------------------------------------
plot_se_check <- function(res, show_naive = FALSE) {
  d <- res$diag
  yl <- range(c(d$sim_sd, d$theo_se, if (show_naive) d$naive_se), na.rm = TRUE)
  plot(d$n, d$sim_sd, log = "xy", pch = 19, ylim = yl,
       xlab = "Sample size n (log scale)", ylab = "SD of sample means (log scale)",
       main = paste0(res$pop$label, ": spread shrinks like 1/sqrt(n)"))
  lines(d$n, d$theo_se, col = PAL$normal, lwd = 2)
  leg <- c("Simulated SD", "Correct theoretical SE")
  cols <- c("black", PAL$normal); lt <- c(NA, 1); pc <- c(19, NA)
  if (show_naive) {
    lines(d$n, d$naive_se, col = "gray40", lwd = 2, lty = 2)
    leg <- c(leg, "Naive i.i.d. sigma/sqrt(n)")
    cols <- c(cols, "gray40"); lt <- c(lt, 2); pc <- c(pc, NA)
  }
  legend("topright", legend = leg, col = cols, lty = lt, pch = pc, bty = "n")
}

# All populations together: |skewness| of the sample mean vs n ----------------
plot_convergence_all <- function(results, criteria = CRITERIA) {
  cols <- c("#1B6CA8", "#5AA469", "#E07B39", "#8E5EA2", "#B22222",
            "#C9A227", "#2CA6A4", "#6B4F3A")
  finite <- Filter(function(r) r$pop$id != "cauchy", results)
  ys <- unlist(lapply(finite, function(r) abs(r$diag$skewness)))
  plot(NULL, log = "xy", xlim = range(N_GRID), ylim = c(min(ys[ys > 0]) * 0.8, max(ys) * 1.2),
       xlab = "Sample size n (log scale)", ylab = "|Skewness| of sample means (log scale)",
       main = "How fast does each population lose its skew?")
  abline(h = criteria$skew, lty = 2, col = "gray30")
  text(max(N_GRID), criteria$skew, paste("threshold", criteria$skew), pos = 3, cex = 0.75, offset = 0.3)
  abline(v = 30, lty = 3, col = "gray60")
  text(30, max(ys), "n = 30", pos = 4, cex = 0.75, col = "gray40")
  for (i in seq_along(finite)) {
    d <- finite[[i]]$diag
    y <- pmax(abs(d$skewness), 1e-3)
    lines(d$n, y, col = cols[i], lwd = 2)
    points(d$n, y, col = cols[i], pch = 19, cex = 0.7)
  }
  legend("bottomleft", bg = "white", box.col = "gray85", legend = vapply(finite, function(r) r$pop$label, ""),
         col = cols[seq_along(finite)], lwd = 2, bty = "n", cex = 0.75)
}

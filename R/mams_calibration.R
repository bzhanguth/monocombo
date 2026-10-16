# Internal helpers for MAMS stage validation and finite-sample calibration.

.normalize_mams_t <- function(j = NULL, t = NULL) {
  if (is.null(j)) {
    j <- if (is.null(t)) 3L else length(t)
  }
  j <- as.integer(j)
  if (length(j) != 1L || is.na(j) || !j %in% c(2L, 3L)) {
    stop("For MAMS designs, `j` must be 2 or 3.", call. = FALSE)
  }

  if (is.null(t)) {
    t <- seq(1 / j, 1, length.out = j)
  }
  t <- as.numeric(t)
  if (length(t) != j || anyNA(t) || any(!is.finite(t))) {
    stop("For MAMS designs, `t` must contain exactly `j` finite cumulative stage fractions.", call. = FALSE)
  }
  if (any(t <= 0) || any(t > 1) || any(diff(t) <= 0) || abs(t[j] - 1) > 1e-8) {
    stop("MAMS stage fractions must be strictly increasing in (0, 1] and end at 1.", call. = FALSE)
  }

  list(j = j, t = t)
}

.with_preserved_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

.mams_analytical_fixed <- function(...) {
  .with_preserved_seed(104729L, mc_pro_analytical(...))
}

.proportion_z_vector <- function(success_treatment, success_control, n_treatment, n_control) {
  p_treatment <- success_treatment / n_treatment
  p_control <- success_control / n_control
  se <- sqrt(
    p_treatment * (1 - p_treatment) / n_treatment +
      p_control * (1 - p_control) / n_control
  )
  z <- (p_treatment - p_control) / se
  fallback <- !is.finite(z) | !is.finite(se) | se <= 0
  if (any(fallback)) {
    n_treatment_adjusted <- n_treatment + 2
    n_control_adjusted <- n_control + 2
    p_treatment_adjusted <- (success_treatment + 1) / n_treatment_adjusted
    p_control_adjusted <- (success_control + 1) / n_control_adjusted
    se_adjusted <- sqrt(
      p_treatment_adjusted * (1 - p_treatment_adjusted) / n_treatment_adjusted +
        p_control_adjusted * (1 - p_control_adjusted) / n_control_adjusted
    )
    z_adjusted <- (p_treatment_adjusted - p_control_adjusted) / se_adjusted
    z[fallback] <- z_adjusted[fallback]
  }
  z[!is.finite(z)] <- 0
  z
}

.make_mams_cache <- function(n, j, t, p0, p1, p2, nsim, seed) {
  stage <- .normalize_mams_t(j, t)
  j <- stage$j
  t <- stage$t

  if (length(n) != 1L || !is.finite(n) || n <= 0 || abs(n - round(n)) > 1e-8) {
    stop("`n` must be a positive integer for finite-sample MAMS calibration.", call. = FALSE)
  }
  n <- as.integer(round(n))
  nsim <- as.integer(nsim)
  if (length(nsim) != 1L || is.na(nsim) || nsim < 1000L) {
    stop("`calibration_nsim` must be an integer of at least 1000.", call. = FALSE)
  }
  if (length(p0) != 1L || !is.finite(p0) || p0 <= 0 || p0 >= 1) {
    stop("`p0` must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (length(p1) != 1L || length(p2) != 1L ||
      !is.finite(p1) || !is.finite(p2) ||
      p1 < 0 || p1 > 1 || p2 < 0 || p2 > 1) {
    stop("`p1` and `p2` must lie between 0 and 1.", call. = FALSE)
  }
  if (length(seed) != 1L || !is.finite(seed)) {
    stop("The calibration seed must be a finite scalar.", call. = FALSE)
  }

  n_cumulative <- as.integer(round(n * t))
  n_increment <- c(n_cumulative[1L], diff(n_cumulative))
  if (any(n_increment <= 0L)) {
    stop("The selected sample size produces an empty MAMS stage.", call. = FALSE)
  }

  draws <- .with_preserved_seed(seed, {
    list(
      control = matrix(
        rbinom(nsim * j, rep(n_increment, each = nsim), p0),
        nrow = nsim, ncol = j
      ),
      treatment1 = matrix(
        rbinom(nsim * j, rep(n_increment, each = nsim), p1),
        nrow = nsim, ncol = j
      ),
      treatment2 = matrix(
        rbinom(nsim * j, rep(n_increment, each = nsim), p2),
        nrow = nsim, ncol = j
      )
    )
  })

  cumulative <- function(x) {
    if (j > 1L) {
      for (stage_index in 2:j) {
        x[, stage_index] <- x[, stage_index - 1L] + x[, stage_index]
      }
    }
    x
  }
  y0 <- cumulative(draws$control)
  y1 <- cumulative(draws$treatment1)
  y2 <- cumulative(draws$treatment2)

  z1 <- z2 <- matrix(0, nrow = nsim, ncol = j)
  for (stage_index in seq_len(j)) {
    z1[, stage_index] <- .proportion_z_vector(
      y1[, stage_index], y0[, stage_index],
      n_cumulative[stage_index], n_cumulative[stage_index]
    )
    z2[, stage_index] <- .proportion_z_vector(
      y2[, stage_index], y0[, stage_index],
      n_cumulative[stage_index], n_cumulative[stage_index]
    )
  }

  list(
    z1 = z1,
    z2 = z2,
    j = j,
    t = t,
    n = n,
    nsim = nsim,
    seed = as.integer(seed)
  )
}

.make_mams_null_cache <- function(n, j, t, p0, nsim, seed) {
  .make_mams_cache(
    n = n, j = j, t = t,
    p0 = p0, p1 = p0, p2 = p0,
    nsim = nsim, seed = seed
  )
}

.mams_reject_indicator <- function(z, upper) {
  nsim <- nrow(z)
  j <- ncol(z)
  alive <- rep(TRUE, nsim)
  rejected <- rep(FALSE, nsim)

  for (stage_index in seq_len(j)) {
    efficacy <- alive & z[, stage_index] >= upper[stage_index]
    rejected <- rejected | efficacy
    if (stage_index < j) {
      alive <- alive & !efficacy & z[, stage_index] >= 0
    }
  }
  rejected
}

.mams_fwer_from_cache <- function(cache, cv, method) {
  .mams_rejection_summary_from_cache(cache, cv, method)$reject_at_least_one
}

.mams_rejection_summary_from_cache <- function(cache, cv, method) {
  upper <- if (method == "OBF") cv / sqrt(cache$t) else rep(cv, cache$j)
  reject1 <- .mams_reject_indicator(cache$z1, upper)
  reject2 <- .mams_reject_indicator(cache$z2, upper)
  list(
    reject1 = mean(reject1),
    reject2 = mean(reject2),
    reject2not1 = mean(reject2 & !reject1),
    reject_at_least_one = mean(reject1 | reject2)
  )
}

.calibrate_mams_cv <- function(cache, c_base, alpha, method,
                               upper_multiplier = 1.15,
                               c_tolerance = 1e-6) {
  base_fwer <- .mams_fwer_from_cache(cache, c_base, method)
  if (base_fwer <= alpha) {
    return(list(cv = c_base, simulated_fwer = base_fwer, base_fwer = base_fwer))
  }

  lower <- c_base
  upper <- c_base * upper_multiplier
  upper_fwer <- .mams_fwer_from_cache(cache, upper, method)
  expansion_count <- 0L
  while (upper_fwer > alpha && expansion_count < 10L) {
    upper <- upper * 1.10
    upper_fwer <- .mams_fwer_from_cache(cache, upper, method)
    expansion_count <- expansion_count + 1L
  }
  if (upper_fwer > alpha) {
    stop("Unable to bracket a simulation-calibrated MAMS critical value.", call. = FALSE)
  }

  while ((upper - lower) > c_tolerance) {
    middle <- (lower + upper) / 2
    middle_fwer <- .mams_fwer_from_cache(cache, middle, method)
    if (middle_fwer <= alpha) {
      upper <- middle
      upper_fwer <- middle_fwer
    } else {
      lower <- middle
    }
  }

  list(cv = upper, simulated_fwer = upper_fwer, base_fwer = base_fwer)
}

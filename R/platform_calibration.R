# Internal helpers for planned-entry Platform finite-sample calibration.

.normalize_platform_t <- function(t) {
  if (length(t) != 1L || !is.finite(t) || t <= 0 || t >= 1) {
    stop("For Platform designs, `t` must be a finite scalar strictly between 0 and 1.", call. = FALSE)
  }
  as.numeric(t)
}

.next_platform_compatible_n <- function(n, t) {
  t <- .normalize_platform_t(t)
  candidate <- as.integer(ceiling(n))
  while (.platform_stage_sizes(candidate, t)$n1 != candidate) {
    candidate <- candidate + 1L
  }
  candidate
}

.platform_boundaries <- function(cv, method, t, t_arrive = NULL, sensitive = NULL) {
  t <- .normalize_platform_t(t)
  if (!is.null(t_arrive)) t_arrive <- .normalize_platform_t(t_arrive)
  if (identical(method, "POC")) {
    return(c(u11 = cv, u21 = cv, u22 = cv, u32 = cv))
  }
  if (!identical(method, "OBF")) {
    stop("`method` must be either \"OBF\" or \"POC\".", call. = FALSE)
  }
  use_planned_t <- identical(sensitive, "method1") || is.null(t_arrive)
  if (use_planned_t) {
    return(c(
      u11 = cv,
      u21 = cv * sqrt(t),
      u22 = cv * sqrt(t) / sqrt(1 - t),
      u32 = cv * sqrt(t)
    ))
  }
  u_final <- cv * sqrt(t_arrive)
  c(
    u11 = cv,
    u21 = u_final,
    u22 = u_final * sqrt(
      ((1 - t_arrive) + 0.5 * (t + t_arrive)) / (1 - t_arrive)
    ),
    u32 = u_final
  )
}

.platform_stage_sizes <- function(n, t, t_arrive = NULL) {
  t <- .normalize_platform_t(t)
  if (!is.null(t_arrive)) t_arrive <- .normalize_platform_t(t_arrive)
  if (length(n) != 1L || !is.finite(n) || n <= 0 || abs(n - round(n)) > 1e-8) {
    stop("`n` must be a positive integer for finite-sample Platform calibration.", call. = FALSE)
  }
  n <- as.integer(round(n))

  entry_fraction <- if (is.null(t_arrive)) t else t_arrive
  final_fraction <- if (is.null(t_arrive)) t else 0.5 * (t + t_arrive)
  n_stage1 <- as.integer(round(n * entry_fraction))
  n_stage2 <- as.integer(round(n * (1 - entry_fraction)))
  n_stage3 <- as.integer(round(n * final_fraction))
  list(
    stage1 = n_stage1,
    stage2 = n_stage2,
    stage3 = n_stage3,
    n0 = n_stage1 + n_stage2 + n_stage3,
    n1 = n_stage1 + n_stage2,
    n2 = n_stage2 + n_stage3,
    N_total = 2L * n_stage1 + 3L * n_stage2 + 2L * n_stage3,
    r1 = entry_fraction,
    r2 = 1 - entry_fraction,
    r3 = final_fraction
  )
}

.make_platform_cache <- function(n, t, p0, p1, p2, nsim, seed, t_arrive = NULL) {
  t <- .normalize_platform_t(t)
  if (!is.null(t_arrive)) t_arrive <- .normalize_platform_t(t_arrive)
  stage_sizes <- .platform_stage_sizes(n, t, t_arrive)
  n <- as.integer(round(n))
  if (length(nsim) != 1L || !is.finite(nsim) || nsim < 1000 || abs(nsim - round(nsim)) > 0) {
    stop("`calibration_nsim` must be an integer of at least 1000.", call. = FALSE)
  }
  nsim <- as.integer(nsim)
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

  r10 <- r11 <- r12 <- stage_sizes$r1
  r20 <- r21 <- r22 <- stage_sizes$r2
  r30 <- r32 <- stage_sizes$r3
  sizes <- c(
    stage_sizes$stage1, stage_sizes$stage1, stage_sizes$stage1,
    stage_sizes$stage2, stage_sizes$stage2, stage_sizes$stage2,
    stage_sizes$stage3, stage_sizes$stage3
  )
  probabilities <- c(p0, p1, p2, p0, p1, p2, p0, p2)

  draws <- .with_preserved_seed(seed, {
    matrix(
      rbinom(
        nsim * length(sizes),
        size = rep(sizes, times = nsim),
        prob = rep(probabilities, times = nsim)
      ),
      nrow = nsim,
      byrow = TRUE
    )
  })
  y0s1 <- draws[, 1L]
  y1s1 <- draws[, 2L]
  y2s1 <- draws[, 3L]
  y0s2 <- draws[, 4L]
  y1s2 <- draws[, 5L]
  y2s2 <- draws[, 6L]
  y0s3 <- draws[, 7L]
  y2s3 <- draws[, 8L]

  z11 <- .proportion_z_vector(y1s1, y0s1, sizes[2L], sizes[1L])
  z21 <- .proportion_z_vector(
    y1s1 + y1s2, y0s1 + y0s2,
    round(n * (r11 + r21)), round(n * (r10 + r20))
  )
  z22 <- .proportion_z_vector(y2s2, y0s2, sizes[6L], sizes[4L])
  z32 <- .proportion_z_vector(
    y2s2 + y2s3, y0s2 + y0s3,
    round(n * (r22 + r32)), round(n * (r20 + r30))
  )

  yc <- y0s1 + 2 * y0s2 + y0s3
  yt <- y1s1 + y1s2 + y2s2 + y2s3
  n_e <- (r11 + r21 + r22 + r32) * n
  n_c <- (r10 + 2 * r20 + r30) * n
  p_e <- yt / n_e
  p_c <- yc / n_c
  var_e <- p_e * (1 - p_e) / n_e
  var_c <- (
    n * r10 * p_c * (1 - p_c) +
      4 * n * r20 * p_c * (1 - p_c) +
      n * r30 * p_c * (1 - p_c)
  ) / n_c^2
  se_ct <- sqrt(var_e + var_c)
  z_ct <- (p_e - p_c) / se_ct
  fallback <- !is.finite(z_ct) | !is.finite(se_ct) | se_ct <= 0
  if (any(fallback)) {
    n_e_adjusted <- n_e + 2
    n_c_adjusted <- n_c + 2
    p_e_adjusted <- (yt + 1) / n_e_adjusted
    p_c_adjusted <- (yc + 1) / n_c_adjusted
    var_e_adjusted <- p_e_adjusted * (1 - p_e_adjusted) / n_e_adjusted
    var_c_adjusted <- (
      n * r10 * p_c_adjusted * (1 - p_c_adjusted) +
        4 * n * r20 * p_c_adjusted * (1 - p_c_adjusted) +
        n * r30 * p_c_adjusted * (1 - p_c_adjusted) +
        2 * p_c_adjusted * (1 - p_c_adjusted)
    ) / n_c_adjusted^2
    se_adjusted <- sqrt(var_e_adjusted + var_c_adjusted)
    z_adjusted <- (p_e_adjusted - p_c_adjusted) / se_adjusted
    z_ct[fallback] <- z_adjusted[fallback]
  }
  z_ct[!is.finite(z_ct)] <- 0

  list(
    z11 = z11, z21 = z21, z22 = z22, z32 = z32, z_ct = z_ct,
    n = n, t = t, t_arrive = t_arrive,
    stage_sizes = stage_sizes, nsim = nsim, seed = as.integer(seed)
  )
}

.make_platform_null_cache <- function(n, t, p0, nsim, seed, t_arrive = NULL) {
  .make_platform_cache(
    n = n, t = t, p0 = p0, p1 = p0, p2 = p0,
    nsim = nsim, seed = seed, t_arrive = t_arrive
  )
}

.platform_rejection_summary_from_cache <- function(cache, cv, method, order, sensitive = NULL) {
  upper <- .platform_boundaries(
    cv, method, cache$t, cache$t_arrive, sensitive
  )
  z11 <- cache$z11
  z21 <- cache$z21
  z22 <- cache$z22
  z32 <- cache$z32

  early_t1 <- z11 >= upper["u11"]
  continue_t1 <- !early_t1 & z11 >= 0
  combination_open <- !early_t1

  if (identical(order, "true")) {
    t1_selected_late <- continue_t1 & z21 >= upper["u21"]
    reject2_with_t1 <- t1_selected_late & (
      z22 >= upper["u22"] |
        z32 >= upper["u32"] |
        (z32 < upper["u32"] & cache$z_ct > upper["u32"])
    )
    t1_not_selected <- combination_open & !t1_selected_late
    reject2_without_t1 <- t1_not_selected & (
      z22 >= upper["u22"] |
        (z22 >= 0 & z22 < upper["u22"] & z32 >= upper["u32"])
    )
    reject2 <- reject2_with_t1 | reject2_without_t1
    # In the order-restricted decision tree, a late monotherapy crossing is
    # retained only when the combination is also selected; otherwise the
    # contradictory-evidence path is classified as no selection.
    reject1 <- early_t1 | reject2_with_t1
  } else if (identical(order, "false")) {
    reject1 <- early_t1 | (continue_t1 & z21 >= upper["u21"])
    reject2 <- combination_open & (
      z22 >= upper["u22"] |
        (z22 >= 0 & z22 < upper["u22"] & z32 >= upper["u32"])
    )
  } else {
    stop("`order` must be either \"true\" or \"false\".", call. = FALSE)
  }

  list(
    reject1 = mean(reject1),
    reject2 = mean(reject2),
    reject2not1 = mean(reject2 & !reject1),
    reject1not2 = mean(reject1 & !reject2),
    reject_at_least_one = mean(reject1 | reject2)
  )
}

.platform_fwer_from_cache <- function(cache, cv, method, order, sensitive = NULL) {
  .platform_rejection_summary_from_cache(
    cache, cv, method, order, sensitive
  )$reject_at_least_one
}

.calibrate_platform_cv <- function(cache, c_base, alpha, method, order,
                                   sensitive = NULL,
                                   upper_multiplier = 1.15,
                                   c_tolerance = 1e-6) {
  score <- function(cv) {
    .platform_fwer_from_cache(cache, cv, method, order, sensitive)
  }
  base_fwer <- score(c_base)
  if (base_fwer <= alpha) {
    return(list(cv = c_base, simulated_fwer = base_fwer, base_fwer = base_fwer))
  }

  lower <- c_base
  upper <- c_base * upper_multiplier
  upper_fwer <- score(upper)
  expansion_count <- 0L
  while (upper_fwer > alpha && expansion_count < 10L) {
    upper <- upper * 1.10
    upper_fwer <- score(upper)
    expansion_count <- expansion_count + 1L
  }
  if (upper_fwer > alpha) {
    stop("Unable to bracket a simulation-calibrated Platform critical value.", call. = FALSE)
  }

  while ((upper - lower) > c_tolerance) {
    middle <- (lower + upper) / 2
    middle_fwer <- score(middle)
    if (middle_fwer <= alpha) {
      upper <- middle
      upper_fwer <- middle_fwer
    } else {
      lower <- middle
    }
  }
  list(cv = upper, simulated_fwer = upper_fwer, base_fwer = base_fwer)
}

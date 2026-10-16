test_that("vectorized planned-entry Platform scorer agrees with mc_pro", {
  for (method in c("OBF", "POC")) {
    for (order in c("true", "false")) {
      cache <- monocombo:::.make_platform_cache(
        n = 36, t = 0.5,
        p0 = 0.15, p1 = 0.15, p2 = 0.15,
        nsim = 3000, seed = 8675309
      )
      vectorized <- monocombo:::.platform_rejection_summary_from_cache(
        cache, cv = 3.1, method = method, order = order
      )
      loop <- monocombo::mc_pro(
        n = 36, method = method,
        p0 = 0.15, p1 = 0.15, p2 = 0.15, cv = 3.1,
        t = 0.5, t_arrive = NULL, platform = "true",
        seed = 8675309, order = order, nsim = 3000
      )

      expect_equal(vectorized$reject1, unname(loop[1]), tolerance = 1e-12)
      expect_equal(vectorized$reject2not1, unname(loop[2]), tolerance = 1e-12)
      expect_equal(vectorized$reject2, unname(loop[3]), tolerance = 1e-12)
      expect_equal(vectorized$reject_at_least_one, unname(loop[4]), tolerance = 1e-12)
      if (order == "false") {
        expect_equal(vectorized$reject1not2, unname(loop[5]), tolerance = 1e-12)
      }
    }
  }
})

test_that("vectorized delayed-entry Platform scorer agrees with mc_pro", {
  for (method in c("OBF", "POC")) {
    for (order in c("true", "false")) {
      cache <- monocombo:::.make_platform_cache(
        n = 36, t = 0.5, t_arrive = 0.6,
        p0 = 0.15, p1 = 0.15, p2 = 0.15,
        nsim = 3000, seed = 246810
      )
      vectorized <- monocombo:::.platform_rejection_summary_from_cache(
        cache, cv = 2.9, method = method, order = order
      )
      loop <- monocombo::mc_pro(
        n = 36, method = method,
        p0 = 0.15, p1 = 0.15, p2 = 0.15, cv = 2.9,
        t = 0.5, t_arrive = 0.6, platform = "true",
        seed = 246810, order = order, nsim = 3000
      )

      expect_equal(vectorized$reject1, unname(loop[1]), tolerance = 1e-12)
      expect_equal(vectorized$reject2not1, unname(loop[2]), tolerance = 1e-12)
      expect_equal(vectorized$reject2, unname(loop[3]), tolerance = 1e-12)
      expect_equal(vectorized$reject_at_least_one, unname(loop[4]), tolerance = 1e-12)
      if (order == "false") {
        expect_equal(vectorized$reject1not2, unname(loop[5]), tolerance = 1e-12)
      }
    }
  }
})

test_that("planned-entry Platform calibration is reproducible", {
  expect_equal(monocombo:::.next_platform_compatible_n(35, 0.5), 36L)
  expect_equal(monocombo:::.next_platform_compatible_n(37, 0.6), 37L)

  cache1 <- monocombo:::.make_platform_null_cache(
    n = 36, t = 0.5, p0 = 0.15, nsim = 5000, seed = 314159
  )
  cache2 <- monocombo:::.make_platform_null_cache(
    n = 36, t = 0.5, p0 = 0.15, nsim = 5000, seed = 314159
  )
  expect_equal(cache1, cache2)

  calibrated <- monocombo:::.calibrate_platform_cv(
    cache1, c_base = 3.0, alpha = 0.025,
    method = "OBF", order = "true"
  )
  expect_lte(calibrated$simulated_fwer, 0.025)
})

test_that("mc_Design calibrates planned and delayed entry", {
  arguments <- list(
    alpha = 0.025, method = "OBF",
    p0 = 0.15, p1 = 0.4, p2 = 0.4, power = 0.8,
    t = 0.5, platform = "true", order = "true",
    adjust = "true", calibration_nsim = 5000,
    calibration_seed = 271828
  )
  design1 <- do.call(monocombo::mc_Design, arguments)
  design2 <- do.call(monocombo::mc_Design, arguments)

  expect_equal(design1, design2)
  expect_null(attr(design1, "calibration"))
  planned_null <- monocombo:::.make_platform_null_cache(
    n = design1$n1, t = 0.5, p0 = 0.15,
    nsim = 5000, seed = 271828
  )
  expect_lte(
    monocombo:::.platform_fwer_from_cache(
      planned_null, design1$u11, "OBF", "true"
    ),
    0.025
  )

  delayed <- do.call(
    monocombo::mc_Design,
    c(arguments, list(t_arrive = 0.6))
  )
  expect_null(attr(delayed, "calibration"))
  delayed_null <- monocombo:::.make_platform_null_cache(
    n = delayed$n1, t = 0.5, t_arrive = 0.6, p0 = 0.15,
    nsim = 5000, seed = 271828
  )
  expect_lte(
    monocombo:::.platform_fwer_from_cache(
      delayed_null, delayed$u11, "OBF", "true"
    ),
    0.025
  )
})

test_that("fixed-N calibration matches simulated integer stage sizes", {
  design <- monocombo::mc_Design(
    alpha = 0.025, method = "OBF",
    p0 = 0.15, p1 = 0.4, p2 = 0.4, power = 0.8,
    t = 0.5, t_arrive = 0.6, N_max = "true",
    platform = "true", order = "true", adjust = "true",
    calibration_nsim = 5000, calibration_seed = 271828
  )
  expect_equal(unname(unlist(design[c("n0", "n1", "n2", "N_total")])), c(56, 36, 34, 126))
  expect_null(attr(design, "calibration"))
  fixed_null <- monocombo:::.make_platform_null_cache(
    n = design$n1, t = 0.5, t_arrive = 0.6, p0 = 0.15,
    nsim = 5000, seed = 271828
  )
  expect_lte(
    monocombo:::.platform_fwer_from_cache(
      fixed_null, design$u11, "OBF", "true"
    ),
    0.025
  )
})

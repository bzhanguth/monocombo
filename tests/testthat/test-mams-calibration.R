test_that("MAMS stage fractions are validated consistently", {
  expect_equal(
    monocombo:::.normalize_mams_t(2, c(0.5, 1))$t,
    c(0.5, 1)
  )
  expect_error(
    monocombo:::.normalize_mams_t(2, c(0.5, 0.75, 1)),
    "exactly `j`"
  )
  expect_error(
    monocombo:::.normalize_mams_t(3, c(1 / 3, 0.8, 0.9)),
    "end at 1"
  )
})

test_that("vectorized MAMS FWER agrees with the path simulator", {
  cases <- list(
    list(j = 2, t = c(0.5, 1), n = 40),
    list(j = 3, t = c(1 / 3, 2 / 3, 1), n = 42)
  )

  for (case in cases) {
    for (method in c("OBF", "POC")) {
      cache <- monocombo:::.make_mams_null_cache(
        n = case$n, j = case$j, t = case$t,
        p0 = 0.15, nsim = 2000, seed = 777 + case$j
      )
      vectorized <- monocombo:::.mams_fwer_from_cache(cache, 2.3, method)
      loop <- unname(mc_pro(
        n = case$n, method = method,
        p0 = 0.15, p1 = 0.15, p2 = 0.15, cv = 2.3,
        j = case$j, t = case$t, platform = "false",
        seed = 777 + case$j, order = "true", nsim = 2000
      )[4])
      expect_equal(vectorized, loop, tolerance = 1e-12)
    }
  }
})

test_that("MAMS calibration is reproducible and returns only the final design", {
  design_args <- list(
    alpha = 0.025, method = "OBF",
    p0 = 0.15, p1 = 0.4, p2 = 0.4, power = 0.8,
    j = 2, t = c(0.5, 1), platform = "false",
    adjust = "true", calibration_nsim = 5000,
    calibration_seed = 271828
  )
  design1 <- do.call(mc_Design, design_args)
  design2 <- do.call(mc_Design, design_args)

  expect_equal(design1, design2)
  expect_setequal(
    names(design1),
    c("u1", "l1", "u2", "n0", "n1", "n2", "N_total")
  )
  expect_s3_class(design1, "data.frame")
  expect_equal(unname(unlist(design1[c("n0", "n1", "n2")])), rep(design1$n1, 3))

  null_cache <- monocombo:::.make_mams_null_cache(
    n = design1$n1, j = 2, t = c(0.5, 1),
    p0 = 0.15, nsim = 5000, seed = 271828
  )
  expect_lte(
    monocombo:::.mams_fwer_from_cache(null_cache, design1$u2, "OBF"),
    0.025
  )

  analytical_power <- monocombo:::.mams_analytical_fixed(
    design1$n1, "OBF", 0.15, 0.4, 0.4, design1$u2,
    j = 2, t = c(0.5, 1), platform = "false"
  )
  expect_gte(unname(analytical_power[4]), 0.8)

  power_cache <- monocombo:::.make_mams_cache(
    n = design1$n1, j = 2, t = c(0.5, 1),
    p0 = 0.15, p1 = 0.4, p2 = 0.4,
    nsim = 5000, seed = 271829
  )
  simulated_power <- monocombo:::.mams_rejection_summary_from_cache(
    power_cache, design1$u2, "OBF"
  )$reject_at_least_one
  expect_gte(simulated_power, 0.8)
})

test_that("three-stage MAMS uses the same one-row format", {
  design <- mc_Design(
    alpha = 0.025, method = "OBF",
    p0 = 0.15, p1 = 0.4, p2 = 0.4, power = 0.8,
    j = 3, t = c(1 / 3, 2 / 3, 1), platform = "false",
    adjust = "false"
  )

  expect_s3_class(design, "data.frame")
  expect_setequal(
    names(design),
    c("u1", "l1", "u2", "l2", "u3", "n0", "n1", "n2", "N_total")
  )
  expect_equal(design$N_total, design$n0 + design$n1 + design$n2)
})

test_that("one-row MAMS output feeds directly into mc_OC", {
  design <- mc_Design(
    alpha = 0.025, method = "OBF",
    p0 = 0.15, p1 = 0.4, p2 = 0.4, power = 0.8,
    j = 2, t = c(0.5, 1), platform = "false",
    adjust = "false"
  )

  result <- mc_OC(
    n = design$n1, cv = design$u2,
    p0 = 0.15, p1 = 0.4, p2 = 0.4,
    t = c(0.5, 1), platform = "false", j = 2,
    method = "OBF", seed = 123, nsim = 1000
  )

  expect_equal(nrow(result), 1L)
  expect_equal(ncol(result), 16L)
})

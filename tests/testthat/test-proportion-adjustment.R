test_that("proportion_test uses the ordinary Wald statistic when SE is positive", {
  expected <- (8 / 18 - 4 / 18) /
    sqrt((8 / 18) * (10 / 18) / 18 + (4 / 18) * (14 / 18) / 18)

  expect_equal(
    proportion_test(success1 = 8, success2 = 4, n1 = 18, n2 = 18),
    expected
  )
})

test_that("SE-zero cases use the add-one boundary safeguard", {
  expected_opposite_boundaries <- (19 / 20 - 1 / 20) /
    sqrt((19 / 20) * (1 / 20) / 20 + (1 / 20) * (19 / 20) / 20)

  expect_equal(proportion_test(0, 0, 18, 18), 0)
  expect_equal(proportion_test(18, 18, 18, 18), 0)
  expect_equal(
    proportion_test(18, 0, 18, 18),
    expected_opposite_boundaries
  )
})

test_that("one boundary group does not trigger the safeguard when SE is positive", {
  expected <- (0 / 18 - 4 / 18) /
    sqrt(0 + (4 / 18) * (14 / 18) / 18)

  expect_equal(proportion_test(0, 4, 18, 18), expected)
})

test_that("vectorized calibration statistic uses the same conditional rule", {
  observed <- monocombo:::.proportion_z_vector(
    success_treatment = c(8, 0, 18, 18),
    success_control = c(4, 0, 18, 0),
    n_treatment = 18,
    n_control = 18
  )
  expected <- c(
    proportion_test(8, 4, 18, 18),
    proportion_test(0, 0, 18, 18),
    proportion_test(18, 18, 18, 18),
    proportion_test(18, 0, 18, 18)
  )

  expect_equal(observed, expected)
})

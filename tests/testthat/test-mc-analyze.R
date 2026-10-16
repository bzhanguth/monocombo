test_that("mc_Analyze uses the weighted Stage-2 control variance for z312", {
  observed <- monocombo::mc_Analyze(
    data = list(
      y0s1 = 3, y0s2 = 4, y0s3 = 5,
      y1s1 = 7, y1s2 = 9,
      y2s2 = 7, y2s3 = 7,
      n = 36
    ),
    u = 2.132,
    t = 0.5,
    method = "OBF"
  )

  n_s1 <- n_s2 <- n_s3 <- 18
  n_combined <- n_s1 + 2 * n_s2 + n_s3
  p_t <- 30 / n_combined
  p_c <- 16 / n_combined
  expected <- round(
    (p_t - p_c) /
      sqrt(
        p_t * (1 - p_t) / n_combined +
          (n_s1 + 4 * n_s2 + n_s3) *
            p_c * (1 - p_c) / n_combined^2
      ),
    3
  )

  expect_equal(expected, 2.328)
  expect_equal(as.numeric(observed["z312"]), expected)
  expect_equal(unname(observed["Stage_3_Decision"]), "Select T1 and T2")
})

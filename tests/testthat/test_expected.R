library(purgeR)
context("Expected inbreeding coefficients")

testthat::test_that("Standard inbreeding", {
  testthat::expect_equal(purgeR::exp_F(Ne = 10, t = 0), 0.0, tolerance = 1e-5)
  testthat::expect_equal(purgeR::exp_F(Ne = 10, t = 50), 0.923055, tolerance = 1e-5)
  testthat::expect_equal(purgeR::exp_F(Ne = 50, t = 50), 0.3949939, tolerance = 1e-5)
  testthat::expect_equal(purgeR::exp_Fa(Ne = 10, t = 0), 0.0, tolerance = 1e-5)
  testthat::expect_equal(purgeR::exp_Fa(Ne = 10, t = 10), 0.9005597, tolerance = 1e-5)
  testthat::expect_equal(purgeR::exp_Fa(Ne = 50, t = 10), 0.3638145, tolerance = 1e-5)
  testthat::expect_equal(purgeR::exp_F(Ne = 10, t = 3), purgeR::exp_Fa(Ne = 10, t = 3), tolerance = 1e-5)
  testthat::expect_equal(purgeR::exp_F(Ne = 50, t = 3), purgeR::exp_Fa(Ne = 50, t = 3), tolerance = 1e-5)
  testthat::expect_equal(purgeR::exp_g(Ne = 10, t = 50, d = 0.0), purgeR::exp_F(Ne = 10, t = 50), tolerance = 1e-5)
  testthat::expect_equal(purgeR::exp_g(Ne = 10, t = 50, d = 0.2), 0.08010975, tolerance = 1e-5)
  testthat::expect_equal(purgeR::exp_g(Ne = 10, t = 50, d = 0.5), 0.004408481, tolerance = 1e-5)
  testthat::expect_equal(purgeR::exp_g(Ne = 50, t = 50, d = 0.5), 0.01616987, tolerance = 1e-5)
  testthat::expect_equal(purgeR::exp_g(Ne = 10, t = 100000, d = 0.1), 0.16666667, tolerance = 1e-5)
  testthat::expect_equal(purgeR::exp_g(Ne = 10, t = 100000, d = 0.2), 0.06976744, tolerance = 1e-5)
})

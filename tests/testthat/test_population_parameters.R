library(purgeR)
context("Population parameters")

data(atlas)
testthat::test_that("Number of founders and ancestors", {
  Ntarget <- sum(atlas$target)
  testthat::expect_equal(Ntarget, 176)
  Ntest <- purgeR::pop_Nancestors(atlas, reference = "target", seed = 1234)
  testthat::expect_equal(Ntarget, Ntest$Nr)
  testthat::expect_error(purgeR::pop_Nancestors(atlas, reference = "pom"), "Failed to coerce 'reference' values: All NAs.")
  atlas_badref <- atlas
  atlas_badref$target <- 0
  testthat::expect_error(purgeR::pop_Nancestors(atlas_badref, reference = "target"), "At least one 'reference' value should be TRUE.")
  atlas_badref$target <- atlas$target
  atlas_badref[943:948, ]$target <- NA
  testthat::expect_equal(purgeR::pop_Nancestors(atlas_badref, reference = "target")$Nr, 170)
  testthat::expect_equal(nrow(atlas[which(atlas$dam == 0 & atlas$sire == 0), ]), 5)
  testthat::expect_equal(Ntest$Nf, 4)
  testthat::expect_equal(Ntest$Na, 249)
  testthat::expect_equal(Ntest$Nfe, 3.583086, tolerance = 1e-5)
  testthat::expect_equal(Ntest$Nae, 3.325158, tolerance = 1e-5)
  testthat::expect_equal(Ntest$Ng, 2.002087, tolerance = 1e-5)
  testthat::expect_equal(Ntest$se_Ng, 0.4601456, tolerance = 1e-5)
  testthat::expect_equal(purgeR::pop_Nfe(atlas, reference = "target"), Ntest$Nfe, tolerance = 1e-5)
  testthat::expect_equal(purgeR::pop_Nae(atlas, reference = "target"), Ntest$Nae, tolerance = 1e-5)
  testthat::expect_equal(purgeR::pop_Ng(atlas, reference = "target", seed = 1234)$Ng, Ntest$Ng, tolerance = 1e-5)
  testthat::expect_equal(purgeR::pop_Ng(atlas, reference = "target", seed = 1234)$Ng, Ntest$Ng, tolerance = 1e-5)
  darwin_test <- darwin %>%
    purgeR::ped_rename(id = "Individual", dam = "Mother", sire = "Father") %>%
    dplyr::mutate(ref = ifelse(id > 60, 1, 0)) %>%
    purgeR::pop_Nancestors(reference = "ref", seed = 1234, skip_Ng = TRUE)
  testthat::expect_equal(darwin_test$Nae, 4.8, tolerance = 1e-5)
  darwin_test <- darwin %>%
    purgeR::ped_rename(id = "Individual", dam = "Mother", sire = "Father") %>%
    dplyr::mutate(ref = ifelse(id > 50, 1, 0)) %>%
    purgeR::pop_Nancestors(reference = "ref", seed = 1234, skip_Ng = TRUE)
  testthat::expect_equal(darwin_test$Nae, 6.145455, tolerance = 1e-5)
  darwin_test <- darwin %>%
    purgeR::ped_rename(id = "Individual", dam = "Mother", sire = "Father") %>%
    dplyr::mutate(ref = 1) %>%
    purgeR::pop_Nancestors(reference = "ref", seed = 1234, skip_Ng = TRUE)
  testthat::expect_equal(darwin_test$Nae, 22.24503, tolerance = 1e-5)
})

testthat::test_that("Number equivalent to complete generations", {
  testthat::expect_equal(base::ncol(atlas), 10)
  atlas_t <- purgeR::pop_t(atlas)
  testthat::expect_equal(base::ncol(atlas_t), 11)
  testthat::expect_equal(base::colnames(atlas_t)[length(atlas_t)], "t")
  testthat::expect_equal(base::min(atlas_t$t), 0.0, tolerance = 1e-5)
  testthat::expect_equal(base::max(atlas_t$t), 10.0639, tolerance = 1e-5)
  testthat::expect_equal(utils::tail(atlas_t$t, n = 1), 9.309570, tolerance = 1e-5)
  testthat::expect_equal(base::sum(atlas_t$t), 5500.627, tolerance = 1e-5)
  testthat::expect_warning(purgeR::pop_t(atlas_t))
  testthat::expect_error(atlas_t %>% ip_F() %>% dplyr::mutate(t = -t) %>% purgeR::pop_Ne(Fcol = "Fi", tcol = "t"), "Generations cannot take negative values.")
})

testthat::test_that("Effective population size", {
  testthat::expect_error(purgeR::pop_Ne(atlas))
  atlas_Ne <- atlas %>%
    purgeR::ip_F() %>%
    purgeR::pop_t() %>%
    purgeR::pop_Ne(Fcol = "Fi", tcol = "t")
  testthat::expect_equal(atlas_Ne$Ne, 8.184803, tolerance = 1e-5)
  testthat::expect_equal(atlas_Ne$se_Ne, 0.2498695, tolerance = 1e-5)
  atlas_RP_Ne <- atlas %>%
    purgeR::ip_F() %>%
    purgeR::pop_t() %>%
    dplyr::filter(target == 1) %>%
    purgeR::pop_Ne(Fcol = "Fi", tcol = "t")
  testthat::expect_equal(atlas_RP_Ne$Ne, 14.01041, tolerance = 1e-5)
  testthat::expect_equal(atlas_RP_Ne$se_Ne, 0.1707653, tolerance = 1e-5)
})

testthat::test_that("Hardy-Weinberg disequilibrium", {
  testthat::expect_equal(purgeR::pop_hwd(atlas), -0.009726449, tolerance = 1e-5)
  testthat::expect_equal(purgeR::pop_hwd(atlas, reference = "target"), -0.02753781, tolerance = 1e-5)
})

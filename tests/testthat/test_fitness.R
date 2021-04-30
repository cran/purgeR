library(purgeR)
context("Fitness values")

data(dorcas)
testthat::test_that("Offspring productivity", {
  dorcas_prod <- purgeR::w_offspring(dorcas, name_to = 'P')
  testthat::expect_equal(base::sum(dorcas_prod$P), 2469)
  testthat::expect_equal(base::max(dorcas_prod$P), 75)
  dorcas_prod <- purgeR::w_offspring(dorcas, name_to = 'P', sire_offspring = FALSE)
  testthat::expect_equal(base::sum(dorcas_prod$P), 1242)
  testthat::expect_equal(base::max(dorcas_prod$P), 18)
})

testthat::test_that("Grandoffspring productivity", {
  dorcas_prod <- purgeR::w_grandoffspring(dorcas, name_to = 'P')
  testthat::expect_equal(base::sum(dorcas_prod$P), 4302)
  testthat::expect_equal(base::max(dorcas_prod$P), 124)
})

testthat::test_that("Reproductive value", {
  # Pedigree from Hunter et al (2019)
  id <- c("A1", "A2", "A3", "A4", "A5", "A6",
          "B1", "B2", "B3", "B4",
          "C1", "C2", "C3", "C4")
  dam <- c("0", "0", "0", "0", "0", "0",
           "A2", "A2", "A2", "A4",
           "B2", "B2", "A4", "A6")
  sire <- c("0", "0", "0", "0", "0", "0",
            "A1", "A1", "A1", "A5",
            "B1", "B3", "B3", "A5")
  t <- c(0, 0, 0, 0, 0, 0,
         1, 1, 1, 1,
         2, 2, 2, 2)
  ped <- tibble::tibble(id, dam, sire, t)
  ped <- purgeR::ped_rename(ped, keep_names = TRUE) %>% dplyr::mutate(reference = ifelse(t == 1, TRUE, FALSE))
  ped_rv <- purgeR::w_reproductive_value(ped, reference = "reference", name_to = "R", enable_correction = FALSE)
  #testthat::expect_equal(ped$R, c(1, 1, 1, 1.5, 1.5, 1.5, 1.5, 2, 2, 1, 1.25, 1.25, 1, 1), tolerance = 1e-3)
  testthat::expect_equal(ped_rv$R, c(0, 0, 0, 0, 0, 0, 1.5, 2, 2, 1, 0, 0, 0, 0), tolerance = 1e-3)
  ped_rv <- purgeR::w_reproductive_value(ped, reference = "reference", name_to = "R", enable_correction = TRUE)
  testthat::expect_equal(ped_rv$R, c(0, 0, 0, 0, 0, 0, 0.11538462, 0.15384615, 0.15384615, 0.07692308, 0, 0, 0, 0), tolerance = 1e-5)

  ped <- ped %>% dplyr::mutate(reference = ifelse(id < 11, TRUE, FALSE), target = ifelse(id %in% c(11, 12), TRUE, FALSE))
  ped_rv <- purgeR::w_reproductive_value(ped, reference = "reference", name_to = "R", enable_correction = TRUE, target = "target")
  testthat::expect_equal(ped_rv$R, c(0.08333333, 0.08333333, 0.08333333, 0.08333333, 0.08333333, 0.08333333, 0.125, 0.16666667, 0.125, 0.08333333, 0, 0, 0, 0), tolerance = 1e-5)

  testthat::expect_error(ped %>%
                           purgeR::w_reproductive_value(reference = "reference", name_to = "R", target = "reference"),
                         "Cannot use reference individuals as target at the same time.")
  testthat::expect_error(ped %>%
                           dplyr::mutate(target = ifelse(id == 12, TRUE, FALSE)) %>% 
                           purgeR::w_reproductive_value(reference = "reference", name_to = "R", target = "target", generation_wise = TRUE),
                         "Cannot define 'target' individuals under generation wise mode.")
  # testthat::expect_error(ped %>%
  #                          dplyr::mutate(target = ifelse(id == 8, TRUE, FALSE),
  #                                        reference = ifelse(id == 12,TRUE,FALSE)) %>%
  #                          purgeR::w_reproductive_value(reference = "reference", name_to = "R", target = "target"),
  #                        "Found reference individual who is a descendant of the target population.")
  # testthat::expect_warning(ped %>%
  #                            dplyr::mutate(target = ifelse(id == 6, TRUE, FALSE)) %>%
  #                            purgeR::w_reproductive_value(reference = "reference", name_to = "R", target = "target"),
  #                          "Target individuals should always have 'id' lower than reference individuals.")

  arrui_rv <- arrui %>%
    purgeR::pop_t() %>%
    dplyr::mutate(t = plyr::round_any(t, 1), t = as.integer(t)) %>%
    w_reproductive_value(reference = "t", name_to = "R", generation_wise = TRUE)
  testthat::expect_equal(base::mean(arrui_rv$R), 0.01190894, tolerance = 1e-5)
  testthat::expect_equal(base::max(arrui_rv$R), 0.7060084, tolerance = 1e-5)
  testthat::expect_equal(base::sum(arrui_rv$R), 4.525398, tolerance = 1e-5)
})
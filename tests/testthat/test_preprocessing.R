library(purgeR)
context("Pedigree preprocessing")

data(arrui)
testthat::test_that("Mandatory columns are properly processed", {
  arrui_error <- arrui
  colnames(arrui_error)[3] <- c("not_sire")
  testthat::expect_error(purgeR::ped_rename(arrui_error))
  colnames(arrui_error)[3] <- c("sire")
  colnames(arrui_error)[2] <- c("not_dam")
  testthat::expect_error(purgeR::ped_rename(arrui_error))
  colnames(arrui_error)[2] <- c("dam")
  colnames(arrui_error)[1] <- c("not_id")
  testthat::expect_error(purgeR::ped_rename(arrui_error))
  testthat::expect_success(testthat::expect_s3_class(purgeR::ped_rename(arrui), "data.frame"))
  testthat::expect_success(testthat::expect_s3_class(purgeR::ped_clean(arrui, "survival15"), "data.frame"))
  testthat::expect_error(testthat::expect_s3_class(purgeR::ped_clean(arrui, "survival30"), "data.frame"), "Column not found: survival30")
})

testthat::test_that("Arrui pedigree is renamed and cleaned", {
  arrui_rename <- purgeR::ped_rename(arrui)
  testthat::expect_equal(arrui_rename$id, arrui$id)
  testthat::expect_error(purgeR::ped_clean(arrui_rename), "argument \"value_from\" is missing, with no default")
  arrui_clean <- purgeR::ped_clean(arrui_rename, "survival15")
  testthat::expect_equal(base::nrow(arrui_clean), 374)
  testthat::expect_equal(arrui_clean$id, 1:base::nrow(arrui_clean)) # clean also renames
  testthat::expect_true(base::all(c(309, 310, 330, 331, 332, 333) %in% arrui$id))
  # testthat::expect_false(any(c(309, 310, 330, 331, 332, 333)) %in% arrui_clean$id)
  testthat::expect_true(base::all(c(309, 310, 330, 331, 332, 333) %in% arrui_clean$id))
  arrui_clean <- purgeR::ped_clean(arrui_rename, "prod")
  testthat::expect_equal(base::nrow(arrui_clean), 134)
  testthat::expect_equal(arrui_clean$id, 1:base::nrow(arrui_clean))
  arrui_final <- purgeR::ped_rename(arrui_clean)
  testthat::expect_equal(arrui_final$id, 1:base::nrow(arrui_final))
  testthat::expect_equal(base::nrow(arrui_final), base::nrow(arrui_clean))
})

testthat::test_that("Rename special arguments work", {
  data(darwin)
  testthat::expect_error(purgeR::ped_rename(darwin), "Mandatory column 'id' not found")
  testthat::expect_error(purgeR::ped_rename(darwin, id = "Individual"), "Mandatory column 'dam' not found")
  testthat::expect_error(purgeR::ped_rename(darwin, id = "Individual", dam = "Mother"), "Mandatory column 'sire' not found")
  darwin_rename <- purgeR::ped_rename(darwin, id = "Individual", dam = "Mother", sire = "Father", keep_names = TRUE)
  testthat::expect_equal(colnames(darwin), c("Individual", "Mother", "Father"))
  testthat::expect_equal(colnames(darwin_rename), c("id", "dam", "sire", "names"))
  testthat::expect_equal(darwin$Individual, darwin_rename$names)
  testthat::expect_equal(class(darwin_rename$id), "integer")
  testthat::expect_equal(class(darwin_rename$dam), "integer")
  testthat::expect_equal(class(darwin_rename$sire), "integer")
  testthat::expect_equal(darwin_rename$id, 1:nrow(darwin))
  testthat::expect_equal(darwin_rename[52, ]$id, 52)
  testthat::expect_equal(darwin_rename[52, ]$dam, 44)
  testthat::expect_equal(darwin_rename[52, ]$sire, 43)
  testthat::expect_warning(purgeR::ped_rename(darwin_rename, keep_names = TRUE))
})

testthat::test_that("Input types are handled", {
  arrui_rename <- purgeR::ped_rename(arrui)
  arrui_str <- arrui_rename
  arrui_str$id <- base::as.character(arrui_rename$id)
  testthat::expect_success(testthat::expect_s3_class(purgeR::ped_rename(arrui_str), "data.frame"))
  testthat::expect_error(purgeR::ped_clean(arrui_str, "survival15"), "Mandatory 'id', 'dam' and 'sire' columns need to be of type integer")
  arrui_str <- arrui_rename
  arrui_str$dam <- base::as.character(arrui_rename$dam)
  testthat::expect_error(purgeR::ped_clean(arrui_str, "survival15"), "Mandatory 'id', 'dam' and 'sire' columns need to be of type integer")
  arrui_str <- arrui_rename
  arrui_str$sire <- base::as.character(arrui_rename$sire)
  testthat::expect_error(purgeR::ped_clean(arrui_str, "survival15"), "Mandatory 'id', 'dam' and 'sire' columns need to be of type integer")
  arrui_str <- arrui_rename
  arrui_str$prod <- base::as.character(arrui_rename$prod)
  testthat::expect_equal(nrow(purgeR::ped_clean(arrui_str, "prod")), nrow(purgeR::ped_clean(arrui, "prod")))
})

testthat::test_that("Impossible pedigrees should always return error", {
  # Selfing
  #ped_error <- data.frame(id = c(1, 2, 3), dam = c(0, 0, 2), sire = c(0, 0, 2))
  #testthat::expect_error(purgeR::rename(ped_error))
  #ped_error <- data.frame(id = c("A", "B", "C"), dam = c("X", "X", "B"), sire = c("X", "X", "B"))
  #testthat::expect_error(purgeR::rename(ped_error))
  # Circular kinship
  ped_error <- data.frame(id =c(1L, 2L, 3L), dam = c(1L, 0L, 1L), sire = c(1L, 0L, 2L))
  testthat::expect_error(purgeR::ped_rename(ped_error))
  ped_error <- data.frame(id = c(1L, 2, 3), dam = c(1, 0, 1), sire = c(0L, 0L, 2L))
  testthat::expect_error(purgeR::ped_rename(ped_error))
  ped_error <- data.frame(id = c(1L, 2L, 3L), dam = c(0L, 0L, 1L), sire = c(0L, 2L, 2L))
  testthat::expect_error(purgeR::ped_rename(ped_error))
  ped_error <- data.frame(id = c("A", "B", "C"), dam = c("A", "X", "A"), sire = c("A", "X", "B"))
  testthat::expect_error(purgeR::ped_rename(ped_error))
  ped_error <- data.frame(id = c("A", "B", "C"), dam = c("A", "X", "A"), sire = c("X", "X", "B"))
  testthat::expect_error(purgeR::ped_rename(ped_error))
  ped_error <- data.frame(id = c("A", "B", "C"), dam = c("X", "X", "A"), sire = c("X", "B", "B"))
  testthat::expect_error(purgeR::ped_rename(ped_error))
  # Unknown individuals
  ped_error <- data.frame(id = c(0L, 1L, 2L), dam = c(0L, 0L, 0L), sire = c(0L, 0L, 2L))
  testthat::expect_error(purgeR::ped_rename(ped_error))
  # Repeated individuals
  ped_error <- data.frame(id = c(1L, 2L, 1L), dam = c(0L, 0L, 1L), sire = c(0L, 0L, 2L))
  testthat::expect_error(purgeR::ped_rename(ped_error))
  ped_error <- data.frame(id = c("A", "B", "A"), dam = c("X", "X", "A"), sire = c("X", "X", "B"))
  testthat::expect_error(purgeR::ped_rename(ped_error))
  # Space-time continuum ruptures
  ped_error <- data.frame(id = c(1L, 2L, 3L), dam = c(2L, 0L, 0L), sire = c(3L, 0L, 0L))
  testthat::expect_error(purgeR::ped_rename(ped_error))
  ped_error <- data.frame(id = c("A", "B", "C"), dam = c("B", "X", "X"), sire = c("C", "X", "X"))
  testthat::expect_error(purgeR::ped_rename(ped_error))
  ped_no_error <- data.frame(id = c(1L, 2L, 3L), dam = c(4L, 5L, 2L), sire = c(6L, 7L, 1L))
  testthat::expect_success(testthat::expect_s3_class(purgeR::ped_rename(ped_no_error), "data.frame"))
  # Individuals out of order
  ped_error <- data.frame(id = c(1L, 4L, 3L), dam = c(0L, 0L, 1L), sire = c(0L, 0L, 4L))
  testthat::expect_error(purgeR::ped_clean(ped_error, value_from =  "prod"))
})

testthat::test_that("Sort pedigrees", {
  # Sort pedigree with IDs as characters
  set.seed(1234)
  darwin_unsort <- darwin[base::sample(1:nrow(darwin)),]
  testthat::expect_error(purgeR::ip_F(darwin_unsort))
  darwin_sort <- purgeR::ped_sort(darwin_unsort, id = "Individual", dam = "Mother", sire = "Father")
  darwin_rename <- purgeR::ped_rename(darwin, id = "Individual", dam = "Mother", sire = "Father", keep_names = TRUE)
  testthat::expect_equal(mean(purgeR::ip_F(darwin_rename)$Fi), mean(purgeR::ip_F(darwin_sort)$Fi))
  # Sort pedigree with IDs as integers
  arrui_unsort <- arrui[base::sample(1:nrow(arrui)),]
  testthat::expect_error(purgeR::ip_F(arrui_unsort))
  arrui_sort <- purgeR::ped_sort(arrui_unsort)
  testthat::expect_equal(mean(purgeR::ip_F(arrui)$Fi), mean(purgeR::ip_F(arrui_sort)$Fi))
})
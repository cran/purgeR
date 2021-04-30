#' Offspring
#'
#' Counts the number of offspring for individuals in the pedigree.
#'
#' @template ped-arg
#' @template nameto-arg
#' @param dam_offspring Compute dam offspring (TRUE by default).
#' @param sire_offspring Compute sire offspring (TRUE by default).
#' @return The input dataframe, plus an additional column indicating the total number of offspring.
#' @examples
#' data(dama)
#' w_offspring(dama, name_to = "P")
#' @export
w_offspring <- function(ped, name_to, dam_offspring = TRUE, sire_offspring = TRUE) {

  # Check input errors
  check_basic(ped, "id", "dam", "sire")
  check_not_col(base::colnames(ped), name_to)
  check_bool(dam_offspring)
  check_bool(sire_offspring)
  if (!dam_offspring & !sire_offspring) stop("At least one of dam or sire productivity must be set TRUE.")

  # Return productivity
  id_names <- base::as.character(ped$id)
  dam_freq <- base::table(ped$dam)
  sire_freq <- base::table(ped$sire)
  x <- ifelse((id_names %in% base::names(dam_freq)) & dam_offspring,
    dam_freq[id_names],
    ifelse((id_names %in% base::names(sire_freq)) & sire_offspring,
      sire_freq[id_names],
      0.0
    )
  )
  ped[, name_to] <- x
  ped
}

#' Grandoffspring
#'
#' Counts the number of grandoffspring for individuals in the pedigree.
#'
#' @template ped-arg
#' @template nameto-arg
#' @return The input dataframe, plus an additional column indicating the total number of grandoffspring.
#' @examples
#' data(dama)
#' w_grandoffspring(dama, name_to = "G")
#' @export
w_grandoffspring <- function(ped, name_to) {
  
  # Check input errors
  check_basic(ped, "id", "dam", "sire")
  check_not_col(base::colnames(ped), name_to)
  
  # Return (grandchildren) productivity
  N <- base::nrow(ped)
  productivity <- base::numeric(N)
  for (i in 1:N) {
    offspring <- ped[base::which(ped$dam == i | ped$sire == i),]$id
    if (base::length(offspring) > 0) {
      grandchildren <- base::nrow(ped[base::which(ped$dam %in% offspring | ped$sire %in% offspring),])
      productivity[i] <- grandchildren
    }
  }
  ped[, name_to] <- productivity
  ped
}

#' Reproductive value
#'
#' Computes the reproductive value following the method by Hunter et al. (2019).
#' This is a measure of how well a gene originated in a set of 'reference' individuals
#' is represented in a different set of 'target' individuals. By default, fitness is
#' computed for individuals in the reference population, using all of their descendants
#' as target. A generation-wise mode can also be enabled, to compute fitness
#' contributions consecutively from one generation to the next.
#' 
#' A reference population must be defined, which represents a set of individuals
#' whose reproductive value is to be calculated.
#' By default, genetic contributions to the rest of individuals in the pedigree
#' is assumed, but a target population can also be defined, restricting the set
#' of individuals accounted when computing the reproductive value. 
#' This could represent for example a cohort of alive individuals.
#'
#' @template ped-arg
#' @template reference-arg
#' @template nameto-arg
#' @template target-arg
#' @param enable_correction Correct reproductive values (enabled by default).
#' @param generation_wise Assume that the reference population is a vector of integers indicating generation numbers. Reproductive values will be computed generation by generation independently (except for the last one).
#' @return The input dataframe, plus an additional column with reproductive values for the reference and target populations assumed.
#' @references
#' \itemize{
#'   \item{Hunter DC et al. 2019. Pedigree-based estimation of reproductive value. Journal of Heredity 10(4): 433-444.}
#' }
#' @examples
#' library(dplyr)
#' library(magrittr)
#' # Pedigree used in Hunter et al. (2019)
#' id <- c("A1", "A2", "A3", "A4", "A5", "A6",
#'         "B1", "B2", "B3", "B4",
#'         "C1", "C2", "C3", "C4")
#' dam <- c("0", "0", "0", "0", "0", "0",
#'          "A2", "A2", "A2", "A4",
#'          "B2", "B2", "A4", "A6")
#' sire <- c("0", "0", "0", "0", "0", "0",
#'           "A1", "A1", "A1", "A5",
#'           "B1", "B3", "B3", "A5")
#' t <- c(0, 0, 0, 0, 0, 0,
#'        1, 1, 1, 1,
#'        2, 2, 2, 2)
#'ped <- tibble::tibble(id, dam, sire, t)
#'ped <- purgeR::ped_rename(ped, keep_names = TRUE) %>%
#'  dplyr::mutate(reference = ifelse(t == 1, TRUE, FALSE))
#'purgeR::w_reproductive_value(ped, reference = "reference", name_to = "R")
#' @export
w_reproductive_value <- function(ped, reference, name_to, target = NULL, enable_correction = TRUE, generation_wise = FALSE) {

  # Check input errors
  check_basic(ped, "id", "dam", "sire")
  check_col(base::colnames(ped), reference)
  check_not_col(base::colnames(ped), name_to)
  ref <- check_reference(ped, reference, "reference")
  tgt <- check_target(ped, reference, target, "target")
  check_bool(enable_correction)

  # Return reproductive value
  check_bool(generation_wise)
  if (!generation_wise) .Call(`_purgeR_reproductive_value`, ped, ref, name_to, tgt, enable_correction)
  # Generation-wise mode
  else {
    if (!base::is.null(target)) stop ("Cannot define 'target' individuals under generation wise mode.")
    check_tcol(ped, tcol = reference, compute = FALSE, force_int = TRUE)
    tmp_ped <- ped
    tmp_reference <- base::paste("tmp_", reference, sep = "")
    tmp_name_to <- base::paste("tmp_", name_to, sep = "")
    rv <- base::numeric(base::nrow(ped))
    generations <- base::sort(base::unique(ped[, reference]))
    generations <- generations[1:(base::length(generations)-1)]
    for (i in generations) {
      tmp_ped[, tmp_reference] <- ifelse(tmp_ped[, reference] == i, TRUE, FALSE)
      tmp_ref <- check_reference(tmp_ped, tmp_reference, "tmp_reference")
      tmp_rv <- .Call(`_purgeR_reproductive_value`, tmp_ped, tmp_ref, tmp_name_to, target, enable_correction)[, tmp_name_to]
      rv <- ifelse(tmp_rv != 0.0, tmp_rv, rv)
    }
    ped[, name_to] <- rv
    ped
  }
}
#' Remove individuals not used for purging analyses
#'
#' Remove individuals that are not necessary for purging analyses involving fitness.
#' This will reduce the size of the pedigree, and speed up the computation of inbreeding
#' parameters.
#' Individuals removed include those with unknown (NA)
#' values of a given parameter and that do not have any descendant in the
#' pedigree with known values of that parameter.
#' Cleaned pedigrees will automatically have individual identities
#' renamed (see \code{\link[purgeR]{ped_rename}}).
#'
#' @template ped-arg
#' @template valuefrom-arg
#' @return A dataframe with the pedigree cleaned for the specificed parameter (column) provided.
#' @seealso \code{\link[purgeR]{ped_rename}}
#' @examples
#' data(arrui)
#' nrow(arrui)
#' arrui <- ped_clean(arrui, "survival15")
#' nrow(arrui)
#' @export
ped_clean <- function(ped, value_from) {
  check_basic(ped)
  check_col(base::colnames(ped), value_from)
  idx <- .Call(`_purgeR_evaluate`, ped, value_from)
  ped <- ped[idx, ]
  purgeR::ped_rename(ped)
}

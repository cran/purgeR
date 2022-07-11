#' Rename individuals in a pedigree from 1 to N
#'
#' Functions in \bold{purgeR} require individuals to be named with integers from 1 to N.
#' This takes a dataframe containing a pedigree, and rename individuals having
#' names in any format to that required by other functions in \pkg{purgeR}. The
#' process will also check that the pedigree format is suitable for other functions
#' in the package.
#' 
#' @template ped-arg
#' @template id-arg
#' @template dam-arg
#' @template sire-arg
#' @template keepnames-arg
#' @return A dataframe with the pedigree's identities renamed.
#' @seealso \code{\link{ped_clean}}
#' @examples
#' data(darwin)
#' darwin <- ped_rename(darwin, id = "Individual", dam = "Mother", sire = "Father", keep_names = TRUE)
#' head(darwin)
#' @export
ped_rename <- function(ped, id = "id", dam = "dam", sire = "sire", keep_names = FALSE) {

  check_basic(ped, id_name = id, dam_name = dam, sire_name = sire, when_rename = TRUE)
  check_bool(keep_names)
  if (keep_names) check_not_col(base::colnames(ped), "names")

  ped_renamed <- ped
  ped_renamed <- dplyr::rename(ped_renamed, "id" = id, "dam" = dam, "sire" = sire)
  ped_renamed <- as.data.frame(ped_renamed)
  names <- base::as.character(ped_renamed[["id"]])
  ped_renamed["id"] <- 1:base::nrow(ped)
  ped_renamed["dam"] <- base::as.character(ped_renamed[["dam"]])
  ped_renamed["sire"] <- base::as.character(ped_renamed[["sire"]])
  ped_renamed["dam"] <- base::match(ped_renamed[["dam"]], names)
  ped_renamed["sire"] <- base::match(ped_renamed[["sire"]], names)
  ped_renamed["dam"] <- ifelse(is.na(ped_renamed[["dam"]]), 0, ped_renamed[["dam"]])
  ped_renamed["sire"] <- ifelse(is.na(ped_renamed[["sire"]]), 0, ped_renamed[["sire"]])
  ped_renamed["dam"] <- base::as.integer(ped_renamed[["dam"]])
  ped_renamed["sire"] <- base::as.integer(ped_renamed[["sire"]])

  ped_renamed <- base::as.data.frame(ped_renamed)
  if (keep_names) {
    ped_renamed["names"] <- names
  }
  ped_renamed
}
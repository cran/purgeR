#' Deviation from Hardy-Weinberg equilibrium
#'
#' Computes the deviation from Hardy-Weinberg equilibrium following Caballero and Toro (2000).
#' 
#' @template ped-arg
#' @template reference-arg
#' @return A numeric value indicating the deviation from Hardy-Weinberg equilibrium.
#' @seealso \code{\link{pop_Ne}}
#' @references
#' \itemize{
#'   \item{Caballero A, Toro M. 2000. Interrelations between effective population size and other pedigree tools for the management of conserved populations. Genet. Res. 75: 331-343.}
#' }
#' @examples
#' data(atlas)
#' pop_hwd(dama)
#' @export
pop_hwd <- function(ped, reference = NULL) {
  check_basic(ped, "id", "dam", "sire")
  if (!base::is.null(reference)) reference <- check_reference(ped, reference)
  .Call(`_purgeR_hwd`, ped, reference)
}
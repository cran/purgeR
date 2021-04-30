#' Number of equivalent complete generations
#'
#' Computes the number of equivalent complete generations (\emph{t}), as defined by Boichard et al (1997).
#' 
#' @template ped-arg
#' @template nameto-arg
#' @return The input dataframe, plus an additional column corresponding to the number of equivalent complete generations of every individual (named "t" by default).
#' @seealso \code{\link{pop_Ne}}
#' @references
#' \itemize{
#'   \item{Boichard D, Maignel L, Verrier E. 1997. The value of using probabilities of gene origin to measure genetic variability in a population. Genet. Sel. Evol., 29: 5-23.}
#' }
#' @examples
#' data(dama)
#' pop_t(dama)
#' @export
pop_t <-function(ped, name_to = "t") {
  check_basic(ped)
  check_not_col(base::colnames(ped), name_to)
  N <- base::nrow(ped)
  t <- base::numeric(N)
  for (i in 1:N) {
    if (ped$dam[i]) t[i] <- t[i] + 0.5 + 0.5 * t[ped$dam[i]]
    if (ped$sire[i]) t[i] <- t[i] + 0.5 + 0.5 * t[ped$sire[i]]
  }
  ped[name_to] <- t;
  ped
}
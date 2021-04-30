#' Effective population size
#'
#' Estimate the effective population size (\emph{Ne}). This is computed from the increase in individual inbreeding, following the method described by Gutiérrez et al (2008, 2009).
#'
#' @template ped-arg
#' @param Fcol Name of column with inbreeding coefficient values.
#' @param tcol Name of column with generation numbers.
#' @return A list with the effective population size (Ne) and its standard error (se_Ne).
#' @seealso \code{\link{ip_F}}, \code{\link{pop_t}}
#' @encoding UTF-8
#' @references
#' \itemize{
#'   \item{Gutiérrez JP, Cervantes I, Molina A, Valera M, Goyache F. 2008. Individual increase in inbreeding allows estimating effective sizes from pedigrees. Genet. Sel. Evol. 40: 359-378.}
#'   \item{Gutiérrez JP, Cervantes I, Goyache F. 2009. Improving the estimation of realized effective population sizes in farm animals. J. Anim. Breed. Genet. 126: 327-332.}
#' }
#' @examples
#' data(atlas)
#' atlas <- ip_F(atlas) # compute inbreeding, appending column "F"
#' atlas <- pop_t(atlas) # compute generations, appending column "t"
#' pop_Ne(atlas, Fcol = "Fi", tcol = "t")
#' @export
pop_Ne <- function(ped, Fcol, tcol) {
  #check_basic(ped, "id", "dam", "sire")
  F_ <- check_Fcol(ped, Fcol, compute = FALSE)
  t_ <- check_tcol(ped, tcol, compute = FALSE)
  delta_ <- delta_Fi(F_, t_)
  Ne_ <- Ne_delta(delta_)
  se_Ne_ <- se_Ne_delta(delta_)
  list(
    Ne = Ne_,
    se_Ne = se_Ne_
  )
}

#' Individual inbreeding variation
#'
#' Computes the increase in inbreeding coefficient for a given individual
#'
#' @name delta_Fi
#' @param Fi Individual inbreeding coefficient.
#' @param t Individual generation number.
#' @return Individual variation in inbreeding.
delta_Fi <- function(Fi, t) {
  1.0 - (1.0 - Fi)^(1.0 / (t - 1))
}

#' Realized effective population size (mean)
#'
#' Computes the mean realized effective population size.
#' Note this function expected a mean delta_F value for all individuals in the reference population
#'
#' @name Ne_delta
#' @param delta Vector of individual variations in inbreeding.
#' @return Mean effective population size.
Ne_delta <- function(delta) {
  delta_F <- base::mean(delta, na.rm = TRUE)
  1.0 / (2.0 * delta_F)
}

#' Realized effective population size (standard error)
#'
#' Computes the standard error of the realized effective population size.
#' Note this function expects the mean and standard deviation of delta F, as well as the total number of individuals in the reference population
#'
#' @name Ne_delta
#' @param delta Vector of individual variations in inbreeding.
#' @return Standard error of the effective population size.
se_Ne_delta <- function(delta) {
  N <- length(delta)
  delta_F <- base::mean(delta, na.rm = TRUE)
  sd_delta_F <- stats::sd(delta, na.rm = TRUE)
  Ne_ <- Ne_delta(delta_F)
  (2.0 / sqrt(N)) * Ne_ * Ne_ * sd_delta_F
}

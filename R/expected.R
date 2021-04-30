#' Expected inbreeding coefficient
#'
#' Estimates the expected inbreeding coefficient (F) as a function of the effective population size and generation number
#'
#' Computation of the inbreeding coefficient uses the classical formula:
#' 
#' F(t) = 1 - (1 - 1/2N) ^ t
#'
#' @template Ne-arg
#' @template t-arg
#' @return The inbreeding coefficient
#' @seealso \code{\link{ip_F}}
#' @encoding UTF-8
#' @references
#' \itemize{
#'   \item{Falconer DS, Mackay TFC. 1996. Introduction to Quantitative Genetics. 4th edition. Longman, Essex, U.K.}
#' }
#' @examples
#' exp_F(Ne = 50, t = 0)
#' exp_F(Ne = 50, t = 50)
#' exp_F(Ne = 10, t = 50)
#' @export
exp_F <- function(Ne, t) {
  check_Ne(Ne)
  check_int(t)
  if (t == 0) return (0.0)
  1 - (1 - 1/(2*Ne))^t
}

#' Expected ancestral inbreeding coefficient
#'
#' Estimates the expected ancestral inbreeding coefficient (Fa) as a function of the effective population size and generation number
#'
#' Computation of the ancestral inbreeding coefficient uses the adaptation from Ballou's (1997) formula, as in López-Cortegano et al. (2018):
#' 
#' Fa(t) = 1 - (1 - 1/2N) ^ (1/2 (t-1)t)
#'
#' @template Ne-arg
#' @template t-arg
#' @return The ancestral inbreeding coefficient
#' @seealso \code{\link{ip_Fa}}
#' @encoding UTF-8
#' @references
#' \itemize{
#'   \item{Ballou JD. 1997. Ancestral inbreeding only minimally affects inbreeding depression in mammalian populations. J Hered. 88:169–178.}
#'   \item{López-Cortegano E et al. 2018. Detection of genetic purging and predictive value of purging parameters estimated in pedigreed populations. Heredity 121(1): 38-51.}
#' }
#' @examples
#' exp_Fa(Ne = 50, t = 0)
#' exp_Fa(Ne = 50, t = 50)
#' exp_Fa(Ne = 10, t = 50)
#' @export
exp_Fa <- function(Ne, t) {
  check_Ne(Ne)
  check_int(t)
  if (t == 0) return (0.0)
  1 - (1 - 1/(2*Ne))^(0.5 * t * (t-1))
}

#' Expected purged inbreeding coefficient
#'
#' Estimates the expected purged inbreeding coefficient (Fa) as a function of the effective population size, generation number, and purging coefficient
#'
#' Computation of the purged inbreeding coefficient is calculated as in García-Dorado (2012):
#' 
#' g(t) = [ (1 - 1/2N) g(t-1)  + 1/2N] * [1 - 2d F(t-1)]
#'
#' @template Ne-arg
#' @template t-arg
#' @template d-arg
#' @return The purged inbreeding coefficient
#' @seealso \code{\link{ip_g}}
#' @encoding UTF-8
#' @references
#' \itemize{
#'   \item{García-Dorado. 2012. Understanding and predicting the fitness decline of shrunk populations: Inbreeding, purging, mutation, and standard selection. Genetics 190: 1-16.}
#' }
#' @examples
#' exp_g(Ne = 50, t = 0, d = 0.15)
#' exp_g(Ne = 50, t = 50, d = 0.15)
#' exp_g(Ne = 10, t = 50, d = 0.15)
#' @export
exp_g <- function(Ne, t, d) {
  check_Ne(Ne)
  check_int(t)
  check_d(d)
  if (t == 0) return (0.0)
  k <- 1 / (2*Ne)
  ((1 - k)*purgeR::exp_g(Ne, t-1, d) + k) * (1 - 2*d*purgeR::exp_F(Ne, t-1))
}
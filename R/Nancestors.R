#' Population founders and ancestors
#'
#' Estimate the total and effective number of founders and ancestors
#' in a pedigree, as well as the number of founder genome equivalents
#' (see details on these parameters below).
#' Note that a reference population (RP) must be defined, so that founders
#' and ancestors are referred to the set of individuals belonging to
#' that RP. This is set by means of a boolean vector passed as argument.
#'
#' The total number of founders (\emph{Nf}) and ancestors (\emph{Na})
#' are calculated simply as the count of founders and ancestors of
#' individuals belonging to the reference population (RP).
#' Founders here are defined as individuals with both parentals unknown.
#'
#' The effective number of founders (\emph{Nfe}) is the number of equally
#' contributing founders, that would account for observed genetic diversity
#' in the RP, while the effective number of ancestors (\emph{Nae})
#' is defined as the minimum number of ancestors, founders or not, required
#' to account for the genetic diversity observed in the RP.
#' These parameters are computed from the probability of gene origin,
#' following methods in Tahmoorespur and Sheikhloo (2011).
#'
#' While the previous parameters account for diversity loss due to bottlenecks
#' at the level of founders or ancestors, other sources of random loss of
#' alleles (such as drift) can be accounted by means of the number of
#' founder genome equivalents (\emph{Ng}, Caballero and Toro 2000).
#' This parameter is estimated via Monte Carlo simulation of allele segregation,
#' following Boichard et al. (1997).
#'
#' @name pop_Nancestors
#' @aliases pop_Nf pop_Nfe pop_Na pop_Nae pop_Ng
#' @template ped-arg
#' @template reference-arg
#' @param nboot Number of bootstrap iterations (for computing \emph{Ng}).
#' @template seed-arg
#' @param skip_Ng Skip \emph{Ng} computation or not (FALSE by default).
#' @return A dataframe containing population size estimates for founders and ancestors:
#' \itemize{
#'   \item{\emph{Nr} - } Total number of individuals in the RP
#'   \item{\emph{Nf} - } Total number of founders
#'   \item{\emph{Nfe} - } Effective number of founders
#'   \item{\emph{Na} - } Total number of ancestors
#'   \item{\emph{Nae} - } Effective number of ancestors
#'   \item{\emph{Ng} - } Number of founder genome equivalents
#'   \item{\emph{se_Ng} - } Standard error of Ng
#' }
#' If some of the auxiliary functions is used (e.g. \emph{pop_Nr}), only the corresponding numerical estimate will be returned.
#' In the case of \emph{pop_Ng}, a list object is returned, with the number of founder genome equivalents (Ng) and its standard error (se_Ng).
#' @references
#' \itemize{
#'   \item{Boichard D, Maignel L, Verrier E. 1997. The value of using probabilities of gene origin to measure genetic variability in a population. Genet. Sel. Evol. 29: 5-23.}
#'   \item{Caballero A, Toro M. 2000. Interrelations between effective population size and other pedigree tools for the management of conserved populations. Genet. Res. 75: 331-343.}
#'   \item{Tahmoorespur M, Sheikhloo M. 2011. Pedigree analysis of the closed nucleus of Iranian Baluchi sheep. Small Rumin. Res. 99: 1-6.}
#'}
#' @examples
#' data(arrui)
#' pop_Nancestors(arrui, reference = "target", skip_Ng = TRUE)
#' @export
pop_Nancestors <- function(ped, reference, nboot = 10000L, seed = NULL, skip_Ng = FALSE) {
  check_basic(ped, "id", "dam", "sire")
  ref <- check_reference(ped, reference, "reference")
  check_int(nboot)
  if (nboot < 100 || nboot > 100000000) stop("Number of bootstrap iterations expected in the range[100, 100000000]")
  if (!base::is.null(seed)) check_int(seed)
  check_bool(skip_Ng)
  rp <- ped[which(ped[reference] == TRUE & !is.na(ped[reference])), ]$id
  .Call(`_purgeR_ancestors`, ped, ref, rp - 1, nboot, seed, skip_Ng)
}

#' @export
pop_Nf <- function(ped, reference) {
  purgeR::pop_Nancestors(ped, reference, skip_Ng = TRUE)$Nf
}

#' @export
pop_Nfe <- function(ped, reference) {
  purgeR::pop_Nancestors(ped, reference, skip_Ng = TRUE)$Nfe
}

#' @export
pop_Na <- function(ped, reference) {
  purgeR::pop_Nancestors(ped, reference, skip_Ng = TRUE)$Na
}

#' @export
pop_Nae <- function(ped, reference) {
  purgeR::pop_Nancestors(ped, reference, skip_Ng = TRUE)$Nae
}

#' @export
pop_Ng <- function(ped, reference, nboot = 10000L, seed = NULL) {
  Ng_ <- purgeR::pop_Nancestors(ped, reference, nboot, seed, skip_Ng = FALSE)
  list(
    Ng = Ng_$Ng,
    se_Ng = Ng_$se_Ng
  )
}

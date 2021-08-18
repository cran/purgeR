#' Opportunity of purging
#'
#' The potential reduction in individual inbreeding load can be
#' estimated by means of the opportunity of purging (\emph{O}) and expressed
#' opportunity of purging (\emph{Oe}) parameters described by Gulisija
#' and Crow (2007). Whereas \emph{O} relates to the total potential reduction
#' of the inbreeding load in an individual, as a consequence of it having
#' inbred ancestors, \emph{Oe} relates to the expressed potential reduction of the
#' inbreeding load. Only \emph{Oe} is computed by default. Estimates of \emph{O}
#' and \emph{Oe} need to be corrected in complex pedigrees (see Details below).
#' Both corrected (named "O" and "Oe" by default), and non-corrected (suffixed
#' with "_raw") are returned.
#' 
#' Model used here assume fully recessive, high effect size alleles (Gulisija and Crow, 2007).
#' 
#' In simple pedigrees, the opportunity of purging (\emph{O}) and the expressed
#' opportunity of purging (\emph{Oe}) are estimated as in Gulisija and Crow (2007).
#' For complex pedigrees involving more than one autozygous individual per
#' path from an individual to an ancestor, \emph{O} and \emph{Oe} in the closer ancestors
#' need to be discounted for what was already accounted for in their predecessors.
#' To solve this problem, Gulisija and Crow (2007) provide expression to 
#' correct \emph{O} and \emph{Oe} (see equations 21 and 22 in the manuscript). 
#' 
#' Here, an heuristic approach is used to prevent the inflation of \emph{O} and \emph{Oe},
#' and avoid the use of additional looped expressions that may result in an
#' excessive computational cost. To do so, only the contribution of the most recent
#' ancestors in a path will be considered. Specifically, the method skips contributions
#' from "far" ancestors \emph{k}, such that \emph{Fj(k) > 0}, where \emph{j} is an intermediate ancestor,
#' both referred to an individual \emph{i} of interest. \emph{Fj(k)} refers to the partial
#' inbreeding of \emph{j} for alleles derived from \emph{k} (see \code{\link{ip_Fij}}).
#' This may not provide exact values of \emph{O} and \emph{Oe}, but we expect little bias, since
#' more distant ancestors also contribute lesser to \emph{O} and \emph{Oe}.
#' 
#' Both types of estimates (corrected and non-corrected) are returned (non-corrected
#' estimates, prefixed with "_raw").
#' 
#' @template ped-arg
#' @param name_Oe A string naming the new output column for the expressed opportunity of purging (defaults to "Oe")
#' @param compute_O Enable computation of total opportunity of purging (disabled by default)
#' @param name_O A string naming the new output column for total opportunity of purging (defaults to "O")
#' @template Fcol-arg
#' @template ncores-arg
#' @param genedrop Number of genedrop iterations run to compute partial inbreedng. If set to zero (as default), exact coefficients are computed.
#' @param seed Sets a seed for the random number generator (only if genedrop is enabled).
#' @param complex Enable correction for complex pedigrees (deprecated in v1.3, both raw and corrected measures of "Oe" are returned now).
#' @return The input dataframe, plus an additional column containing \emph{Oe} and \emph{Oe_raw} estimates (additional columns for \emph{O} can appended by enabling \code{compute_O = TRUE}).
#' @seealso \code{\link{ip_Fij}}
#' @encoding UTF-8
#' @references
#' \itemize{
#'   \item{Gulisija D, Crow JF. 2007. Inferring purging from pedigree data. Evolution 61(5): 1043-1051.}
#' }
#' @examples
#' # Original pedigree file in Gulisija & Crow (2007)
#' pedigree <- tibble::tibble(
#'   id = c("M", "K", "J", "a", "c", "b", "e", "d", "I"),
#'   dam = c("0", "0", "0", "K", "M", "a", "c", "c", "e"),
#'   sire = c("0", "0", "0", "J", "a", "J", "b", "b", "d")
#' )
#' pedigree <- purgeR::ped_rename(pedigree, keep_names = TRUE)
#' ip_op(pedigree, compute_O = TRUE)
#' @export
ip_op <- function(ped, name_Oe = "Oe", compute_O = FALSE, name_O = "O", Fcol = NULL, ncores = 1L, genedrop = 0, seed = NULL, complex = NULL) {
  check_basic(ped, "id", "dam", "sire")
  F_ <- check_Fcol(ped, Fcol)
  check_not_col(base::colnames(ped), name_O)
  check_not_col(base::colnames(ped), name_Oe)
  check_bool(compute_O)
  if (!base::is.null(complex)) base::warning("Argument 'complex' was deprecated in v1.3. It will be ignored.")

  pkm <- ip_Fij(ped, mode = "all", Fcol = Fcol, ncores = ncores, genedrop = genedrop, seed = seed)
  .Call(`_purgeR_op`, ped, pkm, Fcol = F_, name_O, name_Oe, "_raw", compute_O)
}

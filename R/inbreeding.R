#' Inbreeding coefficient
#'
#' Computes the standard inbreeding coefficient (\emph{F}).
#' This is the probability that two alleles on a locus are identical by descent (Falconer and Mackay 1996, Wright 1922), calculated from the genealogical coancestry matrix (Malécot 1948).
#' 
#' @template ped-arg
#' @template nameto-arg
#' @return The input dataframe, plus an additional column with individual inbreeding coefficient values (named "Fi" by default).
#' @encoding UTF-8
#' @references
#' \itemize{
#'   \item{Falconer DS, Mackay TFC. 1996. Introduction to Quantitative Genetics. 4th edition. Longman, Essex, U.K.}
#'   \item{Malécot G, 1948. Les Mathématiques de l’hérédité. Masson & Cie., Paris.}
#'   \item{Wright S. 1922. Coefficients of inbreeding and relationship. The American Naturalist 56: 330-338.}
#'}
#' @seealso \code{\link{exp_F}}
#' @examples
#' data(dama)
#' ip_F(dama)
#' @export
ip_F <- function(ped, name_to = "Fi") {
  check_basic(ped)
  check_not_col(base::colnames(ped), name_to)
  .Call(`_purgeR_F`, ped, name_to)
}

#' Ancestral inbreeding coefficient
#'
#' Computes the ancestral inbreeding coefficient (\emph{Fa}).
#' This is the probability that an allele has been in homozygosity in at least one ancestor (Ballou 1997).
#' A genedrop approach is included to compute unbiased estimates of \emph{Fa} (Baumung et al. 2015).
#' 
#' @template ped-arg
#' @template nameto-arg
#' @param genedrop Number of genedrop iterations to run. If set to zero (as default), Ballou's Fa is computed.
#' @template seed-arg
#' @template Fcol-arg
#' @return The input dataframe, plus an additional column with individual ancestral inbreeding coefficient values (named "Fa" by default).
#' @references
#' \itemize{
#'   \item{Ballou JD. 1997. Ancestral inbreeding only minimally affects inbreeding depression in mammalian populations. J Hered. 88:169–178.}
#'   \item{Baumung et al. 2015. GRAIN: A computer program to calculate ancestral and partial inbreeding coefficients using a gene dropping approach. Journal of Animal Breeding and Genetics 132: 100-108.}
#' }
#' @seealso \code{\link{ip_F}}, \code{\link{exp_Fa}}
#' @examples
#' data(dama)
#' ip_Fa(dama) # Compute F on the go (won't be kept in the pedigree).
#' dama <- ip_F(dama)
#' dama <- ip_Fa(dama, Fcol = 'Fi') # If F is computed in advance.
#' @export
ip_Fa <- function(ped, name_to = "Fa", genedrop = 0, seed = NULL, Fcol = NULL) {
  check_basic(ped)
  check_int(genedrop)
  if (!base::is.null(seed)) check_int(seed)
  F_ <- check_Fcol(ped, Fcol)
  check_not_col(base::colnames(ped), name_to)
  .Call(`_purgeR_Fa`, ped, F_, name_to, genedrop, seed)
}

#' Purged inbreeding coefficient
#'
#' Computes the purged inbreeding coefficient (\emph{g}).
#' This is the probability that two alleles on a locus are identical by descent,
#' but relative to deleterious recessive alleles (García-Dorado 2012). The reduction
#' in \emph{g} relative to standard inbreeding (\emph{F}) is given by an effective purging
#' coefficient (\emph{d}), that measures the strength of the deleterious recessive
#' component in the genome. The coefficient \emph{g} is computed following the methods
#' for pedigrees in García-Dorado (2012) and García-Dorado et al. (2016).
#' 
#' @template ped-arg
#' @template d-arg
#' @template nameto-arg
#' @template Fcol-arg
#' @return The input dataframe, plus an additional column containing purged inbreeding coefficient values (named "g" and followed by the purging coefficient value by default).
#' @seealso \code{\link{ip_F}} \code{\link{exp_g}}
#' @encoding UTF-8
#' @references
#' \itemize{
#'   \item{García-Dorado. 2012. Understanding and predicting the fitness decline of shrunk populations: Inbreeding, purging, mutation, and standard selection. Genetics 190: 1-16.}
#'   \item{García-Dorado et al. 2016. Predictive model and software for inbreeding-purging analysis of pedigreed populations. G3 6: 3593-3601.}
#' }
#' @examples
#' data(dama)
#' ip_g(dama, d = 0.23)
#' @export
ip_g <- function(ped, d, name_to = "g<d>", Fcol = NULL) {
  check_basic(ped, "id", "dam", "sire")
  check_d(d)
  F_ <- check_Fcol(ped, Fcol)
  d_label <- base::format(d, scientific = FALSE)
  d_label <- paste("g", d_label, sep = "")
  if (name_to != "g<d>") d_label <- name_to
  check_not_col(base::colnames(ped), d_label)
  .Call(`_purgeR_g`, ped, d, F_, d_label)
}

#' Sort individuals (with ancestors on top of descendants)
#'
#' Individuals can be sorted according to the pedigree structure, without need of birth dates.
#' In the sorted pedigree, descendants will always be placed in rows with higher index number
#' than that of their ancestors. This way, individuals born first will tend to be in the top
#' of the pedigree. Younger individuals, and individuals with no descendants will tend to be
#' placed at the bottom.
#' This function uses the sorting algorithm developed by Zhang et al (2009).
#' After sorting, individuals will be renamed from 1 to N using \code{\link{ped_rename}}.
#' 
#' @template ped-arg
#' @template id-arg
#' @template dam-arg
#' @template sire-arg
#' @template keepnames-arg
#' @return A sorted pedigree dataframe (with ancestors on top of descendants).
#' @references
#' \itemize{
#'   \item{Zhang Z, Li C, Todhunter RJ, Lust G, Goonewardene L, Wang Z. 2009. An algorithm to sort complex pedigrees chronologically without birthdates. J Anim Vet Adv. 8 (1): 177-182.}
#' }
#' @seealso \code{\link{ped_rename}}
#' @examples
#' data(darwin)
#' # Here we reshuffle rows in the pedigree. It won't be usable for other functions in the package
#' darwin <- darwin[sample(1:nrow(darwin)), ]
#' # Below, we sort the pedigree again. The order might not be the same as before.
#' # But ancestors will always be placed on top of descendants,
#' # making the pedigree usable for other functions in the package.
#' darwin <- ped_sort(darwin, id = "Individual", dam = "Mother", sire = "Father", keep_names = TRUE)
#' @export
ped_sort <- function(ped, id = "id", dam = "dam", sire = "sire", keep_names = FALSE) {
  check_basic(ped, id_name = id, dam_name = dam, sire_name = sire, when_rename = TRUE, when_sort = TRUE)
  p <- ped
  N <- base::nrow(p)
  S <- base::rep(x = TRUE, times = N)
  G <- base::rep(x = 0, times = N)
  t <- ped[NULL]
  t_G <- base::integer()
  ped <- sort_step(p, id, dam, sire, t, S, G, t_G)
  check_nrows(ped, N)
  purgeR::ped_rename(ped, id, dam, sire, keep_names = keep_names)
}

#' Sorting steps
#'
#' Recursive function that computes steps for sorting algorithm described by Zhang et al (2009).
#'
#' @name ped_sort_i
#' @param p Pedigree to sort (used as template)
#' @param id A string naming the column with individual identities. It will be renamed to its default value 'id'.
#' @param dam A string naming the column with maternal identities. It will be renamed to its default value 'dam'.
#' @param sire A string naming the column with paternal identities. It will be renamed to its default value 'sire'.
#' @param t Template for the new sorted pedigree
#' @param S Vector of assumed parent individuals
#' @param G Vector of generation numbers (0 identifies the youngest)
#' @param t_G Vector G for the new sorted pedigree
#' @template check-return
#' @seealso \code{\link{ped_sort}}
#' @references
#' \itemize{
#'   \item{Zhang Z, Li C, Todhunter RJ, Lust G, Goonewardene L, Wang Z. 2009. An algorithm to sort complex pedigrees chronologically without birthdates. J Anim Vet Adv. 8 (1): 177-182.}
#' }
#' @return Filled template for the sorted pedigree. Once recursion ends, it returns the sorted pedigree
sort_step <- function(p, id, dam, sire, t, S, G, t_G) {
  in_dam <- base::ifelse(p[[id]] %in% p[[dam]], TRUE, FALSE)
  in_sire <- base::ifelse(p[[id]] %in% p[[sire]], TRUE, FALSE)
  S <- base::ifelse(!in_dam & !in_sire, FALSE, S)
  G <- base::ifelse(!in_dam & !in_sire, G, G+1)
  pns <- p[!S, ]
  p <- p[S, ]
  t <- base::rbind(t, pns)
  t_G <- base::c(t_G, G[!S])
  if (base::nrow(p) == 0) {
    t[base::order(t_G, decreasing = TRUE), ]
  } else {
    G <- G[S]
    S <- S[S]
    sort_step (p, id, dam, sire, t, S, G, t_G)
  }
}
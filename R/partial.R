#' Partial inbreeding coefficient
#'
#' Computes partial inbreeding coefficients, \emph{Fi(j)}.
#' A coefficient \emph{Fi(j)} can be read as the probability of individual \emph{i} being
#' homozygous for alleles derived from ancestor \emph{j}.
#' It is calculated following the tabular method described by Gulisija & Crow (2007).
#' Optionally, it can be estimated via genedrop simulation.
#' 
#' @template ped-arg
#' @param mode Defines the set of ancestors considered when computing partial inbreeding. It can be set as:
#' "founder" for inbreeding conditional to founders only (default),
#' "all" for all individuals in the pedigree (it may take long to compute in large pedigrees),
#' and "custom" for individuals identities given in a integer vector (see 'ancestors' argument).
#' @param ancestors Under the "custom" run mode, it defines a vector of ancestors that will be considered
#' when computing partial inbreeding values.
#' @template Fcol-arg
#' @param genedrop Number of genedrop iterations to run. If set to zero (as default), exact coefficients are computed.
#' @param seed Sets a seed for the random number generator (only if genedrop is enabled).
#' @template ncores-arg
#' @return A matrix of partial inbreeding coefficients. \emph{Fi(j)} values can thus be read from row i and column j.
#' In the resultant matrix, there are as many rows as individuals in the pedigree, and as many columns as ancestors used.
#' Columns will be named and sorted by ancestor identity.
#' @references
#' \itemize{
#'   \item{Gulisija D, Crow JF. 2007. Inferring purging from pedigree data. Evolution 61(5): 1043-1051.}
#' }
#' @seealso \code{\link{ip_F}}
#' @examples
#' # Original pedigree file in Gulisija & Crow (2007)
#' pedigree <- tibble::tibble(
#'   id = c("M", "K", "J", "a", "c", "b", "e", "d", "I"),
#'   dam = c("0", "0", "0", "K", "M", "a", "c", "c", "e"),
#'   sire = c("0", "0", "0", "J", "a", "J", "b", "b", "d")
#' )
#' pedigree <- purgeR::ped_rename(pedigree, keep_names = TRUE)
#'
#' # Partial inbreeding relative to founder ancestors
#' m <- ip_Fij(pedigree)
#' # Note that in the example above, the sum of the values in
#' # rows will equal the vector of inbreeding coefficients
#' # i.e. base::rowSums(m) equals purgeR::ip_F(pedigree)$Fi
#' 
#' # Compute partial inbreeding relative to an arbitrary ancestor
#' # with id = 3 (i.e. individual named "J")
#' anc <- as.integer(c(3))
#' m <- ip_Fij(pedigree, mode = "custom", ancestors = anc)
#' @export
ip_Fij <- function(ped, mode = "founders", ancestors = NULL, Fcol = NULL, genedrop = 0,  seed = NULL, ncores = 1L) {
  
  # Checks and Initialize mode
  check_basic(ped, "id", "dam", "sire")
  check_length(mode)
  if (!mode %in% c("all", "custom", "founders")) stop("Unknown mode value. Select one of 'founders', 'all', 'custom'", call. = FALSE)
  if (mode == "custom") {
    if (base::is.null(ancestors)) stop("mode='custom' requires a vector of ancestors ('ancestors' argument)", call. = FALSE)
    check_ancestors(ped[["id"]], ancestors)
  }
  check_int(ncores)
  check_int(genedrop)
  if (!base::is.null(seed)) check_int(seed)
  
  # Set inbreeding
  F_ <- check_Fcol(ped, Fcol)
  
  # Set ancestors
  N <- base::nrow(ped)
  if (mode == "founders") {
    ancestors_ids <- ped[ped[["dam"]] == 0 & ped[["sire"]] == 0, ][["id"]]
    if (!base::is.null(ancestors)) stop("Vector of ancestors should only be given with argument mode='custom'", call. = FALSE)
  } else if (mode == "all") {
    ancestors_ids <- ped[["id"]];
    if (!base::is.null(ancestors)) stop("Vector of ancestors should only be given with argument mode='custom'", call. = FALSE)
  } else if (mode == "custom") {
    if (!base::is.null(ancestors)) {
      ancestors_ids <- ancestors
    } else stop("mode='custom' requires a vector of ancestors ('ancestors' argument)", call. = FALSE)
  } else stop("Unknown mode value. Select one of 'founders', 'all', 'custom'", call. = FALSE)

  # Index and map ancestors
  ancestors_idx <- idx_ancestors(ancestors_ids, N)
  mapa <- map_ancestors(ped, ancestors_idx)

  # Compute partial inbreeding
  pi <- Fij_core(ped, ancestors_ids, ancestors_idx, F_, mapa, ncores, genedrop, seed)
  pi
}

#' Index ancestors
#'
#' Creates a vector of length N (the number of individuals)
#' Only coordinates for valid ancestors will be given
#' 
#' @param ids Ancestor identities
#' @param N Total number of individuals
#' @return A logical matrix.
idx_ancestors <- function(ids, N) {
  M <- base::length(ids)
  idx <- base::integer(N)
  idx[ids] <- 1:M
  base::ifelse(idx == 0, NA, idx)
}

#' Map ancestors
#'
#' Creates a logical matrix that indicates whether an individual i (in columns) is ancestor of other j (in rows)
#' For example, matrix[, 1] will indicate descendants of id = 1
#' And matrix[1, ] indicates ancestors of id = 1
#' 
#' @template ped-arg
#' @param idx Index of ancestors to map
#' @return A logical matrix.
map_ancestors <- function (ped, idx) {
  N <- base::nrow(ped)
  M <- base::sum(!is.na(idx))
  map <- base::matrix(data = FALSE, nrow = N, ncol = N)
  for (i in 1:N) {
    if (ped[["dam"]][i] > 0) {
      map[i, ped[["dam"]][i]] <- TRUE
      map[i, ][1:ped[["dam"]][i]] <- base::ifelse(map[ped[["dam"]][i], ][1:ped[["dam"]][i]], TRUE, map[i, ][1:ped[["dam"]][i]])
    }
    if (ped[["sire"]][i] > 0) {
      map[i, ped[["sire"]][i]] <- TRUE
      map[i, ][1:ped[["sire"]][i]] <- base::ifelse(map[ped[["sire"]][i], ][1:ped[["sire"]][i]], TRUE, map[i, ][1:ped[["sire"]][i]])
    }
  }
  base::as.matrix(map[, !base::is.na(idx)])
}

#' Partial inbreeding coefficient (core function)
#'
#' Computes partial inbreeding coefficients, Fi(j).
#' A coefficient Fi(j) can be read as the probability of individual i being
#' homozygous for alleles derived from ancestor j
#' 
#' @name Fij_core
#' @template ped-arg
#' @param ancestors Vector of the identities to be assumed as founder ancestors.
#' @param ancestors_idx Index of ancestors.
#' @param Fi Vector of inbreeding coefficients.
#' @param mapa Map of ancestors
#' @template ncores-arg
#' @param genedrop Enable genedrop simulation
#' @template seed-arg
#' @return A matrix of partial inbreeding coefficients. Fi(j) values can thus be read from row i and column j.
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @import progress
Fij_core <- function(ped, ancestors, ancestors_idx, Fi, mapa, ncores = 1, genedrop, seed) {

  N <- base::nrow(ped)
  M <- base::length(ancestors)

  # Save partial kinship matrix in a hashmap
  pi <- base::matrix(data = 0.0, nrow = N, ncol = M)
  pb <- progress::progress_bar$new(format = "[:bar] :percent",
                                   total = M,
                                   width = 56,
                                   complete = '*',
                                   current = '*',
                                   clear = FALSE,
                                   show_after = 0)

  max_ncores <- parallel::detectCores() - 1
  pb$message("Computing partial kinship matrix. This may take a while.")
  pb$tick(0)
  if (ncores < 1) stop ("The minimum number of cores is 1", call. = FALSE)
  else if (ncores > max_ncores) stop(paste("It seems that your system has access to ", max_ncores+1, " cores. Use is limited to ", max_ncores , sep = ""), call. = FALSE)
  else if (ncores == 1) {
    for (i in 1:M) {
      pb$tick(0)
      map_i <- mapa[, i]
      pi[, i] <- Fij_core_i_cpp(ped[["dam"]], ped[["sire"]], ancestors[i] - 1, map_i, Fi, genedrop, seed)
      pb$tick()
    }
  } else {
    cl <- parallel::makeCluster(ncores)
    doSNOW::registerDoSNOW(cl)
    progress <- function(n) {
      pb$tick()
      #pb$tick(tokens = list(ancestor = ancestors[n]))
    } 
    opts <- list(progress = progress)
    pi <- foreach::foreach (i = 1:M,
                            .combine = "cbind",
                            .packages = c("purgeR"),
                            .export = c("Fij_core_i_cpp"),
                            .noexport = c("_purgeR_Fij_core_i_cpp"),
                            .options.snow = opts
                            ) %dopar% {
      map_i <- mapa[, i]
      pi[, i] <- Fij_core_i_cpp(ped[["dam"]], ped[["sire"]], ancestors[i] - 1, map_i, Fi, genedrop, seed)
    }
    parallel::stopCluster(cl)
  }
  base::colnames(pi) <- ancestors
  pi
}
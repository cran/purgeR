#' Check basic
#'
#' This function will group some other checking functions, that should be run everytime when using
#' functions in this package, to avoid unexpected errors.
#' 
#' @name check_basic
#' @template ped-arg
#' @inheritParams check_names
#' @param when_rename True when called from ped_rename function. It softs checks on individual ID column name and types
#' @param when_sort True when called from ped_sort function. It softs checks on pedigree sorting
#' @template check-return
check_basic <- function(ped, id_name = "id", dam_name = "dam", sire_name = "sire", when_rename = FALSE, when_sort = FALSE) {
  check_df(ped)
  check_names(ped, id_name = id_name, dam_name = dam_name, sire_name = sire_name); # if id, dam and sire columns exist
  id_vector <- base::as.vector(base::unlist(ped[[id_name]]))
  dam_vector <- base::as.vector(base::unlist(ped[[dam_name]]))
  sire_vector <- base::as.vector(base::unlist(ped[[sire_name]]))
  check_zero_id(id_vector); # ids named 0
  check_repeat_id(id_vector); # repeated ids
  check_order(id = id_vector, dam = dam_vector, sire = sire_vector, when_sort); # parents declared before descendants
  #check_selfing(copy_ped, id, dam, sire); // selfing and id==(dam||sire)
  if (!when_rename) {
    check_types(id_vector,  dam_vector, sire_vector);
    check_index(id_vector);
  }
}

#' Check that mandatory column names are included
#'
#' Columns for id, dam and sire are mandatory. This function checks that they are named in the pedigree.
#' The function works with arbitrary column names (not 'id', 'dam' and 'sire') to work with ped_rename()
#' 
#' @name check_names
#' @template ped-arg
#' @param id_name Column name for individual id.
#' @param dam_name Column name for dam.
#' @param sire_name Column name for sire.
#' @template check-return
check_names <- function(ped, id_name = "id", dam_name = "dam", sire_name = "sire") {
  check_length(id_name)
  check_length(dam_name)
  check_length(sire_name)
  names <- base::colnames(ped)
  if (!id_name %in% names) stop("Mandatory column 'id' not found", call. = FALSE)
  if (!dam_name %in% names) stop("Mandatory column 'dam' not found", call. = FALSE)
  if (!sire_name %in% names) stop("Mandatory column 'sire' not found", call. = FALSE)
}

#' Check that mandatory column names are of type int
#'
#' Columns for id, dam and sire are mandatory, and required to be of type integer
#' 
#' @name check_types
#' @param id Vector of individual ids.
#' @param dam Vector of dam ids.
#' @param sire Vector of sire ids.
#' @template check-return
check_types <- function(id, dam, sire) {
  if (base::class(id) != "integer") stop("Mandatory 'id', 'dam' and 'sire' columns need to be of type integer", call. = FALSE)
  else if (base::class(dam) != "integer") stop("Mandatory 'id', 'dam' and 'sire' columns need to be of type integer", call. = FALSE)
  else if (base::class(sire) != "integer") stop("Mandatory 'id', 'dam' and 'sire' columns need to be of type integer", call. = FALSE)
}

#' Check individuals named zero
#'
#' Individual id cannot equal zero (0). This is reserved to dams and sires.
#' 
#' @name check_zero_id
#' @param id Vector of individual ids.
#' @template check-return
check_zero_id <- function(id) {
  if (base::any(id == 0)) stop("Individual IDs cannot be set to 0 (it is reserved to unknown parents)", call. = FALSE)
}

#' Check repeated ids
#'
#' Individual id are unique and cannot be repeated
#' 
#' @name check_repeat_id
#' @param id Vector of individual ids.
#' @template check-return
check_repeat_id <- function(id) {
  if(base::anyDuplicated(id)) stop("There are repeated individual IDs", call. = FALSE)
}

#' Check individual order
#'
#' Individuals must be sorted from older to younger
#'  
#' @name check_order
#' @param id Vector of individual ids.
#' @param dam Vector of dam ids.
#' @param sire Vector of sire ids.
#' @param soft_sorting If TRUE checking is relaxed, allowing descendants to be declared before ancestors
#' @template check-return
check_order <- function(id, dam, sire, soft_sorting = FALSE) {
  if (base::any(base::is.na(id))) {
    stop("Individual ids cannot contain missing values (NA)", call. = FALSE)
  }
  if (base::any(base::ifelse(id[!is.na(dam)] == dam[!is.na(dam)], TRUE, FALSE)) | base::any(base::ifelse(id[!is.na(sire)] == sire[!is.na(sire)], TRUE, FALSE))) {
    stop("Individuals cannot be born from themselves!", call. = FALSE)
  }
  if (!soft_sorting) {
    N <- base::length(id)
    idx <- 1:N
    idx_dam <- base::match(id, dam)
    idx_sire <- base::match(id, sire)
    if (base::any(base::ifelse(base::is.na(idx_dam), FALSE, idx_dam <= idx)) | base::any(base::ifelse(base::is.na(idx_sire), FALSE, idx_sire <= idx))) {
      stop("Dams and sires must be declared before their offspring!", call. = FALSE)
    }
  }
}

#' Check individual index
#'
#' Renamed individuals must be named by their index (from 1 to N)
#' 
#' @name check_index
#' @param id Column of individual ids.
#' @template check-return
check_index <- function(id) {
  N <- base::length(id)
  if(!base::identical(id, 1:N)) stop("Individuals must be named from 1 to N", call. = FALSE)
}

#' Check pedigree class
#'
#' The pedigree must be of object class 'data.frame'.
#' 
#' @name check_df
#' @param obj Object to test
#' @template check-return
check_df <- function(obj) {
  if (!base::is.data.frame(obj)) stop("Pedigree must be of class 'data.frame'", call. = FALSE)
}

#' Check that optional column is included
#'
#' Some functions require additional columns. Check that they are named in the pedigree.
#' 
#' @name check_col
#' @param names Column names (all)
#' @param name Column name to check.
#' @template check-return
check_col <- function(names, name) {
  check_length(name, "Expected one column name, but more detected")
  if (!name %in% names) stop(paste("Column not found:", name, sep = " "), call. = FALSE)
}

#' Check if optional column is included
#'
#' Some functions require additional columns. Check if they are already named in the pedigree.
#' 
#' @name check_not_col
#' @inheritParams check_col
#' @template check-return
check_not_col <- function(names, name) {
  check_length(name, "Expected one column name, but more detected")
  if (name %in% names) warning(paste("Column already exists (it will be overwritten)", name, sep = " "), call. = FALSE)
}

#' Check if a variable is boolean or not
#'
#' Can be used to test arguments that need to be of logical (boolean) type
#' 
#' @name check_bool
#' @param variable Variable to test
#' @template check-return
check_bool <- function(variable) {
  check_length(variable, "Expected boolean of length 1")
  if (!is.logical(variable)) stop("Expected boolean (TRUE/FALSE) argument.", call. = FALSE)
}

#' Check if a variable is a positive integer or not
#'
#' Can be used to test arguments that need to be integers
#' 
#' @name check_int
#' @param variable Variable to test
#' @template check-return
check_int <- function(variable) {
  check_length(variable, "Expected single integer value")
  if ((!class(variable) %in% c("integer", "numeric")) | (variable < 0)) stop("Expected a positive integer argument value.", call. = FALSE)
}

#' Check if a vector contains NA values
#'
#' Return warning when NA values are present
#' 
#' @name check_na
#' @param variable Variable to test
#' @template check-return
check_na <- function(variable) {
  if (base::all(base::is.na(variable))) stop("Cannot input a vector of NAs", call. = FALSE)
  if (base::any(base::is.na(variable))) warning("NAs can cause unexpected behavior. Remove them", call. = FALSE)
}

#' Check if a variable has length >1
#'
#' Used to test arguments that need to be of length 1
#' 
#' @name check_length
#' @param variable Variable to test
#' @param message Error message to display
#' @template check-return
check_length <- function(variable, message = "Expected value of length 1") {
  if (base::length(variable) > 1) stop(message, call. = FALSE)
}

#' Check observed and expected number of rows
#'
#' Expected and observed number of rows must be equal.
#' 
#' @name check_nrows
#' @param df Dataframe to test
#' @param exp Expected number of rows
#' @param message Error message to display
#' @template check-return
check_nrows <- function(df, exp, message = "Expected value of length 1") {
  check_df(df)
  check_int(exp)
  if (base::nrow(df) != exp) stop("Unexpected number of rows returned. Please contact the developer.", call. = FALSE)
}

#' Check columns with inbreeding values
#'
#' Takes a column name, and checks its use as inbreeding coefficient.
#' It should name a numeric vector, with values in the range [0,1]
#'
#' @name check_Fcol
#' @template ped-arg
#' @template Fcol-arg
#' @param compute Compute inbreeding if Fcol is NULL
#' @return Vector of inbreeding values (if checks are successful)
check_Fcol <- function(ped, Fcol, compute = TRUE) {
  if (base::is.null(Fcol) & compute) {
    F_ <- purgeR::ip_F(ped[, c("id", "dam", "sire")])[["Fi"]]
    return (F_)
  } else if (base::is.null(Fcol)) stop("check_Fcol is designed to return a vector of F values.", call. = FALSE)

  check_col(base::colnames(ped), Fcol)
  if (!class(ped[[Fcol]]) %in% c("numeric")) stop ("Inbreeding needs to be of numeric type.", call. = FALSE)
  else F_ <- ped[[Fcol]]
  if (max(F_, na.rm = TRUE) > 1.0 | min(F_, na.rm = TRUE)  < 0.0)  stop ("Inbreeding needs to be in the range [0,1].", call. = FALSE)
  Fcol <- F_
  check_na(Fcol)
  F_
}

#' Check columns with generation numbers
#'
#' Takes a column name, and checks its use as generation numbers.
#' It should name a numeric vector, with values >= 0.
#'
#' @name check_tcol
#' @template ped-arg
#' @template tcol-arg
#' @param compute Compute generation numbers if tcol is NULL
#' @param force_int Generation numbers must be integers (disabled by default)
#' @return Vector of generation numbers (if checks are successful)
check_tcol <- function(ped, tcol, compute = TRUE, force_int = FALSE) {
  if (is.null(tcol) & compute) {
    t_ <- purgeR::pop_t(ped)[["t"]]
    return (t_)
  }
  check_col(base::colnames(ped), tcol)
  if (!force_int & !class(ped[[tcol]]) %in% c("integer", "numeric")) stop ("Generations need to be of numeric type.", call. = FALSE)
  else if (force_int & !class(ped[[tcol]]) %in% c("integer")) stop ("Generations need to be of integer type.", call. = FALSE)
  else t_ <- ped[[tcol]]
  if (min(t_, na.rm = TRUE)  < 0)  stop ("Generations cannot take negative values.", call. = FALSE)
  tcol <- t_
  check_na(tcol)
  t_
}

#' Check columns with reference individuals
#'
#' Takes a column name, and checks its use as reference.
#' It should name a boolean vector (or coercible to it),
#' with at least one TRUE value.
#'
#' @name check_reference
#' @template ped-arg
#' @template reference-arg
#' @return Vector of reference numbers (if checks are successful)
check_reference <- function(ped, reference) {
  check_length(reference, "Expected one column name, but more detected")
  reference <- base::as.logical(ped[[reference]])
  if (base::all(base::is.na(reference))) stop(paste("Failed to coerce 'reference' values: All NAs.", sep = ""), call. = FALSE)
  else if (base::sum(reference, na.rm = TRUE) == 0) stop(paste("At least one 'reference' value should be TRUE.", sep = ""), call. = FALSE)
  check_na(reference)
  reference
}

#' Check columns with target individuals
#'
#' Takes a column name, and checks its use as target.
#' It should name a boolean vector (or coercible to it),
#' with at least one TRUE value.
#'
#' @name check_reference
#' @template ped-arg
#' @template reference-arg
#' @param target Target column
#' @param variable To be used in printed messages
#' @return Vector of target numbers (if checks are successful)
check_target <- function(ped, reference, target, variable) {
  if (base::is.null(target)) return (target)
  id_ref <- ped[["id"]][ped[[reference]]]
  id_tgt <- ped[["id"]][ped[[target]]]
  target <- base::as.logical(ped[[target]])
  if (base::all(base::is.na(target))) stop(paste("Failed to coerce '", variable, "' values: All NAs.", sep = ""), call. = FALSE)
  else if (base::sum(target, na.rm = TRUE) == 0) stop(paste("At least one '", variable, "' value should be TRUE.", sep = ""), call. = FALSE)
  else if (base::any(id_ref %in% id_tgt)) stop("Cannot use reference individuals as target at the same time.", call. = FALSE)
  if (base::any(id_tgt < base::max(id_ref, na.rm = TRUE))) base::warning("Target individuals should always have 'id' lower than reference individuals.", call. = FALSE)
  check_na(target)
  target
}

#' Check purging coefficient
#'
#' The purging coefficient must be a number between 0 and 0.5
#'
#' @name check_d
#' @param d Purging coefficient (taking values between 0.0 and 0.5).
#' @template check-return
check_d <- function(d) {
  check_length(d, "Expected single numeric value between 0 and 0.5")
  if ((!class(d) %in% c("integer","numeric")) | (d < 0.0) | (d > 0.5)) stop("Expected a numeric value in the range [0, 0.5]", call. = FALSE)
}

#' Check Ne
#'
#' The effective population size (Ne) must be a number higher than 0
#'
#' @name check_Ne
#' @param Ne Effective population size
#' @template check-return
check_Ne <- function(Ne) {
  check_length(Ne)
  if ((!class(Ne) %in% c("integer","numeric")) | (Ne < 1.0)) stop("Expected a numeric value higher or equal than 1", call. = FALSE)
}

#' Check ancestor individuals
#'
#' Takes a column name, and checks its use as target.
#' It should name a boolean vector (or coercible to it),
#' with at least one TRUE value.
#'
#' @name check_ancestors
#' @param id Vector of individual ids.
#' @param ancestors Vector of ancestor ids.
#' @template check-return
check_ancestors <- function(id, ancestors) {
  if (base::class(ancestors) != "integer") {
    if (base::class(ancestors) == "numeric") stop("'ancestors' column needs to be of type integer (not 'numeric', use as.integer())", call. = FALSE)
    else stop("'ancestors' column needs to be of type integer", call. = FALSE)
  }
  check_zero_id(ancestors)
  check_repeat_id(ancestors)
  if (!base::all(ancestors %in% id)) stop("All ancestors must have a valid id", call. = FALSE)
  sorted_ancestors <- base::sort(ancestors)
  if(!base::identical(ancestors, sorted_ancestors)) stop("Ancestor ids must be sorted", call. = FALSE)
  check_na(ancestors)
}
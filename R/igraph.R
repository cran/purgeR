#' Input for igraph
#'
#' Processes a pedigree into a list with two objects, one dataframe of edges,
#' and a dataframe of vertices, which can be used as input for functions of the
#' igraph package.
#' 
#' @template ped-arg
#' @return A list with one dataframe 'edges' and another 'vertices', each following igraph format.
#' 
#' The 'edges' dataframe will contain two columns in addition to the defaults "from" and "to":
#' 1) 'from_parent' indicates whether the vertex from which the edge originates represents a mother ("dam") or a father ("sire").
#' 2) 'to_parent' indicates whether the vertex to which the edge is directed represents a mother ("dam"), father ("sire") or none ("NA").  
#' @seealso \code{\link{ped_rename}}, \code{\link[igraph]{graph_from_data_frame}}
#' @examples
#' data(atlas)
#' atlas_graph <- ped_graph(atlas)
#' G <- igraph::graph_from_data_frame(d = atlas_graph$edges,
#'                                    vertices = atlas_graph$vertices,
#'                                    directed = TRUE)
#' @export
ped_graph <- function(ped) {
  check_basic(ped, "id", "dam", "sire")

  vertices <- ped[, !(base::names(ped) %in% c("dam", "sire"))]
  base::rownames(vertices) <- NULL

  # edges <- ped[ped[["dam"]] != 0 & ped[["sire"]] != 0, ]
  # edges <- dplyr::bind_rows(dplyr::transmute(edges, from = dam, to = id),
  #                           dplyr::transmute(edges, from = sire, to = id))
  filtered_rows_from_dam <- ped[ped[["dam"]] != 0, c("dam", "id")]
  filtered_rows_from_sire <- ped[ped[["sire"]] != 0, c("sire", "id")]
  colnames(filtered_rows_from_dam) <- c("from", "to")
  colnames(filtered_rows_from_sire) <- c("from", "to")
  filtered_rows_from_dam["from_parent"] <- "dam"
  filtered_rows_from_sire["from_parent"] <- "sire"
  edges <- base::rbind(filtered_rows_from_dam, filtered_rows_from_sire)

  to_parent <- base::ifelse(edges[["to"]] %in% ped[["dam"]], "dam",
               base::ifelse(edges[["to"]] %in% ped[["sire"]], "sire", NA))
  edges["to_parent"] <- to_parent

  edges <- edges[base::order(edges[["from"]], edges[["to"]]), ]
  base::rownames(edges) <- NULL
  list(edges = edges, vertices = vertices)
}
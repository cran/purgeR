% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partial.R
\name{map_ancestors}
\alias{map_ancestors}
\title{Map ancestors}
\usage{
map_ancestors(ped, idx)
}
\arguments{
\item{ped}{A dataframe containing the pedigree. Individual (id), maternal (dam), and paternal (sire) identities are mandatory columns.}

\item{idx}{Index of ancestors to map}
}
\value{
A logical matrix.
}
\description{
Creates a logical matrix that indicates whether an individual i (in columns) is ancestor of other j (in rows)
For example, matrix[, 1] will indicate descendants of id = 1
And matrix[1, ] indicates ancestors of id = 1
}

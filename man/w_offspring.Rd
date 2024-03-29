% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitness.R
\name{w_offspring}
\alias{w_offspring}
\title{Offspring}
\usage{
w_offspring(ped, name_to, dam_offspring = TRUE, sire_offspring = TRUE)
}
\arguments{
\item{ped}{A dataframe containing the pedigree. Individual (id), maternal (dam), and paternal (sire) identities are mandatory columns.}

\item{name_to}{A string naming the new output column.}

\item{dam_offspring}{Compute dam offspring (TRUE by default).}

\item{sire_offspring}{Compute sire offspring (TRUE by default).}
}
\value{
The input dataframe, plus an additional column indicating the total number of offspring.
}
\description{
Counts the number of offspring for individuals in the pedigree.
}
\examples{
data(arrui)
dama <- w_offspring(arrui, name_to = "P")
head(arrui)
}

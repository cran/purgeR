% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ne.R
\name{Ne_delta}
\alias{Ne_delta}
\alias{se_Ne_delta}
\title{Realized effective population size (mean)}
\usage{
Ne_delta(delta)

se_Ne_delta(delta)
}
\arguments{
\item{delta}{Vector of individual variations in inbreeding.}
}
\value{
Mean effective population size.

Standard error of the effective population size.
}
\description{
Computes the mean realized effective population size.
Note this function expected a mean delta_F value for all individuals in the reference population

Computes the standard error of the realized effective population size.
Note this function expects the mean and standard deviation of delta F, as well as the total number of individuals in the reference population
}

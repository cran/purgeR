% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check.R
\name{check_nrows}
\alias{check_nrows}
\title{Check observed and expected number of rows}
\usage{
check_nrows(df, exp, message = "Expected value of length 1")
}
\arguments{
\item{df}{Dataframe to test}

\item{exp}{Expected number of rows}

\item{message}{Error message to display}
}
\value{
No return value. Will print an error message if checking fail.
}
\description{
Expected and observed number of rows must be equal.
}

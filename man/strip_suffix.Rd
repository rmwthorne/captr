% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import.R
\name{strip_suffix}
\alias{strip_suffix}
\title{Removes _# if at the end of a string}
\usage{
strip_suffix(x)
}
\arguments{
\item{x}{A character object}
}
\value{
The cleaned version of \code{x}.
}
\description{
Removes _# if at the end of a string
}
\examples{
strip_suffix("MYC_1")
}

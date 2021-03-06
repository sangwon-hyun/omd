% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{fill_in_empty}
\alias{fill_in_empty}
\title{Fill in empty by taking the average of the 4 x 4 window around it.}
\usage{
fill_in_empty(mat, wind = 4)
}
\arguments{
\item{mat}{Matrix.}

\item{wind}{Window size, defaults to 4.}
}
\value{
Matrix of same size as \code{mat}.
}
\description{
Fill in empty by taking the average of the 4 x 4 window around it.
}

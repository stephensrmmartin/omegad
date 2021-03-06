% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/miscHelpers.R
\name{.array_extract}
\alias{.array_extract}
\title{Extract subarray, while maintaining shape.}
\usage{
.array_extract(array, slice)
}
\arguments{
\item{array}{Array to slice up.}

\item{slice}{Which slice to take.}
}
\value{
Array of one less dimensionality. If an array of matrices, then a matrix. If an array of vectors, then a vector. If an array of row vectors, then a row vector. And so on.
}
\description{
Given some array of any dimension, return a slice off the end, in the desired shape.
}
\details{
This function largely solves the one thing in R that has made me ragequit.
When given an array of some dimensions ([N, P, S]), this function will extract out the given slice (e.g., [N,P,2]), and return an [N,P] matrix.
This is true whether there is only ONE row or ONE column.
Please, for the love of all that is holy, just return a slice of an array, so we can do some matrix algebra over several matrices, vectors, or row vectors.
}
\author{
Stephen R. Martin
}
\keyword{internal}

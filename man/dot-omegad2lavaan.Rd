% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reliability.R
\name{.omegad2lavaan}
\alias{.omegad2lavaan}
\title{lavaanify omegad models.}
\usage{
.omegad2lavaan(object)
}
\arguments{
\item{object}{omegad object.}
}
\value{
list containing the Character string of lavaan model, equivalent options, and data frame for use in do.call.
}
\description{
Converts omegad model specification to a similar lavaan model.
}
\author{
Stephen R. Martin
}
\keyword{internal}

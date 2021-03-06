% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rstantools.R
\name{log_lik.omegad}
\alias{log_lik.omegad}
\alias{log_lik}
\title{Compute log likelihoods for each observation (and item).}
\usage{
\method{log_lik}{omegad}(object, items = FALSE, ...)
}
\arguments{
\item{object}{omegad object.}

\item{items}{Logical (Default: FALSE). Whether to compute log_lik per item, per observation (TRUE), or just per observation (FALSE).}

\item{...}{Not used.}
}
\value{
See \code{\link[rstantools]{log_lik}} and \code{\link{loo.omegad}}.
}
\description{
log_lik method for omegad objects.
}
\details{
log_lik method for omegad objects.
}
\seealso{
\code{\link[rstantools]{log_lik}}
}
\author{
Stephen R. Martin
}

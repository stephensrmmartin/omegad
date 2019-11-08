##' Extracts the number of posterior samples tored in a fitted omegad model.
##'
##' @title Extract number of posterior samples.
##' @param object omegad object.
##' @return Integer. How many posterior post-warmup samples were drawn.
##' @author Stephen R. Martin
##' @importFrom rstantools nsamples
##' @export nsamples
##' @export
nsamples.omegad <- function(object){
   n.samps <- (object$fit@sim$iter - object$fit@sim$warmup)*object$fit@sim$chains
   return(n.samps)
}
##' log_lik method for omegad objects.
##'
##' log_lik method for omegad objects.
##' @title Compute log likelihoods for each observation (and item).
##' @param object omegad object.
##' @param items Logical (Default: FALSE). Whether to compute log_lik per item, per observation (TRUE), or just per observation (FALSE).
##' @return See \code{\link[rstantools]{log_lik}} and \code{\link{loo.omegad}}.
##' @author Stephen R. Martin
##' @importFrom rstantools log_lik
##' @export log_lik
##' @export
##' @seealso \code{\link[rstantools]{log_lik}}
log_lik.omegad <- function(object, items = FALSE) {
    loo.omegad(object, items, log_lik = TRUE)
}

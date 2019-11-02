##' Extracts the number of posterior samples tored in a fitted omegad model.
##'
##' @title Extract number of posterior samples.
##' @param object omegad object.
##' @return Integer. How many posterior post-warmup samples were drawn.
##' @author Stephen R. Martin
##' @importFrom rstantools nsamples
##' @keywords internal
##' @export
nsamples.omegad <- function(object){
   n.samps <- (object$fit@sim$iter - object$fit@sim$warmup)*object$fit@sim$chains
   return(n.samps)
}

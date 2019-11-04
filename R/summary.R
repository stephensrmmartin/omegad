##' Print method for omegad object.
##'
##' Prints metadata for omegad object.
##' @title Print method for omegad objects.
##' @param x omegad object.
##' @param ... Not used.
##' @return x (invisible).
##' @author Stephen R. Martin
##' @export
print.omegad <- function(x, ...) {
    cat("Formula: \n")
    lapply(x$formula, function(x){
        cat("\t", deparse(x), "\n")
    })
    cat("Number of Observations: ", x$meta$N, "\n")
    cat("Number of Indicators: ", x$meta$J, "\n")
    cat("Number of Factors: ", x$meta$F, "\n")
    cat("Dependency Model: ", .get_model_description(x), "\n")
    cat("Chains: ", x$fit@sim$chains, "\n")
    cat("Time: ", max(rowSums(rstan::get_elapsed_time(x$fit))), " seconds \n")
    cat("Finished:", x$fit@date, "\n")

    return(invisible(x))
    
}

summary.omegad <- function(object, prob = .95, ...) {
    probs <- .prob_to_probs(prob)
    
}
##' Prints summary.omegad objects.
##'
##' Prints summary.omegad objects.
##' @title Print method for summary.omegad object.
##' @param x summary.omegad object.
##' @param ... Not used.
##' @return x (invisible)
##' @author Stephen R. Martin
##' @export
print.summary.omegad <- function(x, ...) {
    
    return(invisible(x))
}

.get_model_description <- function(object) {
    if (object$meta$gp) {
        out <- c("\n \t [Factor -> Error Factor]",
                 paste0("\t Univariate gaussian process: Additive linear and exponential quadratic kernels (", object$meta$M, " basis functions) ")
                 )
        if (object$meta$exo) {
            out <- c(out[1], "\t [Exogenous -> Error Factor]", out[2])
        }
        out <- paste0(out, collapse = "\n")
    } else {
        out <- c("\n \t [Factors -> Error Factors]", "\t Covariance only")
        if (object$meta$exo) {
            out <- c(out, "\t [Exogenous -> Factor]", "\t Linear model")
        }
        out <- paste0(out, collapse = "\n")
    }

    return(out)
}

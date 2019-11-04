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
        cat(as.character(x), "\n")
    })
    cat("Number of Observations: ", x$meta$N, "\n")
    cat("Number of Indicators: ", x$meta$J, "\n")
    cat("Number of Factors: ", x$meta$F, "\n")
    cat("Dependency Model: ", .get_model_description(x), "\n")
    cat("Chains: ")
    cat("Cores: ")
    cat("Time: ")

    return(invisible(x))
    
}

summary.omegad <- function(object, prob = .95, ...) {
    
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
        out <- c("[Factor -> Error Factor]",
                 paste0("Univariate gaussian process: Additive linear and exponential quadratic kernels (", object$meta$M, " basis functions) ")
                 )
        if (object$meta$exo) {
            out <- c(out[1], "[Exogenous -> Error Factor]", out[2])
        }
        out <- paste0(out, collapse = "\n")
    } else {
        out <- c("[Factors -> Error Factors]", "Covariance only")
        if (object$meta$exo) {
            out <- c(out, "[Exogenous -> Factor]", "Linear model")
        }
        out <- paste0(out, collapse = "\n")
    }

    return(out)
}

##' @title Operator for testing NULL and returning expr if NULL
##' @param object Object to test NULL on.
##' @param expr Expression to evaluate.
##' @return object if not NULL, results of expression if object is NULL
##' @author Stephen R Martin
##' @keywords internal
%IfNull% <- function(object, expr) {
    if(is.null(object)) {
        return(eval(expr))
    } else {
        return(object)
    }
    
}
##' @title Construct named list from arguments.
##' @param ... Objects for list.
##' @return Named list.
##' @author Stephen R Martin
##' @keywords internal
nlist <- function(...) {
    mc <- match.call()
    out <- list(...)

    not_named <- is.null(names(out))
    is_named <- if(not_named) {
                    FALSE
                } else {
                    nzchar(names(out))
                }

    args <- as.character(mc)[-1] # Not the fn name.
    if(not_named) {
        names(out) <- args
    } else {
        names(out)[!is_named] <- args[!is_named]
    }

    out
}

##' Extracts samples from stan fit, converts to R-friendly vectors, matrices, arrays, for each sample.
##'
##' @title Extract and transform samples to R vectors, matrices, arrays.
##' @param object stanfit.
##' @param par Character. Name of parameter to extract and transform.
##' @param nsamples Numeric (Optional). Number of samples to use.
##' @param transpose Logical (Default: FALSE). Only applicable if two dimensional.
##' Intended for use when an array of vectors is used.
##' When TRUE, will transpose such that each row is a vector element, and each column is a vector.
##' @param stan Logical (Default: FALSE). If true, will not repermute; will be in same order as stan provides.
##' @param ... Not used.
##' @return If par is a vector, then a PxS array.
##' If par is matrix, then a RxCxS array.
##' If par is array of matrices, then RxCxA1xA2x...xS array.
##' If par is array of vectors, then a AxRxS array (unless transpose = TRUE).
##' @author Stephen R. Martin
##' @keywords internal
.extract_transform <- function(object, par, nsamples = NULL, transpose = FALSE, stan = FALSE, ...) {
    samps <- as.matrix(object, pars = par)
    if (!is.null(nsamples)) {
        samps <- samps[1:nsamples, , drop = FALSE]
    }
    nsamples <- nrow(samps)

    samps.names <- colnames(samps)
    rex_innerBrackets <- ".*\\[(\\d*(?:,\\d*)*)\\]"
    innerBrackets <- gsub(rex_innerBrackets, "\\1", samps.names)

    innerSplit <- strsplit(innerBrackets, split = ",")
    innerSplit.num <- lapply(innerSplit, as.numeric)
    innerSplit.num <- do.call(rbind, innerSplit.num)
    maxDims <- apply(innerSplit.num, 2, max)
    lengthDims <- length(maxDims)

    structured_samps <- array(t(samps), dim=c(maxDims, nsamples))
    if (lengthDims == 2 & transpose & !stan) {
        structured_samps <- aperm(structured_samps, c(2, 1, 3))
    }
    if (lengthDims > 2 & !stan) {
        structured_samps <- aperm(structured_samps, c((lengthDims - 1):(lengthDims), 1:(lengthDims - 2), lengthDims + 1))
    }
    return(structured_samps)
}

##' Given some array of any dimension, return a slice off the end, in the desired shape.
##'
##' This function largely solves the one thing in R that has made me ragequit.
##' When given an array of some dimensions ([N, P, S]), this function will extract out the given slice (e.g., [N,P,2]), and return an [N,P] matrix.
##' This is true whether there is only ONE row or ONE column.
##' Please, for the love of all that is holy, just return a slice of an array, so we can do some matrix algebra over several matrices, vectors, or row vectors.
##' @title Extract subarray, while maintaining shape.
##' @param array Array to slice up.
##' @param slice Which slice to take.
##' @return Array of one less dimensionality. If an array of matrices, then a matrix. If an array of vectors, then a vector. If an array of row vectors, then a row vector. And so on.
##' @author Stephen R. Martin
##' @keywords internal
.array_extract <- function(array, slice) {
    dims <- dim(array)
    lastDem <- length(dims)
    dimExtract <- lapply(dims[1:(lastDem - 1)], function(x){
        1:x
    })
    args <- c(list(array), dimExtract, slice, drop = FALSE)
    arraySub <- do.call(`[`, args)
    out <- array(arraySub, dim = dims[1:(lastDem - 1)])
    return(out)
}

.prob_to_probs <- function(prob) {
    c((1 - prob)/2, 1 - (1 - prob)/2)
}

.summarize_chains <- function(object, chains, par) {
    samps <- rstan::extract(object$fit, par, permuted = FALSE)
    samps <- samps[,chains,, drop = FALSE]
    rstan::monitor(samps, warmup = 0, digits_summary = 4)
}

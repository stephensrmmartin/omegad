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


.array_extract <- function(array, slice) {
    dims <- dim(array)
    lastDem <- length(dims)
    
}

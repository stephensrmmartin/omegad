##' Generic function for reliability. We include this function for compatibility with lavaan and semTools.
##'
##' semTools includes a lavaan-specific reliability function. When omegad is loaded, reliability is masked.
##' For convenience, we define a generic \code{reliability}, with lavaan methods.
##' @title Generic function for reliability.
##' @param object omegad, or lavaan object.
##' @param ... Passed onto methods.
##' @return See methods.
##' @author Stephen R. Martin
##' @export
##' @seealso \code{\link[semTools]{reliability}}
reliability <- function(object, ...) {
    UseMethod("reliability")
}

##' @describeIn reliability lavaan method.
##' @export
reliability.lavaan <- function(object, ...) {
    semTools::reliability(object, ...)
}

##' @describeIn reliability lavaan.mi method.
##' @export
reliability.lavaan.mi <- function(object, ...) {
    semTools::reliability(object, ...)
 }

##' Computes summary statistics for omega reliabilities.
##'
##' \code{omegad} produces observation-specific omega coefficients (See \code{\link{predict.omegad}} and \code{\link{fitted.omegad}} for extracting these).
##' Users may want to compute summary statistics for the omega coefficients, and to compare these to traditional estimates.
##'
##' @param prob Numeric (Default: .95). The amount of probability mass to include within the credible interval. Default values provide a 95\% credible interval.
##' @return Matrix containing omega summaries for each factor.
##' @inheritParams predict.omegad
##' @describeIn reliability omegad method.
reliability.omegad <- function(object, prob = .95, ...) {
    probs <- .prob_to_probs(prob)

    omgs <- fitted(object, summary = FALSE)[c("omega1","omega2")]
    omg.1.0 <- .extract_transform(object$fit, "omega1_expected")
    omg.2.0 <- .extract_transform(object$fit, "omega2_expected")

    .fun <- function(x){
        m <- mean(x)
        sd <- sd(x)
        ci <- quantile(x, probs)
        L <- ci[1]
        U <- ci[2]
        out <- c(m, sd, L, U)
        names(out) <- c("Mean", "SD", paste0("Q", probs * 100, "%"))
        return(out)
    }

    omgs.mean <- lapply(omgs, function(x) {
        out <- apply(x, 2:3, mean)
        out <- apply(out, 1, .fun)
        return(t(out))
    })
    omgs.sd <- lapply(omgs, function(x) {
        out <- apply(x, 2:3, sd)
        out <- apply(out, 1, .fun)
        return(t(out))
    })

    omg.1.0.sum <- apply(omg.1.0, 1, .fun)
    omg.2.0.sum <- apply(omg.2.0, 1, .fun)
    omg.1.0.sum <- t(omg.1.0.sum)
    omg.2.0.sum <- t(omg.2.0.sum)
    rownames(omg.1.0.sum) <- rownames(omg.2.0.sum) <- unlist(object$meta$fnames$factor)

    out <- list(
            omega1_cond = omg.1.0.sum,
            omega2_cond = omg.2.0.sum,
            omega1_exp = omgs.mean[["omega1"]],
            omega2_exp = omgs.mean[["omega2"]],
            omega1_sd = omgs.sd[["omega1"]],
            omega2_sd = omgs.sd[["omega2"]]
        )
    class(out) <- "reliability.omegad"
    return(out)
}
##' @param x reliability.omegad object.
##' @return x (Invisible).
##' @export
##' @describeIn reliability Print method for reliability.omegad objects.
print.reliability.omegad <- function(x, ...) {
    cat("Omega 1 \n")
    
    cat("[E(Omega | Error_f = 0)] \n")
    print(x$omega1_cond)
    cat("\n")

    cat("[E(Omega_i)] \n")
    print(x$omega1_exp)
    cat("\n")

    cat("[SD(Omega_i)] \n")
    print(x$omega1_sd)
    cat("\n")

    cat("Omega 2 \n")

    cat("[E(Omega | Error_f = 0)] \n")
    print(x$omega2_cond)
    cat("\n")

    cat("[E(Omega_i)] \n")
    print(x$omega2_exp)
    cat("\n")

    cat("[SD(Omega_i)] \n")
    print(x$omega1_sd)

    return(invisible(x))
}
##' Converts omegad model specification to a similar lavaan model.
##'
##' 
##' @title lavaanify omegad models.
##' @param object omegad object.
##' @return list containing the Character string of lavaan model, equivalent options, and data frame for use in do.call.
##' @author Stephen R. Martin
##' @keywords internal
.omegad2lavaan <- function(object) {
    formPieces <- object$meta$fnames
    F <- object$meta$F
    formChars <- lapply(1:F, function(x) {
        paste0(formPieces$factor[[x]], " =~ ", paste0(formPieces$indicator[[x]], collapse = " +"), " \n")
    })
    modelChar <- paste0(unlist(formChars))

    out <- list(model = modelChar, data = object$data)
    return(out)
}
##' @title Compute omegas using lavaan.
##' @param object omegad object.
##' @return Reliability output from \code{\link[semTools]{reliability}}.
##' @author Stephen R. Martin
##' @keywords internal
.lavOmega <- function(object) {
    fun <- get("cfa", asNamespace("lavaan"))
    lavOut <- do.call('fun', .omegad2lavaan(object))
    reliability(lavOut)
}

##' @title Fit similar model in lavaan.
##' @param object omegad object.
##' @return lavaan fit
##' @author Stephen R. Martin
##' @keywords inetrnal
.lavFit <- function(object) {
    fun <- get("cfa", asNamespace("lavaan"))
    lavOut <- do.call('fun', .omegad2lavaan(object))
    return(lavOut)
}

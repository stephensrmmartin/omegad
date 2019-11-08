##' LOO method for omegad objects.
##'
##' LOO method for omegad objects.
##' @title omegad method for computing leave-one-out CV scores.
##' @param x omegad object.
##' @param items Logical (Default: FALSE). Whether to compute LOO per item, per observation (TRUE), or just per observation (FALSE).
##' @param ... 
##' @return LOO object (if items = FALSE). List of LOO objects, one for each item (if items = TRUE).
##' @author Stephen R. Martin
##' @import loo
##' @export loo
##' @export
loo.omegad <- function(x, items = FALSE, ...) {
    dots <- list(...)
    r.LL <- dots$log_lik
    if (is.null(r.LL)) {
        r.LL <- FALSE
    }

    
    S <- nsamples(x)
    N <- x$meta$N
    J <- x$meta$J
    F <- x$meta$F
    xobs <- x$stan_data$x

    theta <- .extract_transform(x$fit, "theta")
    theta_loc <- theta[, 1:F, , drop = FALSE]
    theta_sca <- theta[, (F + 1):(2 * F), , drop = FALSE]
    lambda_loc <- .extract_transform(x$fit, "lambda_loc_mat")
    lambda_sca <- .extract_transform(x$fit, "lambda_sca_mat")
    nu_loc <- t(.extract_transform(x$fit, "nu_loc"))
    nu_sca <- t(.extract_transform(x$fit, "nu_sca"))

    xhat <- array(0, dim = c(N, J, S))
    shat <- array(0, dim = c(N, J, S))
    LL <- array(0, dim = c(N, J, S))

    for (s in 1:S) {
        xhat[,, s] <- matrix(1, nrow = N, ncol = 1) %*% .array_extract(t(nu_loc), s) + .array_extract(theta_loc, s) %*% .array_extract(lambda_loc, s)
        shat[,, s] <- exp(matrix(1, nrow = N, ncol = 1) %*% .array_extract(t(nu_sca), s) + .array_extract(theta_sca, s) %*% .array_extract(lambda_sca, s))
        LL[,, s] <- dnorm(xobs, xhat[,,s], shat[,,s], log = TRUE)
    }

    if (items) {
        log_lik <- lapply(seq_len(J), function(j) {
            out <- t(LL[, J, ])
            out
        })
        ## r_eff <- lapply(log_lik, function(j) {
        ##     out <- relative_eff(exp(j))
        ##     out
        ## })
        looOut <- lapply(seq_len(J), function(j) {
            ## loo(log_lik[[j]], r_eff = r_eff[[j]])
            loo(log_lik[[j]])
        })
        names(looOut) <- colnames(xobs)
    } else {
        log_lik <- apply(LL, c(1, 3), function(x) {
            sum(x)
        })
        log_lik <- t(log_lik)
        ## r_eff <- relative_eff(exp(log_lik))
        ## looOut <- loo(log_lik, r_eff = r_eff)
        looOut <- loo(log_lik)
    }

    if (r.LL) {
        return(log_lik)
    } else {
        return(looOut)
    }
}

##' Predicts omega reliability from provided data.
##'
##' Empty for now.
##' @title Predict method for omegad object.
##' @param object omegad object.
##' @param newdata Data.frame. Must contain the exogenous (if used) and latent factor names.
##' @param summary Logical (Default: TRUE). Whether to return summary or array of samples.
##' @param prob Numeric (Default: .95). 
##' @param nsamples Numeric (Default: All). Number of posterior samples to use.
##' @param error Logical (Default: TRUE). Whether the predicted latent variables are sampled with stochastic error or not.
##' @param ... Unused.
##' @param samples Numeric (Default: All). Number of posterior samples to use.
##' @return Data.frame (summary=TRUE) or array (summary=FALSE).
##' @author Stephen R. Martin
##' @export
predict.omegad <- function(object, newdata, summary, prob = .95, nsamples = NULL, error = TRUE, ...) {
    if (is.null(nsamples)) {
        nsamples <- nsamples(omegad)
    }


    d <- .parse_formula.predict(object, newdata)
    if (object$meta$gp) {
        predOut <- .predict_gp_posterior(object, d, samples, error)
    }
}

##' Converts object and newdata into appropriate design matrices.
##'
##' @title Create data structures needed for prediction.
##' @inheritParams predict.omegad
##' @return List of data structures for prediction.
##' @author Stephen R. Martin
##' @keywords internal
.parse_formula.predict <- function(object, newdata) {
    if (any(is.na(newdata))) {
        stop("Predict requires complete data.")
    }

    # Exogenous MM #
    exoForm <- object$meta$modelForms$exo
    mm.exo <- model.matrix(exoForm, newdata)

    # Latent MM #
    latNames <- unlist(object$meta$fnames$factor)
    if (any(!(latNames %in% colnames(newdata)))) {
        stop("Latent factor scores must be included, and named.")
    }
    mm.lat <- as.matrix(newdata[,latNames, drop = FALSE])

    # Misc #
    N <- nrow(mm.lat)
    M <- object$meta$M
    F <- object$meta$F
    P <- ncol(mm.exo)

    out <- list(N = N,
                F = F,
                P = P,
                M = M,
                theta_loc = mm.lat,
                exo_x = mm.exo)
    return(out)
}

##############
## GP STUFF ##
##############
# yhat = xhat %*% (spds * gp_z(xs))
# You need to save the gp_z values, and estimate across the posterior of them.
# f(x) ~ MVN(0, K(x,x)) ## Standard GP
# f(x) = L(K)*gp_z, gp_z ~ N(0,1) ## Standard GP
# f(x) = phi*(D[spds]*gp_z), gp_z ~ N(0,1) ## GPA

.predict_gp_posterior <- function(object, data, nsamples, error){
    theta_loc <- data$theta_loc
    exo_x <- data$exo_x
    N <- data$N
    F <- data$F
    P <- object$meta$P
    M <- data$M
    exo <- object$meta$exo
    L <- 3*5/2
    lambdas <- .lambdas(L, M)
    F_inds <- object$stan_data$F_inds

    # Create Basis Functions
    gp_theta_phis <- lapply(seq_len(F), function(f) {
        .basis_phis(L, M, theta_loc[,f])
    })
    if (exo) {
        gp_exo_phis <- lapply(seq_len(P - 1), function(p) {
            .basis_phis(L, M, exo_x[, p + 1])
        })
    }

    # Extract posterior samples of needed params.
    gp_z <- .extract_transform(object$fit, "gp_z", nsamples = nsamples)
    gp_alpha <- .extract_transform(object$fit, "gp_alpha", nsamples = nsamples)
    gp_rho <- .extract_transform(object$fit, "gp_rho", nsamples = nsamples)
    lambda_loc_mat <- .extract_transform(object$fit, "lambda_loc_mat", nsamples = nsamples)
    lambda_sca_mat <- .extract_transform(object$fit, "lambda_sca_mat", nsamples = nsamples)
    nu_sca <- .extract_transform(object$fit, "nu_sca", nsamples = nsamples)
}

.predict_gp <- function(x, gp_linear_beta, gp_alpha, gp_rho, gp_z, M) {
    L <- 3*5/2
    lambdas <- .lambdas(L, M)
    x_phi <- .basis_phis(L, M, x)
    preds <- .spd_gp_fast(x_phi, gp_alpha, gp_rho, lambdas, gp_z)
    preds
}

.lambdas <- function(L, M) {
    ((pi * 1:M) / (2 * L))^2
}

.basis_phi <- function(L, m, x) {
    (1 / sqrt(L)) * sin((pi * m) / (2 * L) * (x + L))
}

.basis_phis <- function(L, M, x) {
    phis <- matrix(0, nrow = length(x), ncol = M)
    for (m in 1:M) {
        phis[, m] <- .basis_phi(L, m, x)
    }
    phis
}

.spd <- function(alpha, rho, lambda) {
    (alpha^2) * sqrt(2 * pi) * rho * exp(-.5 * (rho^2) * lambda^2)
}

.spds <- function(alpha, rho, lambdas) {
    .spd(alpha, rho, lambdas)
}

.spd_gp_fast <- function(x_phi, alpha, rho, lambdas, gp_z) {
    spds <- .spds(alpha, rho, lambdas)
    x_phi %*% (spds * gp_z)
}

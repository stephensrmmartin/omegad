##' Predicts omega reliability from provided data.
##'
##' Empty for now.
##' @title Predict method for omegad object.
##' @param object omegad object.
##' @param newdata Data.frame. Must contain the exogenous (if used) and latent factor names.
##' @param summary Logical (Default: TRUE). Whether to return summary or array of samples.
##' @param prob Numeric (Default: .95). 
##' @param samples Numeric (Default: All). Number of posterior samples to use.
##' @param error Logical (Default: TRUE). Whether the predicted latent variables are sampled with stochastic error or not.
##' @param ... Unused.
##' @return Data.frame (summary=TRUE) or array (summary=FALSE).
##' @author Stephen R. Martin
##' @export
predict.omegad <- function(object, newdata, summary, prob = .95, samples, error = TRUE, ...) {
    
}

.parse_formula.predict <- function(object, newdata) {
    
}

##############
## GP STUFF ##
##############
# yhat = xhat %*% (spds * gp_z(xs))
# You need to save the gp_z values, and estimate across the posterior of them.
# f(x) ~ MVN(0, K(x,x)) ## Standard GP
# f(x) = L(K)*gp_z, gp_z ~ N(0,1) ## Standard GP
# f(x) = phi*(D[spds]*gp_z), gp_z ~ N(0,1) ## GPA

.predict_gp_posterior <- function(){
    F <- ncol(theta_loc)
    L <- 3*5/2
    lambdas <- .lambdas(L, M)

    x_phis <- lapply(1:F, function(f) {
        .basis_phis(L, M, theta_loc[,f])
    })
    x_exo_phis <- lapply(1:P, function(p) {
        .basis_phis(L, M, exo_x[,p])
    })

    
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

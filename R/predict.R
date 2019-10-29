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
predict.omegad <- function(object, newdata, summary = TRUE, prob = .95, nsamples = NULL, error = TRUE, ...) {
    if (is.null(nsamples)) {
        nsamples <- nsamples(object)
    }
    probs <- .prob_to_probs(prob)


    d <- .parse_formula.predict(object, newdata)
    if (object$meta$gp) {
        predOut <- .predict_gp_posterior(object, d, nsamples, error)
    }

    if (summary) {
        predOut <- lapply(predOut, function(x) {
            out <- apply(x, c(1,2), function(y) {
                m <- mean(y)
                sd <- sd(y)
                ci <- quantile(y, probs)
                L <- ci[1]
                U <- ci[2]
                c(mean = m, sd = sd, L = L, U = U)
            })
            out <- aperm(out, c(2, 1, 3))
            colnames(out) <- c("Mean","SD",paste0("Q", probs[1]*100), paste0("Q", probs[2]*100))
            dimnames(out)[3] <- unlist(object$meta$fnames$factor)
            return(out)
        })
        ## predOut <- lapply(predOut, function(x){aperm(x, c(2,1,3))})
    } else {
        predOut <- lapply(predOut, function(x) {
            colnames(x) <- unlist(object$meta$fnames$factor)
            return(x)
        })
    }

    return(predOut)
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

#STEPHEN:: You need to implement .array_extract to ensure column and row vectors remain as such, rather than downgraded to just a vector.

.predict_gp_posterior <- function(object, data, nsamples, error){
    theta_loc <- data$theta_loc
    exo_x <- data$exo_x
    N <- data$N
    F <- data$F
    P <- object$meta$P
    M <- data$M
    exo <- object$meta$exo
    ## F_inds <- object$stan_data$F_inds
    ## F_inds_num <- sapply(1:F, function(x){
    ##     sum(F_inds[x,] != 0)
    ## })
    F_inds <- lapply(1:F, function(x){object$stan_data$F_inds[x,]})
    F_inds_num <- sapply(F_inds, function(x){sum(x != 0)})

    L <- 3*5/2
    lambdas <- .lambdas(L, M)

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
    gp_linear <- .extract_transform(object$fit, "gp_linear", nsamples = nsamples)
    gp_alpha <- .extract_transform(object$fit, "gp_alpha", nsamples = nsamples)
    gp_rho <- .extract_transform(object$fit, "gp_rho", nsamples = nsamples)
    lambda_loc_mat <- .extract_transform(object$fit, "lambda_loc_mat", nsamples = nsamples)
    lambda_sca_mat <- .extract_transform(object$fit, "lambda_sca_mat", nsamples = nsamples)
    nu_sca <- .extract_transform(object$fit, "nu_sca", nsamples = nsamples)
    theta_cor_L <- .extract_transform(object$fit, "theta_cor", nsamples = nsamples)
    theta_cor_L <- array(apply(theta_cor_L, 3, function(x){t(chol(x))}), dim=c(F,F,nsamples))
    if (exo) {
        exo_gp_z <- .extract_transform(object$fit, "exo_gp_z", nsamples = nsamples)
        exo_gp_linear <- .extract_transform(object$fit, "exo_gp_linear", nsamples = nsamples)
        exo_gp_alpha <- .extract_transform(object$fit, "exo_gp_alpha", nsamples = nsamples)
        exo_gp_rho <- .extract_transform(object$fit, "exo_gp_rho", nsamples = nsamples)
    }

    theta_sca <- array(0, dim=c(N, F, nsamples))
    omega1 <- array(0, dim = c(N, F, nsamples))
    omega2 <- array(0, dim = c(N, F, nsamples))
    if(error) {
        theta_sca[,,] <- rnorm(N*F*nsamples, 0, 1)
    }

    for (s in 1:nsamples) {
        for(f in 1:F){
            # GP part
            theta_sca[,f ,s] <- theta_sca[,f ,s] +
                .spd_gp_fast(gp_theta_phis[[f]], gp_alpha[f,s], gp_rho[f,s], lambdas, gp_z[,f,s]) +
                theta_loc[,f] * gp_linear[f, s]
            # Exogenous part
            if (exo) {
                for (p in 1:(P-1)) {
                    theta_sca[,f , s] <- theta_sca[,f , s] +
                        .spd_gp_fast(exo_gp_phis[[f]], exo_gp_alpha[p,f,s], gp_rho[p,f,s], lambdas, exo_gp_z[,f, p,s]) +
                        exo_x[,(p + 1)]*exo_gp_linear[p,f,s]
                }
            }
        }
        # Compute omegas
        shat <- exp(matrix(1,ncol=1,nrow=N)%*%t(.array_extract(nu_sca,s)) + .array_extract(theta_sca, s) %*% .array_extract(lambda_sca_mat, s))
        omega1[,,s] <- omega_one(.array_extract(lambda_loc_mat, s), F_inds, F_inds_num, shat)
        omega2[,,s] <- omega_two(.array_extract(lambda_loc_mat, s), F_inds, F_inds_num, .array_extract(theta_cor_L, s), shat)
    }
    out <- list(theta_sca = theta_sca, omega1 = omega1, omega2 = omega2)
    return(out)
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
    (alpha^2) * sqrt(2 * pi) * rho * exp(-.5 * (rho^2) * lambda)
}

.spds <- function(alpha, rho, lambdas) {
    .spd(alpha, rho, lambdas)
}

.spd_gp_fast <- function(x_phi, alpha, rho, lambdas, gp_z) {
    spds <- .spds(alpha, rho, lambdas)
    x_phi %*% (spds * gp_z)
}

# Omega Functions #

.omega_one <- function(lambda_loc_mat, F_inds, F_inds_num, shat) {
    N <- nrow(shat)
    F <- nrow(lambda_loc_mat)
    
}

.omega_two <- function(lambda_loc_mat, F_inds, F_inds_num, theta_cor, shat) {
    
}

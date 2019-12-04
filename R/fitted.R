##' Extracts estimated Error factor scores and reliability estimates for the fitted data.
##'
##' @title Extract fitted error factor scores, omega1, and omega2 values.
##' @inheritParams predict.omegad
##' @inherit predict.omegad return 
##' @author Stephen R. Martin
##' @export
fitted.omegad <- function(object, summary = TRUE, prob = .95, ...) {
    N <- object$meta$N
    S <- nsamples(object)
    F <- object$meta$F
    probs <- .prob_to_probs(prob)

    ## theta <- .extract_transform(object$fit, "theta")
    ## theta_sca <- theta[,(F+1):(F*2),,drop=FALSE]
    ## theta_loc <- theta[,1:(F),,drop=FALSE]
    ## omega1 <- .extract_transform(object$fit, "omega1")
    ## omega2 <- .extract_transform(object$fit, "omega2")

    theta <- .extract_transform(object$fit, "theta")
    theta_sca <- theta[,(F+1):(F*2),,drop=FALSE]
    theta_loc <- theta[,1:(F),,drop=FALSE]
    lambda_loc_mat <- .extract_transform(object$fit, "lambda_loc_mat")
    lambda_sca_mat <- .extract_transform(object$fit, "lambda_sca_mat")
    nu_sca <- .extract_transform(object$fit, "nu_sca")
    theta_cor <- .extract_transform(object$fit, "theta_cor")
    theta_cor_L <- array(apply(theta_cor, 3, function(s) {
        x <- s[1:F, 1:F]
        t(chol(x))
    }), dim = c(F, F, S))

    F_inds <- lapply(1:F, function(x) {object$stan_data$F_inds[x,]})
    F_inds_num <- sapply(F_inds, function(x) {sum(x != 0)})

    omega1 <- array(0, dim = c(N, F, S))
    omega2 <- array(0, dim = c(N, F, S))
    omega_total <- array(0, dim = c(N, 1, S))

    for(s in 1:S) {
        shat <- exp(matrix(1, ncol = 1, nrow = N) %*% t(.array_extract(nu_sca, s)) + .array_extract(theta_sca, s) %*% .array_extract(lambda_sca_mat, s))
        lambda_loc_mat_s <- .array_extract(lambda_loc_mat, s)
        theta_cor_L_s <- .array_extract(theta_cor_L, s)
        omega1[,, s] <- omega_one(lambda_loc_mat_s, F_inds, F_inds_num, shat)
        omega2[,, s] <- omega_two(lambda_loc_mat_s, F_inds, F_inds_num, theta_cor_L_s, shat)
        omega_total[,,s] <- omega_total(lambda_loc_mat_s, theta_cor_L_s, shat)
    }

    fitOut <- list(theta_sca = theta_sca, omega1 = omega1, omega2 = omega2, theta_loc = theta_loc, omega_total = omega_total)

    if (summary) {
        fitOut <- lapply(fitOut, function(x) {
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
            return(out)
        })
        fitOut[c("theta_sca","omega1","omega2","theta_loc")] <- lapply(fitOut[c("theta_sca","omega1","omega2","theta_loc")], function(x){
            dimnames(x)[[3]] <- unlist(object$meta$fnames$factor)
            return(x)
        })
        fitOut$omega_total <- fitOut$omega_total[,,1]
    } else {
        fitOut[c("theta_sca","omega1","omega2","theta_loc")] <- lapply(fitOut[c("theta_sca","omega1","omega2","theta_loc")], function(x){
            colnames(x) <- unlist(object$meta$fnames$factor)
            return(x)
        })
    }

    return(fitOut)
    
}

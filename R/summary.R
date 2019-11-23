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
##' Summarizes omegad object.
##'
##' Summarizes omegad object.
##' @title Summary method for omegad objects.
##' @param object omegad object.
##' @param prob Numeric (Default: .95). The amount of probability mass to include within the credible interval. Default values provide a 95\% credible interval.
##' @param std.lv Logical (Default: FALSE). Whether to compute loadings with standardized latents (TRUE) or not (FALSE).
##' @param ... Not used.
##' @return List containing "summary", "meta" (meta-data), and "diagnostics" (BFMI, Rhats, n_eff, max treedepth, divergences). "summary" is a list containing summaries (Mean, SD, intervals). Dimensions provided in brackets. J = number of items, N = number of observations, F = number of factors, P = number of exogenous predictors. Items and factors are named according to the model formula:
##' \describe{
##' \item{nu_loc}{[J, 4]. Expected values (intercepts) of the indicators. (I.e., when factors are at zero.)}
##' \item{nu_sca}{[J, 4]. Expected residual standard deviation of the indicators. (I.e., when Error factors are at zero.)}
##' \item{lambda_loc_mat}{[J, 4, F]. The loadings relating the latent factor to the indicators.}
##' \item{lambda_sca_mat}{[J, 4, F]. The loadings relating the latent Error factor to the residual SD of the indicators.}
##' \item{theta_cor}{[F*2, 4, F*2] if not using a GP. [F, 4, F] if using a GP. The correlation between all factors and Error factors (when not using a GP), or between factors (when using a GP).}
##' \item{gp_linear}{[F, 4]. The linear coefficient between factors and their corresponding error factors.}
##' \item{gp_alpha}{[F, 4]. The GP "alpha" parameter between factors and their corresponding error factors. This describes the marginal SD in the Error factor due to the non-linear GP.}
##' \item{gp_rho}{[F, 4]. The GP "rho" parameter (i.e., length scale). This controls the "wiggliness" of the function between factors and their corresponding Error factors. The smaller it is, the more influence nearby inputs have on the function output.}
##' \item{exo_beta}{[P, 4, F]. The linear coefficient between exogenous predictors and the Error factors.}
##' \item{exo_gp_linear}{[P, 4, F]. Exogenous GP model parameter. See gp_linear.}
##' \item{exo_gp_alpha}{[P, 4, F]. Exogenous GP model parameter. See gp_alpha.}
##' \item{exo_gp_rho}{[P, 4, F]. Exogenous GP model parameter. See gp_rho.}
##' }
##' @author Stephen R. Martin
##' @export
summary.omegad <- function(object, prob = .95, std.lv = TRUE, ...) {
    probs <- .prob_to_probs(prob)
    F <- object$meta$F
    F_inds <- object$stan_data$F_inds
    F_inds_num <- rowSums(F_inds != 0)
    fnames <- unlist(object$meta$fnames$factor)
    inames.all <- colnames(object$stan_data$x)
    inames <- lapply(seq_len(F), function(f) {
        inames.all[F_inds[f, seq_len(F_inds_num[f])]]
    })
    enames <- object$meta$enames
    exo <- object$meta$exo
    gp <- object$meta$gp

    .summary <- function(x) {
        m <- mean(x)
        sd <- sd(x)
        ci <- quantile(x, probs)
        L <- ci[1]
        U <- ci[2]
        out <- c(m, sd, L, U)
        names(out) <- c("Mean","SD",paste0("Q", probs*100))
        return(out)
    }

    latent_sd <- .get_latent_vars(object, prob, SD = TRUE, summary = FALSE)
    latent_mean <- .get_latent_means(object, prob, summary = FALSE)

    nu_loc <- .extract_transform(object$fit, "nu_loc")
    nu_sca <- .extract_transform(object$fit, "nu_sca")
    lambda_loc_mat <- .extract_transform(object$fit, "lambda_loc_mat")
    lambda_sca_mat <- .extract_transform(object$fit, "lambda_sca_mat")

    if (std.lv) {
        for (s in 1:nsamples(object)) {
            nu_loc[, s] <- nu_loc[, s] + t(latent_mean[s, 1:F, drop=FALSE] %*% .array_extract(lambda_loc_mat, s))
            nu_sca[, s] <- nu_sca[, s] + t(latent_mean[s, (F+1):(F*2), drop=FALSE] %*% .array_extract(lambda_sca_mat, s))
            for (f in 1:F) {
                lambda_loc_mat[f, , s] <- lambda_loc_mat[f, , s] * latent_sd[s, f]
                lambda_sca_mat[f, , s] <- lambda_sca_mat[f, , s] * latent_sd[s, (F + f)]
            }
        }
    }

    nu_loc_sum <- aperm(apply(nu_loc, 1, .summary), c(2,1))
    nu_sca_sum <- aperm(apply(nu_sca, 1, .summary), c(2,1))
    dimnames(nu_loc_sum)[[1]] <- dimnames(nu_sca_sum)[[1]] <- inames.all


    lambda_loc_mat_sum <- aperm(apply(lambda_loc_mat, c(1,2), .summary), c(3,1,2))
    lambda_sca_mat_sum <- aperm(apply(lambda_sca_mat, c(1,2), .summary), c(3,1,2))
    dimnames(lambda_loc_mat_sum)[[3]] <- dimnames(lambda_sca_mat_sum)[[3]] <- fnames
    dimnames(lambda_loc_mat_sum)[[1]] <- dimnames(lambda_sca_mat_sum)[[1]] <- inames.all

    theta_cor <- .extract_transform(object$fit, "theta_cor")
    theta_cor_sum <- aperm(apply(theta_cor, c(1, 2), .summary), c(2, 1, 3))

    # Collects names of objects to return via mget.
    outNames <- c("nu_loc_sum", "nu_sca_sum", "lambda_loc_mat_sum", "lambda_sca_mat_sum", "theta_cor_sum")

    if (gp) {
        gp_linear <- .extract_transform(object$fit, "gp_linear")
        gp_alpha <- .extract_transform(object$fit, "gp_alpha")
        gp_rho <- .extract_transform(object$fit, "gp_rho")

        if (std.lv) {
            for(s in seq_len(nsamples(object))) {
                gp_linear[, s] <- gp_linear[, s] / latent_sd[s, (F+1):(2*F)]
                gp_alpha[, s] <- gp_alpha[, s] / latent_sd[s, (F+1):(2*F)]
            }
        }

        gp_linear_sum <- aperm(apply(gp_linear, 1, .summary), c(2, 1))
        gp_alpha_sum <- aperm(apply(gp_alpha, 1, .summary), c(2, 1))
        gp_rho_sum <- aperm(apply(gp_rho, 1, .summary), c(2, 1))

        dimnames(gp_linear_sum)[[1]] <- dimnames(gp_alpha_sum)[[1]] <- dimnames(gp_rho_sum)[[1]] <- fnames
        dimnames(theta_cor_sum)[[1]] <- dimnames(theta_cor_sum)[[3]] <- fnames

        outNames <- c(outNames, "gp_linear_sum", "gp_alpha_sum", "gp_rho_sum")

        if (exo) {
            exo_gp_linear <- .extract_transform(object$fit, "exo_gp_linear")
            exo_gp_alpha <- .extract_transform(object$fit, "exo_gp_alpha")
            exo_gp_rho <- .extract_transform(object$fit, "exo_gp_rho")

            if (std.lv) {
                for(s in seq_len(nsamples(object))) {
                    for (f in seq_len(F)) {
                        exo_gp_linear[, f, s] <- exo_gp_linear[, f, s]/latent_sd[s, (F + f)]
                        exo_gp_alpha[, f, s] <- exo_gp_alpha[, f, s]/latent_sd[s, (F + f)]
                    }
                }
            }

            exo_gp_linear_sum <- aperm(apply(exo_gp_linear, c(1, 2), .summary), c(2, 1, 3))
            exo_gp_alpha_sum <- aperm(apply(exo_gp_alpha, c(1, 2), .summary), c(2, 1, 3))
            exo_gp_rho_sum <- aperm(apply(exo_gp_rho, c(1, 2), .summary), c(2, 1, 3))
            dimnames(exo_gp_linear_sum)[[3]] <- dimnames(exo_gp_alpha_sum)[[3]] <- dimnames(exo_gp_rho_sum)[[3]] <- paste0(fnames, "_Error")
            dimnames(exo_gp_linear_sum)[[1]] <- dimnames(exo_gp_alpha_sum)[[1]] <- dimnames(exo_gp_rho_sum)[[1]] <- enames$terms

            outNames <- c(outNames, "exo_gp_linear_sum", "exo_gp_alpha_sum", "exo_gp_rho_sum")
        }
    } else { # Not GP
       dimnames(theta_cor_sum)[[1]] <- dimnames(theta_cor_sum)[[3]] <- c(fnames, paste0(fnames, "_Error"))

       if (exo) {
           exo_beta <- .extract_transform(object$fit, "exo_beta")
           if (std.lv) {
               for (s in seq_len(nsamples(object))) {
                   for (f in seq_len(F)) {
                       exo_beta[, f, s] <- exo_beta[, f, s] / latent_sd[s, (F + f)]
                   }
               }
           }
           exo_beta_sum <- aperm(apply(exo_beta, c(1, 2), .summary), c(2, 1, 3))
           dimnames(exo_beta_sum)[[3]] <- paste0(fnames, "_Error")
           dimnames(exo_beta_sum)[[1]] <- enames$terms

           outNames <- c(outNames, "exo_beta_sum")
       }
    }
    out <- mget(outNames)
    names(out) <- gsub("_sum", "", outNames)
    out <- list(summary = out)

    out$meta <- object$meta
    out$meta$std.lv <- std.lv
    out$diagnostics <- object$diagnostics

    dots <- list(...)
    if (is.null(dots$digits)) {
        dots$digits <- 3
    }
    out$meta$digits <- dots$digits

    class(out) <- "summary.omegad"

    return(out)
}
##' Prints summary.omegad objects.
##'
##' Prints summary.omegad objects.
##' @title Print method for omegad summaries.
##' @param x summary.omegad object.
##' @param ... Not used.
##' @return x (invisible)
##' @author Stephen R. Martin
##' @export
print.summary.omegad <- function(x, ...) {
    digits <- x$meta$digits
    dots <- list(...)
    if (!is.null(dots$digits)) {
       digits <- dots$digits 
    }

    .cat_line()
    cat("Diagnostics: \n")
    .print_diagnostics(x$diagnostics)
    .cat_line()
    cat("Model Description: \n")
    cat(.get_model_description(x), "\n")
    .cat_line()

    # Loadings
    cat("Loadings: \n")
    for (f in 1:x$meta$F) {
        fname <- x$meta$fnames$factor[[f]]
        inames <- x$meta$fnames$indicator[[f]]
        cat(fname, "[Location] \n")
        print(x$summary$lambda_loc_mat[inames,,fname], digits = digits)
        cat("\n")

        cat(fname, "[Error] \n" )
        print(x$summary$lambda_sca_mat[inames,,fname], digits = digits)
        cat("\n")
    }

    # Intercepts
    cat("Intercepts: \n")
    cat("[Location] \n")
    print(x$summary$nu_loc, digits = digits)
    cat("\n")
    cat("[Log SD] \n")
    print(x$summary$nu_sca, digits = digits)
    cat("\n")

    # Covariance
    if (x$meta$F > 1 | !x$meta$gp) {
        cat("Latent Correlations: \n")
        print(x$summary$theta_cor[,"Mean",], digits = digits)
        cat("\n")
    }

    # GP
    if (x$meta$gp) {
        cat("Approximate Gaussian Process: \n")
        cat("[Linear] \n")
        print(x$summary$gp_linear, digits = digits)
        cat("[Alpha] \n")
        print(x$summary$gp_alpha, digits = digits)
        cat("[Length scale] \n")
        print(x$summary$gp_rho, digits = digits)
        cat("\n")
    }
    # Exogenous
    if (x$meta$exo) {
        cat("Exogenous Variables: \n")
        if (x$meta$gp) {
            for (f in 1:x$meta$F) {
                fname <- paste0(x$meta$fnames$factor[[f]],"_Error")
                cat(fname, "\n")
                print(x$summary$exo_gp_linear[,, fname], digits = digits)
                print(x$summary$exo_gp_alpha[,, fname], digits = digits)
                print(x$summary$exo_gp_rho[,, fname], digits = digits)
                cat("\n")
            }
        } else {
            for (f in 1:x$meta$F) {
                fname <- paste0(x$meta$fnames$factor[[f]],"_Error")
                cat(fname, "\n")
                print(x$summary$exo_beta[,, fname], digits = digits)
                cat("\n")
            }
        }

    }
    return(invisible(x))
}

.get_model_description <- function(object) {
    if (object$meta$gp) {
        out <- c("\n\t[Factor -> Error Factor]",
                 paste0("\tUnivariate gaussian process: Additive linear and exponential quadratic kernels (", object$meta$M, " basis functions) ")
                 )
        if (object$meta$exo) {
            out <- c(out[1], "\t[Exogenous -> Error Factor]", out[2])
        }
        out <- paste0(out, collapse = "\n")
    } else {
        out <- c("\n \t[Factors <-> Error Factors]", "\tCovariance only")
        if (object$meta$exo) {
            out <- c(out, "\t[Exogenous -> Error Factors]", "\tLinear model")
        }
        out <- paste0(out, collapse = "\n")
    }

    return(out)
}

.cat_line <- function(n = 40) {
    str <- c(paste0(rep("-", n), collapse=""), "\n")
    cat(str)
}

.print_diagnostics <- function(diagnostics) {
    d <- diagnostics
    rh <- d$rhats[!is.na(d$rhats)]
    if (sum(rh > 1.1)) {
        cat("\t Rhats: Failed \n")
        cat("\t Some Rhats > 1.1. Do not interpret results. The largest 10 are:\n")
        print(head(sort(rh, decreasing = TRUE), 10))
    } else {
        cat("\t Rhats: Passed \n")
    }

    if (d$div > 0) {
        cat("\t Divergent Transitions: Failed - ", d$div, "divergent transitions detected. \n")
    } else {
        cat("\t Divergent Transitions: Passed \n")
    }

    if (d$tree.max > 0) {
        cat("\t Max treedepth hit:", d$tree.max, "\n")
    }
    if (any(d$bfmi < .2)) {
       cat("\t Low E-BFMI detected in chains", which(d$bfmi < .2), "\n") 
    }
}
##' Computes latent variances from fitted.
##'
##' Standard decomp rules could be followed.
##' GPs add alpha^2 marginal variance, and the linear functions follow standard path tracing rules.
##' However, a more straightforward approach is to just directly compute the variance from the fitted scores, which uses fewer assumptions and inherits all uncertainty from the params anyway.
##' @title Compute latent variances from fitted.
##' @param object omegad object.
##' @param prob Probability mass for interval.
##' @param SD Logical (Default: FALSE). Whether to compute SDs (TRUE) or Vars (FALSE).
##' @param summary Logical (Default: TRUE). Whether to summarize (TRUE) or return samples (FALSE).
##' @return Named vector of latent variances.
##' @author Stephen R. Martin
##' @keywords internal
.get_latent_vars <- function(object, prob, SD = FALSE, summary = TRUE) {
    probs <- .prob_to_probs(prob)
    fnames <- unlist(object$meta$fnames$factor)
    fit <- fitted(object, summary = FALSE)
    theta_loc <- fit$theta_loc
    theta_sca <- fit$theta_sca
    if(SD) {
        f <- sd
    } else {
        f <- var
    }
    theta_loc_vars <- t(apply(theta_loc, c(2,3), f))
    theta_sca_vars <- t(apply(theta_sca, c(2,3), f))
    theta_vars <- cbind(theta_loc_vars, theta_sca_vars)
    colnames(theta_vars) <- c(fnames, paste0(fnames,"_Error"))
    out <- theta_vars
    if (!summary) {
        return(out)
    }

    out <- apply(theta_vars, 2, function(x) {
        M <- mean(x)
        S <- sd(x)
        ci <- quantile(x, probs)
        L <- ci[1]
        U <- ci[2]
        out <- c(M, S, L, U)
        names(out) <- c("Mean","SD",paste0("Q",probs*100))
        return(out)
    })
    out <- t(out)
    return(out)
}

.get_latent_means <- function(object, prob, summary = TRUE) {
    probs <- .prob_to_probs(prob)
    fnames <- unlist(object$meta$fnames$factor)
    fit <- fitted(object, summary = FALSE)
    theta_loc <- fit$theta_loc
    theta_sca <- fit$theta_sca
    theta_loc_mean <- t(apply(theta_loc, c(2,3), mean))
    theta_sca_mean <- t(apply(theta_sca, c(2,3), mean))
    theta_mean <- cbind(theta_loc_mean, theta_sca_mean)
    colnames(theta_mean) <- c(fnames, paste0(fnames,"_Error"))
    out <- theta_mean
    if (!summary) {
        return(out)
    }

    out <- apply(theta_mean, 2, function(x) {
        M <- mean(x)
        S <- sd(x)
        ci <- quantile(x, probs)
        L <- ci[1]
        U <- ci[2]
        out <- c(M, S, L, U)
        names(out) <- c("Mean","SD",paste0("Q",probs*100))
        return(out)
    })
    out <- t(out)
    return(out)
}

##' Recompute parameter outputs as standardized values.
##'
##' Loadings and intercepts can be recomputed assuming standardized values via linear transformations.
##' For loadings:
##' \eqn{\beta_{a}SD_a = \beta_{b}SD_b}
##' \eqn{\beta_{un}SD_{un} = \beta_{st}}
##' \eqn{\beta_{0c} = \beta_{0} + \beta_1\bar x}
##' \eqn{\nu_c = \nu + \alpha\Lambda}
##' \eqn{\Lambda^c = D\Lambda}
##' If a standardized endogenous variable (Var(theta) = 1), then exogenous predictors have:
##' \eqn{\beta_{st} = \beta_{un}/SD_{latent}}
##' \eqn{\sigma^2_\epsilon = 1/Var(latent)}
##'
##' For accuracy, I think we need to compute these over each posterior sample.
##' On second thought, maybe not. Scaling the lambdas by E(SD) seemed to work fine (?). Need to doublecheck posterior SDs.
##' To do so, should transform theta_sca on each s by the sd at s; compare that to using E(theta_sca[i])/E(SD).
##' This actually does work - The Posterior SDs are a touch high, but it's extremely close. Expected values are correct.
##' @title Standardize summary output.
##' @param x summary.omegad object.
##' @return Summary object (with standardized values).
##' @author Stephen R. Martin
##' @keywords internal
.standardize_output <- function(x) {
    
}

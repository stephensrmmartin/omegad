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
summary.omegad <- function(object, prob = .95, ...) {
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

    nu_loc <- .extract_transform(object$fit, "nu_loc")
    nu_sca <- .extract_transform(object$fit, "nu_sca")
    lambda_loc_mat <- .extract_transform(object$fit, "lambda_loc_mat")
    lambda_sca_mat <- .extract_transform(object$fit, "lambda_sca_mat")

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
    out$diagnostics <- object$diagnostics

    class(out) <- "summary.omegad"

    return(out)
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
    .cat_line()
    cat("Diagnostics: \n")
    # .print_diagnostics(x)
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
        print(x$summary$lambda_loc_mat[inames,,fname])
        cat("\n")

        cat(fname, "[Error] \n" )
        print(x$summary$lambda_sca_mat[inames,,fname])
        cat("\n")
    }

    # Intercepts
    cat("Intercepts: \n")
    cat("[Location] \n")
    print(x$summary$nu_loc)
    cat("\n")
    cat("[Log SD] \n")
    print(x$summary$nu_sca)
    cat("\n")

    # Covariance
    if (x$meta$F > 1 | !x$meta$gp) {
        cat("Latent Correlations: \n")
        print(x$summary$theta_cor[,"Mean",])
        cat("\n")
    }

    # GP
    if (x$meta$gp) {
        cat("Approximate Gaussian Process: \n")
        cat("[Linear] \n")
        print(x$summary$gp_linear)
        cat("[Alpha] \n")
        print(x$summary$gp_alpha)
        cat("[Length scale] \n")
        print(x$summary$gp_rho)
        cat("\n")
    }
    # Exogenous
    if (x$meta$exo) {
        cat("Exogenous Variables: \n")
        if (x$meta$gp) {
            for (f in 1:x$meta$F) {
                fname <- x$meta$fnames$factor[[f]]
                cat(fname, "\n")
                print(x$summary$exo_gp_linear[,, fname])
                print(x$summary$exo_gp_alpha[,, fname])
                print(x$summary$exo_gp_rho[,, fname])
                cat("\n")
            }
        } else {
            for (f in 1:x$meta$F) {
                fname <- x$meta$fnames$factor[[f]]
                cat(fname, "\n")
                print(x$summary$exo_beta[,, fname])
                cat("\n")
            }
        }

    }
    return(invisible(x))
}

.get_model_description <- function(object) {
    if (object$meta$gp) {
        out <- c("\n \t [Factor -> Error Factor]",
                 paste0("\t Univariate gaussian process: Additive linear and exponential quadratic kernels (", object$meta$M, " basis functions) ")
                 )
        if (object$meta$exo) {
            out <- c(out[1], "\t [Exogenous -> Error Factor]", out[2])
        }
        out <- paste0(out, collapse = "\n")
    } else {
        out <- c("\n \t [Factors -> Error Factors]", "\t Covariance only")
        if (object$meta$exo) {
            out <- c(out, "\t [Exogenous -> Factor]", "\t Linear model")
        }
        out <- paste0(out, collapse = "\n")
    }

    return(out)
}

.cat_line <- function(n = 40) {
    str <- c(paste0(rep("-", n), collapse=""), "\n")
    cat(str)
}

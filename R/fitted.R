##' Extracts estimated Error factor scores and reliability estimates for the fitted data.
##'
##' @title Extract fitted error factor scores, omega1, and omega2 values.
##' @inheritParams predict.omegad
##' @inherit predict.omegad return 
##' @author Stephen R. Martin
##' @export
fitted.omegad <- function(object, summary = TRUE, prob = .95, ...) {
    F <- object$meta$F
    probs <- .prob_to_probs(prob)

    theta <- .extract_transform(object$fit, "theta")
    theta_sca <- theta[,(F+1):(F*2),,drop=FALSE]
    theta_loc <- theta[,1:(F),,drop=FALSE]
    omega1 <- .extract_transform(object$fit, "omega1")
    omega2 <- .extract_transform(object$fit, "omega2")

    fitOut <- list(theta_sca = theta_sca, omega1 = omega1, omega2 = omega2, theta_loc = theta_loc)

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
            dimnames(out)[[3]] <- unlist(object$meta$fnames$factor)
            return(out)
        })
        ## predOut <- lapply(predOut, function(x){aperm(x, c(2,1,3))})
    } else {
        fitOut <- lapply(fitOut, function(x) {
            colnames(x) <- unlist(object$meta$fnames$factor)
            return(x)
        })
    }

    return(fitOut)
    
}

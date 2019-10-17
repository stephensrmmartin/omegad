##' Generates data for univariate GP model.
##'
##' Generates dataset for GP model.
##' The (univariate) GP model is such that theta_sca = fun(theta_loc).
##'
##' @title Datagen for GP.
##' @param N Number of persons.
##' @param lambda_loc (Standardized) Loadings for theta_loc.
##' @param lambda_sca_std Standardized loadings for theta_sca.
##' @param nu_loc Intercepts for location models.
##' @param nu_sca Intercepts for scale models.
##' @param fun True function for theta_sca = fun(theta_loc)
##' @return List of params, meta-data (meta), and stan data (data).
##' @author Stephen R. Martin
##' @keywords internal
datagen.gp.uni <- function(N,lambda_loc,lambda_sca_std,nu_loc,nu_sca,fun){
  N <- N
  J <- length(lambda_loc)
  theta_loc <- rnorm(N)
  theta_sca <- fun(theta_loc)
  theta_sca_var <- var(theta_sca)
  lambda_sca <- sqrt((lambda_sca_std)^2/theta_sca_var)
  xhat <- theta_loc %*% matrix(lambda_loc,nrow=1) + rep(1,N)%*%matrix(nu_loc,nrow=1)
  shat <- exp(theta_sca%*%matrix(lambda_sca,nrow=1) + rep(1,N)%*%matrix(nu_sca,nrow=1))
  x <- matrix(rnorm(N*J,as.vector(xhat),as.vector(shat)),ncol=J)
  colnames(x) <- paste0('x',1:J)

  params <- mget(c('theta_loc','theta_sca','lambda_loc','lambda_sca','nu_loc','nu_sca'))
  meta <- mget(c('N','J'))
  data <- list(N=N,J=J,`F`=1,F_inds=matrix(1:J,nrow=1),x=x)
  return(mget(c('params','meta','data')))
}

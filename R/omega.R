
#' Fit Omega model
#'
#' Fits the "Reliability factor" model for varying omega coefficients.
#'
#' \code{omegad} takes a formula, or list of formulas, describing the factor structure.
#' If a unidimensional structure is desired (i.e., one latent variable), then the formula
#' should be as follows:
#'
#' \code{factorName ~ indicator1 + indicator2 + indicatorJ}
#'
#' When a multidimensional structure is desired, then the formula argument should be a list of formulas.
#' For example, if one has 5 items for agreeableness and openness:
#'
#' \code{list(agree ~ A1 + A2 + A3 + A4 + A5, open ~ O1 + O2 + O3 + O4 + O5)}
#'
#' If \code{\link[lavaan]{lavaan}} is familiar to you, then a lavaan formula such as:
#' \code{agree =~ A1 + A2} is identical to an \code{omegad} formula of \code{agree ~ A1 + A2}.
#' 
#' In addition, the ``Error'' factor(s) can be predicted using \code{Error ~ x1 + x2 + ...}.
#' Currently, the provided covariates predict all Error factor scores.
#' Moreover, exogenous variables must be observed (as of now).
#'
#' Note that omegad imposes a positivity constraint on all loadings.
#' This means that *all* reverse-scored items should be reflected so that all items are in the same direction.
#'
#' @param formula Formula specifying the factor structure. See details.
#' @param data Data frame containing indicators as columns.
#' @param gp Logical (Default: FALSE). Whether to fit an approximate gaussian process model between latent factors (and exogenous variables) and their corresponding Error factors.
#' @param M Integer (Default: 10). Only applicable if \code{gp = TRUE}. The number of basis functions to use in the approximate gaussian processes. Larger M provides a better approximation at a greater computational cost.
#' @param std.ov Logical (Default: TRUE). Whether to standardize the indicators predictors. The implemented priors assume approximately standardized indicators.
#' @param ... Options passed onto \code{\link[rstan]{sampling}}.
#'
#' @importFrom rstan sampling
#' @importFrom parallel detectCores
#'
#' @return omegad object.
#' @export
#'
omegad <- function(formula, data, ...) {
  dots <- list(...)
  if (is.null(dots$gp)) {
    gp <- FALSE
  } else {
    gp <- dots$gp
    dots$gp <- NULL
  }
  if (is.null(dots$M)) {
      M <- 10
  } else {
      M <- dots$M
      dots$M <- NULL
  }
  if (is.null(dots$std.ov)) {
      std.ov <- TRUE
  } else {
      std.ov <- dots$std.ov
      dots$std.ov <- NULL
  }
  if (is.null(dots$cores)) {
    dots$cores <- getOption("mc.cores")
    if (is.null(dots$cores)) {
      dots$cores <- detectCores()
    }
  }
  if (is.null(dots$control)) {
    dots$control <- list(adapt_delta = .95)
  }
  if (is.null(dots$control$adapt_delta)) {
    dots$control$adapt_delta <- .95
  }
  if (is.null(dots$chains)) {
    dots$chains <- 4
  }
  if (is.null(dots$init)) {
      dots$init <- 0
  }

  d <- .parse_formula(formula, data, std.ov)
  pars <- c("lambda_loc_mat",
            "lambda_sca_mat",
            "nu_loc",
            "nu_sca",
            "theta_cor",
            "theta",
            ## "omega1",
            ## "omega2",
            "omega1_expected",
            "omega2_expected",
            "omega_total_expected")
  if (gp) {
    d$stan_data$M <- M
    model <- stanmodels$relFactorGeneralGPBP
    pars <- c(pars,
              "gp_alpha",
              "gp_rho",
              "gp_linear",
              "exo_gp_alpha",
              "exo_gp_rho",
              "exo_gp_linear",
              "gp_z",
              "exo_gp_z")
  } else {
    pars <- c(pars, "exo_beta")
    model <- stanmodels$relFactorGeneral
  }
  args <- c(list(object = model, data = d$stan_data, pars = pars), dots)
  stanOut <- do.call("sampling", args = args)

  meta <- list(gp = gp,
               M = M,
               exo = d$exo,
               P = d$stan_data$P,
               N = d$stan_data$N,
               J = d$stan_data$J,
               F = d$stan_data$F,
               fnames = d$fnames,
               enames = d$enames,
               modelForms = d$modelForms,
               std.ov = std.ov,
               scaling = d$scaling
               )
  if (!is.list(formula)) {
      formula <- list(formula)
  }
  out <- list(formula = formula,
              data = d$model.frame,
              stan_data = d$stan_data,
              fit = stanOut,
              meta = meta)
  out$diagnostics <- .get_diagnostics(out)
  class(out) <- "omegad"

  return(out)
}

#' Parses formula and data into stan_data structure
#'
#' @param formula Formula or list of formulas with factor name on LHS, indicators on RHS.
#' @param data data.frame containing indicators.
#' @param std.ov Logical. Whether to standardize indicators.
#' @param ... Currently not used.
#'
#' @import Formula
#' @return stan_data list.
#' @keywords internal
#'
.parse_formula <- function(formula, data, std.ov = TRUE, ...) {
  if (!is.list(formula)) {
    forms <- list(formula)
  } else {
    forms <- formula
  }
  # Convert to Formula (easier)
  forms <- lapply(forms, function(f) {
    as.Formula(f)
  })

  # Check structure
  for (f in seq_len(length(forms))) {
    if (length(forms[[f]])[1] == 0) {
      stop("Factor name should be provided on LHS of formula.")
    }
  }

  # Separate formulas
  ## forms.rhs <- lapply(forms, FUN = function(f) {
  ##   formula(f, lhs = 0, rhs = 1)
  ## })

  # Get names
  fnames <- .get_names_formulaList(forms, formula = TRUE)

  # Combine all rhs formulas for model frame
  form.rhs <- do.call(c, fnames$indicator)
  form.rhs <- as.Formula(paste0("~ ", paste0(unique(form.rhs), collapse = " + ")))

  # Model frame
  mf <- model.frame(form.rhs, data = data, na.action = na.pass)
  if (any(is.na(mf))) {
    n_before <- nrow(mf)
    mf <- mf[complete.cases(mf), ]
    n_after <- nrow(mf)
    warning("Removing ", n_before - n_after, " incomplete cases.")
  }

  # Rescale indicators and non-factor exogenous variables?
  ## if (std.ov) {
  ##     which.nonfactor <- sapply(mf, is.numeric)
  ##     mf[,which.nonfactor] <- scale(mf[,which.nonfactor])
  ## }

  # Extract out Exogenous variables
  whichExoForm <- which(fnames$factor == "Error")
  if (any(whichExoForm)) {
      exoForm.rhs <- formula(forms[[whichExoForm]], lhs = 0, rhs = 1)
      forms[[whichExoForm]] <- NULL
      exo <- TRUE
      ## forms.rhs[[whichExoForm]] <- NULL
  } else {
      exoForm.rhs <- as.Formula(~ 1)
      exo <- FALSE
  }
  mm.exo <- model.matrix(exoForm.rhs, mf)
  P <- ncol(mm.exo)

  ## Indicator matrix
  fnames <- .get_names_formulaList(forms, formula = TRUE)
  form.rhs <- do.call(c, fnames$indicator)
  form.rhs <- as.Formula(paste0("~ ", paste0(unique(form.rhs), collapse = " + ")))
  mm <- model.matrix(form.rhs, mf)[, -1] ## No intercept

  scaling <- list(mean = apply(mm, 2, mean), sd = apply(mm, 2, sd))
  if (std.ov) {
      mm <- scale(mm)
  }

  # Loading indicator matrix
  F_inds <- .get_loadings_matrix(forms, mm)

  # Misc Stan data
  N <- nrow(mm)
  J <- ncol(mm)
  `F` <- nrow(F_inds)

  # Misc
  enames <- list(terms = attr(terms(exoForm.rhs), "term.labels"), vars = all.vars(exoForm.rhs))
  modelForms <- list(latent = forms, exo = exoForm.rhs)


  out <- list(stan_data = list(N = N, J = J, `F` = `F`, F_inds = F_inds, x = mm, exo_x = mm.exo, P = P),
              model.frame = mf,
              exo = exo,
              enames = enames,
              fnames = fnames,
              modelForms = modelForms,
              scaling = scaling)
  return(out)
}

#' Takes formula and returns names
#'
#' @param formList List of 2-sided Formulas
#' @param formula Logical (Default: FALSE). If TRUE, returns the rhs formula terms (e.g., I(x^2), rather than x)
#'
#' @return List containing factor names (factor) and indicator names (indicator).
#' @keywords internal
.get_names_formulaList <- function(formList, formula=FALSE) {
  fNames <- lapply(formList, function(f) {
    all.vars(f)[1]
  })
  iNames <- lapply(formList, function(f) {
    if (formula) {
      attr(terms(f), "term.labels")
    } else {
      all.vars(f)[-1]
    }
  })
  list(factor = fNames, indicator = iNames)
}

#' Create FxJ loadings matrix, padded with zeroes
#'
#' Takes (RHS) formula list and model *matrix*
#'
#' @param formList Formula list (one entry per factor)
#' @param mm Model matrix
#'
#' @return FxJ matrix with the indicator indices for each factor. Padded with zeroes.
#' @keywords internal
#'
.get_loadings_matrix <- function(formList, mm) {
  fnames <- .get_names_formulaList(formList, formula = TRUE)
  inames <- colnames(mm)
  `F` <- length(formList)
  J <- ncol(mm)

  lmat <- matrix(0, nrow = `F`, ncol = J)
  for (f in seq_len(F)) {
    f.num <- length(fnames$indicator[[f]])
    lmat[f, seq_len(f.num)] <- match(fnames$indicator[[f]], inames)
  }
  row.names(lmat) <- do.call(c, fnames$factor)
  lmat

}

.get_diagnostics <- function(object) {
    params <- c("lambda_loc_mat", "lambda_sca_mat", "nu_loc", "nu_sca", "theta_cor")
    if (object$meta$gp) {
        params <- c(params, "gp_linear", "gp_alpha", "gp_rho")
        if (object$meta$exo) {
            params <- c(params, "exo_gp_linear", "exo_gp_alpha", "exo_gp_rho")
        }
    } else {
        if(object$meta$exo) {
            params <- c(params, "exo_beta")
        }
    }
  rhats <- rstan::summary(object$fit, pars = params)$summary[, "Rhat"]

  n_effs <- rstan::summary(object$fit, pars = params)$summary[, "n_eff"]

  div <- rstan::get_num_divergent(object$fit)

  tree.max <- rstan::get_num_max_treedepth(object$fit)

  bfmi <- rstan::get_bfmi(object$fit)

  return(mget(c("rhats", "n_effs", "div", "tree.max", "bfmi")))
}

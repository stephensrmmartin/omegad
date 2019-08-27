
#' Fit Omega model
#'
#' @param formula Formula specifying the factor structure. See details.
#' @param data Data frame containing indicators as columns.
#' @param ... Options passed onto sampling.
#'
#' @importFrom rstan sampling
#'
#' @return omegad object.
#' @export
#'
#' @examples
omega <- function(formula,data,...){

}

#' Parses formula and data into stan_data structure
#'
#' @param formula ...
#' @param data ...
#' @param ... ...
#'
#' @import Formula
#' @return stan_data list.
#' @keywords internal
#'
.parse_formula <- function(formula,data,...){
  if(!is.list(formula)) {
    forms <- list(formula)
  } else {
    forms <- formula
  }
  # Convert to Formula (easier)
  forms <- lapply(forms,function(f){
    as.Formula(f)
  })

  # Check structure
  for(f in 1:length(forms)){
    if(length(forms[[f]])[1] == 0){
      stop('Factor name should be provided on LHS of formula.')
    }
  }

  # Separate formulas
  forms.rhs <- lapply(forms,FUN=function(f){
    formula(f,lhs=0,rhs=1)
  })

  # Get names
  fnames <- .get_names_formulaList(forms,formula=TRUE)

  # Combine all rhs formulas for model frame
  form.rhs <- do.call(c,fnames$indicator)
  form.rhs <- as.Formula(paste0('~ ',paste0(unique(form.rhs),collapse = ' + ')))

  # Model frame
  mf <- model.frame(form.rhs,data=data,na.action = na.pass)
  if(any(is.na(mf))){
    n_before <- nrow(mf)
    mf <- mf[complete.cases(mf),]
    n_after <- nrow(mf)
    warning('Removing ',n_before - n_after, ' incomplete cases.')
  }
  ## Remove incompletes
  mm <- model.matrix(form.rhs,mf)[,-1] ## No intercept

  # Loading indicator matrix
  F_inds <- .get_loadings_matrix(forms,mm)

  # Misc Stan data
  N <- nrow(mm)
  J <- ncol(mm)
  `F` <- length(forms.rhs)

  list(N=N,J=J,`F`=`F`,F_inds=F_inds,x=mm)

}

#' Takes formula and returns names
#'
#' @param formList List of 2-sided Formulas
#' @param formula Logical (Default: FALSE). If TRUE, returns the rhs formula terms (e.g., I(x^2), rather than x)
#'
#' @return
#' @keywords internal
#'
.get_names_formulaList <- function(formList,formula=FALSE){
  fNames <- lapply(formList,function(f){
    all.vars(f)[1]
  })
  iNames <- lapply(formList,function(f){
    if(formula){
      attr(terms(f),'term.labels')
    } else {
      all.vars(f)[-1]
    }
  })
  list(factor=fNames,indicator=iNames)
}

#' Create FxJ loadings matrix, padded with zeroes
#'
#' Takes (RHS) formula list and model *matrix*
#'
#' @param names Output of .get_names_formulaList
#' @param formList Formula list (one entry per factor)
#' @param mf Model frame
#'
#' @return
#' @keywords internal
#'
.get_loadings_matrix <- function(formList,mm){
  fnames <- .get_names_formulaList(formList,formula=TRUE)
  inames <- colnames(mm)
  `F` <- length(formList)
  J <- ncol(mm)

  lmat <- matrix(0,nrow=`F`,ncol=J)
  for(f in 1:`F`){
    f.num <- length(fnames$indicator[[f]])
    lmat[f,1:f.num] <- match(fnames$indicator[[f]],inames)
  }
  row.names(lmat) <- do.call(c,fnames$factor)
  lmat

}

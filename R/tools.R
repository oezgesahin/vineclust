#' internal function
#' @noRd
fit_info <- function(x,...){
  cat("logLik =", round(x$output$loglik, 2), "  ")
  cat("BIC =", round(x$output$bic, 2), "  ")
  cat("ICL =", round(x$output$icl, 2), "  ")
  cat("npars =", x$output$total_pars, "  ")
  cat("initial clustering =", x$output$init_clustering, "  ")
  cat("ECM iterations =", x$output$iteration, "  ")
  cat("\n")
}

#' internal function
#' @noRd
initial_df_check <- function(data){
  data <- data.frame(data)
  if(is.null(dim(data)[1])){
    stop("data must contain at least 2 variables")
  }
  else{
    nrow_df <- dim(data)[1]
    ncol_df <- dim(data)[2]
  }
  if(nrow_df == 0)
    stop("data does not contain any observations")

  if(any(sapply(data, is.na)))
    stop("data must not contain any missing values")

  if(!all(sapply(data, is.numeric)))
    stop("data must be numeric")
}

#' internal function
#' @noRd
initial_args_check <- function(data, total_comp, is_cvine, init_vinestr, init_trunclevel, init_mar, init_bicop,
                               methods, threshold, maxit){
  ncol_df <- dim(data)[2]
  if(is.null(total_comp))
    stop("number of components to be fitted must be specified")

  if(!all(is.element(is_cvine,c(NA,0,1))))
    stop("C-vine specification must be NA, 0, or 1")

  if(is.matrix(init_vinestr)){
    if(VineCopula::RVineMatrixCheck(init_vinestr) != 1)
      stop("Vine structure must be a valid R-vine matrix, see VineCopula::RVineMatrix")
  }

  if(!is.na(init_trunclevel)){
    if(!is.numeric(init_trunclevel))
      stop("truncation level must be numeric")

    if(init_trunclevel > ncol_df-1)
      stop("truncation level must be smaller than (number of variables (p) - 1)")
  }

  if(!all(is.element(init_mar, c(NA, 'std', 'norm', 'logis', 'gamma', 'lnorm', 'llogis', 'snorm', 'sstd', 'cauchy', 'weibull'))))
    stop("Marginal distribution set has not yet been implemented")

  if(!all(is.element(init_bicop, c(NA,1,2,3,4,5,6,7,8,10,13,14,16,17,18,20,23,24,26,27,28,30,33,34,36,37,38,40))))
    stop("Bivariate copula family set must be correctly specified, see vcmm")

  if(!all(is.element(methods,c(NA,'kmeans', 'gmm', 'hcVVV'))))
    stop("Initial clustering method set has not yet been implemented")

  if(!is.numeric(threshold))
    stop("threshold must be numeric")

  if(threshold < 0)
    stop("threshold must be larger than or equal to 0")

  if(!is.numeric(maxit))
    stop("maxit must be numeric")

  if(maxit < 1)
    stop("maxit must be larger than or equal to 1")
}

#' internal function
#' @noRd
sim_args_check <- function(dims, obs, margin, margin_pars, RVMs){
  if(!is.numeric(dims))
    stop("dims must be numeric")

  if(dims < 2)
    stop("data must contain at least 2 variables")

  total_comp <- length(obs)
  if(total_comp < 2)
    stop("data must contain at least 2 components")

  var_mar <- dim(margin)[1]
  comp_mar <- dim(margin)[2]

  if(comp_mar != total_comp)
    stop("margin matrix does not match with total number of components")

  for(j in 1:total_comp){
    df <- data.frame(margin_pars[,,j])
    if(!all(sapply(df, is.double)))
      stop("marginal parameters must be numeric")
  }

  npar_par <- dim(margin_pars)[1]
  var_par <- dim(margin_pars)[2]
  comp_par <- dim(margin_pars)[3]

  if(comp_par != total_comp)
    stop("marginal parameters array does not match with total number of components")
}

#' internal function
#' @noRd
mar_RVM_check <- function(margin, margin_pars, RVMs){
  if(is.null(dim(margin)[1])){
    stop("margin matrix must contain at least 2 components")
  }
  else{
    var_mar <- dim(margin)[1]
    comp_mar <- dim(margin)[2]
  }

  if(!all(is.element(c(margin), c('Cauchy', 'Student-t', 'Normal', 'Logistic', 'Gamma',
                                  'Lognormal', 'Loglogistic', 'Skew Normal', 'Skew Student-t', 'Weibull'))))
    stop("marginal distribution has not yet been implemented")

  if(is.null(dim(margin_pars)[1])){
    stop("marginal parameters array must contain at least 2 components")
  }
  else{
    npar_par <- dim(margin_pars)[1]
    var_par <- dim(margin_pars)[2]
    comp_par <- dim(margin_pars)[3]
  }

  if(comp_par != comp_mar | var_par != var_mar)
    stop("marginal parameters array does not match with margin matrix")

  if(npar_par > 4)
    stop("marginal parameters must not be (yet) larger than 4")

  total_comp_rvm <- length(RVMs)

  if(comp_par != total_comp_rvm)
    stop("All components' R-vine matrix must be specified")
}

#' internal function
#' @noRd
dens_args_check <- function(x, mix_probs, total_comp){
  if(is.null(dim(x)[1])){
    if(!is.numeric(x))
      stop("x must be numeric")
  }
  else{
    x <- data.frame(x)
    if(any(sapply(x, is.na)))
      stop("x must not contain any missing values")

    if(!all(sapply(x, is.numeric)))
      stop("x must be numeric")
  }
  if(!is.null(dim(mix_probs)[1]))
    stop("mix_probs must be vector")

  if(!is.numeric(mix_probs))
    stop("mix_probs must be numeric")

  if(length(mix_probs) != total_comp)
    stop("mix_probs must match with total components")

  if(sum(mix_probs) != 1)
    stop("sum of mix_probs must be 1")
}

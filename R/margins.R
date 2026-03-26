#' internal function
#' @noRd
normalize_margin_data <- function(data){
  if(is.list(data)){
    data <- unlist(data, use.names = FALSE)
  }
  if(!is.numeric(data)){
    stop("margin data must be numeric")
  }
  if(length(data) == 0){
    stop("margin data must contain at least one observation")
  }
  if(anyNA(data) || any(!is.finite(data))){
    stop("margin data must be finite and non-missing")
  }
  as.numeric(data)
}

#' internal function
#' @noRd
margin_scale_floor <- function(data){
  scale_ref <- max(1, abs(data), na.rm = TRUE)
  max(sqrt(.Machine$double.eps), scale_ref * 1e-8)
}

#' internal function
#' @noRd
stabilize_positive_density <- function(values){
  values[is.na(values) | values <= 0] <- .Machine$double.xmin
  values[is.infinite(values)] <- .Machine$double.xmax
  values
}

#' internal function
#' @noRd
regularize_margin_params <- function(fam, params, scale_floor){
  params <- as.numeric(params)
  if(fam %in% c("Normal", "Lognormal", "Logistic", "Cauchy", "Skew Normal", "Student-t", "Skew Student-t")){
    params[2] <- max(params[2], scale_floor)
  }
  if(fam %in% c("Gamma", "Loglogistic", "Weibull")){
    params[1] <- max(params[1], scale_floor)
    params[2] <- max(params[2], scale_floor)
  }
  if(fam %in% c("Student-t", "Skew Student-t")){
    params[3] <- max(params[3], 2.0001)
  }
  if(fam == "Skew Normal"){
    params[3] <- max(params[3], scale_floor)
  }
  if(fam == "Skew Student-t"){
    params[4] <- max(params[4], scale_floor)
  }
  params
}

#' internal function
#' @noRd
pdf_cdf_quant_margin <- function(data, dist_name, params, value_name){
  data <- normalize_margin_data(data)
  value_name <- match.arg(value_name, c("cdf", "pdf", "quant"))
  required_pars <- switch(dist_name,
                          "Cauchy" = 2,
                          "Gamma" = 2,
                          "Logistic" = 2,
                          "Loglogistic" = 2,
                          "Lognormal" = 2,
                          "Normal" = 2,
                          "Skew Normal" = 3,
                          "Student-t" = 3,
                          "Skew Student-t" = 4,
                          "Weibull" = 2,
                          NULL)
  if(is.null(required_pars)){
    stop("Unsupported marginal distribution: ", dist_name)
  }
  if(length(params) < required_pars){
    stop("Incorrect number of marginal parameters supplied for ", dist_name)
  }
  params <- as.numeric(params[seq_len(required_pars)])
  if(anyNA(params) || any(!is.finite(params))){
    stop("Marginal parameters must be finite")
  }
  if(value_name == "quant" && any(data < 0 | data > 1)){
    stop("Quantile inputs must lie in [0, 1]")
  }
  if(dist_name %in% c("Cauchy", "Logistic", "Lognormal", "Normal", "Skew Normal", "Student-t", "Skew Student-t") &&
     params[2] <= 0){
    stop("Scale parameters must be strictly positive")
  }
  if(dist_name %in% c("Gamma", "Loglogistic", "Weibull") && any(params[1:2] <= 0)){
    stop("Shape and rate/scale parameters must be strictly positive")
  }
  if(dist_name %in% c("Student-t", "Skew Student-t") && params[3] <= 0){
    stop("Degrees of freedom must be strictly positive")
  }
  if(dist_name == "Skew Normal" && params[3] <= 0){
    stop("Skew parameter must be strictly positive")
  }
  if(dist_name == "Skew Student-t" && params[4] <= 0){
    stop("Skew parameter must be strictly positive")
  }
  values <- switch(dist_name,
                   "Cauchy" = switch(value_name,
                                     cdf = pcauchy(data, location = params[1], scale = params[2]),
                                     pdf = dcauchy(data, location = params[1], scale = params[2]),
                                     quant = qcauchy(data, location = params[1], scale = params[2])),
                   "Gamma" = switch(value_name,
                                    cdf = pgamma(data, shape = params[1], rate = params[2]),
                                    pdf = dgamma(data, shape = params[1], rate = params[2]),
                                    quant = qgamma(data, shape = params[1], rate = params[2])),
                   "Logistic" = switch(value_name,
                                       cdf = plogis(data, location = params[1], scale = params[2]),
                                       pdf = dlogis(data, location = params[1], scale = params[2]),
                                       quant = qlogis(data, location = params[1], scale = params[2])),
                   "Loglogistic" = switch(value_name,
                                          cdf = {
                                            values <- numeric(length(data))
                                            positive_idx <- data > 0
                                            values[positive_idx] <- 1 - 1 / (1 + (data[positive_idx] * params[2])^params[1])
                                            values
                                          },
                                          pdf = {
                                            values <- numeric(length(data))
                                            positive_idx <- data > 0
                                            values[positive_idx] <- (params[1] * params[2] * (data[positive_idx] * params[2])^(params[1] - 1)) /
                                              (1 + (data[positive_idx] * params[2])^params[1])^2
                                            values
                                          },
                                          quant = 1 / params[2] * (data / (1 - data))^(1 / params[1])),
                   "Lognormal" = switch(value_name,
                                        cdf = plnorm(data, meanlog = params[1], sdlog = params[2]),
                                        pdf = dlnorm(data, meanlog = params[1], sdlog = params[2]),
                                        quant = qlnorm(data, meanlog = params[1], sdlog = params[2])),
                   "Normal" = switch(value_name,
                                     cdf = pnorm(data, mean = params[1], sd = params[2]),
                                     pdf = dnorm(data, mean = params[1], sd = params[2]),
                                     quant = qnorm(data, mean = params[1], sd = params[2])),
                   "Skew Normal" = switch(value_name,
                                          cdf = fGarch::psnorm(data, mean = params[1], sd = params[2], xi = params[3]),
                                          pdf = fGarch::dsnorm(data, mean = params[1], sd = params[2], xi = params[3]),
                                          quant = fGarch::qsnorm(data, mean = params[1], sd = params[2], xi = params[3])),
                   "Student-t" = switch(value_name,
                                        cdf = fGarch::pstd(data, mean = params[1], sd = params[2], nu = params[3]),
                                        pdf = fGarch::dstd(data, mean = params[1], sd = params[2], nu = params[3]),
                                        quant = fGarch::qstd(data, mean = params[1], sd = params[2], nu = params[3])),
                   "Skew Student-t" = switch(value_name,
                                             cdf = fGarch::psstd(data, mean = params[1], sd = params[2], nu = params[3], xi = params[4]),
                                             pdf = fGarch::dsstd(data, mean = params[1], sd = params[2], nu = params[3], xi = params[4]),
                                             quant = fGarch::qsstd(data, mean = params[1], sd = params[2], nu = params[3], xi = params[4])),
                   "Weibull" = switch(value_name,
                                      cdf = stats::pweibull(data, shape = params[1], scale = params[2]),
                                      pdf = stats::dweibull(data, shape = params[1], scale = params[2]),
                                      quant = stats::qweibull(data, shape = params[1], scale = params[2])))
  values
}

#' internal function
#' @noRd
fit_margin <- function(data, min_value, margin_fam){
  data <- normalize_margin_data(data)
  min_value <- min(data)
  all_fams <- c('norm', 'logis', 'gamma', 'lnorm', 'llogis', 'cauchy', 'weibull')
  positive_support_fams <- c('gamma', 'lnorm', 'llogis', 'weibull')
  use_all_fams <- length(margin_fam) == 1 && is.na(margin_fam)
  requested_fams <- character(0)
  if(!use_all_fams){
    requested_fams <- unique(as.character(margin_fam[!is.na(margin_fam)]))
  }
  if(use_all_fams){
    fam_set <- if(min_value > 0) all_fams else setdiff(all_fams, positive_support_fams)
  } else{
    fam_set <- intersect(setdiff(requested_fams, c('std', 'snorm', 'sstd')), all_fams)
    if(min_value <= 0){
      fam_set <- setdiff(fam_set, positive_support_fams)
    }
  }
  bic_win <- Inf
  fam <- NULL
  par_mar <- NULL
  scale_floor <- margin_scale_floor(data)
  data_sd <- max(sd(data), scale_floor)
  lower_loc <- min(data)
  upper_loc <- max(data)
  if(lower_loc == upper_loc){
    lower_loc <- lower_loc - data_sd
    upper_loc <- upper_loc + data_sd
  }
  update_best_fit <- function(opt_fit, total_pars, fam_name){
    if(!is.null(opt_fit) && is.finite(opt_fit$value)){
      bic_fit <- 2 * opt_fit$value + log(length(data)) * total_pars
      if(bic_fit < bic_win){
        bic_win <<- bic_fit
        fam <<- fam_name
        par_mar <<- regularize_margin_params(fam_name, opt_fit$par, scale_floor)
      }
    }
  }
  if(length(fam_set) != 0){
    uML_fit <- tryCatch(univariateML::model_select(data, models = fam_set, criterion = 'bic'),
                        error = function(e) NULL)
    if(!is.null(uML_fit) && is.finite(attr(uML_fit, 'logLik'))){
      fam <- attr(uML_fit, 'model')
      par_mar <- regularize_margin_params(fam, as.numeric(uML_fit), scale_floor)
      bic_win <- -2 * attr(uML_fit, 'logLik') + log(length(data)) * length(par_mar)
    }
  }
  if(use_all_fams || any(requested_fams == 'std')){
    par_t <- c(mean(data), data_sd, 10)
    fn_t <- function(par){
      mdt <- stabilize_positive_density(fGarch::dstd(data, par[1], par[2], par[3]))
      -sum(log(mdt))
    }
    opt_t <- tryCatch(optim(par_t, fn_t,
                            lower = c(lower_loc, scale_floor, 2.0001),
                            upper = c(upper_loc, max(100 * data_sd, 10 * scale_floor), 100),
                            method = "L-BFGS-B"),
                      error = function(e) NULL)
    update_best_fit(opt_t, 3, 'Student-t')
  }
  if(use_all_fams || any(requested_fams == 'snorm')){
    par_snorm <- c(mean(data), data_sd, 5)
    fn_snorm <- function(par){
      mdsnorm <- stabilize_positive_density(fGarch::dsnorm(data, par[1], par[2], par[3]))
      -sum(log(mdsnorm))
    }
    opt_snorm <- tryCatch(optim(par_snorm, fn_snorm,
                                lower = c(lower_loc, scale_floor, scale_floor),
                                upper = c(upper_loc, max(100 * data_sd, 10 * scale_floor), 100),
                                method = "L-BFGS-B"),
                          error = function(e) NULL)
    update_best_fit(opt_snorm, 3, 'Skew Normal')
  }
  if(use_all_fams || any(requested_fams == 'sstd')){
    par_sstd <- c(mean(data), data_sd, 10, 5)
    fn_sstd <- function(par){
      mdsstd <- stabilize_positive_density(fGarch::dsstd(data, par[1], par[2], par[3], par[4]))
      -sum(log(mdsstd))
    }
    opt_sstd <- tryCatch(optim(par_sstd, fn_sstd,
                               lower = c(lower_loc, scale_floor, 2.0001, scale_floor),
                               upper = c(upper_loc, max(100 * data_sd, 10 * scale_floor), 100, 100),
                               method = "L-BFGS-B"),
                         error = function(e) NULL)
    update_best_fit(opt_sstd, 4, 'Skew Student-t')
  }
  if(is.null(fam) || is.null(par_mar)){
    stop("No admissible marginal family could be fit to the supplied data")
  }
  result_mar <- list('par_mar'=par_mar, 'fam'=fam)
  result_mar
}

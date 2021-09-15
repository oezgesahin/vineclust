#' internal function
#' @noRd
pdf_cdf_quant_margin <- function(data, dist_name, params, value_name){
  if(typeof(data)=='list'){data <- unlist(data, use.names = FALSE)}
  if(dist_name == 'Cauchy'&& value_name == 'cdf'){
    values <- pcauchy(data, location=params[1], scale=params[2])
  }
  if(dist_name == 'Cauchy'&& value_name == 'pdf'){
    values <- dcauchy(data, location=params[1], scale=params[2])
  }
  if(dist_name == 'Cauchy'&& value_name == 'quant'){
    values <- qcauchy(data, location=params[1], scale=params[2])
  }
  if(dist_name == 'Gamma'&& value_name == 'cdf'){
    values <- pgamma(data, shape=params[1], rate=params[2])
  }
  if(dist_name == 'Gamma'&& value_name == 'pdf'){
    values <- dgamma(data, shape=params[1], rate=params[2])
  }
  if(dist_name == 'Gamma'&& value_name == 'quant'){
    values <- qgamma(data, shape=params[1], rate=params[2])
  }
  if(dist_name == 'Logistic'&& value_name == 'cdf'){
    values <- plogis(data, location=params[1], scale=params[2])
  }
  if(dist_name == 'Logistic'&& value_name == 'pdf'){
    values <- dlogis(data, location=params[1], scale=params[2])
  }
  if(dist_name == 'Logistic'&& value_name == 'quant'){
    values <- qlogis(data, location=params[1], scale=params[2])
  }
  if(dist_name == 'Loglogistic'&& value_name == 'cdf'){
    values <- actuar::pllogis(data, shape=params[1], rate=params[2])
  }
  if(dist_name == 'Loglogistic'&& value_name == 'pdf'){
    values <- actuar::dllogis(data, shape=params[1], rate=params[2])
  }
  if(dist_name == 'Loglogistic'&& value_name == 'quant'){
    values <- actuar::qllogis(data, shape=params[1], rate=params[2])
  }
  if(dist_name == 'Lognormal'&& value_name == 'cdf'){
    values <- plnorm(data, meanlog=params[1], sdlog=params[2])
  }
  if(dist_name == 'Lognormal'&& value_name == 'pdf'){
    values <- dlnorm(data, meanlog=params[1], sdlog=params[2])
  }
  if(dist_name == 'Lognormal'&& value_name == 'quant'){
    values <- qlnorm(data, meanlog=params[1], sdlog=params[2])
  }
  if(dist_name == 'Normal'&& value_name == 'cdf'){
    values <- pnorm(data, mean=params[1], sd=params[2])
  }
  if(dist_name == 'Normal'&& value_name == 'pdf'){
    values <- dnorm(data, mean=params[1], sd=params[2])
  }
  if(dist_name == 'Normal'&& value_name == 'quant'){
    values <- qnorm(data, mean=params[1], sd=params[2])
  }
  if(dist_name == 'Skew Normal'&& value_name == 'cdf'){
    values <- fGarch::psnorm(data, mean=params[1], sd=params[2], xi=params[3])
  }
  if(dist_name == 'Skew Normal'&& value_name == 'pdf'){
    values <- fGarch::dsnorm(data, mean=params[1], sd=params[2], xi=params[3])
  }
  if(dist_name == 'Skew Normal'&& value_name == 'quant'){
    values <- fGarch::qsnorm(data, mean=params[1], sd=params[2], xi=params[3])
  }
  if(dist_name == 'Student-t'&& value_name == 'cdf'){
    values <- fGarch::pstd(data, mean=params[1], sd=params[2], nu=params[3])
  }
  if(dist_name == 'Student-t'&& value_name == 'pdf'){
    values <- fGarch::dstd(data, mean=params[1], sd=params[2], nu=params[3])
  }
  if(dist_name == 'Student-t'&& value_name == 'quant'){
    values <- fGarch::qstd(data, mean=params[1], sd=params[2], nu=params[3])
  }
  if(dist_name == 'Skew Student-t'&& value_name == 'cdf'){
    values <- fGarch::psstd(data, mean=params[1], sd=params[2], nu=params[3], xi=params[4])
  }
  if(dist_name == 'Skew Student-t'&& value_name == 'pdf'){
    values <- fGarch::dsstd(data, mean=params[1], sd=params[2], nu=params[3], xi=params[4])
  }
  if(dist_name == 'Skew Student-t'&& value_name == 'quant'){
    values <- fGarch::qsstd(data, mean=params[1], sd=params[2], nu=params[3], xi=params[4])
  }
  values
}

#' internal function
#' @noRd
fit_margin <- function(data, min_value, margin_fam){
  bic_win <- 10000000
  all_fams <- c('norm', 'logis', 'gamma', 'lnorm', 'llogis', 'cauchy')
  if(min_value>0 && is.na(margin_fam)){fam_set <- all_fams}
  else if (min_value<=0 && is.na(margin_fam)){fam_set <- c('norm', 'logis','cauchy')}
  else{fam_set <- setdiff(margin_fam, c('std', 'snorm', 'sstd'))}
  if(length(fam_set)!=0){
    uML_fit <- univariateML::model_select(data, models = fam_set, criterion = 'bic')
    bic_win <- -2*attr(uML_fit, 'logLik') + log(length(data))*2
    fam <- attr(uML_fit, 'model')
    par_mar <- as.numeric(uML_fit)
  }
  if(any(margin_fam=='std') || is.na(margin_fam)){
    par_t <- c(mean(data), sd(data), 10)
    fn_t <- function(par){
      mdt <- fGarch::dstd(data, par[1], par[2], par[3])
      mdt[which(mdt == 0)] <- 1e-100
      -sum(log(mdt))
    }
    opt_t <- optim(par_t, fn_t,
                   lower = c(min(data), 0.01 * sd(data), 2.0001),
                   upper = c(max(data), 100 * sd(data), 100),
                   method = "L-BFGS-B")
    bic_t <- 2*opt_t$value + log(length(data))*3
    if(bic_t < bic_win){
      bic_win <- bic_t
      fam <- 'Student-t'
      par_mar <- opt_t$par
    }
  }
  if(any(margin_fam=='snorm') || is.na(margin_fam)){
    par_snorm <- c(mean(data), sd(data), 5)
    fn_snorm <- function(par){
      mdsnorm <- fGarch::dsnorm(data, par[1], par[2], par[3])
      mdsnorm[which(mdsnorm == 0)] <- 1e-100
      -sum(log(mdsnorm))
    }
    opt_snorm <- optim(par_snorm, fn_snorm,
                       lower = c(min(data), 0.01 * sd(data), 0.0001),
                       upper = c(max(data), 100 * sd(data), 100),
                       method = "L-BFGS-B")
    bic_snorm <- 2*opt_snorm$value + log(length(data))*3
    if(bic_snorm < bic_win){
      bic_win <- bic_snorm
      fam <- 'Skew Normal'
      par_mar <- opt_snorm$par
    }
  }
  if(any(margin_fam=='sstd') || is.na(margin_fam)){
    par_sstd <- c(mean(data), sd(data), 10, 5)
    fn_sstd <- function(par){
      mdsstd <- fGarch::dsstd(data, par[1], par[2], par[3], par[4])
      mdsstd[which(mdsstd == 0)] <- 1e-100
      -sum(log(mdsstd))
    }
    opt_sstd <- optim(par_sstd, fn_sstd,
                      lower = c(min(data), 0.01 * sd(data), 2.0001, 0.0001),
                      upper = c(max(data), 100 * sd(data), 100, 100),
                      method = "L-BFGS-B")
    bic_sstd <- 2*opt_sstd$value + log(length(data))*4
    if(bic_sstd < bic_win){
      fam <- 'Skew Student-t'
      par_mar <- opt_sstd$par
    }
  }
  result_mar <- list('par_mar'=par_mar, 'fam'=fam)
  result_mar
}

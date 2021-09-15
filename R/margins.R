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
    values <- actuar::pllogis(data, shape=params[1], scale=params[2])
  }
  if(dist_name == 'Loglogistic'&& value_name == 'pdf'){
    values <- actuar::dllogis(data, shape=params[1], scale=params[2])
  }
  if(dist_name == 'Loglogistic'&& value_name == 'quant'){
    values <- actuar::qllogis(data, shape=params[1], scale=params[2])
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
  all_fams <- c('std', 'norm', 'logis', 'gamma', 'lnorm', 'llogis', 'snorm', 'sstd', 'cauchy')
  if(min_value>0 && is.na(margin_fam)){fam_set <- all_fams}
  else if (min_value<=0 && is.na(margin_fam)){fam_set <- c('std', 'norm', 'logis', 'snorm', 'sstd', 'cauchy')}
  else{fam_set <- margin_fam}
  uML_fit <- univariateML::model_select(data, models = fam_set, criterion = 'bic')
  fam <- attr(uML_fit, 'model')
  par_mar <- as.numeric(uML_fit)
  result_mar <- list('par_mar'=par_mar, 'fam'=fam)
  result_mar
}

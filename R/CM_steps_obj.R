#' internal function
#' @noRd
margin_param_count <- function(fam){
  if(fam == "Skew Student-t") return(4L)
  if(fam %in% c("Skew Normal", "Student-t")) return(3L)
  2L
}

#' internal function
#' @noRd
trim_u_data <- function(u_data){
  pmin(pmax(u_data, .Machine$double.eps), 1 - .Machine$double.eps)
}

#' internal function
#' @noRd
encode_margin_params <- function(params, fam, scale_floor){
  params <- regularize_margin_params(fam, params, scale_floor)
  if(fam %in% c("Normal", "Lognormal", "Logistic", "Cauchy")){
    return(c(params[1], log(params[2])))
  }
  if(fam %in% c("Gamma", "Loglogistic")){
    return(log(params[1:2]))
  }
  if(fam == "Skew Normal"){
    return(c(params[1], log(params[2]), log(params[3])))
  }
  if(fam == "Student-t"){
    return(c(params[1], log(params[2]), log(params[3] - 2)))
  }
  if(fam == "Skew Student-t"){
    return(c(params[1], log(params[2]), log(params[3] - 2), log(params[4])))
  }
  stop("Unsupported marginal distribution: ", fam)
}

#' internal function
#' @noRd
decode_margin_params <- function(pars, fam){
  pars <- as.numeric(pars[seq_len(margin_param_count(fam))])
  if(fam %in% c("Normal", "Lognormal", "Logistic", "Cauchy")){
    return(c(pars[1], exp(pars[2])))
  }
  if(fam %in% c("Gamma", "Loglogistic")){
    return(exp(pars[1:2]))
  }
  if(fam == "Skew Normal"){
    return(c(pars[1], exp(pars[2]), exp(pars[3])))
  }
  if(fam == "Student-t"){
    return(c(pars[1], exp(pars[2]), 2 + exp(pars[3])))
  }
  if(fam == "Skew Student-t"){
    return(c(pars[1], exp(pars[2]), 2 + exp(pars[3]), exp(pars[4])))
  }
  stop("Unsupported marginal distribution: ", fam)
}

#' internal function
#' @noRd
assign_margin_params <- function(marginal_par, p, fam, params){
  count <- margin_param_count(fam)
  marginal_par[, p] <- 0
  marginal_par[seq_len(count), p] <- params[seq_len(count)]
  marginal_par
}

#' internal function
#' @noRd
bound_margin_params <- function(fam, params, data, scale_floor){
  data_sd <- max(sd(data), scale_floor)
  max_scale <- max(100 * data_sd, 10 * scale_floor)
  if(fam %in% c("Normal", "Lognormal", "Logistic", "Cauchy", "Skew Normal", "Student-t", "Skew Student-t") &&
     length(params) >= 2){
    params[2] <- min(max(params[2], scale_floor), max_scale)
  }
  if(fam %in% c("Student-t", "Skew Student-t") && length(params) >= 3){
    params[3] <- min(max(params[3], 2.0001), 100)
  }
  if(fam == "Skew Normal" && length(params) >= 3){
    params[3] <- min(max(params[3], scale_floor), 100)
  }
  if(fam == "Skew Student-t" && length(params) >= 4){
    params[4] <- min(max(params[4], scale_floor), 100)
  }
  params
}

#' internal function
#' @noRd
rvine_density <- function(data, vine_structures, family_sets,
                          cop_params_1, cop_params_2){
  RVM <- VineCopula::RVineMatrix(Matrix=vine_structures,family=family_sets,
                                 par=cop_params_1,par2=cop_params_2)
  rvine_pdf <- VineCopula::RVinePDF(data, RVM)
  rvine_pdf
}

#' internal function
#' @noRd
CM_step_mixture_probs <- function(z_values) apply(z_values, 2, mean)

#' internal function
#' @noRd
CM_step_margin_params <- function(pars, data, marginal_par, p, marginal_families,
                                  z_values, vine_structures, family_sets,
                                  cop_params_1, cop_params_2){
  total_features <- dim(data)[2]
  penalty_value <- 1e300
  fam <- marginal_families[p]
  scale_floor <- margin_scale_floor(data[, p])
  decoded_pars <- decode_margin_params(pars, fam)
  if(any(!is.finite(decoded_pars))){
    return(penalty_value)
  }
  decoded_pars <- regularize_margin_params(fam, decoded_pars, scale_floor)
  decoded_pars <- bound_margin_params(fam, decoded_pars, data[, p], scale_floor)
  marginal_par <- assign_margin_params(marginal_par, p, fam, decoded_pars)
  u_data <- tryCatch(
    sapply(1:total_features, function(x) pdf_cdf_quant_margin(data[,x],marginal_families[x],
                                                              marginal_par[,x], 'cdf')),
    error = function(e) NULL
  )
  if(is.null(u_data) || any(!is.finite(u_data))){
    return(penalty_value)
  }
  u_data <- trim_u_data(u_data)
  rvine_densities <- tryCatch(
    rvine_density(u_data, vine_structures, family_sets, cop_params_1, cop_params_2),
    error = function(e) NULL
  )
  if(is.null(rvine_densities) || any(!is.finite(rvine_densities))){
    return(penalty_value)
  }
  rvine_densities <- stabilize_positive_density(rvine_densities)
  margin_densities <- tryCatch(
    sapply(1:total_features, function(x) pdf_cdf_quant_margin(data[,x],marginal_families[x],
                                                              marginal_par[,x], 'pdf')),
    error = function(e) NULL
  )
  if(is.null(margin_densities) || any(!is.finite(margin_densities))){
    return(penalty_value)
  }
  margin_density <- rep(1, dim(margin_densities)[1])
  for(t in 1:total_features){
    margin_density <- margin_density * margin_densities[,t]
  }
  density <- stabilize_positive_density(rvine_densities * margin_density)
  objective_value <- -sum(z_values*log(density))
  if(!is.finite(objective_value)){
    return(penalty_value)
  }
  objective_value
}

#' internal function
#' @noRd
CM_steps <- function(data, vine_structure, family_set, cop_params_j, cop_params_2_j, z_value,
                     marginal_fam, marginal_par, maxit){
  total_features <- dim(data)[2]
  cop_param <- cop_params_j
  cop_param_2 <- cop_params_2_j
  #CM-step 2
  for(p in 1:total_features){
    fam <- marginal_fam[p]
    scale_floor <- margin_scale_floor(data[, p])
    current_pars <- regularize_margin_params(fam, marginal_par[, p], scale_floor)
    current_pars <- bound_margin_params(fam, current_pars, data[, p], scale_floor)
    marginal_par <- assign_margin_params(marginal_par, p, fam, current_pars)
    pars <- encode_margin_params(current_pars, fam, scale_floor)
    opt_margins <- optim(par=pars, CM_step_margin_params, data=data, marginal_par=marginal_par, p=p,
                         marginal_families=marginal_fam, z_values=z_value, vine_structures=vine_structure,
                         family_sets=family_set, cop_params_1=cop_param, cop_params_2=cop_param_2, method = "BFGS",
                         control = list(maxit=maxit))
    optimized_par <- opt_margins$par
    if(any(!is.finite(optimized_par))){
      stop("Non-finite optimized marginal parameters for variable ", p,
           " (", fam, ")")
    }
  optimized_margin_par <- decode_margin_params(optimized_par, fam)
    if(any(!is.finite(optimized_margin_par))){
      stop("Non-finite transformed marginal parameters for variable ", p,
           " (", fam, ")")
    }
    optimized_margin_par <- regularize_margin_params(fam, optimized_margin_par, scale_floor)
    optimized_margin_par <- bound_margin_params(fam, optimized_margin_par, data[, p], scale_floor)
    marginal_par <- assign_margin_params(marginal_par, p, fam, optimized_margin_par)
  }
  #CM-step 3
  udata <- sapply(1:total_features, function(x) pdf_cdf_quant_margin(data[,x],marginal_fam[x],
                                                               marginal_par[,x], 'cdf'))
  udata <- trim_u_data(udata)
  RVM <- VineCopula::RVineMatrix(Matrix=vine_structure,family=family_set, par=cop_param, par2=cop_param_2)
  seq_RVM <-VineCopula::RVineSeqEst(udata, RVM, weights=z_value,  progress=FALSE)
  cop_param <- seq_RVM$par
  cop_param_2 <- seq_RVM$par2
  result <- list("marginal_par"=marginal_par, "cop_param"=cop_param, "cop_param_2"=cop_param_2, "u_data"=udata)
  result
}

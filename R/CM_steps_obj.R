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
  pars <- as.numeric(pars[seq_len(margin_param_count(marginal_families[p]))])
  if(marginal_families[p]=='Normal' || marginal_families[p]=='Lognormal' || marginal_families[p]=='Logistic'){
    marginal_par[1,p] <- pars[1]
    marginal_par[2,p] <- exp(pars[2])
  }
  else if(marginal_families[p]=='Cauchy'){
    marginal_par[1,p] <- pars[1]
    marginal_par[2,p] <- exp(pars[2])
  }
  else if(marginal_families[p]=='Skew Normal'  || marginal_families[p]=='Student-t'){
    marginal_par[1,p] <- pars[1]
    marginal_par[2,p] <- pars[2]
    marginal_par[3,p] <- pars[3]
  }
  else if(marginal_families[p]=='Skew Student-t'){
    marginal_par[1,p] <- pars[1]
    marginal_par[2,p] <- pars[2]
    marginal_par[3,p] <- pars[3]
    marginal_par[4,p] <- pars[4]
  }
  else{
    marginal_par[1,p] <- exp(pars[1])
    marginal_par[2,p] <- exp(pars[2])
  }
  u_data <- sapply(1:total_features, function(x) pdf_cdf_quant_margin(data[,x],marginal_families[x],
                                                                marginal_par[,x], 'cdf'))
  u_data <- trim_u_data(u_data)
  rvine_densities <- rvine_density(u_data, vine_structures, family_sets, cop_params_1, cop_params_2)
  margin_densities <- sapply(1:total_features, function(x) pdf_cdf_quant_margin(data[,x],marginal_families[x],
                                                                         marginal_par[,x], 'pdf'))
  margin_density <- rep(1, dim(margin_densities)[1])
  for(t in 1:total_features){
    margin_density <- margin_density * margin_densities[,t]
  }
  density <- rvine_densities*margin_density
  density[which(density == 0)] <- 1e-100
  -sum(z_values*log(density))
}

#' internal function
#' @noRd
CM_steps <- function(data, vine_structure, family_set, cop_params_j, cop_params_2_j, z_value,
                     marginal_fam, marginal_par, maxit){
  total_features <- dim(data)[2]
  cop_param <- cop_params_j
  cop_param_2 <- cop_params_2_j
  data_min <- apply(data, 2, min)
  data_max <- apply(data, 2, max)
  data_sd <- pmax(apply(data, 2, sd), sqrt(.Machine$double.eps))
  #CM-step 2
  for(p in 1:total_features){
    if(marginal_fam[p]=='Normal' || marginal_fam[p]=='Lognormal' || marginal_fam[p]=='Logistic'){
      marginal_par[2,p] <- log(marginal_par[2,p])
    }
    else if(marginal_fam[p]=='Cauchy'){
      marginal_par[2,p] <- log(marginal_par[2,p])
    }
    else if(marginal_fam[p]!='Skew Normal'  && marginal_fam[p]!='Student-t' && marginal_fam[p]!='Skew Student-t'){
      marginal_par[1,p] <- log(marginal_par[1,p])
      marginal_par[2,p] <- log(marginal_par[2,p])
    }
    pars <- marginal_par[seq_len(margin_param_count(marginal_fam[p])),p]
    if(marginal_fam[p]=='Student-t'){
      opt_margins <- optim(par=pars, CM_step_margin_params, lower = c(data_min[p], 0.01 * data_sd[p], 2.0001),
                           upper = c(data_max[p], 100 * data_sd[p], 100), data=data, marginal_par=marginal_par, p=p,
                           marginal_families=marginal_fam, z_values=z_value, vine_structures=vine_structure,
                           family_sets=family_set, cop_params_1=cop_param, cop_params_2=cop_param_2, method = "L-BFGS-B",
                           control = list(maxit=maxit))
    }
    else if(marginal_fam[p]=='Cauchy'){
      opt_margins <- optim(par=pars, CM_step_margin_params, data=data, marginal_par=marginal_par, p=p,
                           marginal_families=marginal_fam, z_values=z_value, vine_structures=vine_structure,
                           family_sets=family_set, cop_params_1=cop_param, cop_params_2=cop_param_2, method = "BFGS",
                           control = list(maxit=maxit))
    }
    else if(marginal_fam[p]=='Skew Normal'){
      opt_margins <- optim(par=pars, CM_step_margin_params, lower = c(data_min[p], 0.01 * data_sd[p], 0.0001),
                           upper = c(data_max[p], 100 * data_sd[p], 100), data=data, marginal_par=marginal_par, p=p,
                           marginal_families=marginal_fam, z_values=z_value, vine_structures=vine_structure,
                           family_sets=family_set, cop_params_1=cop_param, cop_params_2=cop_param_2, method = "L-BFGS-B",
                           control = list(maxit=maxit))
    }
    else if(marginal_fam[p]=='Skew Student-t'){
      opt_margins <- optim(par=pars, CM_step_margin_params, lower = c(data_min[p], 0.01 * data_sd[p], 2.0001, 0.0001),
                           upper = c(data_max[p], 100 * data_sd[p], 100, 100), data=data, marginal_par=marginal_par, p=p,
                           marginal_families=marginal_fam, z_values=z_value, vine_structures=vine_structure,
                           family_sets=family_set, cop_params_1=cop_param, cop_params_2=cop_param_2, method = "L-BFGS-B",
                           control = list(maxit=maxit))
    }
    else{
      opt_margins <- optim(par=pars, CM_step_margin_params, data=data, marginal_par=marginal_par, p=p,
                           marginal_families=marginal_fam, z_values=z_value, vine_structures=vine_structure,
                           family_sets=family_set, cop_params_1=cop_param, cop_params_2=cop_param_2, method = "BFGS",
                           control = list(maxit=maxit))

    }
    optimized_par <- opt_margins$par
    if(any(!is.finite(optimized_par))){
      stop("Non-finite optimized marginal parameters for variable ", p,
           " (", marginal_fam[p], ")")
    }
    if(marginal_fam[p]=='Normal' || marginal_fam[p]=='Lognormal' || marginal_fam[p]=='Logistic'){
      marginal_par[2,p]  <- exp(optimized_par[2])
      marginal_par[1,p]  <- optimized_par[1]
    }
    else if(marginal_fam[p]=='Cauchy'){
      marginal_par[2,p]  <- exp(optimized_par[2])
      marginal_par[1,p]  <- optimized_par[1]
    }
    else if(marginal_fam[p]=='Skew Normal'  || marginal_fam[p]=='Student-t'){
      marginal_par[1,p]  <- optimized_par[1]
      marginal_par[2,p]  <- optimized_par[2]
      marginal_par[3,p]  <- optimized_par[3]
    }
    else if(marginal_fam[p]=='Skew Student-t'){
      marginal_par[1,p]  <- optimized_par[1]
      marginal_par[2,p]  <- optimized_par[2]
      marginal_par[3,p]  <- optimized_par[3]
      marginal_par[4,p]  <- optimized_par[4]
    }
    else{
      marginal_par[2,p]  <- exp(optimized_par[2])
      marginal_par[1,p]  <- exp(optimized_par[1])
    }
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

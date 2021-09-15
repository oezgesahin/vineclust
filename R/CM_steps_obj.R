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
  pars <- matrix(pars, 4, 1)
  if(marginal_families[p]=='Normal' || marginal_families[p]=='Lognormal' || marginal_families[p]=='Logistic'){
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
  #CM-step 2
  for(p in 1:total_features){
    if(marginal_fam[p]=='Normal' || marginal_fam[p]=='Lognormal' || marginal_fam[p]=='Logistic'){
      marginal_par[2,p] <- log(marginal_par[2,p])
    }
    else if(marginal_fam[p]=='Skew Normal'  || marginal_fam[p]=='Student-t' || marginal_fam[p]=='Skew Student-t'){
      marginal_par <- marginal_par
    }
    else{
      marginal_par[1,p] <- log(marginal_par[1,p])
      marginal_par[2,p] <- log(marginal_par[2,p])
    }
    pars <- marginal_par[,p]
    if(marginal_fam[p]=='Student-t'){
      opt_margins <- optim(par=pars, CM_step_margin_params, lower = c(min(data[,p]), 0.01*sd(data[,p]), 2.0001),
                           upper = c(max(data[,p]), 100*sd(data[,p]), 100), data=data, marginal_par=marginal_par, p=p,
                           marginal_families=marginal_fam, z_values=z_value, vine_structures=vine_structure,
                           family_sets=family_set, cop_params=cop_param, cop_params_2=cop_param_2, method = "L-BFGS-B",
                           control = list(maxit=maxit))
    }
    else if(marginal_fam[p]=='Skew Normal'){
      opt_margins <- optim(par=pars, CM_step_margin_params, lower = c(min(data[,p]), 0.01*sd(data[,p]), 0.0001),
                           upper = c(max(data[,p]), 100*sd(data[,p]), 100), data=data, marginal_par=marginal_par, p=p,
                           marginal_families=marginal_fam, z_values=z_value, vine_structures=vine_structure,
                           family_sets=family_set, cop_params=cop_param, cop_params_2=cop_param_2, method = "L-BFGS-B",
                           control = list(maxit=maxit))
    }
    else if(marginal_fam[p]=='Skew Student-t'){
      opt_margins <- optim(par=pars, CM_step_margin_params, lower = c(min(data[,p]), 0.01*sd(data[,p]), 2.0001, 0.0001),
                           upper = c(max(data[,p]), 100*sd(data[,p]), 100, 100), data=data, marginal_par=marginal_par, p=p,
                           marginal_families=marginal_fam, z_values=z_value, vine_structures=vine_structure,
                           family_sets=family_set, cop_params=cop_param, cop_params_2=cop_param_2, method = "L-BFGS-B",
                           control = list(maxit=maxit))
    }
    else{
      opt_margins <- optim(par=pars, CM_step_margin_params, data=data, marginal_par=marginal_par, p=p,
                           marginal_families=marginal_fam, z_values=z_value, vine_structures=vine_structure,
                           family_sets=family_set, cop_params=cop_param, cop_params_2=cop_param_2, method = "BFGS",
                           control = list(maxit=maxit))

    }
    optimized_par <- opt_margins$par
    if(marginal_fam[p]=='Normal' || marginal_fam[p]=='Lognormal' || marginal_fam[p]=='Logistic'){
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
  RVM <- VineCopula::RVineMatrix(Matrix=vine_structure,family=family_set, par=cop_param, par2=cop_param_2)
  seq_RVM <-VineCopula::RVineSeqEst(udata, RVM, weights=z_value,  progress=FALSE)
  cop_param <- seq_RVM$par
  cop_param_2 <- seq_RVM$par2
  result <- list("marginal_par"=marginal_par, "cop_param"=cop_param, "cop_param_2"=cop_param_2, "u_data"=udata)
  result
}

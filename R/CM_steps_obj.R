#' internal function
#' @noRd
rvine_density <- function(data, vine_structures, family_sets, cop_params_1, cop_params_2){
  RVM <- VineCopula::RVineMatrix(Matrix=vine_structures,family=family_sets,par=cop_params_1,par2=cop_params_2)
  rvine_pdf <- VineCopula::RVinePDF(data, RVM)
  rvine_pdf
}

#' internal function
#' @noRd
CM_step_mixture_probs <- function(z_values) apply(z_values, 2, mean)

#' internal function
#' @noRd
CM_step_margin_params <- function(pars, data, marginal_par, p, marginal_families, z_values,
                                  vine_structures, family_sets, cop_params_1, cop_params_2){
  total_features <- dim(data)[2]
  pars <- matrix(pars, 4, 1)
  if(marginal_families[p]=='Normal' || marginal_families[p]=='Lognormal' || marginal_families[p]=='Logistic'){
    marginal_par[1,p] <- pars[1]
    marginal_par[2,p] <- exp(pars[2])
  }
  else if(marginal_families[p]=='Skew Normal'  || marginal_families[p]=='Student-t'){
    marginal_par[1,p] <- pars[1]
    marginal_par[2,p] <- exp(pars[2])
    marginal_par[3,p] <- exp(pars[3])
  }
  else if(marginal_families[p]=='Skew Student-t'){
    marginal_par[1,p] <- pars[1]
    marginal_par[2,p] <- exp(pars[2])
    marginal_par[3,p] <- exp(pars[3])
    marginal_par[4,p] <- exp(pars[4])
  }
  else{
    marginal_par[1,p] <- exp(pars[1])
    marginal_par[2,p] <- exp(pars[2])
  }
  u_data <- sapply(1:total_features, function(x) pdf_cdf_quant_margin(data[,x],marginal_families[x],
                                                                marginal_par[,x], 'cdf'))
  rvine_densities <- rvine_density(u_data, vine_structures, family_sets, cop_params_1, cop_params_2)
  margin_densities <-sapply(1:total_features, function(x) pdf_cdf_quant_margin(data[,x],marginal_families[x],
                                                                         marginal_par[,x], 'pdf'))
  margin_density <- rep(1, dim(margin_densities)[1])
  for(t in 1:total_features){
    margin_density <- margin_density * margin_densities[,t]
  }
  density <- rvine_densities*margin_density
  density[which(density == 0)] <- 1e-100
  -sum(z_values*log(density))
}

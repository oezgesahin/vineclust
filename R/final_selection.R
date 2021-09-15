#' internal function
#' @noRd
final_selection <- function(data, total_cluster, final_cvine, final_vinestr, final_trunclevel, mix_probs, p_probs,
                            iteration, init_method, final_mar, final_bicop){
  if(is.na(final_cvine)) final_cvine <- 0
  if(is.na(final_trunclevel)) final_trunclevel <- ncol(data) - 1
  if(all(is.na(final_bicop))) final_bicop <- c(1,2,3,4,5,6,7,8,10,13,14,16,17,18,
                                               20,23,24,26,27,28,30,33,34,36,37,38,40)
  data <- cbind(data, apply(p_probs,1,function(x) which(x==max(x))))
  total_obs <- dim(data)[1]
  total_features <- dim(data)[2]-1
  data_cluster <- list()
  one_par <- vector()
  two_par <- vector()
  u_data <- array(0, dim=c(total_obs, total_features, total_cluster))
  cop_params <- array(0, dim=c(total_features, total_features, total_cluster))
  cop_params_2 <- array(0, dim=c(total_features, total_features, total_cluster))
  family_sets <- array(0, dim=c(total_features, total_features, total_cluster))
  vine_structures <- array(0, dim=c(total_features, total_features, total_cluster))
  marginal_fams <- matrix(0,total_features, total_cluster)
  marginal_params <- array(0, dim=c(4, total_features, total_cluster))
  rvine_densities <- matrix(0, total_obs, total_cluster)
  total_margin_dens <- matrix(0, total_obs, total_cluster)
  margin_densities <- array(0, dim=c(total_obs, total_features, total_cluster))
  lik_points <- matrix(0,dim(data)[1],total_cluster)
  for(j in 1:total_cluster){
    data_cluster[[j]] <- data[data[,(total_features+1)] == j,1:total_features]
    for(i in 1:total_features){
      min_value <- min(data_cluster[[j]][,i])
      model_margin <- fit_margin(data_cluster[[j]][,i], min_value, final_mar)
      marginal_fams[i,j] <- model_margin$fam
      marginal_params[1,i,j] <- model_margin$par_mar[1]
      marginal_params[2,i,j] <- model_margin$par_mar[2]
      marginal_params[3,i,j] <- model_margin$par_mar[3]
      marginal_params[4,i,j] <- model_margin$par_mar[4]
    }
    u_data[,,j] <- sapply(1:total_features, function(x) pdf_cdf_quant_margin(data[,x],marginal_fams[x,j],
                                                                       marginal_params[,x,j], 'cdf'))
    if(is.matrix(final_vinestr)){
      fit_rvine <- VineCopula::RVineCopSelect(u_data[data[,(total_features+1)] == j,,j], familyset=final_bicop,
                                              Matrix=final_vinestr, trunclevel=final_trunclevel)
    }else{
      fit_rvine <- VineCopula::RVineStructureSelect(u_data[data[,(total_features+1)] == j,,j], familyset=final_bicop,
                                                    type=final_cvine, trunclevel=final_trunclevel)
    }
    vine_structures[,,j] <- fit_rvine$Matrix
    family_sets[,,j] <- fit_rvine$family
    cop_params[,,j] <- fit_rvine$par
    cop_params_2[,,j] <- fit_rvine$par2
    has_one_par_0 <- apply(family_sets[,,j], 1, function(row) sum(row == 1))
    has_one_par_1 <- apply(family_sets[,,j], 1, function(row) sum(row < 7 & row > 2))
    has_one_par_2 <- apply(family_sets[,,j], 1, function(row) sum(row < 17 & row > 12))
    has_one_par_3 <- apply(family_sets[,,j], 1, function(row) sum(row < 27 & row > 22))
    has_one_par_4 <- apply(family_sets[,,j], 1, function(row) sum(row < 37 & row > 32))
    has_two_par_0 <- apply(family_sets[,,j], 1, function(row) sum(row == 2))
    has_two_par_1 <- apply(family_sets[,,j], 1, function(row) sum(row < 11 & row > 6))
    has_two_par_2 <- apply(family_sets[,,j], 1, function(row) sum(row < 21 & row > 16))
    has_two_par_3 <- apply(family_sets[,,j], 1, function(row) sum(row < 31 & row > 26))
    has_two_par_4 <- apply(family_sets[,,j], 1, function(row) sum(row < 41 & row > 36))
    one_par[j] <- sum(has_one_par_0) + sum(has_one_par_1) + sum(has_one_par_2) + sum(has_one_par_3) + sum(has_one_par_4)
    two_par[j] <- sum(has_two_par_0) + sum(has_two_par_1) + sum(has_two_par_2) + sum(has_two_par_3) + sum(has_two_par_4)
  }
  data <- data[,1:total_features]
  rvine_densities <-sapply(1:total_cluster, function(j) rvine_density(u_data[,,j], vine_structures[,,j], family_sets[,,j],
                                                                      cop_params[,,j], cop_params_2[,,j]))
  for(j in 1:total_cluster){
    density <- rep(1, dim(margin_densities)[1])
    margin_densities[,,j]<-sapply(1:total_features, function(x) pdf_cdf_quant_margin(data[,x],marginal_fams[x,j],
                                                                               marginal_params[,x,j], 'pdf'))
    for(t in 1:total_features){
      density <- density * margin_densities[,t,j]
    }
    total_margin_dens[,j] <- density
  }
  lik_points <-sapply(1:total_cluster, function(j) mix_probs[j]*total_margin_dens[,j]*rvine_densities[,j])
  lik_per_obs <- apply(lik_points, 1, sum)
  z_values <- lik_points/rep(lik_per_obs, total_cluster)
  loglik <- sum(log(lik_per_obs))
  total_mar_pars <- 0
  for(j in 1:total_cluster){
    for(i in 1:total_features){
      if(marginal_fams[i,j]=='Skew Normal'  || marginal_fams[i,j]=='Student-t'){total_mar_pars <- total_mar_pars + 3}
      else if(marginal_fams[i,j]=='Skew Student-t'){total_mar_pars <- total_mar_pars + 4}
      else{total_mar_pars <- total_mar_pars + 2}
    }
  }
  total_cop_pars <- sum(one_par)+2*sum(two_par)
  total_mix_pars <- total_cluster-1
  total_pars <- total_mar_pars + total_cop_pars + total_mix_pars
  bic_cop <- (-2)*loglik + log(total_obs)*total_pars
  class <- apply(z_values,1,function(x) which(x==max(x)))
  const <- 0
  if(length(unique(class)) == total_cluster){
    for(i in 1:total_obs){
      cl <- class[i]
      const <- const + log(z_values[i, cl])
    }
    icl <- bic_cop-2*const
  }
  else{icl <- bic_cop}
  output <- list("loglik"=loglik, "bic"=bic_cop, "icl"=icl, "init_clustering"=init_method, "iteration"=iteration,
                 "total_pars"=total_pars,"mixture_prob"=mix_probs, "margin" = marginal_fams,
                 "marginal_param"=marginal_params, "vine_structure"=vine_structures, "bicop_familyset"=family_sets,
                 "bicop_param"=cop_params, "bicop_param2"=cop_params_2, "z_values"=z_values)
  output
}

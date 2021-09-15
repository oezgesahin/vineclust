#' internal function
#' @noRd
initial_clustering <- function(data, total_cluster, is_cvine, init_vinestr, init_trunclevel, init_mar,
                               init_bicop, clustering_method){
  if(is.na(is_cvine)) is_cvine <- 0
  if(all(is.na(init_bicop))) {init_bicop <- c(1,2,3,4,5,6,7,8,10,13,14,16,17,18,20,23,24,26,27,28,30,33,34,36,37,38,40)}
  total_features <- dim(data)[2]
  total_obs <- dim(data)[1]
  u_data_to_cluster <- list()
  data_to_cluster <- list()
  u_data <- array(0, dim=c(total_obs, total_features, total_cluster))
  cop_params <- array(0, dim=c(total_features, total_features, total_cluster))
  cop_params_2 <- array(0, dim=c(total_features, total_features, total_cluster))
  family_sets <- array(0, dim=c(total_features, total_features, total_cluster))
  vine_structures <- array(0, dim=c(total_features, total_features, total_cluster))
  marginal_fams <- matrix(0,total_features, total_cluster)
  marginal_params <- array(0, dim=c(4, total_features, total_cluster))
  if(clustering_method == 'gmm'){
    gmm_fit <- mclust::Mclust(data, total_cluster, verbose=FALSE)
  }
  if(clustering_method == 'kmeans'){
    kmeans_fit <- kmeans(scale(data), total_cluster)
  }
  if(clustering_method == 'hcVVV'){
    hcVVV_fit <- mclust::hcVVV(data=scale(data), alpha = 1)
    hcVVV_cl <- mclust::hclass(hcVVV_fit, total_cluster)
  }
  mix_probs <- vector()
  for(j in 1:total_cluster){
    if(clustering_method == 'gmm'){
      data_to_cluster[[j]] <- data[gmm_fit$classification == j,]
    }
    if(clustering_method == 'kmeans'){
      data_to_cluster[[j]] <- data[kmeans_fit$cluster == j,]
    }
    if(clustering_method == 'hcVVV'){
      data_to_cluster[[j]] <- data[hcVVV_cl == j,]
    }
    for(i in 1:total_features){
      min_value <- min(data_to_cluster[[j]][,i])
      model_margin <- fit_margin(data_to_cluster[[j]][,i], min_value, init_mar)
      marginal_fams[i,j] <- model_margin$fam
      marginal_params[1,i,j] <- model_margin$par_mar[1]
      marginal_params[2,i,j] <- model_margin$par_mar[2]
      if(model_margin$fam=='Skew Normal' || model_margin$fam=='Student-t'){
        marginal_params[3,i,j] <- model_margin$par_mar[3]
        }
      if(model_margin$fam=='Skew Student-t'){
        marginal_params[3,i,j] <- model_margin$par_mar[3]
        marginal_params[4,i,j] <- model_margin$par_mar[4]
      }
    }
    u_data[,,j] <- sapply(1:total_features, function(x) pdf_cdf_quant_margin(data[,x],marginal_fams[x,j],
                                                                       marginal_params[,x,j], 'cdf'))
    if(clustering_method == 'gmm'){
      if(is.matrix(init_vinestr)){
        fit_rvm <- VineCopula::RVineCopSelect(u_data[gmm_fit$classification == j,,j], familyset=init_bicop, Matrix=init_vinestr, trunclevel=init_trunclevel)
        vine_structures[,,j] <- init_vinestr
      }
      else{
        fit_rvm <- VineCopula::RVineStructureSelect(u_data[gmm_fit$classification == j,,j], familyset=init_bicop, type=is_cvine, trunclevel=init_trunclevel)
        vine_structures[,,j] <- fit_rvm$Matrix
      }
    }
    if(clustering_method == 'kmeans'){
      if(is.matrix(init_vinestr)){
        fit_rvm <- VineCopula::RVineCopSelect(u_data[kmeans_fit$cluster == j,,j], familyset=init_bicop, Matrix=init_vinestr, trunclevel=init_trunclevel)
        vine_structures[,,j] <- init_vinestr
      }
      else{
        fit_rvm <- VineCopula::RVineStructureSelect(u_data[kmeans_fit$cluster == j,,j], familyset=init_bicop, type=is_cvine, trunclevel=init_trunclevel)
        vine_structures[,,j] <- fit_rvm$Matrix
      }
    }
    if(clustering_method == 'hcVVV'){
      if(is.matrix(init_vinestr)){
        fit_rvm <- VineCopula::RVineCopSelect(u_data[hcVVV_cl == j,,j], familyset=init_bicop, Matrix=init_vinestr, trunclevel=init_trunclevel)
        vine_structures[,,j] <- init_vinestr
      }
      else{
        fit_rvm <- VineCopula::RVineStructureSelect(u_data[hcVVV_cl == j,,j], familyset=init_bicop, type=is_cvine, trunclevel=init_trunclevel)
        vine_structures[,,j] <- fit_rvm$Matrix
      }
    }
    family_sets[,,j] <- fit_rvm$family
    cop_params[,,j] <- fit_rvm$par
    cop_params_2[,,j] <- fit_rvm$par2
    mix_probs[j] <- length(data_to_cluster[[j]][,1])/total_obs
  }
  result <- list("u_data"=u_data, "mix_probs"=mix_probs, "marginal_fams" = marginal_fams, "marginal_params"=marginal_params,
                 "family_sets"=family_sets,"vine_structures"=vine_structures, "cop_params"=cop_params, "cop_params_2"=cop_params_2)
  result
}

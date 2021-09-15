#' Vine copula based  mixture model distributions
#'
#' Density and random generation for the vine copula based mixture model distributions.
#' @name vcmm_methods
#' @aliases dvcmm pvcmm rvcmm
#' @param x  A vector of length d or a d-column matrix for evaluation points, where d is the number of variables (p)
#' @param margin A matrix, containing the name of univariate marginal distributions in each component as given in \code{\link{vcmm}}.
#' Rows correspond to variable, columns correspond to component.
#' @param margin_pars An array, specifiying the univariate marginal distributions' parameters in each component. First, second, third dimensions
#' specify the parameter, variable, component, respectively.
#' @param RVMs A list, containing the R-vine copula model of each component.
#' \link[VineCopula]{RVineMatrix} describes the encoding of a R-vine copula model.
#' @param mix_probs A vector of length k (number of components), containing mixture proportion of each component
#'
#' @return `dvcmm()` returns the density, and `rvcmm()` returns random deviates.
#'
#' @examples
#' # Generate data with 3 variables from a vine copula based mixture model
#' # with 2 components, the first/second component has 300/600 observations.
#' dims <- 3
#' obs <- c(300,600)
#' RVMs <- list()
#' RVMs[[1]] <- VineCopula::RVineMatrix(Matrix=matrix(c(1,3,2,0,3,2,0,0,2),dims,dims),
#'                         family=matrix(c(0,3,4,0,0,14,0,0,0),dims,dims),
#'                         par=matrix(c(0,0.5,2.5,0,0,5,0,0,0),dims,dims),
#'                         par2=matrix(sample(0, dims*dims, replace=TRUE),dims,dims))
#' RVMs[[2]] <- VineCopula::RVineMatrix(Matrix=matrix(c(1,3,2,0,3,2,0,0,2), dims,dims),
#'                          family=matrix(c(0,6,5,0,0,13,0,0,0), dims,dims),
#'                          par=matrix(c(0,2,14,0,0,1,0,0,0),dims,dims),
#'                          par2=matrix(sample(0, dims*dims, replace=TRUE),dims,dims))
#' margin <- matrix(c('Normal', 'Gamma', 'Lognormal', 'Lognormal', 'Normal', 'Student-t'), 3, 2)
#' margin_pars <- array(0, dim=c(4, 3, 2))
#' margin_pars[,1,1] <- c(1, 2, 0, 0)
#' margin_pars[,1,2] <- c(1.5, 0.4, 0, 0)
#' margin_pars[,2,1] <- c(1, 0.2, 0, 0)
#' margin_pars[,2,2] <- c(18, 5, 0, 0)
#' margin_pars[,3,1] <- c(0.8, 0.8, 0, 0)
#' margin_pars[,3,2] <- c(4, 2, 5, 0)
#' x_data <- rvcmm(dims, obs, margin, margin_pars, RVMs)
#'
#' @importFrom stats qgamma qlnorm qlogis qnorm qcauchy
#' @importFrom fGarch qstd qsnorm qsstd
#' @importFrom actuar qllogis
#' @rdname vcmm_methods
#' @export

dvcmm <- function(x, margin, margin_pars, RVMs, mix_probs){
  total_features <- dim(margin)[1]
  total_comp <- dim(margin)[2]
  mar_RVM_check(margin, margin_pars, RVMs)
  dens_args_check(x, mix_probs, total_comp)
  if(is.null(dim(x)[1])){
    total_obs <- 1
    x <- matrix(x, 1, total_features)
  }
  else{total_obs <- dim(x)[1]}
  rvine_densities <- matrix(0, total_obs, total_comp)
  total_margin_dens <- matrix(0, total_obs, total_comp)
  margin_densities <- array(0, dim=c(total_obs, total_features, total_comp))
  lik_points <- matrix(0,total_obs,total_comp)
  u_data <- array(0, dim=c(total_obs, total_features, total_comp))
  for(j in 1:total_comp){
    density <- rep(1, total_obs)
    u_data[,,j] <- sapply(1:total_features, function(i) pdf_cdf_quant_margin(x[,i],margin[i,j],
                                                                             margin_pars[,i,j], 'cdf'))
    margin_densities[,,j] <- sapply(1:total_features, function(i) pdf_cdf_quant_margin(x[,i],margin[i,j],
                                                                                       margin_pars[,i,j], 'pdf'))
    for(t in 1:total_features){
      density <- density * margin_densities[,t,j]
    }
    total_margin_dens[,j] <- density
    rvine_densities[,j] <- VineCopula::RVinePDF(u_data[,,j], RVMs[[j]])
  }
  lik_points <- sapply(1:total_comp, function(j) mix_probs[j]*total_margin_dens[,j]*rvine_densities[,j])
  if(dim(x)[1] == 1){lik_per_obs <- sum(lik_points)}
  else{lik_per_obs <- apply(lik_points, 1, sum)}
  lik_per_obs
}

#' @rdname vcmm_methods
#' @param dims An integer, specifying the number of variables
#' @param obs A vector, containing number of observations of each component (sum(obs) = i)
#' @export

rvcmm <- function(dims, obs, margin, margin_pars, RVMs){
  mar_RVM_check(margin, margin_pars, RVMs)
  sim_args_check(dims, obs, margin, margin_pars, RVMs)
  total_comp <- length(obs)
  total_obs <- sum(obs)
  data_to_cluster <- matrix(0, total_obs, dims+1)
  row <- 1
  for(component in 1:total_comp){
    RVM <- RVMs[[component]]
    u_data <- VineCopula::RVineSim(obs[component], RVM)
    x_data_mtr <- matrix(0, obs[component], dims)
    x_data_mtr <- sapply(1:dims, function(x) pdf_cdf_quant_margin(u_data[,x],margin[x,component],
                                                                            margin_pars[,x,component], 'quant'))
    data_to_cluster[row:(obs[component]+row-1),1:dims] <- x_data_mtr
    data_to_cluster[row:(obs[component]+row-1),dims+1] <- rep(component, obs[component])
    row <- row + obs[component]
  }
  x_data <- data.frame(data_to_cluster)
  names(x_data)[ncol(x_data)] <- "comp_id"
  x_data
}





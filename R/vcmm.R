#' Model-Based Clustering with Vine Copulas
#'
#' It fits vine copula based mixture model distributions to the continuous data
#' for a given number of components as described in Sahin and Czado (2021)
#' and use its results for clustering.
#'
#' @param data 	A matrix or data frame of observations. Categorical/discrete variables not (yet) allowed.
#' If a matrix or data frame, rows correspond to observations (i) and columns correspond to variables (p).
#' @param total_comp An integer specifying the numbers of mixture components (clusters)
#' @param is_cvine An integer specifying if the type of components' vine tree structure is C-vine
#' before/after the ECM phase of clustering.
#' * 0 = R-vine (default)
#' * 1 = C-vine
#' @param vinestr A matrix specifying vine tree structures before/after the ECM phase of clustering.
#' The default is automatic selection.
#' \link[VineCopula]{RVineMatrixCheck} checks for a valid R-vine matrix.
#' @param trunclevel An integer showing the level of truncation for vine tree structures before the ECM phase of clustering.
#' The default is 1.
#' @param mar A vector of character strings indicating the parametric univariate marginal distributions
#' to be fitted before/after the ECM phase of clustering.
#' The default is c('cauchy','gamma','llogis','lnorm','logis','norm','snorm','std', 'sstd').
#' Other distributions not (yet) allowed.
#' @param bicop A vector of integers denoting the parametric bivariate copula families to be fitted
#' before/after the ECM phase of clustering.
#' The default is c(1,2,3,4,5,6,7,8,10,13,14,16,17,18,20,23,24,26,27,28,30,33,34,36,37,38,40).
#' \link[VineCopula]{BiCop} describes the available families with their specifications.
#' @param methods A vector of character strings indicating initial clustering method(s) to have a partition
#' for model selection before the ECM phase of clustering. Current options:
#' * 'kmeans' (default)
#' * c('kmeans', 'gmm', 'hcVVV')
#' @param threshold A numeric, stopping the ECM phase of clustering. The default is 1e-4.
#' @param maxit An integer, specifying the maximum number of iterations in the CM-step 2 optimization. The default is 10.
#'
#' @return An object of class vcmm result. It contains the elements
#' \describe{
#' \item{cluster}{the vector with the classification of observations}
#' \item{output}{a list containing the fitted VCMM.}
#' } Use `print.vcmm_res()` to obtain log-likelihood, BIC, ICL, number of estimated parameters, initial clustering method used
#'  and total number of ECM iterations for the fitted VCMM. `summary.vcmm_res()` shows the fitted vine tree structures and
#'  univariate marginal distributions, bivariate copula families with the estimated parameters, as well as
#'  mixture proportions of each component.
#'
#' @references
#' Sahin and Czado (2021), Vine copula mixture models and clustering for non-Gaussian data, Econometrics and Statistics.
#' doi: 10.1016/j.ecosta.2021.08.011
#'
#' @seealso [dvcmm()], [rvcmm()]
#'
#' @examples
#' # Simulate 3-dimensianal data from parametric vine copula based mixture model
#' # with 2 components on x-scale. Each component has 500 observations.
#' dims <- 3
#' obs <- c(500,500)
#' RVMs <- list()
#' RVMs[[1]] <- VineCopula::RVineMatrix(Matrix=matrix(c(1,3,2,0,3,2,0,0,2),dims,dims),
#'                         family=matrix(c(0,3,4,0,0,14,0,0,0),dims,dims),
#'                         par=matrix(c(0,0.8571429,2.5,0,0,5,0,0,0),dims,dims),
#'                         par2=matrix(sample(0, dims*dims, replace=TRUE),dims,dims))
#' RVMs[[2]] <- VineCopula::RVineMatrix(Matrix=matrix(c(1,3,2,0,3,2,0,0,2), dims,dims),
#'                         family=matrix(c(0,6,5,0,0,13,0,0,0), dims,dims),
#'                         par=matrix(c(0,1.443813,11.43621,0,0,2,0,0,0),dims,dims),
#'                         par2=matrix(sample(0, dims*dims, replace=TRUE),dims,dims))
#' margin <- matrix(c('Normal', 'Gamma', 'Lognormal', 'Lognormal', 'Normal', 'Gamma'), 3, 2)
#' margin_pars <- array(0, dim=c(2, 3, 2))
#' margin_pars[,1,1] <- c(1, 2)
#' margin_pars[,1,2] <- c(1.5, 0.4)
#' margin_pars[,2,1] <- c(1, 0.2)
#' margin_pars[,2,2] <- c(18, 5)
#' margin_pars[,3,1] <- c(0.8, 0.8)
#' margin_pars[,3,2] <- c(1, 0.2)
#' x_data <- rvcmm(dims, obs, margin, margin_pars, RVMs)
#'
#' # Example-1: fit parametric 3 dimensional R-vine copula based mixture model with 2 components,
#' # using all univariate marginal distributions and bivariate copula families allowed,
#' # initial partition of k-means
#' \dontrun{
#' fit <- vcmm(x_data[,1:3], total_comp = 2)
#' print(fit)
#' summary(fit)
#' table(x_data$comp_id, fit$cluster)
#' # evaluate the density at (X1, X2, X3) = (1,2,3) for the fitted vcmm
#' # after encoding fitted vine copula
#' RVMs_fitted <- list()
#' RVMs_fitted[[1]] <- VineCopula::RVineMatrix(Matrix=fit$output$vine_structure[,,1],
#'                         family=fit$output$bicop_familyset[,,1],
#'                         par=fit$output$bicop_param[,,1],
#'                         par2=fit$output$bicop_param2[,,1])
#' RVMs_fitted[[2]] <- VineCopula::RVineMatrix(Matrix=fit$output$vine_structure[,,2],
#'                         family=fit$output$bicop_familyset[,,2],
#'                         par=fit$output$bicop_param[,,2],
#'                         par2=fit$output$bicop_param2[,,2])
#' dvcmm(c(1, 2, 3), fit$output$margin, fit$output$marginal_param,
#' RVMs_fitted, fit$output$mixture_prob)
#' }
#'
#' @export
#'
#' @import VineCopula
#' @import mclust
#' @import univariateML
#' @importFrom parallel mclapply
#' @importFrom fGarch psnorm dsnorm pstd dstd psstd dsstd
#' @importFrom actuar dllogis pllogis
#' @importFrom stats dgamma dlnorm dlogis dnorm dcauchy kmeans optim pgamma plnorm plogis pnorm pcauchy sd

vcmm <- function(data, total_comp, is_cvine=NA, vinestr=NA, trunclevel=1, mar=NA, bicop=NA,
                 methods=c('kmeans'),  threshold=0.0001, maxit=10){
  initial_df_check(data)
  initial_args_check(data, total_comp, is_cvine, vinestr, trunclevel, mar, bicop,
                     methods, threshold, maxit)
  final_cvine <- is_cvine
  final_vinestr <- vinestr
  final_trunclevel <- NA
  final_mar <- mar
  final_bicop <- bicop
  winner_bic <- 1000000
  for(method in methods){
    initial_out <- initial_clustering(data, total_comp, is_cvine, vinestr, trunclevel, mar, bicop, method)
    marginal_params <- initial_out$marginal_params
    marginal_fams <- initial_out$marginal_fams
    u_data <- initial_out$u_data
    vine_structures <- initial_out$vine_structures
    family_sets <- initial_out$family_sets
    cop_params <- initial_out$cop_params
    cop_params_2 <- initial_out$cop_params_2
    mix_probs <- initial_out$mix_probs
    total_obs <- dim(data)[1]
    total_features <- dim(data)[2]
    iteration <- 1
    loglik_res <- vector()
    cond <- TRUE
    while(cond==TRUE){
      rvine_densities <- matrix(0, total_obs, total_comp)
      total_margin_dens <- matrix(0, total_obs, total_comp)
      margin_densities <- array(0, dim=c(total_obs, total_features, total_comp))
      lik_points <- matrix(0,dim(data)[1],total_comp)
      rvine_densities <-sapply(1:total_comp, function(j) rvine_density(u_data[,,j], vine_structures[,,j], family_sets[,,j],
                                                                          cop_params[,,j], cop_params_2[,,j]))
      for(j in 1:total_comp){
        density <- rep(1, dim(margin_densities)[1])
        margin_densities[,,j]<-sapply(1:total_features, function(x) pdf_cdf_quant_margin(data[,x],marginal_fams[x,j],
                                                                                         marginal_params[,x,j], 'pdf'))
        for(t in 1:total_features){
          density <- density * margin_densities[,t,j]
        }
        total_margin_dens[,j] <- density
      }
      lik_points <-sapply(1:total_comp, function(j) mix_probs[j]*total_margin_dens[,j]*rvine_densities[,j])
      lik_per_obs <- apply(lik_points, 1, sum)
      lik_per_obs[which(lik_per_obs == 0)] <- 1e-100
      loglik <- sum(log(lik_per_obs))
      loglik_res[iteration] <- loglik
      if(iteration > 2){
        if ((abs(loglik_res[iteration]-loglik_res[iteration-1])/abs(loglik_res[iteration-1]))
            <= threshold){
          cond <- FALSE
          break
        }
      }
      #if(((iteration-1) - floor((iteration-1)/10)*10)==0 & (iteration-1) > 0) cat(iteration-1, "ECM iterations are complete", "\n")
      #E-step
      z_values <- lik_points/rep(lik_per_obs, total_comp)
      #CM-steps:
      #CM-step 1
      mix_probs <- CM_step_mixture_probs(z_values)
      #CM-step 2 and 3
      CMS <- parallel::mclapply(1:total_comp, function(x) CM_steps(data, vine_structures[,,x], family_sets[,,x], cop_params[,,x],
                                                         cop_params_2[,,x], z_values[,x], marginal_fams[,x], marginal_params[,,x],
                                                         maxit))
      for(j in 1:total_comp){
        marginal_params[,,j] <- CMS[[j]]$marginal_par
        cop_params[,,j] <-  CMS[[j]]$cop_param
        cop_params_2[,,j] <- CMS[[j]]$cop_param_2
        u_data[,,j] <- CMS[[j]]$u_data
      }
      iteration <- iteration + 1
    }
    iteration <- iteration - 1
    mix_probs <- CM_step_mixture_probs(z_values)
    final_out <- final_selection(data, total_comp, final_cvine, final_vinestr, final_trunclevel, mix_probs, z_values,
                                 iteration, method, final_mar, final_bicop)
    vcmm_bic <- final_out$bic
    if(vcmm_bic < winner_bic){
      winner_bic <- vcmm_bic
      out <- final_out
      vcmm_class <- apply(out$z_values,1,function(x) which(x==max(x)))
    }
  }
  winner_bic <- round(winner_bic, 0)
  out_list <- list("output"=out, "cluster"=vcmm_class)
  class(out_list) <- "vcmm_res"
  out_list
}

#' @export
print.vcmm_res <- function(x, ...) {
  fit_info(x)
  invisible(x)
}

#' @export
summary.vcmm_res <- function(object, ...) {
  list(
    margins = object$output$margin,
    marginal_pars =  object$output$marginal_param,
    copula = object$output$bicop_familyset,
    copula_first_par = object$output$bicop_param,
    copula_second_par = object$output$bicop_param2,
    vine_structure =  object$output$vine_structure,
    mixture_probs = object$output$mixture_prob
  )
}




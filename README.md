
# vineclust: Model-based clustering with vine copulas

An R package that fits vine copula based mixture model distributions to
the continuous data for a given number of components as proposed in
[VCMM algorithm](https://arxiv.org/pdf/2102.03257.pdf) and use its
results for clustering.

It depends on [VineCopula](https://github.com/tnagler/VineCopula),
[actuar](https://github.com/cran/actuar),
[fGarch](https://github.com/cran/fGarch),
[mclust](https://github.com/cran/mclust), and
[univariateML](https://github.com/JonasMoss/univariateML).

## Installation

You can install the development version from
[GitHub](https://github.com/oezgesahin) with:

``` r
# install.packages("devtools")
devtools::install_github("oezgesahin/vineclust")
```

## Package overview

Below is an overview of some functions and features.

  - `vcmm()`: fits vine copula based mixture model distributions to the
    continuous data for a given number of components. Returns an object
    of class `vcmm_res()`. The class has the following methods:
      - `print`: a brief overview of the model statistics.
      - `summary`: list of fitted model components, including selected
        vine tree structures, bivariate copula families, univariate
        marginal distributions, and estimated parameters.
  - `dvcmm(), rvcmm()`: density and random generation for the vine
    copula based mixture model distributions.

### Bivariate copula families

This package works with a wide range of parametric bivariate copula
families for bivarite or multivariate clustering using vine copulas.
Specifically, it allows fitting elliptical (Gaussian, Student-t) and
Archimedean (Clayton, Gumbel, Frank, Joe, BB1, BB6, and BB8) copulas
with their possible 90, 180, 270 degrees rotations to cover a large
range of dependence patterns. Their encoding is detailed on
[VineCopula](https://github.com/tnagler/VineCopula).

### Univariate marginal distributions

This package currently includes following unimodal univariate marginal
distributions.

  - `cauchy(a,b)`: Cauchy distribution with location parameter a and
    scale parameter b,
  - `gamma(a,b)`: gamma distribution with shape parameter a and rate
    parameter b,
  - `llogis(a,b)`: log-logistic distribution with shape parameter a and
    scale parameter b,
  - `lnorm(a,b)`: log-normal distribution with mean parameter a and
    standard deviation parameter b on the logarithmic scale,
  - `logis(a,b)`: logistic distribution with location parameter a and
    scale parameter b,
  - `norm(a,b)`: normal distribution with mean parameter a and standard
    deviation parameter b,
  - `snorm(a,b,c)`: skew normal distribution with location parameter a,
    scale parameter b, and skewness parameter c.
  - `std(a,b,c)`: Student’s t distribution with location parameter a,
    scale parameter b, and shape parameter c,
  - `sstd(a,b,c,d)`: skew Student’s t distribution with location
    parameter a, scale parameter b, shape parameter c, and skewness
    parameter d.

### Initial partition methods

This package currently implements following partition approaches to have
starting values.

  - `kmeans`: performs k-means clustering (Hartigan-Wong) on given data
    after scaling,
  - `hcVVV`: performs model-based hierarchical clustering on given data
    after scaling,
  - `gmm`: performs model-based clustering with Gaussian mixture models
    on given data.

## Usage

``` r
library(vineclust)
# data from UCI Machine Learning Repository 
data_wisc <- read.csv("http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data", header = FALSE)
```

### Model fitting

``` r
 # R-vine copula based mixture model with total components 2
fit <- vcmm(data=data_wisc[,c(15,27,29,30)], total_comp=2)
 # Model statistics
print(fit)
#> logLik = 2156.03   BIC = -4210.55   ICL = -4174.56   npars = 16   initial clustering = kmeans   ECM iterations = 5
# Fitted vine copula distributions
summary(fit) 
#> $margins
#>      [,1]          [,2]         
#> [1,] "Lognormal"   "Lognormal"  
#> [2,] "Lognormal"   "Gamma"      
#> [3,] "Lognormal"   "Skew Normal"
#> [4,] "Skew Normal" "Normal"     
#> 
#> $marginal_pars
#> , , 1
#> 
#>           [,1]       [,2]       [,3]      [,4]
#> [1,] 1.3147399 -1.9226732 -0.7996311 0.1869493
#> [2,] 0.5059076  0.1380584  0.3534658 0.0393310
#> [3,]        NA         NA         NA 1.8611051
#> [4,]        NA         NA         NA        NA
#> 
#> , , 2
#> 
#>           [,1]      [,2]        [,3]       [,4]
#> [1,] 0.6408718  42.03979   0.1447363 0.07210234
#> [2,] 0.3795854 340.54956   0.1093514 0.03237977
#> [3,]        NA        NA 306.7345684         NA
#> [4,]        NA        NA          NA         NA
#> 
#> 
#> $copula
#> , , 1
#> 
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0    0    0
#> [2,]    1    0    0    0
#> [3,]   16    5    0    0
#> [4,]    5    1    5    0
#> 
#> , , 2
#> 
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0    0    0
#> [2,]   33    0    0    0
#> [3,]   33    2    0    0
#> [4,]   10    3   20    0
#> 
#> 
#> $copula_first_par
#> , , 1
#> 
#>           [,1]       [,2]     [,3] [,4]
#> [1,]  0.000000  0.0000000 0.000000    0
#> [2,] -0.376350  0.0000000 0.000000    0
#> [3,]  1.107328 -1.4521446 0.000000    0
#> [4,]  1.990249  0.3762734 4.296215    0
#> 
#> , , 2
#> 
#>             [,1]        [,2]     [,3] [,4]
#> [1,]  0.00000000  0.00000000 0.000000    0
#> [2,] -0.06998957  0.00000000 0.000000    0
#> [3,] -0.19001573 -0.06769775 0.000000    0
#> [4,]  1.27858180  0.42239359 4.071188    0
#> 
#> 
#> $copula_second_par
#> , , 1
#> 
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0    0    0
#> [2,]    0    0    0    0
#> [3,]    0    0    0    0
#> [4,]    0    0    0    0
#> 
#> , , 2
#> 
#>           [,1]     [,2]      [,3] [,4]
#> [1,] 0.0000000  0.00000 0.0000000    0
#> [2,] 0.0000000  0.00000 0.0000000    0
#> [3,] 0.0000000 20.07937 0.0000000    0
#> [4,] 0.8923816  0.00000 0.9621593    0
#> 
#> 
#> $vine_structure
#> , , 1
#> 
#>      [,1] [,2] [,3] [,4]
#> [1,]    2    0    0    0
#> [2,]    1    1    0    0
#> [3,]    4    3    3    0
#> [4,]    3    4    4    4
#> 
#> , , 2
#> 
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    0    0
#> [2,]    3    2    0    0
#> [3,]    2    3    3    0
#> [4,]    4    4    4    4
#> 
#> 
#> $mixture_probs
#> [1] 0.3677205 0.6322795
# Evaluate the density of the fitted model at (2.747, 0.1467, 0.13, 0.05334)
RVMs_fitted <- list()
RVMs_fitted[[1]] <- VineCopula::RVineMatrix(Matrix=fit$output$vine_structure[,,1],
                        family=fit$output$bicop_familyset[,,1],
                        par=fit$output$bicop_param[,,1],
                        par2=fit$output$bicop_param2[,,1])
RVMs_fitted[[2]] <- VineCopula::RVineMatrix(Matrix=fit$output$vine_structure[,,2],
                        family=fit$output$bicop_familyset[,,2],
                        par=fit$output$bicop_param[,,2],
                        par2=fit$output$bicop_param2[,,2])
dvcmm(c(2.747, 0.1467, 0.13, 0.05334), fit$output$margin, fit$output$marginal_param, RVMs_fitted, fit$output$mixture_prob)
#> [1] 21.67364
```

``` r
# C-vine copula based mixture model
fit_cvine <- vcmm(data=data_wisc[,c(15,27,29,30)], total_comp=2, is_cvine=1)  
#> 10 ECM iterations are complete
# Confusion matrix w.r.t. true classification
table(fit_cvine$cluster, data_wisc$V2) 
#>    
#>       B   M
#>   1 324  18
#>   2  33 194
```

``` r
# Fit only bivariate Clayton copula for pairs of variables in both components
fit_clayton <- vcmm(data=data_wisc[,c(15,27,29,30)], total_comp=2, bicop=c(3))  
```

``` r
# Fit only univariate normal distribution for margins in both components
fit_norm <- vcmm(data=data_wisc[,c(15,27,29,30)], total_comp=2, mar=c("norm"))  
```

``` r
# Fix vine tree structures of both components 
fit_fix_vinestr <- vcmm(data=data_wisc[,c(15,27,29,30)], total_comp=2, vinestr=matrix(c(1,2,3,4,0,2,4,3,0,0,4,3,0,0,0,3),4,4)) 
#> 10 ECM iterations are complete
```

``` r
# Run ECM iterations shorter with a smaller threshold than the default threshold
fit_sthr <- vcmm(data=data_wisc[,c(15,27,29,30)], total_comp=2, threshold=0.001)  
```

``` r
# Use a different initial partition approach from k-means
fit_best_init <- vcmm(data=data_wisc[,c(15,27,29,30)], total_comp=2, methods=c("gmm"))  
```

### Simulation

``` r
# Simulation setup given in Section 5.2 of the paper at https://arxiv.org/pdf/2102.03257.pdf
dims <- 3
obs <- c(500,500) 
RVMs <- list()
RVMs[[1]] <- VineCopula::RVineMatrix(Matrix=matrix(c(1,3,2,0,3,2,0,0,2),dims,dims),
                        family=matrix(c(0,3,4,0,0,14,0,0,0),dims,dims),
                        par=matrix(c(0,0.8571429,2.5,0,0,5,0,0,0),dims,dims),
                        par2=matrix(sample(0, dims*dims, replace=TRUE),dims,dims)) 
RVMs[[2]] <- VineCopula::RVineMatrix(Matrix=matrix(c(1,3,2,0,3,2,0,0,2), dims,dims),
                        family=matrix(c(0,6,5,0,0,13,0,0,0), dims,dims),
                        par=matrix(c(0,1.443813,11.43621,0,0,2,0,0,0),dims,dims),
                        par2=matrix(sample(0, dims*dims, replace=TRUE),dims,dims))
margin <- matrix(c('Normal', 'Gamma', 'Lognormal', 'Lognormal', 'Normal', 'Gamma'), 3, 2) 
margin_pars <- array(0, dim=c(2, 3, 2))
margin_pars[,1,1] <- c(1, 2)
margin_pars[,1,2] <- c(1.5, 0.4)
margin_pars[,2,1] <- c(1, 0.2)
margin_pars[,2,2] <- c(18, 5)
margin_pars[,3,1] <- c(0.8, 0.8)
margin_pars[,3,2] <- c(1, 0.2)
x_data <- rvcmm(dims, obs, margin, margin_pars, RVMs)
```

## Contact

Please contact <ozge.sahin@tum.de> if you have any questions.

## References - Citation

Please cite our paper if you use the algorithm or code in your own work:

Sahin, {"O}., & Czado, C. (2021). Vine copula mixture models and
clustering for non-gaussian data. arXiv:2102.03257.
[pdf](https://arxiv.org/pdf/2102.03257.pdf)

RVMs <- list()
dims <- 3
RVMs[[1]] <- VineCopula::RVineMatrix(Matrix=matrix(c(1,3,2,0,3,2,0,0,2),dims,dims),
                                     family=matrix(c(0,3,4,0,0,14,0,0,0),dims,dims),
                                     par=matrix(c(0,0.5,2.5,0,0,5,0,0,0),dims,dims),
                                     par2=matrix(sample(0, dims*dims, replace=TRUE),dims,dims))
RVMs[[2]] <- VineCopula::RVineMatrix(Matrix=matrix(c(1,3,2,0,3,2,0,0,2), dims,dims),
                                     family=matrix(c(0,6,5,0,0,13,0,0,0), dims,dims),
                                     par=matrix(c(0,2,14,0,0,1,0,0,0),dims,dims),
                                     par2=matrix(sample(0, dims*dims, replace=TRUE),dims,dims))
margin <- matrix(c('Normal', 'Gamma', 'Lognormal', 'Lognormal', 'Normal', 'Student-t'), 3, 2)
margin_pars <- array(0, dim=c(4, 3, 2))
margin_pars[,1,1] <- c(1, 2, 0, 0)
margin_pars[,1,2] <- c(1.5, 0.4, 0, 0)
margin_pars[,2,1] <- c(1, 0.2, 0, 0)
margin_pars[,2,2] <- c(18, 5, 0, 0)
margin_pars[,3,1] <- c(0.8, 0.8, 0, 0)
margin_pars[,3,2] <- c(4, 2, 5, 0)

set.seed(111)
x_data <- rvcmm(dims=3, obs=c(100,200), margin, margin_pars, RVMs)

test_that("d/r vcmm work", {
  expect_false(all(rvcmm(dims=3, obs=c(100,200), margin, margin_pars, RVMs) ==
                     rvcmm(dims=3, obs=c(100,200), margin, margin_pars, RVMs)))
  set.seed(111)
  expect_true(all(x_data == rvcmm(dims=3, obs=c(100,200), margin=margin, margin_pars=margin_pars, RVMs=RVMs)))
  expect_gte(min(dvcmm(x_data, margin, margin_pars, RVMs, c(1/3,2/3))), 0)
})

test_that("catches incorrect arguments", {
  expect_error(dvcmm(x_data, margin[,1], margin_pars, RVMs, c(1/3, 2/3)))
  expect_error(dvcmm(x_data, margin, margin_pars[,,1], RVMs, c(1/3, 2/3)))
  expect_error(dvcmm(x_data, margin, margin_pars, RVMs[[1]], c(1/3, 2/3)))
  expect_error(dvcmm(x_data, margin, margin_pars, RVMs, c(1/3, 1/3, 1/3)))
  expect_error(dvcmm(x_data, margin, margin_pars, RVMs, c(1/3, 1)))
  expect_error(dvcmm(x_data,  matrix(c('xyz', 'xyz', 'xyz', 'xyz', 'xyz', 'xyz'), 3, 2), margin_pars, RVMs, c(1/3, 2/3)))
  expect_error(dvcmm(matrix(c(1,2,3,1,2,NA),2,3), margin, margin_pars, RVMs, c(1/3, 2/3)))
  expect_error(dvcmm(matrix(c(1,2,3,1,2,'x'),2,3), margin, margin_pars, RVMs, c(1/3, 2/3)))
  expect_error(dvcmm(matrix(c(1,2,3,1,2,3),2,3), margin, margin_pars, RVMs, 1))
  expect_error(rvcmm(dims='x', obs=c(100,200), margin, margin_pars, RVMs))
  expect_error(rvcmm(dims=999, obs=c(100,200), margin, margin_pars, RVMs))
  expect_error(rvcmm(dims=3, obs=c(100,200,300), margin, margin_pars, RVMs))
  expect_error(rvcmm(dims=1, obs=c(100,200), margin, margin_pars, RVMs))
  expect_error(rvcmm(dims=3, obs=c(100,200), margin[,1], margin_pars, RVMs))
  expect_error(rvcmm(dims=3, obs=c(100,200), margin, margin_pars[,,1], RVMs))
  expect_error(rvcmm(dims=3, obs=c(100,200), margin, margin_pars, RVMs[[1]]))
  expect_error(rvcmm(dims=3, obs=c(100,200),  matrix(c('xyz', 'xyz', 'xyz', 'xyz', 'xyz', 'xyz'), 3, 2), margin_pars, RVMs))
})

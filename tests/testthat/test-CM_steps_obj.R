# simulate data
set.seed(111)
cl1 <- sapply(1:2, function(i) rnorm(100, -10, 5))
cl2 <- sapply(1:2, function(i) rnorm(100, 10, 5))
x_data <- data.frame(rbind(cl1, cl2))

test_that("CM-step 1 works", {
  fit <- CM_step_mixture_probs(matrix(c(rep(0.9, 100), rep(0.1, 100), rep(0.1, 100), rep(0.9, 100)), 200, 2))
  expect_type(fit, "double")
  expect_lte(max(fit), 1)
  expect_gte(min(fit), 0)
  expect_equal(sum(fit), 1)
})


test_that("CM-steps 2 and 3 work", {
  fit <- CM_steps(data=x_data, vine_structure=matrix(c(1,2,0,2),2,2), family_set=matrix(c(0,1,0,0),2,2),
                  cop_params_j=matrix(c(0,0.7,0,0),2,2), cop_params_2_j=matrix(0,2,2),
                  z_value=c(rep(0.9, 120), rep(0.1, 80)),
                  marginal_fam=c('Normal', 'Skew Normal'),
                  marginal_par=matrix(c(-9, 4.5, 0, 0, -3, 2, 4, 0), 4, 2), maxit=2)
  expect_identical(
    names(fit),
    c(
      "marginal_par", "cop_param", "cop_param_2", "u_data"
    )
  )
  expect_type(fit$marginal_par, "double")
  expect_type(fit$cop_param, "double")
  expect_type(fit$cop_param_2, "double")
  expect_type(fit$u_data, "double")
})

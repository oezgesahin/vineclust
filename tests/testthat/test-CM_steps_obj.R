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

test_that("log-scale parameterization round-trips across supported families", {
  family_params <- list(
    "Normal" = c(1, 2),
    "Lognormal" = c(0.2, 0.8),
    "Logistic" = c(-1, 1.5),
    "Cauchy" = c(-3, 2.5),
    "Gamma" = c(2, 4),
    "Loglogistic" = c(1.3, 0.7),
    "Skew Normal" = c(0.5, 1.2, 2),
    "Student-t" = c(-0.5, 1.5, 6),
    "Skew Student-t" = c(0.3, 1.1, 7, 1.8)
  )
  for(fam in names(family_params)){
    params <- family_params[[fam]]
    encoded <- encode_margin_params(params, fam, 1e-8)
    decoded <- decode_margin_params(encoded, fam)
    expect_equal(decoded, params, tolerance = 1e-10)
  }
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
  expect_lte(fit$marginal_par[2, 2], max(100 * sd(x_data[[2]]), 10 * margin_scale_floor(x_data[[2]])))
  expect_lte(fit$marginal_par[3, 2], 100)
})

test_that("CM-step keeps log-scale-constrained parameters admissible", {
  set.seed(333)
  mixed_data <- data.frame(
    x1 = rgamma(200, shape = 2, rate = 1),
    x2 = rt(200, df = 6) + 1
  )
  fit <- CM_steps(
    data = mixed_data,
    vine_structure = matrix(c(1, 2, 0, 2), 2, 2),
    family_set = matrix(c(0, 1, 0, 0), 2, 2),
    cop_params_j = matrix(c(0, 0.2, 0, 0), 2, 2),
    cop_params_2_j = matrix(0, 2, 2),
    z_value = rep(0.5, 200),
    marginal_fam = c("Gamma", "Student-t"),
    marginal_par = matrix(c(2, 1, 0, 0, 1, 1.5, 6, 0), 4, 2),
    maxit = 2
  )
  expect_gt(fit$marginal_par[1, 1], 0)
  expect_gt(fit$marginal_par[2, 1], 0)
  expect_gt(fit$marginal_par[2, 2], 0)
  expect_gt(fit$marginal_par[3, 2], 2)
})

test_that("CM-step handles negative-location Cauchy margins", {
  set.seed(222)
  cauchy_data <- data.frame(
    x1 = rcauchy(200, location = -40, scale = 20),
    x2 = rnorm(200, mean = 3, sd = 1.5)
  )
  fit <- CM_steps(
    data = cauchy_data,
    vine_structure = matrix(c(1, 2, 0, 2), 2, 2),
    family_set = matrix(c(0, 1, 0, 0), 2, 2),
    cop_params_j = matrix(c(0, 0.2, 0, 0), 2, 2),
    cop_params_2_j = matrix(0, 2, 2),
    z_value = rep(0.5, 200),
    marginal_fam = c("Cauchy", "Normal"),
    marginal_par = matrix(c(-35, 18, 0, 0, 3, 1.5, 0, 0), 4, 2),
    maxit = 2
  )
  expect_true(all(is.finite(fit$marginal_par[, 1])))
  expect_lt(fit$marginal_par[1, 1], 0)
  expect_gt(fit$marginal_par[2, 1], 0)
})

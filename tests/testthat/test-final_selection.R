# simulate data
set.seed(111)
cl1 <- sapply(1:2, function(i) rnorm(100, -10, 5))
cl2 <- sapply(1:2, function(i) rnorm(100, 10, 5))
x_data <- data.frame(rbind(cl1, cl2))

test_that("final model selection works properly", {
  fit <- final_selection(data=x_data, total_cluster=2, final_cvine=NA, final_vinestr=NA,
                            final_trunclevel=NA, mix_probs=c(0.5, 0.5),
                            p_probs=matrix(c(rep(0.9, 100), rep(0.1, 100), rep(0.1, 100), rep(0.9, 100)), 200, 2),
                            iteration=1, init_method="gmm", final_mar=NA, final_bicop=NA)
  expect_identical(
    names(fit),
    c(
      "loglik", "bic", "icl", "init_clustering", "iteration", "total_pars",
      "mixture_prob", "margin","marginal_param", "vine_structure", "bicop_familyset",
      "bicop_param","bicop_param2", "z_values"
    )
  )
  expect_type(fit$loglik, "double")
  expect_type(fit$bic, "double")
  expect_type(fit$icl, "double")
  expect_type(fit$init_clustering, "character")
  expect_type(fit$iteration, "double")
  expect_type(fit$total_pars, "double")
  expect_type(fit$mixture_prob, "double")
  expect_type(fit$margin, "character")
  expect_type(fit$marginal_param, "double")
  expect_type(fit$vine_structure, "double")
  expect_type(fit$bicop_familyset, "double")
  expect_type(fit$bicop_param, "double")
  expect_type(fit$bicop_param2, "double")
  expect_type(fit$z_values, "double")
})

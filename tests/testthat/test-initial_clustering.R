# simulate data
set.seed(111)
cl1 <- sapply(1:2, function(i) rnorm(100, -10, 5))
cl2 <- sapply(1:2, function(i) rnorm(100, 10, 5))
x_data <- data.frame(rbind(cl1, cl2))

test_that("other initial methods work properly", {
  for (met in c("gmm", "hcVVV")) {
    fit <- initial_clustering(data=x_data, total_cluster=2, is_cvine=NA, init_vinestr=NA,
                              init_trunclevel=1, init_mar=NA,
                              init_bicop=NA, clustering_method=met)
    expect_identical(
      names(fit),
      c(
        "u_data", "mix_probs", "marginal_fams", "marginal_params",
        "family_sets", "vine_structures",  "cop_params",  "cop_params_2"
      )
    )
    expect_type(fit$u_data, "double")
    expect_type(fit$mix_probs, "double")
    expect_type(fit$marginal_fams, "character")
    expect_type(fit$marginal_params, "double")
    expect_type(fit$family_sets, "double")
    expect_type(fit$vine_structures, "double")
    expect_type(fit$cop_params, "double")
    expect_type(fit$cop_params_2, "double")
  }
})




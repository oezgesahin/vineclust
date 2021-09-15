# simulate data
set.seed(111)
cl1 <- sapply(1:2, function(i) rnorm(100, -10, 5))
cl2 <- sapply(1:2, function(i) rnorm(100, 10, 5))
x_data <- data.frame(rbind(cl1, cl2))

test_that("catches incorrect input", {
  expect_error(vcmm(x_data[,1], total_comp=2))
  expect_error(vcmm(matrix(2,2,2,'c'),total_comp=2))
  expect_error(vcmm(matrix(2,2,2,'c'),total_comp=2))
  expect_error(vcmm(matrix(2,2,2,NA),total_comp=2))
  expect_error(vcmm(total_comp=2))
  expect_error(vcmm(x_data))
  expect_error(vcmm(x_data, total_comp=2, is_cvine=2))
  expect_error(vcmm(x_data, total_comp=2, trclevel='xyz'))
  expect_error(vcmm(x_data, total_comp=2, trclevel=9999))
  expect_error(vcmm(x_data, total_comp=2, mar='whatisthis'))
  expect_error(vcmm(x_data, total_comp=2, mar=1))
  expect_error(vcmm(x_data, total_comp=2, bicop='xyz'))
  expect_error(vcmm(x_data, total_comp=2, bicop=9999))
  expect_error(vcmm(x_data, total_comp=2, methods = "xyz"))
  expect_error(vcmm(x_data, total_comp=2, methods = 9999))
  expect_error(vcmm(x_data, total_comp=2, threshold=-1))
  expect_error(vcmm(x_data, total_comp=2, threshold='xyz'))
  expect_error(vcmm(x_data, total_comp=2, maxit=-1))
  expect_error(vcmm(x_data, total_comp=2, maxit='xyz'))
})

fit <- vcmm(x_data, 2)

test_that("returns correct 'vcmm_res' object", {
  expect_s3_class(fit, "vcmm_res")
  expect_identical(
    names(fit),
    c(
      "output", "cluster"
    )
  )
})

test_that("proper print/summary generic functions", {
  expect_output(print(fit))
  svcmm <- summary(fit)
  expect_type(svcmm$margins, "character")
  expect_type(svcmm$marginal_pars, "double")
  expect_type(svcmm$copula, "double")
  expect_type(svcmm$copula_first_par, "double")
  expect_type(svcmm$copula_second_par, "double")
  expect_type(svcmm$vine_structure, "double")
  expect_type(svcmm$mixture_probs, "double")
})

test_that("fixed mar/bicop selection works", {
  expect_type(
    vcmm(x_data, total_comp=2, mar='norm', bicop=1)$output$margin,
    "character"
  )
})

test_that("other threshold works", {
  expect_type(
    vcmm(x_data, total_comp=2, threshold=0.1)$output$margin,
    "character"
  )
})


test_that("other maxit works", {
  expect_type(
    vcmm(x_data, total_comp=2, maxit=2)$output$margin,
    "character"
  )
})


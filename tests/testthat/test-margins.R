# simulate data
set.seed(111)
x_data <- rgamma(100, 1, 2)
fit <- fit_margin(x_data, 1, NA)
pars <- c(1,2)
prob <- 0.5

test_that("marginal fitting returns correct output", {
  expect_identical(
    names(fit),
    c(
      "par_mar", "fam"
    )
  )
})

test_that("marginal fitting work", {
  expect_type(fit$par_mar, "double")
  expect_type(fit$fam, "character")
})


test_that("density functions for margins work", {
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Cauchy", pars, "pdf")),
    0
  )
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Gamma", pars, "pdf")),
    0
  )
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Logistic", pars, "pdf")),
    0
  )
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Loglogistic", pars, "pdf")),
    0
  )
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Lognormal", pars, "pdf")),
    0
  )
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Normal", pars, "pdf")),
    0
  )
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Weibull", pars, "pdf")),
    0
  )
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Skew Normal", c(1,2,1), "pdf")),
    0
  )
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Student-t", c(1,2,3), "pdf")),
    0
  )
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Skew Student-t", c(1,2,3,4), "pdf")),
    0
  )
})

test_that("distribution functions for margins work", {
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Cauchy", pars, "cdf")),
    0
  )
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Gamma", pars, "cdf")),
    0
  )
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Logistic", pars, "cdf")),
    0
  )
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Loglogistic", pars, "cdf")),
    0
  )
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Lognormal", pars, "cdf")),
    0
  )
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Normal", pars, "cdf")),
    0
  )
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Weibull", pars, "cdf")),
    0
  )
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Skew Normal", c(1,2,1), "cdf")),
    0
  )
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Student-t", c(1,2,3), "cdf")),
    0
  )
  expect_gte(
    min(pdf_cdf_quant_margin(x_data, "Skew Student-t", c(1,2,3,4), "cdf")),
    0
  )
  expect_lte(
    max(pdf_cdf_quant_margin(x_data, "Cauchy", pars, "cdf")),
    1
  )
  expect_lte(
    max(pdf_cdf_quant_margin(x_data, "Gamma", pars, "cdf")),
    1
  )
  expect_lte(
    max(pdf_cdf_quant_margin(x_data, "Logistic", pars, "cdf")),
    1
  )
  expect_lte(
    max(pdf_cdf_quant_margin(x_data, "Loglogistic", pars, "cdf")),
    1
  )
  expect_lte(
    max(pdf_cdf_quant_margin(x_data, "Lognormal", pars, "cdf")),
    1
  )
  expect_lte(
    max(pdf_cdf_quant_margin(x_data, "Normal", pars, "cdf")),
    1
  )
  expect_lte(
    max(pdf_cdf_quant_margin(x_data, "Weibull", pars, "cdf")),
    1
  )
  expect_lte(
    max(pdf_cdf_quant_margin(x_data, "Skew Normal", c(1,2,1), "cdf")),
    1
  )
  expect_lte(
    max(pdf_cdf_quant_margin(x_data, "Student-t", c(1,2,3), "cdf")),
    1
  )
  expect_lte(
    max(pdf_cdf_quant_margin(x_data, "Skew Student-t", c(1,2,3,4), "cdf")),
    1
  )
})

test_that("quantile functions for margins work", {
  expect_type(
    pdf_cdf_quant_margin(prob, "Cauchy", pars, "quant"),
    "double"
  )
  expect_type(
    pdf_cdf_quant_margin(prob, "Gamma", pars, "quant"),
    "double"
  )
  expect_type(
    pdf_cdf_quant_margin(prob, "Logistic", pars, "quant"),
    "double"
  )
  expect_type(
    pdf_cdf_quant_margin(prob, "Loglogistic", pars, "quant"),
    "double"
  )
  expect_type(
    pdf_cdf_quant_margin(prob, "Lognormal", pars, "quant"),
    "double"
  )
  expect_type(
    pdf_cdf_quant_margin(prob, "Normal", pars, "quant"),
    "double"
  )
  expect_type(
    pdf_cdf_quant_margin(prob, "Weibull", pars, "quant"),
    "double"
  )
  expect_type(
    pdf_cdf_quant_margin(prob, "Skew Normal", c(1,2,1), "quant"),
    "double"
  )
  expect_type(
    pdf_cdf_quant_margin(prob, "Student-t", c(1,2,3), "quant"),
    "double"
  )
  expect_type(
    pdf_cdf_quant_margin(prob, "Skew Student-t", c(1,2,3,1), "quant"),
    "double"
  )
})

test_that("marginal fitting accepts vector family subsets", {
  fit_subset <- fit_margin(x_data, min(x_data), c("norm", "gamma"))
  expect_true(fit_subset$fam %in% c("Normal", "Gamma"))
})

test_that("marginal fitting regularizes degenerate samples", {
  fit_constant <- fit_margin(rep(1, 5), 1, NA)
  expect_true(is.finite(fit_constant$par_mar[2]))
  expect_gt(fit_constant$par_mar[2], 0)
})

test_that("incompatible family sets fail clearly", {
  expect_error(
    fit_margin(c(-1, 0, 1), -1, c("gamma", "lnorm", "llogis")),
    "No admissible marginal family"
  )
})

test_that("unsupported margins fail explicitly", {
  expect_error(
    pdf_cdf_quant_margin(c(0.2, 0.5), "NotAMargin", c(1, 2), "cdf"),
    "Unsupported marginal distribution"
  )
  expect_error(
    pdf_cdf_quant_margin(1.2, "Normal", c(1, 2), "quant"),
    "Quantile inputs must lie in \\[0, 1\\]"
  )
})

test_that("loglogistic support boundary is handled safely", {
  loglogistic_pdf <- pdf_cdf_quant_margin(c(-1, 0, 1), "Loglogistic", c(2, 3), "pdf")
  loglogistic_cdf <- pdf_cdf_quant_margin(c(-1, 0, 1), "Loglogistic", c(2, 3), "cdf")
  expect_equal(loglogistic_pdf[1:2], c(0, 0))
  expect_equal(loglogistic_cdf[1:2], c(0, 0))
})

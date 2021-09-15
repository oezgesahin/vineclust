# simulate data
set.seed(111)
x_data <- rgamma(100, 1, 2)
fit <- fit_margin(x_data, 1, NA)
pars <- c(1,2)

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
    pdf_cdf_quant_margin(x_data[1], "Cauchy", pars, "quant"),
    "double"
  )
  expect_type(
    pdf_cdf_quant_margin(x_data[1], "Gamma", pars, "quant"),
    "double"
  )
  expect_type(
    pdf_cdf_quant_margin(x_data[1], "Logistic", pars, "quant"),
    "double"
  )
  expect_type(
    pdf_cdf_quant_margin(x_data[1], "Loglogistic", pars, "quant"),
    "double"
  )
  expect_type(
    pdf_cdf_quant_margin(x_data[1], "Lognormal", pars, "quant"),
    "double"
  )
  expect_type(
    pdf_cdf_quant_margin(x_data[1], "Normal", pars, "quant"),
    "double"
  )
  expect_type(
    pdf_cdf_quant_margin(x_data[1], "Skew Normal", c(1,2,1), "quant"),
    "double"
  )
  expect_type(
    pdf_cdf_quant_margin(x_data[1], "Student-t", c(1,2,3), "quant"),
    "double"
  )
  expect_type(
    pdf_cdf_quant_margin(x_data[1], "Skew Student-t", c(1,2,3,1), "quant"),
    "double"
  )
})

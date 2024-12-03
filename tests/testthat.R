library(testthat)
library(deepregression)

if (.Platform$OS.type != "windows" &&
  reticulate::py_available() &&
  reticulate::py_module_available("tensorflow") &&
  reticulate::py_module_available("keras") &&
  reticulate::py_module_available("tensorflow_probability")) {
  test_check("deeptrafo")
}

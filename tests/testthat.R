library(testthat)
library(deepregression)

if (reticulate::py_available() &
  reticulate::py_module_available("tensorflow") &
  reticulate::py_module_available("keras") &
  reticulate::py_module_available("tensorflow_probability") &
  .Platform$OS.type != "windows") {
  test_check("deeptrafo")
}

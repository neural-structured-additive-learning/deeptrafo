library(testthat)
library(deepregression)

if (reticulate::py_module_available("tensorflow") &
    reticulate::py_module_available("keras") &
    reticulate::py_module_available("tensorflow_probability")) {
  test_check("deeptrafo")
}

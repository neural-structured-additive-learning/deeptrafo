# Demo ordinal response
# Lucas Kook
# Feb 2022

set.seed(1234)

# Deps --------------------------------------------------------------------

library(tram)
devtools::load_all("../deepregression/")
devtools::load_all(".")

# Data --------------------------------------------------------------------

data("wine", package = "ordinal")

# Model -------------------------------------------------------------------

tm <- Polr(rating ~ temp + contact, data = wine)

m <- deeptrafo(rating ~ 0 + temp + contact, data = wine,
               optimizer = optimizer_adam(learning_rate = 0.1, decay = 1e-4))
fit(m, epochs = 3e2, validation_split = NULL, batch_size = nrow(wine))

coef(tm)
unlist(coef(m, which = "shifting"))

# Unconditional case ------------------------------------------------------

tm <- Polr(rating ~ 1, data = wine)

m <- deeptrafo(rating ~ 1, data = wine,
               weight_options = weight_control(
                 warmstart_weights = list(list(), list(), list("1" = 0)),
                 specific_weight_options = list(list(), list(), list("1" = list(trainable = FALSE)))),
               optimizer = optimizer_adam(learning_rate = 0.1, decay = 1e-4))
coef(m, "shifting")
fit(m, epochs = 3e2, validation_split = NULL, batch_size = nrow(wine))
coef(m, "shifting")

coef(tm, with_baseline = TRUE)
unlist(coef(m, which = "interacting"))[-5] + unlist(coef(m, which = "shifting"))

# Image data --------------------------------------------------------------

mnist <- dataset_mnist()
c(c(x_train, y_train), c(x_test, y_test)) %<-% mnist
x_train <- array_reshape(x_train, c(60000, 28, 28, 1))
x_test <- array_reshape(x_test, c(10000, 28, 28, 1))
nim <- 1e3
x_train <- x_train[1:nim,,,, drop = FALSE] / 255
x_test <- x_test[1:nim,,,, drop = FALSE] / 255
df <- list(y = ordered(y_train[1:nim]), x = x_train)

mim <- keras_model_sequential() %>%
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "relu",
                input_shape = c(28, 28, 1)) %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_conv_2d(filters = 64, kernel_size = c(3, 3), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_conv_2d(filters = 64, kernel_size = c(3, 3), activation = "relu") %>%
  layer_flatten() %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 32, activation = "relu") %>%
  layer_dense(units = 1, use_bias = FALSE)

# Complex shift model
m <- deeptrafo(y ~ mim(x), data = df, list_of_deep_models = list(mim = mim),
               optimizer = optimizer_adam(learning_rate = 1e-4))

# Complex intercept model
# m <- deeptrafo(y | mim(x) ~ 1, data = df, list_of_deep_models = list(mim = mim),
#                optimizer = optimizer_adam(learning_rate = 1e-4))

fit(m, epochs = 10L)
coef(m, which = "interacting")
coef(m, which = "shifting")

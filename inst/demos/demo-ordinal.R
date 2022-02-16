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
unlist(coef(m, which = "h2"))

# Unconditional case ------------------------------------------------------

tm <- Polr(rating ~ 1, data = wine)

m <- deeptrafo(rating ~ 1, data = wine,
               weight_options = weight_control(
                 warmstart_weights = list(list(), list(), list("1" = 0)),
                 specific_weight_options = list(list(), list(), list("1" = list(trainable = FALSE)))),
               optimizer = optimizer_adam(learning_rate = 0.1, decay = 1e-4))
coef(m, "h2")
fit(m, epochs = 3e2, validation_split = NULL, batch_size = nrow(wine))
coef(m, "h2")

coef(tm, with_baseline = TRUE)
unlist(coef(m, which = "h1"))[-5] + unlist(coef(m, which = "h2"))

# Image data --------------------------------------------------------------

mnist <- dataset_mnist()
c(c(x_train, y_train), c(x_test, y_test)) %<-% mnist
x_train <- array_reshape(x_train, c(60000, 28, 28, 1))
x_test <- array_reshape(x_test, c(10000, 28, 28, 1))
x_train <- x_train / 255
x_test <- x_test / 255
df <- list(y = ordered(y_train), x = x_train)

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

m <- deeptrafo(y ~ mim(x), data = df, list_of_deep_models = list(mim = mim),
               optimizer = optimizer_adam(learning_rate = 1e-4))
fit(m, epochs = 5L)
coef(m, which = "h1")
coef(m, which = "h2")
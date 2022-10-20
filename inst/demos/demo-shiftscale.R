# Shift-scale demo with shared weights
# LK 10-22

set.seed(321)

library(deeptrafo)

n <- 1e3
x <- runif(n, min = -4, max = 4)
y <- sin(x) + (0.5 + cos(x)) * rnorm(n)
dat <- data.frame(y = y, x = x)

plot(x, y)

nn <- keras_model_sequential() %>%
  layer_dense(16, "relu", input_shape = 1L) %>%
  layer_dense(16, "relu") %>%
  layer_dense(16, "relu") %>%
  layer_dense(16, "relu")

nn_shift <- keras_model(nn$input, nn$output %>% layer_dense(1, "linear"))
nn_scale <- keras_model(nn$input, nn$output %>% layer_dense(1, "softplus"))

m <- BoxCoxNN(y | scale(x) ~ 0 + shift(x), data = dat, order = 2,
          list_of_deep_models = list(scale = nn_scale, shift = nn_shift),
          optimizer = optimizer_adam(learning_rate = 1e-2, decay = 1e-5))

get_weights(nn) # Common trunk
get_weights(nn_shift) # Shift specific last layer
get_weights(nn_scale) # Scale specific last layer

fit(m, epochs = 1e4, batch_size = n, validation_split = NULL)

# plot(dat$y, predict(m, type = "trafo"))

nd <- data.frame(x = seq(-4, 4, length.out = 100))
preds <- predict(m, newdata = nd, q = (qy <- seq(-3, 3, length.out = 1e3)),
                 type = "pdf")

matplot(qy, t(do.call("cbind", preds)), type = "l")
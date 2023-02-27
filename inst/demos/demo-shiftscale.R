# Shift-scale demo with shared weights
# LK 10-22

set.seed(321)

library("deeptrafo")

n <- 1e3
x <- runif(n, min = -4, max = 4)
y <- sin(x) + (0.5 + cos(x)) * rnorm(n)
dat <- data.frame(y = y, x = x)

nn <- keras_model_sequential() %>%
  layer_dense(16, "relu", input_shape = 1L) %>%
  layer_dense(16, "relu") %>%
  layer_dense(16, "relu") %>%
  layer_dense(16, "relu")

nn_shift <- keras_model(nn$input, nn$output %>% layer_dense(1, "linear"))
nn_scale <- keras_model(nn$input, nn$output %>% layer_dense(1, "softplus"))

m <- LmNN(y | scale(x) ~ 0 + shift(x), data = dat,
          list_of_deep_models = list(scale = nn_scale, shift = nn_shift),
          optimizer = optimizer_adam(learning_rate = 1e-2, decay = 1e-5))

# tl <- m$model$layers[[tidx <- grep("ia_1", unlist(lapply(m$model$layers,
#                       \(x) x$name)))[2]]]
# tmp <- tl$get_weights()
# tmp[[1]][] <- c(0, 0)
# tl$set_weights(tmp)
# freeze_weights(tl)
# freeze_weights(m$model$layers[[25]])

get_weights(nn) # Common trunk
get_weights(nn_shift) # Shift specific last layer
get_weights(nn_scale) # Scale specific last layer

# coef(m, which = "int")
fit(m, epochs = 3e3, batch_size = n, validation_split = NULL)
# coef(m, which = "int")

# plot(dat$y, predict(m, type = "trafo"))

nd <- data.frame(x = seq(-4, 4, length.out = 100))
preds <- predict(m, newdata = nd, q = (qy <- seq(-3, 3, length.out = 1e3)),
                 type = "pdf")

matplot(qy, mpreds <- t(do.call("cbind", preds)), type = "l")

# plot(x, y)
# nd$mode <- qy[apply(do.call("cbind", preds), 1, which.max)]
# lines(nd$x, nd$mode, lwd = 2, col = 2)

contour(nd$x, qy, t(mpreds), col = 2, nlevels = 110)
points(x, y, pch = 20, col = rgb(.1, .1, .1, .5))
# lapply(seq_along(nd$x), \(pr) lines(nd$x[pr] + sqrt(mpreds[, pr]) / 3,
#        qy, col = rgb(.5, .1, .1, .5)))

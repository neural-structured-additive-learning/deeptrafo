# Shift-scale demo with shared weights
# LK 10-22

set.seed(321)

library("deeptrafo")

n <- 1e3
x <- runif(n, min = -4, max = 4)
y <- sin(x) + (0.5 + cos(x)) * rnorm(n)
dat <- data.frame(y = y, x = x)

m <- LmNN(y | s(x) ~ 0 + s(x), data = dat,
          optimizer = optimizer_adam(learning_rate = 1e-2, decay = 1e-5))

fit(m, epochs = 3e3, batch_size = n, validation_split = NULL)

nd <- data.frame(x = seq(-4, 4, length.out = 100))
preds <- predict(m, newdata = nd, q = (qy <- seq(-3, 3, length.out = 1e3)),
                 type = "pdf")

matplot(qy, mpreds <- t(do.call("cbind", preds)), type = "l")

contour(nd$x, qy, t(mpreds), col = 2, nlevels = 110)
points(x, y, pch = 20, col = rgb(.1, .1, .1, .5))

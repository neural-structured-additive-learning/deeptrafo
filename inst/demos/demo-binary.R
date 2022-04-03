# Demo binary response
# Lucas Kook
# Feb 2022

set.seed(1234)

# Deps --------------------------------------------------------------------

library(tram)
devtools::load_all("../deepregression/")
devtools::load_all(".")

# Data --------------------------------------------------------------------

df <- data.frame(
  x = rnorm(100),
  y = ordered(sample(0:1, 100, replace = TRUE))
)

# Model -------------------------------------------------------------------

tm <- glm(y ~ x, data = df, family = "binomial")

m <- deeptrafo(y ~ 0 + x, data = df,
               optimizer = optimizer_adam(learning_rate = 0.1, decay = 1e-4))
fit(m, epochs = 3e2, validation_split = NULL, batch_size = nrow(df))

coef(tm)
c(unlist(coef(m, which = "interacting")), unlist(coef(m, which = "shifting")))

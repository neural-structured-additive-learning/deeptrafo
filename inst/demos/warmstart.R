# Demo continuous response
# Lucas Kook
# Feb 2022

set.seed(1234)

# Deps --------------------------------------------------------------------

library(tram)
devtools::load_all(".")

# Data --------------------------------------------------------------------

data("BostonHousing2", package = "mlbench")

# Model -------------------------------------------------------------------

tm <- Colr(cmedv ~ nox + age, data = BostonHousing2, order = 6,
           support = range(BostonHousing2$cmedv))
cfb <- coef(tm, with_baseline = TRUE)[1:7]

m <- deeptrafo(cmedv ~ 0 + nox + age, data = BostonHousing2, order = 6)

tmp <- get_weights(m$model)
tmp[[2]][] <- coef(tm)[1]
tmp[[3]][] <- coef(tm)[2]
tmp[[1]][] <- c(cfb[1], log(exp(diff(cfb)) - 1 + 1e-6))
set_weights(m$model, tmp)

# Compare shift coefs
coef(m)
coef(m, "shifting")

# Compare baseline coefs
coef(tm, with_baseline = TRUE)
unlist(coef(m, which = "interacting"))

# Compare logLik
logLik(tm)
logLik(m)

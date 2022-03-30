# Demo continuous response -- nonparametric estimator
# Lucas Kook
# Feb 2022

set.seed(1234)

# Deps --------------------------------------------------------------------

library(tram)
devtools::load_all("../deepregression/")
devtools::load_all(".")

# Data --------------------------------------------------------------------

data("BostonHousing2", package = "mlbench")

# Model -------------------------------------------------------------------

tm <- Colr(cmedv ~ nox + chas, data = BostonHousing2, order = 25)

BostonHousing2$ocmedv <- BostonHousing2$cmedv
BostonHousing2$cmedv <- ordered(BostonHousing2$cmedv)

# otm <- Polr(cmedv ~ nox + chas, data = BostonHousing2) # takes long!
# lvls <- as.numeric(levels(BostonHousing2$cmedv))
# cfb <- coef(otm, with_baseline = TRUE)
# plot(lvls[-length(lvls)], cfb[grepl("cmedv", names(cfb))], type = "s")

## Unconditional
m <- deeptrafo(cmedv ~ 1, data = BostonHousing2,
               optimizer = optimizer_adam(learning_rate = 0.1, decay = 1e-4))
fit(m, epochs = 4e2, validation_split = NULL, batch_size = nrow(BostonHousing2))

plot(ecdf(BostonHousing2$ocmedv))
lines(levels(BostonHousing2$cmedv),
      c(plogis(unlist(coef(m, which = "interacting")) +
                 unlist(coef(m))), 1), type = "s", col = "red")

## Conditional
m <- deeptrafo(cmedv ~ 0 + chas + nox, data = BostonHousing2,
               optimizer = optimizer_adam(learning_rate = 0.1, decay = 1e-4))
fit(m, epochs = 8e2, validation_split = NULL, batch_size = nrow(BostonHousing2))

plot(tm, which = "baseline only")
lines(levels(BostonHousing2$cmedv)[-length(levels(BostonHousing2$cmedv))],
      unlist(coef(m, which = "interacting")), type = "s", col = "red")
legend("topleft", c("Bernstein 25", "Nonparametric"), col = c("black", "red"),
       lwd = 2)


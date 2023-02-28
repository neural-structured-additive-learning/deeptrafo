# Demo survival response
# Lucas Kook
# Feb 2022

set.seed(1234)

# Deps --------------------------------------------------------------------

library(tram)
library(survival)
devtools::load_all(".")

# Data --------------------------------------------------------------------

data("GBSG2", package = "TH.data")
GBSG2$surv <- survival::Surv(GBSG2$time, GBSG2$cens)
GBSG2$time <- as.double(GBSG2$time)

# Model -------------------------------------------------------------------

tm <- Coxph(surv ~ horTh + age, data = GBSG2, order = 6, support = range(GBSG2$time))
m <- CoxphNN(surv ~ 0 + horTh + age, data = GBSG2, order = 6)
cxph <- coxph(surv ~ horTh + age, data = GBSG2)

cfb <- coef(tm)
cfx <- coef(tm, with_baseline = TRUE)
cfx <- cfx[!names(cfx) %in% names(cfb)]
tmp <- get_weights(m$model)
tmp[[1]][] <- c(cfx[1], log(exp(diff(cfx)) - 1 + 1e-6))
tmp[[2]][] <- cfb[1]
tmp[[3]][] <- cfb[2]
set_weights(m$model, tmp)

unlist(coef(m, which = "int"))
coef(tm, with_baseline = TRUE)

coef(tm)
unlist(coef(m, which = "shifting"))

# Martingale residuals
all.equal(residuals(m), unname(residuals(tm)), tolerance = 1e-4)
plot(residuals(m), -residuals(cxph, type = "martingale"))
abline(0, 1)

# Unconditional case ------------------------------------------------------

tm <- Colr(surv ~ 1, data = GBSG2, order = 6, support = range(GBSG2$time))
cfb <- coef(tm, with_baseline = TRUE)

m <- deeptrafo(surv ~ 1, data = GBSG2, order = 6,
               optimizer = optimizer_adam(learning_rate = 0.001))
coef(m)

tmp <- get_weights(m$model)
tmp[[2]][] <- 0
tmp[[1]][] <- c(cfb[1], log(exp(diff(cfb)) - 1 + 1e-8))
set_weights(m$model, tmp)

coef(m); coef(m, "shifting")

logLik(tm)
logLik(m)

coef(tm, with_baseline = TRUE)
unlist(coef(m, which = "interacting"))

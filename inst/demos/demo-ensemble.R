# Demo ensemble
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

# Ensemble ----------------------------------------------------------------

ens <- ensemble(m, epochs = 2e2)

coef(ens)
fitted(ens, type = "cdf")

# cross-validation --------------------------------------------------------

cvs <- cv(m, epochs = 2e2, cv_folds = 3L, batch_size = nrow(wine))

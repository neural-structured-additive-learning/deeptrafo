# Demo count response
# Lucas Kook
# Feb 2022

set.seed(1234)

# Deps --------------------------------------------------------------------

library(cotram)
devtools::load_all("../deepregression/")
devtools::load_all(".")

# Data --------------------------------------------------------------------

data("birds", package = "TH.data")
birds$AOT <- c(scale(birds$AOT))
birds$AFS <- c(scale(birds$AFS))

# Model -------------------------------------------------------------------

tm <- cotram(SG5 ~ AOT + AFS, data = birds, method = "logit", log_first = FALSE)

m <- deeptrafo(SG5 ~ 0 + AOT + AFS, data = birds, order = 6,
               optimizer = optimizer_adam(learning_rate = 0.1, decay = 1e-4))
fit(m, epochs = 3e2, validation_split = NULL, batch_size = nrow(birds))

coef(tm)
unlist(coef(m, which = "shifting"))

# Unconditional case ------------------------------------------------------

tm <- Colr(SG5 ~ 1, data = birds, order = 6, supp = range(birds$SG5))
coef(tm, with_baseline = TRUE)

m <- deeptrafo(SG5 ~ 1, data = birds, order = 6,
               weight_options = weight_control(
                 warmstart_weights = list(list(), list(), list("1" = 0)),
                 specific_weight_options = list(list(), list(), list("1" = list(trainable = FALSE)))),
               optimizer = optimizer_adam(learning_rate = 0.05, decay = 1e-3))
coef(m, "shifting")
fit(m, epochs = 2e3, validation_split = NULL, batch_size = nrow(birds))
coef(m, "shifting")

coef(tm, with_baseline = TRUE)
unlist(coef(m, which = "interacting"))

logLik(tm)
logLik(m)

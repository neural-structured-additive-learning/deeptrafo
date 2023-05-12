rm(list = ls())
reticulate::use_condaenv("/Users/flipst3r/opt/anaconda3/envs/r-reticulate", required = TRUE)

library(mlbench)
data("BostonHousing2")

library(tensorflow)
tf$config$run_functions_eagerly(TRUE) # does not construct a static comp graph to be executed later
tf$data$experimental$enable_debug_mode()

#tensorboard("logs/run_a")

# Set up neural network architecture
nn <- \(x) x |>
  layer_dense(input_shape = 1L, units = 2L, activation = "relu") |>
  layer_dense(1L)

# dev version
devtools::load_all("/Users/flipst3r/RStHomeDir/GitHub/deeptrafo")

# Model formula and definition

BostonHousing2_train <- BostonHousing2[1:400, ]
BostonHousing2_test <- BostonHousing2[401:nrow(BostonHousing2), ]
grid_size <- 15L # for evaluation of the density/cdf inside loss

fml <- cmedv|lstat ~ lstat + rm + ptratio + lon + nn(b)

### CRPS

m_crps <- deeptrafo(fml, BostonHousing2_train, latent_distr = "normal", monitor_metric = NULL,
               crps = TRUE,
               grid_size = grid_size,
               addconst_interaction = 0,
               batch_size = 16*grid_size,
               optimizer = optimizer_adam(learning_rate = 0.025),# clipnorm = 1.0),
               return_data = TRUE, list_of_deep_models = list(nn = nn), 
               trafo_options = trafo_control(support = c(5, 50)))

m_crps %>% fit(epochs = 100L,
              batch_size = 16*grid_size, # grid_size x batch_size since the entire distribution for each observation is pushed through the net
              callbacks = list(callback_early_stopping(patience = 5, 
                                                       monitor = "val_loss",
                                                       restore_best_weights = T),
                               callback_reduce_lr_on_plateau(patience = 2, factor = 0.5)),
              validation_split = 0.2,
              shuffle = FALSE) # shuffle FALSE is crucial

## Predict on test

BostonHousing2_test_copy <- BostonHousing2_test
n_size <- nrow(BostonHousing2_test_copy)
BostonHousing2_test_copy <- BostonHousing2_test_copy %>% tidyr::uncount(grid_size)
BostonHousing2_test_copy$y_grid <-  rep(seq(5, 50, length.out = grid_size), n_size)
BostonHousing2_test_copy$ID <- rep(1:n_size, each = grid_size) 

n_test <- length(unique(BostonHousing2_test_copy$ID))

pp_crps <- predict(m_crps, type = "pdf", newdata = BostonHousing2_test_copy)
pp_crps <- split(pp_crps[,1], rep(1:n_test, each = grid_size))

par(mfrow = c(3,3))
for (idx in 1:9) {
  plot(BostonHousing2_test_copy$y_grid[1:grid_size], pp_crps[[idx]])
  abline(v = BostonHousing2_test$cmedv[idx])
}

logLik(m_crps, newdata = BostonHousing2_test_copy)/n_test # divide crps by n_test to get mean instead of sum

### Regular Trafo

m_trafo <- deeptrafo(fml, BostonHousing2_train, latent_distr = "normal", monitor_metric = NULL,
                crps = FALSE,
                grid_size = grid_size,
                addconst_interaction = 0,
                batch_size = 32,
                optimizer = optimizer_adam(learning_rate = 0.001),# clipnorm = 1.0),
                return_data = TRUE, list_of_deep_models = list(nn = nn), 
                trafo_options = trafo_control(support = c(5, 50)))

m_trafo %>% fit(epochs = 100L,
               batch_size = 32, # grid_size x batch_size since the entire distribution for each observation is pushed through the net
               callbacks = list(callback_early_stopping(patience = 5, 
                                                        monitor = "val_loss",
                                                        restore_best_weights = T),
                                callback_reduce_lr_on_plateau(patience = 2, factor = 0.5)),
               validation_split = 0.2,
               shuffle = FALSE) # shuffle FALSE is crucial

logLik(m_trafo, newdata = BostonHousing2_test)

pp_trafo <- predict(m_trafo, type = "pdf", newdata = BostonHousing2_test[,-6], K = grid_size)

par(mfrow = c(3,3))
pp_trafo <- do.call("cbind", pp_trafo) # rows are IDs
for (idx in 1:9) {
  plot(BostonHousing2_test_copy$y_grid[1:grid_size], pp_trafo[idx,])
  abline(v = BostonHousing2_test$cmedv[idx])
}


### TODO:

# compare run-time
# abandon eager execution and see if works on graph
# add a penalty to the loss, already done at layer level?
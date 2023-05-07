rm(list = ls())
reticulate::use_condaenv("/Users/flipst3r/opt/anaconda3/envs/r-reticulate", required = TRUE)

library(mlbench)
data("BostonHousing2")

library(tensorflow)
tf$config$run_functions_eagerly(TRUE) # does not construct a static comp graph to be executed later
tf$data$experimental$enable_debug_mode()

# Set up neural network architecture
nn <- \(x) x |>
  layer_dense(input_shape = 1L, units = 2L, activation = "relu") |>
  layer_dense(1L)

# dev version
devtools::load_all("/Users/flipst3r/RStHomeDir/GitHub/deeptrafo")

# Model formula and definition

BostonHousing2_train <- BostonHousing2[1:400, ]
BostonHousing2_test <- BostonHousing2[401:nrow(BostonHousing2), ]
grid_size <- 8L # for evaluation of the density/cdf inside loss

fml <- cmedv|lstat + rm + ptratio + lon + b ~ lstat + rm + ptratio + lon + nn(b)
mm <- deeptrafo(fml, BostonHousing2_train, latent_distr = "normal", monitor_metric = NULL,
               crps = TRUE,
               grid_size = grid_size,
               batch_size = 32*grid_size,
               optimizer = optimizer_adam(learning_rate = 0.005),#, clipnorm = 2.0),
               return_data = TRUE, list_of_deep_models = list(nn = nn), 
               trafo_options = trafo_control(support = c(5, 50)))

mm %>% fit(epochs = 20L,
          batch_size = 32*grid_size, # grid_size x batch_size since the entire distribution for each observation is pushed through the net
          callbacks = list(callback_early_stopping(patience = 5, 
                                                   monitor = "val_loss",
                                                   restore_best_weights = T),
                           callback_reduce_lr_on_plateau(patience = 2, factor = 0.5)),
          validation_split = 0.2,
          shuffle = FALSE) # shuffle FALSE is crucial

n_size <- nrow(BostonHousing2_test)
BostonHousing2_test <- BostonHousing2_test %>% tidyr::uncount(grid_size)
BostonHousing2_test$y_grid <-  rep(seq(5, 50, length.out = grid_size), n_size)
BostonHousing2_test$ID <- rep(1:n_size, each = grid_size) 

pp <- predict(mm, type = "pdf", newdata = BostonHousing2_test)

logLik(mm, newdata = BostonHousing2_test)

plot(BostonHousing2_test$y_grid[1:30], pp[1:30])

# compare run-time
# abandon eager execution and see if works on graph
# compare course/evolution of logLik/CRPS for training and validation on Boston Housing
# add a penalty to the loss, already done at layer level?



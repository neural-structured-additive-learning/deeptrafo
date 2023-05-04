rm(list = ls())
reticulate::use_condaenv("/Users/flipst3r/opt/anaconda3/envs/r-reticulate", required = TRUE)

library(mlbench)
data("BostonHousing2")

library(tensorflow)
tf$config$run_functions_eagerly(TRUE) # does not construct a static comp graph to be executed later
#tf$data$experimental$enable_debug_mode()

# Set up neural network architecture
nn <- \(x) x |>
  layer_dense(input_shape = 1L, units = 2L, activation = "relu") |>
  layer_dense(1L)

# dev version
devtools::load_all("/Users/flipst3r/RStHomeDir/GitHub/deeptrafo")

# Model formula and definition

BostonHousing2$y_grid <- BostonHousing2$cmedv

grid_size <- 30L
n_size <- nrow(BostonHousing2)
BostonHousing2 <- BostonHousing2 %>% tidyr::uncount(grid_size)
#BostonHousing2 <- data.table::rbindlist(replicate(n = grid_size, expr = BostonHousing2, simplify = FALSE))
BostonHousing2$y_grid <- rep(make_grid(BostonHousing2[["y_grid"]], n = grid_size)$y, n_size)
BostonHousing2$ID <- rep(1:n_size, each = grid_size)

fml <- y_grid ~ age + crim + lon + lat + ptratio + nn(nox)
m <- deeptrafo(fml, BostonHousing2, latent_distr = "normal", monitor_metric = NULL,
               optimizer = optimizer_adam(learning_rate = 0.1),
               return_data = TRUE, list_of_deep_models = list(nn = nn), crps = TRUE,
               trafo_options = trafo_control(support = c(5, 50)))

m %>% fit(epochs = 200L,
          batch_size = 64*grid_size,
          callbacks = list(callback_early_stopping(patience = 5, 
                                                   monitor = "val_loss",
                                                   restore_best_weights = T),
                           callback_reduce_lr_on_plateau(patience = 2, factor = 0.5)),
          validation_split = 0.2,
          shuffle = FALSE)

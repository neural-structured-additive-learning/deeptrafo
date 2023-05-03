rm(list = ls())
reticulate::use_condaenv("/Users/flipst3r/opt/anaconda3/envs/r-reticulate", required = TRUE)
library(tensorflow)
library(keras)
library(reticulate)
tf$config$run_functions_eagerly(TRUE)
tf$data$experimental$enable_debug_mode()
library(mlbench)
source("crps_helper.R")
devtools::load_all("/Users/flipst3r/RStHomeDir/GitHub/deeptrafo")

loss_tracker <- metric_mean(name = "loss")

crps_model <- new_model_class(
  classname = "test_model",
  train_step = function(data) {

    c(x, y) %<-% data
    
    with(tf$GradientTape() %as% tape, {
      
      ### CRPS
        
      # grid granularity for F^-1_Y|X = x evaluation
      y_grid <- tf$linspace(supp_y[1], supp_y[2], grid_gran)
      
      # F_Z
      bd <- tfd_normal(0,1)
      
      # f_Y|X = x
      res <- eval_density(self, x, y_grid, bd, train = TRUE)
      y_pred <- res$y_pred
      
      # F_Y|X = x
      cumulative_df <- eval_cdf(res$densities, y_grid)
      
      crp_scores <- calc_crps(y_true, cumulative_df, y_grid)
      
      loss <- tf$reduce_mean(crp_scores)
        
      # ### Log-Scores
      # 
      # y_pred <- self(x, training = TRUE)
      # bd <- tfd_normal(0,1)
      # tloss <- nll(bd)
      # loss <- tloss(y, y_pred)
    })

    # Compute gradients
    trainable_vars <- self$trainable_variables
    
    gradients <- tape$gradient(loss, trainable_vars)
    
    # Update weights
    self$optimizer$apply_gradients(zip_lists(gradients, trainable_vars))
    # Update metrics (includes the metric that tracks the loss)
    self$compiled_metrics$update_state(y, y_pred)
    
    # Return a named list mapping metric names to current value
    results <- list()
    for (m in self$metrics)
      results[[m$name]] <- m$result()
    results
    
    # Compute our own metrics
    loss_tracker$update_state(loss)
    list(loss = loss_tracker$result())
    
  }, test_step = function(data) {
    
    c(x, y) %<-% data
    
    # grid granularity for F^-1_Y|X = x evaluation
    y_grid <- tf$linspace(supp_y[1], supp_y[2], grid_gran)
    
    # F_Z
    bd <- tfd_normal(0,1)
    
    # f_Y|X = x
    res <- eval_density(self, x, y_grid, bd, train = FALSE)
    y_pred <- res$y_pred
    
    # F_Y|X = x
    cumulative_df <- eval_cdf(res$densities, y_grid)
    
    crp_scores <- calc_crps(y_true, cumulative_df, y_grid)
    
    loss <- tf$reduce_mean(crp_scores)
    
    # Compute our own metrics
    loss_tracker$update_state(loss)
    list(loss = loss_tracker$result())
    
  },
  metrics = mark_active(function() {
    list(loss_tracker)
  })
)

data("BostonHousing2")
BostonHousing2 <- BostonHousing2[1:128,]# faster test runs
range(BostonHousing2$cmedv)

nn <- \(x) x |>
  layer_dense(input_shape = 1L, units = 2L, activation = "relu") |>
  layer_dense(1L)

# needed in scope of "my_class"
y_true <- BostonHousing2$cmedv 
grid_gran <- 10L
supp_y <- c(0, 60)
M_crps <- 15L

# not needed in scope
fml <- cmedv|lon ~ 0 + age + nn(crim) #+ s(dis, df = 3)

m <- deeptrafo(fml, BostonHousing2, latent_distr = "normal", monitor_metric = NULL,
               return_data = TRUE, list_of_deep_models = list(nn = nn),
               trafo_options = trafo_control(order_bsp = 7, 
                                             support = supp_y),
               optimizer = optimizer_adam(learning_rate = 0.15),
               model_fun = crps_model
)

m %>% fit(epochs = 1L,
          validation_split = 0.2,
          shuffle = FALSE,
          callbacks = list(callback_early_stopping(patience = 5, 
                                                   monitor = "val_loss",
                                                   restore_best_weights = T),
                           callback_reduce_lr_on_plateau(patience = 2, factor = 0.5)),
          batch_size = nrow(BostonHousing2))


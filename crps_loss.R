rm(list = ls())
library(mlbench)
library(data.table)
reticulate::use_condaenv("/Users/flipst3r/opt/anaconda3/envs/r-reticulate", required = TRUE)

library(tensorflow)
tf$config$run_functions_eagerly(TRUE) # does not construct a static comp graph to be executed later
#tf$data$experimental$enable_debug_mode()

# Prepare the data
# data("wine", package = "ordinal")
# wine$z <- rnorm(nrow(wine))
# wine$x <- rnorm(nrow(wine))

data("BostonHousing2")
d <- setDT(BostonHousing2)

# Set up neural network architecture
nn <- \(x) x |>
  layer_dense(input_shape = 1L, units = 2L, activation = "relu") |>
  layer_dense(1L)

devtools::load_all("/Users/flipst3r/RStHomeDir/GitHub/deeptrafo")

# Model formula and definition

BostonHousing2$y_grid <- BostonHousing2$cmedv

grid_size <- 50L
n_size <- nrow(BostonHousing2)
BostonHousing2 <- BostonHousing2 %>% tidyr::uncount(grid_size)
#BostonHousing2 <- data.table::rbindlist(replicate(n = grid_size, expr = BostonHousing2, simplify = FALSE))
BostonHousing2$y_grid <- rep(make_grid(BostonHousing2[["y_grid"]], n = grid_size)$y, n_size)
BostonHousing2$ID <- rep(1:n_size, each = grid_size)

# BostonHousing2 <- list(y_grid = list(rnorm(4), rnorm(4)), lon = list(rnorm(4), rnorm(4)),
#                        age = list(rnorm(4), rnorm(4)), crim = list(rnorm(4), rnorm(4)))

fml <- y_grid|lon ~ 0 + age + nn(crim)
m <- deeptrafo(fml, BostonHousing2, latent_distr = "normal", monitor_metric = NULL,
               return_data = TRUE, list_of_deep_models = list(nn = nn), crps = TRUE,
               trafo_options = trafo_control(support = c(0, 50)))

# n <- nrow(BostonHousing2)
# not_found <- TRUE
# i <- 101
# while(not_found) {
#   if (i %% 50L == 0) {
#     candidate = i / 50L
#     if (n %% candidate == 0) {
#       print(i)
#       not_found <- FALSE
#     }
#   }
#   i <- i + 1
# }
# # Overview
# print(m)


# Fit
m %>% fit(epochs = 1,
          batch_size = 32*grid_size, 
          shuffle = FALSE,
          validation_data = NULL,
          validation_split = list())

# Fit model using custom generator
batch_size <- 10
steps_per_epoch <- ceiling(length(y) / batch_size)
model %>% fit_generator(
  generator = custom_generator(x, y, id, batch_size),
  steps_per_epoch = steps_per_epoch,
  epochs = 10
)

fit_generator()

list(a =  , b =)

data.frame()

# Prepare the data
data("wine", package = "ordinal")
wine$z <- rnorm(nrow(wine))
wine$x <- rnorm(nrow(wine))

# Set up neural network architecture
nn <- \(x) x |>
  layer_dense(input_shape = 1L, units = 2L, activation = "relu") |>
  layer_dense(1L)

options(tensorflow.extract.disallow_out_of_bounds = FALSE)
# Model formula and definition
fml <- rating ~ 0 + temp + contact + s(z, df = 3) + nn(x)
m <- deeptrafo(fml, wine, latent_distr = "logistic", monitor_metric = NULL,
               return_data = TRUE, list_of_deep_models = list(nn = nn))

# Overview
print(m)

# Fit
m %>% fit(epochs = 10, batch_size = nrow(wine))

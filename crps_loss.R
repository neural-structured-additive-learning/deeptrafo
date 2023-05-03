rm(list = ls())
library(mlbench)
library(data.table)

reticulate::use_condaenv("/Users/flipst3r/opt/anaconda3/envs/r-reticulate", required = TRUE)


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
fml <- y_grid|lon ~ 0 + age + nn(crim)
m <- deeptrafo(fml, BostonHousing2, latent_distr = "normal", monitor_metric = NULL,
               return_data = TRUE, list_of_deep_models = list(nn = nn), crps = TRUE,
               trafo_options = trafo_control(support = c(0, 50)))

# Overview
print(m)

# Fit
m %>% fit(epochs = 1, batch_size =32)

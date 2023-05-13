rm(list = ls())
reticulate::use_condaenv("/Users/flipst3r/opt/anaconda3/envs/r-reticulate", required = TRUE)

library(data.table)
library(ggridges)
library(ggplot2)
library(ggrepel)
library(ggjoy)
library(tsdl) # available from GitHub (FinYang/tsdl)
devtools::load_all("/Users/flipst3r/RStHomeDir/GitHub/deeptrafo")

tf$config$run_functions_eagerly(TRUE) # does not construct a static comp graph to be executed later
tf$data$experimental$enable_debug_mode()

## -----------------------------------------------------------------------------

d <- subset(tsdl, "Meteorology")
nm <- "Mean maximum temperature in Melbourne: degrees C. Jan 71 – Dec 90."
temp_idx <- sapply(d, attr, "description") == nm
d_ts <- d[temp_idx][[1]]

## -----------------------------------------------------------------------------

n_gr <- 20
min_supp <- 10
max_supp <- 30
M <- 6L
ep <- 1
p <- 3
gr <- seq(min_supp ,max_supp , length.out = n_gr)
time_var <- seq(as.Date("1971-01-01"), as.Date("1990-12-01"), by = "month")
d_ts <- data.table(time = time_var, y = as.numeric(d_ts))
d_ts[, month := factor(month(time))]

d_crps <- deeptrafo:::prep_atpm_crps(d_ts, min_supp, max_supp, n_gr, c(1:p))
d_atpm <- deeptrafo:::prep_atpm_logLik(d_ts, min_supp, max_supp, c(1:p))

time_train <- seq(as.Date("1971-01-01"), as.Date("1987-12-01"), by = "month")
time_validation <- seq(as.Date("1988-01-01"), as.Date("1989-12-01"), by = "month")
time_test <- seq(as.Date("1990-01-01"), as.Date("1990-12-01"), by = "month")

## ----------------------------FORMULA------------------------------------------

lags <- c(paste0("y_lag_", 1:p, collapse = "+"))
atplags <- paste0("atplag(", paste0("y_lag_", 1:p), ")", collapse = "+")
#(fm_atm <- as.formula(paste0("y |", lags, "~ 0 + month +", atplags)))
(fm_atp <- as.formula(paste0("y ~ 0 + month +", atplags)))
(fm_colr <- as.formula(paste0("y ~ 0 + month + ", lags)))

## ----------------------------FITTING-----------------------------------------

# d_atpm_train <- d_atpm[d_atpm$time %in% time_train,]
# d_atpm_val <- as.data.frame(d_atpm[d_atpm$time %in% time_validation,])
# d_atpm_test <- as.data.frame(d_atpm[d_atpm$time %in% time_test,])
# 
# learn_trafo <- \(fm) ColrNN(fm, data = d_atpm, trafo_options = trafo_control(
#   order_bsp = M, support = c(min_supp, max_supp)), tf_seed = 1,
#   optimizer = optimizer_adam(learning_rate = 0.01))
# 
# mods_loglik <- lapply(list(fm_atp, fm_colr), learn_trafo)
# 
# fit_fun <- \(m) m |> fit(epochs = ep,
#                          callbacks = list(
#                            callback_early_stopping(patience = 5, monitor = "val_loss"),
#                            callback_reduce_lr_on_plateau(patience = 3, factor = 0.5, monitor = "val_loss")),
#                          batch_size = 16,
#                          validation_data = list(d_atpm_val, d_atpm_val$y),
#                          )
# 
# lapply(mods_loglik, \(m) {
#   mhist <- fit_fun(m)
#   # plot(mhist)
# })
# 
# logLik(mods_loglik[[1]], newdata = d_atpm_test, convert_fun = function(x) mean(x, na.rm = T))
# logLik(mods_loglik[[2]], newdata = d_atpm_test, convert_fun = function(x) mean(x, na.rm = T))

## ----fitting-CRPS-------------------------------------------------------------

d_crps_train <- as.data.frame(d_crps[d_crps$time %in% time_train,])
d_crps_val <- as.data.frame(d_crps[d_crps$time %in% time_validation,])
d_crps_test <- as.data.frame(d_crps[d_crps$time %in% time_test,])

learn_crps <- \(fm) deeptrafo(fm, data = d_crps_train, trafo_options = trafo_control(
  order_bsp = M, support = c(min_supp, max_supp)), 
  tf_seed = 1, crps = TRUE, grid_size = n_gr,
  #batch_size = 16*grid_size, # only needed in graph execution but not in eager
  addconst_interaction = 0,
  optimizer = optimizer_adam(learning_rate = 0.01))

mods_cprs <- lapply(list(fm_atp, fm_colr), learn_crps)

fit_fun <- \(m) m |> fit(epochs = ep,
                         callbacks = list(
                           callback_early_stopping(patience = 5, monitor = "val_loss"),
                           callback_reduce_lr_on_plateau(patience = 3, factor = 0.5, monitor = "val_loss")),
                         batch_size = nrow(d_crps_train), 
                         validation_data = list(d_crps_val, d_crps_val$y_grid),
                         shuffle = FALSE) # shuffle FALSE is crucial

lapply(mods_cprs, \(m) {
  mhist <- fit_fun(m)
  # plot(mhist)
})

logLik(mods_cprs[[1]], newdata = d_crps_test, criteria = "crps")

d_crps_logLik <- d_crps_test[, .SD[1], by = ID]
d_crps_logLik$y_grid <- d_crps_logLik$y
logLik(mods_cprs[[1]], newdata = d_crps_logLik, criteria = "logLik")


logLik(mods_cprs[[2]], newdata = d_crps_test, criteria = "crps")

d_crps_logLik <- d_crps_test[, .SD[1], by = ID]
d_crps_logLik$y_grid <- d_crps_logLik$y
logLik(mods_cprs[[2]], newdata = d_crps_logLik, criteria = "logLik")

## ----in_sample_cond_densities_ATMs--------------------------------------------
# In-sample densities
m_atp <- mods[[1]]
#m_atp <- mods[[2]]
#m_colr <- mods[[3]]

nd <- as.data.table(deeptrafo:::create_lags("y", d_ts, "atplag(1:p)")$data)
nd <- nd |> dplyr::mutate(y_true = y, y = list(gr)) |> tidyr::unnest(y)
nd$d_atp <- c(predict(m_atp, newdata = nd, type = "pdf"))

nd$y_grid <- nd$y
nd$d_crps <- c(predict(m_crps, newdata = nd, type = "pdf"))

d_density <- nd |>
  tidyr::gather("method", "y_density", d_atp, d_crps) |>
  dplyr::mutate(y_grid = y) |> as.data.table()

t_idx <- seq(as.Date("1977-06-01"), as.Date("1978-05-01"), by = "month")
d_density <- d_density[d_density$time %in% t_idx,]

Sys.setlocale("LC_ALL", "en_GB.UTF-8")
(g_dens <- ggplot() +
  geom_path(data = d_density, aes(x = y_true, y = time, group = method),
            colour="red", size=1.5, alpha = 0.2) +
  geom_point(data = d_density, aes(x = y_true, y = time, group = method),
             colour="red", size=1, shape=4) +
  geom_joy(data = d_density,
           aes(height = y_density, x = y_grid, y = time,
               group = time, fill = factor(time)),
           stat="identity", alpha = 0.7, colour = rgb(0,0,0,0.5)) +
  scale_y_date(date_breaks = "4 months", date_labels = "%b %Y") +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(nrow = 2)) +
  facet_grid(~ method, labeller = as_labeller(c("d_atp" = paste0("AT(", p, ")"),
                                                "d_crps" = "CRPS"))) +
  theme_bw() +
  labs(color = "month") +
  xlab("") +
  ylab("") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank(),
        text = element_text(size=12),
        legend.position = "none",
        rect = element_rect(fill = "transparent")))

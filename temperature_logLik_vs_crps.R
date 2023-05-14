rm(list = ls())
reticulate::use_condaenv("/Users/flipst3r/opt/anaconda3/envs/r-reticulate", required = TRUE)
devtools::load_all("/Users/flipst3r/RStHomeDir/GitHub/deeptrafo")

library(data.table)
library(ggridges)
library(ggplot2)
library(ggrepel)
library(ggjoy)
library(tsdl) # available from GitHub (FinYang/tsdl)

tf$config$run_functions_eagerly(TRUE) # does not construct a static comp graph to be executed later
tf$data$experimental$enable_debug_mode()

## -----------------------------------------------------------------------------

d <- subset(tsdl, "Meteorology")
nm <- "Mean maximum temperature in Melbourne: degrees C. Jan 71 – Dec 90."
temp_idx <- sapply(d, attr, "description") == nm
d_ts <- d[temp_idx][[1]]
time_var <- seq(as.Date("1971-01-01"), as.Date("1990-12-01"), by = "month")
d_ts <- data.table(time = time_var, y = as.numeric(d_ts))
d_ts[, month := factor(month(time))]

## -----------------------------------------------------------------------------

n_gr <- 50
min_supp <- 10
max_supp <- 30
M <- 6L
ep <- 2L
bs <- 16L
lr <- 1e-1
p <- 3 # lags

### different data prep for logLik vs CRPS

# whole distribution on grid evaluation
d_crps <- prep_atpm_crps(d_ts, min_supp, max_supp, n_gr, c(1:p))

# point evaluation
d_loglik <- prep_atpm_logLik(d_ts, c(1:p))

time_train <- seq(as.Date("1971-01-01"), as.Date("1986-12-01"), by = "month")
time_validation <- seq(as.Date("1987-01-01"), as.Date("1989-12-01"), by = "month")
time_test <- seq(as.Date("1990-01-01"), as.Date("1990-12-01"), by = "month")

## ----------------------------FORMULA------------------------------------------

lags <- c(paste0("y_lag_", 1:p, collapse = "+"))
atplags <- paste0("atplag(", paste0("y_lag_", 1:p), ")", collapse = "+")
fm_atp <- as.formula(paste0("y ~ 0 + month +", atplags))

## ----------------------------LOGLIK------------------------------------------

d_loglik_train <- as.data.frame(d_loglik[d_loglik$time %in% time_train,])
d_loglik_val <- as.data.frame(d_loglik[d_loglik$time %in% time_validation,])
d_loglik_test <- as.data.frame(d_loglik[d_loglik$time %in% time_test,])

## Validation

m_loglik <- ColrNN(fm_atp, data = d_loglik_train,
                   trafo_options = trafo_control(order_bsp = M, support = c(min_supp, max_supp)),
                   tf_seed = 1,
                   optimizer = optimizer_adam(learning_rate = lr))

hist_loglik <- m_loglik |> fit(epochs = ep,
                               callbacks = list(
                                 callback_early_stopping(patience = 5, monitor = "val_loss"),
                                 callback_reduce_lr_on_plateau(patience = 3, factor = 0.5, monitor = "val_loss")),
                               batch_size = bs,
                               validation_data = list(d_loglik_val, d_loglik_val$y))

min_epochs <- which.min(hist_loglik$metrics$val_loss)

## Test

m_loglik <- ColrNN(fm_atp, data = rbind(d_loglik_train, d_loglik_val),
                   trafo_options = trafo_control(order_bsp = M, support = c(min_supp, max_supp)),
                   tf_seed = 1, 
                   optimizer = optimizer_adam(learning_rate = lr))

hist_loglik <- m_loglik |> fit(epochs = ep,
                               callbacks = list(),
                               batch_size = bs,
                               validation_data = NULL,
                               validation_split = NULL)

# pls
pls_loglik <- logLik(m_loglik, newdata = d_loglik_test, convert_fun = \(x) -x) # larger is better

## ----------------------------CRPS---------------------------------------------

d_crps_train <- as.data.frame(d_crps[d_crps$time %in% time_train,])
d_crps_val <- as.data.frame(d_crps[d_crps$time %in% time_validation,])
d_crps_test <- d_crps[d_crps$time %in% time_test,]

## Validation

m_crps <- ColrNN(fm_atp, data = d_crps_train,
                 trafo_options = trafo_control(order_bsp = M, support = c(min_supp, max_supp)),
                 tf_seed = 1,
                 crps = TRUE,
                 grid_size = n_gr,
                 #batch_size = bs*n_gr, # only needed in graph execution but not in eager
                 optimizer = optimizer_adam(learning_rate = lr))

hist_crps <- m_crps |> fit(epochs = ep,
                           callbacks = list(
                             callback_early_stopping(patience = 5, monitor = "val_loss"),
                             callback_reduce_lr_on_plateau(patience = 3, factor = 0.5, monitor = "val_loss")),
                           batch_size = bs*n_gr,
                           validation_data = list(d_crps_val, cbind(d_crps_val$y, d_crps_val$ID, d_crps_val$y_grid)),
                           shuffle = FALSE) # shuffle FALSE is crucial

min_epochs <- which.min(hist_loglik$metrics$val_loss)

## Test

m_crps <- ColrNN(fm_atp, data = rbind(d_crps_train, d_crps_val),
                 trafo_options = trafo_control(order_bsp = M, support = c(min_supp, max_supp)),
                 tf_seed = 1,
                 crps = TRUE,
                 grid_size = n_gr,
                 #batch_size = bs*n_gr, # only needed in graph execution but not in eager
                 optimizer = optimizer_adam(learning_rate = lr))

hist_crps <- m_crps |> fit(epochs = ep,
                           callbacks = list(),
                           batch_size = bs*n_gr,
                           validation_data = NULL,
                           validation_split = NULL,
                           shuffle = FALSE) # shuffle FALSE is crucial

logLik(m_crps, newdata = d_crps_test, criteria = "crps")

d_crps_logLik <- d_crps_test[, .SD[1], by = ID]
d_crps_logLik$y_grid <- d_crps_logLik$y
pls_crps <- logLik(m_crps, newdata = d_crps_logLik, criteria = "logLik") # larger is better

saveRDS(c("crps" = pls_crps, "logLik" = pls_loglik), file = "pls_temp.RDS")

## -------------------- out-of-sample densities --------------------------------

# use d_crps_test since it has grid provided
y_true <- d_crps_test$y
d_crps_test$y <- d_crps_test$y_grid
d_crps_test$logLik_pdf <- c(predict(m_loglik, newdata = d_crps_test, type = "pdf"))

d_crps_test$y <- y_true
d_crps_test$crps_pdf <- c(predict(m_crps, newdata = d_crps_test, type = "pdf"))

d_density <- d_crps_test |>
  tidyr::gather("method", "y_density", logLik_pdf, crps_pdf) |> as.data.table()

Sys.setlocale("LC_ALL", "en_GB.UTF-8")
(g_dens <- ggplot() +
  geom_path(data = d_density, aes(x = y, y = time, group = method),
            colour="red", size=1.5, alpha = 0.2) +
  geom_point(data = d_density, aes(x = y, y = time, group = method),
             colour="red", size=1, shape=4) +
  geom_joy(data = d_density,
           aes(height = y_density, x = y_grid, y = time,
               group = time, fill = factor(time)),
           stat="identity", alpha = 0.7, colour = rgb(0,0,0,0.5)) +
  scale_y_date(date_breaks = "4 months", date_labels = "%b %Y") +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(nrow = 2)) +
  facet_grid(~ method, labeller = as_labeller(c("logLik_pdf" = "LogLik",
                                                "crps_pdf" = "CRPS"))) +
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

ggsave("temp_crps_vs_loglik.pdf")

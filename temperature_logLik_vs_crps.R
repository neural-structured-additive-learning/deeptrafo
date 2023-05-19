rm(list = ls())

if (version$os != "darwin20") {
  
  if (grepl("scstepha", getwd())) {
    
    reticulate::use_virtualenv("/cluster/home/scstepha/nsgaSA", required = TRUE)
  } else {
    
    reticulate::use_virtualenv("/cluster/home/phibauma/nsgaSA", required = TRUE)
  }
  
} else {
  reticulate::use_condaenv("/Users/flipst3r/opt/anaconda3/envs/r-reticulate", required = TRUE)
}

devtools::load_all("/Users/flipst3r/RStHomeDir/GitHub/deeptrafo")

#library(deeptrafo)
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

grid_size <- n_gr <- 30
min_supp <- -30
max_supp <- 90
M <- 5L
ep <- 10L
bs <- 64L
lr <- 3e-1
p <- 3 # lags

### different data prep for logLik vs CRPS

# whole distribution on grid evaluation
d_crps <- deeptrafo:::prep_atpm_crps(d_ts, min_supp, max_supp, n_gr, c(1:p))

# point evaluation
d_loglik <- deeptrafo:::prep_atpm_logLik(d_ts, c(1:p))

time_train <- seq(as.Date("1971-01-01"), as.Date("1986-12-01"), by = "month")
time_validation <- seq(as.Date("1987-01-01"), as.Date("1989-12-01"), by = "month")
time_test <- seq(as.Date("1990-01-01"), as.Date("1990-12-01"), by = "month")

## ----------------------------FORMULA------------------------------------------

lags <- c(paste0("y_lag_", 1:p, collapse = "+"))
atplags <- paste0("atplag(", paste0("y_lag_", 1:p), ")", collapse = "+")
fm_atp <- as.formula(paste0("y ~ 0 + month +", atplags))

d_loglik_train <- as.data.frame(d_loglik[d_loglik$time %in% time_train,])
d_loglik_val <- as.data.frame(d_loglik[d_loglik$time %in% time_validation,])
d_loglik_test <- as.data.frame(d_loglik[d_loglik$time %in% time_test,])

d_crps_train <- as.data.frame(d_crps[d_crps$time %in% time_train,])
d_crps_val <- as.data.frame(d_crps[d_crps$time %in% time_validation,])
d_crps_test <- d_crps[d_crps$time %in% time_test,]

## ----------------------------LOGLIK------------------------------------------

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

# m_loglik <- ColrNN(fm_atp, data = rbind(d_loglik_train, d_loglik_val),
#                    trafo_options = trafo_control(order_bsp = M, support = c(min_supp, max_supp)),
#                    tf_seed = 1, 
#                    optimizer = optimizer_adam(learning_rate = lr))
# 
# hist_loglik <- m_loglik |> fit(epochs = min_epochs,
#                                callbacks = list(),
#                                batch_size = bs,
#                                validation_data = NULL,
#                                validation_split = NULL)

# pls
pls_loglikmodel <- logLik(m_loglik, 
                         newdata = d_loglik_test,
                         convert_fun = \(x) c(-x), # larger is better
                         criteria = "logLik")

crps_loglikmodel <- logLik(m_loglik, 
                           newdata = d_crps_test,
                           convert_fun = \(x) x, # smaller is better
                           criteria = "crps")

loglikmodel <- list("pls" = pls_loglikmodel, "crps" = crps_loglikmodel)

## ----------------------------CRPS---------------------------------------------

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

min_epochs <- which.min(hist_crps$metrics$val_loss)

## Test

# m_crps <- ColrNN(fm_atp, data = rbind(d_crps_train, d_crps_val),
#                  trafo_options = trafo_control(order_bsp = M, support = c(min_supp, max_supp)),
#                  tf_seed = 1,
#                  crps = TRUE,
#                  grid_size = n_gr,
#                  #batch_size = bs*n_gr, # only needed in graph execution but not in eager
#                  optimizer = optimizer_adam(learning_rate = lr))
# 
# hist_crps <- m_crps |> fit(epochs = min_epochs,
#                            callbacks = list(),
#                            batch_size = bs*n_gr,
#                            validation_data = NULL,
#                            validation_split = NULL,
#                            shuffle = FALSE) # shuffle FALSE is crucial

crps_crpsmodel <- logLik(m_crps, newdata = d_crps_test, criteria = "crps")

d_crps_logLik <- d_crps_test[, .SD[1], by = ID]
d_crps_logLik$y_grid <- d_crps_logLik$y
pls_crpsmodel <- logLik(m_crps, newdata = d_crps_logLik, criteria = "logLik") # larger is better
crpsmodel <- list("pls" = pls_crpsmodel, "crps" = crps_crpsmodel)

res <- list("crpsmodel" = crpsmodel, "loglikmodel" = loglikmodel)
saveRDS(res, file = "scores.RDS")

## -------------------- out-of-sample densities --------------------------------

# use d_crps_test since it has grid provided
y_true <- d_crps_test$y
d_crps_test$y <- d_crps_test$y_grid
d_crps_test$logLik_pdf <- c(predict(m_loglik, newdata = d_crps_test, type = "pdf"))

d_crps_test$y <- y_true
d_crps_test$crps_pdf <- c(predict(m_crps, newdata = d_crps_test, type = "pdf"))

d_density <- d_crps_test |>
  tidyr::gather("method", "y_density", logLik_pdf, crps_pdf) |> as.data.table()

check_dens <- function(idd) {
  d_crps <- d_density[ID == idd & method == "crps_pdf"]
  d_ll <- d_density[ID == idd & method == "logLik_pdf"]
  c1 <- cumtrapz(d_crps$y_grid, d_crps$y_density)
  c2 <- cumtrapz(d_ll$y_grid, d_ll$y_density)
  list("crps" = c1[n_gr], "logLik" = c2[n_gr])
}
sapply(unique(d_density$ID), \(idd) check_dens(idd))

Sys.setlocale("LC_ALL", "en_GB.UTF-8")
g_dens <- ggplot() +
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
        rect = element_rect(fill = "transparent"))

ggsave("temp_crps_vs_loglik.pdf", plot = g_dens)

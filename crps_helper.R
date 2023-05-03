comp_cdf <- function(x, y) {
  
  # 1-d integration of a density with trapez rule
  
  # x: quantile vector on grid
  # y: densities values belonging to x
  
  number_points <- as.integer(y$shape[1]) - 1L
  integrals <- tf$reshape(tf$constant(list()), c(0L))
  
  for (k in 1L:number_points) {
    integ <- tfp$math$trapz(tf$gather(y, 0L:k), tf$gather(x, 0L:k))
    integrals <- tf$concat(list(integrals, tf$reshape(integ, 1L)), axis = 0L)
  }
  return(integrals)
}

lin_interpol <- function(x_eval, x, y) {
  
  # linear univariate interpolation of a quantile function
  
  # x_eval: scalar prob to be evaluated
  # x: tensor of CDF values
  # y: tensor of quantiles belonging to x
  
  p_val <- tf$reshape(x_eval, 1L)
  
  tf$assert_equal(x$shape, y$shape)
  tf$assert_less(x, 1 + 1e-3)
  tf$assert_greater(x, 0 - 1e-3)
  
  if (tf$reduce_all(tf$equal(x,tf$unique(x)$y))$numpy()) {
    x <- tf$add(x, tf$linspace(1e-5, 1e-4, x$shape[1]))
  }
  
  F_y_hat_diff <- tf$experimental$numpy$diff(x)
  y_grid_diff <- tf$experimental$numpy$diff(y)
  
  tf$assert_greater(F_y_hat_diff, 0 - 1e-16) # monotone (but not strictly) increasing
  tf$assert_greater(y_grid_diff, 0 - 1e-16)
  
  relevant_interval <- tf$searchsorted(x, p_val, side = 'right')
  
  if(relevant_interval$numpy() == 0) {
    relevant_interval <- tf$reshape(1L, 1L)
  }
  
  if(relevant_interval$numpy() == x$shape[1]) {
    relevant_interval <- tf$subtract(x$shape[1], 1L)$numpy()
  }
  
  slopes <- tf$math$divide(y_grid_diff, F_y_hat_diff)
  q <- tf$add(tf$gather(y, relevant_interval - 1L),
              tf$multiply(tf$gather(slopes, relevant_interval - 1L),
                          p_val - tf$gather(x, relevant_interval - 1L)))
  tf$reshape(q, 1L)
}

eval_density <- function(self, data_x, y_grid, F_Z, train = T) {
  
  # evaluates the density for a given model (self)
  
  # batch size
  n <- as.integer(data_x[[1]]$shape)[1]
  
  grid_gran <- as.integer(y_grid$shape[1])
  collect_f_y <- tf$reshape(tf$constant(list()), c(n, 0L))
  
  # iterate across the entire grid
  for (j in 1:grid_gran) {
    
    # h_hat(v_j, X_1...X_n)
    repeated_grid <- rep(y_grid[j]$numpy(), n)
    eval_data <- data.frame(repeated_grid)
    colnames(eval_data) <- m$init_params$response_varname
    
    y_basis <- m$init_params$parsed_formulas_contents$yterms
    
    data_x[[1]] <- tf$constant(y_basis[[1]]$predict_trafo(eval_data))
    data_x[[2]] <- tf$constant(y_basis[[2]]$predict_trafo(eval_data))
    data_x[[3]] <- tf$constant(y_basis[[3]]$predict_trafo(eval_data))
    
    y_pred <- self(data_x, training = train) # forward pass with grid y
    
    # h_hat = h_1_hat + h_2_hat
    h_hat <- layer_add(list(tf_stride_cols(y_pred, 1L),
                            tf_stride_cols(y_pred, 2L)))
    
    # h'
    h_prime <- tf$math$log(tf$clip_by_value(tf_stride_cols(y_pred, 4L),
                                            1e-8, Inf))
    
    # rescaling needed?
    h_hat <- tf$divide(h_hat, tf$norm(h_hat))
    
    # maybe add penalty to the loss like in the likelihood case
    
    # f_Y|X = x
    f_y_dens <- tf$exp(tfd_log_prob(F_Z, h_hat) + h_prime)
    collect_f_y <- tf$concat(list(collect_f_y, f_y_dens), axis = 1L)
    
  }
  
  return(list("densities" =  collect_f_y, "y_pred" = y_pred))
}

eval_cdf <- function(density_collection, y_grid) {
  
  # returns a collection of CDFs for the corresponding collection of densities
  
  grid_size <- as.integer(y_grid$shape[1])
  F_y <- tf$reshape(tf$constant(list()), c(grid_size, 0L))
  b_size <- as.integer(density_collection$shape[1])
  
  for (i in 1L:b_size) {
    
    f_y <- tf$gather(density_collection, indices = c(i - 1L))
    scale_dens <- tfp$math$trapz(f_y, y_grid)
    f_y <- tf$divide(f_y, scale_dens)
    integrated_f_y <- tf$reshape(tf$concat(list(tf$reshape(0, 1L), 
                                                comp_cdf(y_grid, f_y)), axis = 0L), c(grid_size, 1L))
    F_y <- tf$concat(list(F_y, integrated_f_y), axis = 1L)
    
  }
  
  return(tf$transpose(F_y))
}

calc_crps <- function(y_observed, cdf_collection, y_grid) {
  
  # compute the crps by means of the weighted quantile loss
  
  p_grid <- tf$linspace(0.01, 0.99, M_crps)
  crp_scores <- tf$reshape(tf$constant(list()), c(0L))
  b_size <- as.integer(cdf_collection$shape[1])
  
  # cycle through each individual and calc their crps
  for (i in 1L:b_size) {
    
    #tf$print("i:",i)
    y_obs <- y_observed[i]
    
    # collect over grid
    pin_ball_loss <- tf$reshape(tf$constant(list()), c(0L))
    
    for (k in 1L:M_crps) {
      
      F_y_hat <- tf$gather(cdf_collection, indices = c(i - 1L))
      
      # find quantile through linear interpolation of quantile function
      quant <- lin_interpol(p_grid[k], F_y_hat, y_grid)
      
      # pinball loss
      pin_ball <- tf$cast(tf$math$less(y_obs, quant),tf$float32)
      pin_ball <- tf$subtract(pin_ball, p_grid[k])
      pin_ball <- tf$math$multiply(pin_ball, tf$subtract(quant, y_obs))
      
      pin_ball_loss <- tf$concat(list(pin_ball_loss, 
                                      tf$reshape(pin_ball, 1L)), axis = 0L)
    }
    
    scle <- tf$cast(tf$divide(2L, M_crps), tf$float32)
    crps <- tf$math$multiply(scle, tf$reduce_sum(pin_ball_loss))
    crp_scores <- tf$concat(list(crp_scores, 
                                 tf$reshape(crps, 1L)), axis = 0L)
    
  }
  return(crp_scores)
}
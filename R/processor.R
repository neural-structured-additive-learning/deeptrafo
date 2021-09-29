basis_processor <- function(term, data, output_dim = NULL, param_nr, controls){

  name <- makelayername(term, param_nr)
  bfy <- controls$y_basis_fun(data[[extractvar(term)]])
  
  # define layer  
  layer <- function(x, ...)
    return(
      layer_lambda(x, f = function(x) x)
    )
  
  list(
    data_trafo = function() bfy,
    predict_trafo = function(newdata) controls$y_basis_fun(newdata[extractvar(term)]),
    input_dim = as.integer(ncol(bfy)),
    layer = layer
  )
}

basis_prime_processor <- function(term, data, output_dim = NULL, param_nr, controls){
  
  name <- makelayername(term, param_nr)
  bfy <- controls$y_basis_fun_prime(data[[extractvar(term)]])
  
  # define layer  
  layer <- function(x, ...)
    return(
      layer_lambda(x, f = function(x) x)
    )
  
  list(
    data_trafo = function() bfy,
    predict_trafo = function(newdata) controls$y_basis_fun_prime(newdata[extractvar(term)]),
    input_dim = as.integer(ncol(bfy)),
    layer = layer
  )
}

make_atm_processor <- function(ar_lags_fun){
  
  atm_lag_processor <- function(term, data, output_dim = NULL, param_nr=3, controls=NULL)
  {
    
    # gets a formula like ~ -1 + atmlag(lag1, lag2, lag3) 
    name <- makelayername(term, param_nr)
    layer <- ar_lags_fun
    
    list(
      data_trafo = function() data[extractvar(term)],
      predict_trafo = function(newdata) newdata[extractvar(term)],
      input_dim = as.integer(length(extractvar(term))),
      layer = layer
    )
    
  }
  
  return(atm_lag_processor)
  
}
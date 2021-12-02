basis_processor <- function(term, data, output_dim = NULL, param_nr, controls){

  name <- makelayername(term, param_nr)
  is_prime <- extractval(term, "prime")
  term <- gsub("\\,\\s?prime\\s?=.*[^\\)]","",term)
  terms <- split_interaction_terms(term)

  if(is_prime){
    bfy <- controls$y_basis_fun_prime(data[[extractvar(terms[1])]])
  }else{
    bfy <- controls$y_basis_fun(data[[extractvar(terms[1])]])
  }

  suppy <- range(data[[extractvar(terms[1])]])

  # interacting term
  args <- list(term = terms[2], data = data, output_dim = output_dim,
               param_nr = param_nr, controls = controls)
  spec <- deepregression:::get_special(terms[2], specials = names(controls$procs))

  if(is.null(spec)){
    if(terms[2]=="1"){
      iat <- do.call(deepregression:::int_processor, args)
    }else{
      iat <- do.call(deepregression:::lin_processor, args)
    }
  }else iat <- do.call(controls$procs[[spec]], args)

  dim_iat <- iat$input_dim
  dim_basis <- ncol(bfy)
  penalty_iat <- iat$penalty
  penalty_basis <- controls$basis_penalty
  combined_penalty <- deepregression:::combine_penalties(list(penalty_basis,
                                                              penalty_iat),
                                                         c(dim_basis, dim_iat))

  if(is_prime){
    predict_trafo_bs <- function(newdata)
      controls$y_basis_fun_prime(newdata[[extractvar(terms[1])]],
                                 suppy = suppy)
  }else{
    predict_trafo_bs <-  function(newdata)
      controls$y_basis_fun(newdata[[extractvar(terms[1])]],
                           suppy = suppy)
  }

  if(!is_prime){

    thetas_layer <- layer_mono_multi(
      units = output_dim,
      dim_bsp = dim_basis,
      kernel_regularizer = combined_penalty,
      name = name
    )

    layer <- function(x, ...){

      input_theta_y <- tf_stride_cols(x, 1L, dim_basis)
      interact_pred <- tf_stride_cols(x, dim_basis+1L, dim_basis+dim_iat)
      AoB <- tf_row_tensor(input_theta_y, interact_pred)
      return(AoB %>% thetas_layer())

    }

  }else{

    layer <- tf$identity

  }

  plot_fun <- NULL
  get_org_values <- NULL

  if(!is.null(iat$plot_fun)){

    get_org_values <- function()
      return(list(data[[extractvar(terms[1])]],
                  iat$get_org_values()))

    plot_fun <- h1_plotfun(dim_basis)

  }

  data_trafo <- function() cbind(bfy, iat$data_trafo())
  predict_trafo <- function(newdata) cbind(predict_trafo_bs(newdata),
                                           iat$predict_trafo(newdata))

  list(
    data_trafo = data_trafo,
    predict_trafo = predict_trafo,
    input_dim = as.integer(dim_basis + dim_iat),
    layer = layer,
    coef = function(weights) as.matrix(weights),
    partial_effect = function(weights, newdata=NULL){
      X <- if(is.null(newdata)) data_trafo() else
        predict_trafo(newdata)
      return(row_tensor_by_basis(X, dim_basis) %*% weights)
    },
    plot_fun = plot_fun,
    get_org_values = get_org_values,
    get_bfy = function(newdata) predict_trafo_bs(newdata)
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

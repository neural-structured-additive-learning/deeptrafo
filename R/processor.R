ia_processor <- function(term, data, output_dim = NULL, param_nr, controls){

  name <- makelayername(term, param_nr)
  term <- gsub("ia\\((.*)\\)","\\1",term)

  # interacting term
  args <- list(term = term, data = data, output_dim = output_dim,
               param_nr = param_nr, controls = controls)
  spec <- get_special(term, specials = names(controls$procs))

  if (is.null(spec)) {
    if (term == "1") {
      iat <- do.call(int_processor, args)
    } else {
      iat <- do.call(lin_processor, args)
    }
  } else iat <- do.call(controls$procs[[spec]], args)

  dim_iat <- iat$input_dim
  dim_basis <- controls$order_bsp + 1L
  penalty_iat <- iat$penalty
  penalty_basis <- controls$basis_penalty
  combined_penalty <- combine_penalties(list(penalty_basis,
                                             penalty_iat),
                                        c(dim_basis, dim_iat))


  thetas_layer <- layer_mono_multi(
    units = output_dim,
    dim_bsp = dim_basis,
    kernel_regularizer = combined_penalty,
    name = name,
    trafo = controls$trafo
  )

  if (!is.null(spec) && "dnn" %in% names(environment(controls$procs[[spec]]))) # deep network in ia
  {

    layer <- function(bspy, iaterm, ...){

      iat_dnn_out <- iat$layer(iaterm)
      tf_row_tensor(bspy, iat_dnn_out) %>% thetas_layer()

    }

  } else {

    layer <- function(bspy, iaterm, ...){

      tf_row_tensor(bspy, iaterm) %>% thetas_layer()

    }

  }

  if (!is.null(iat$plot_fun)) {

    get_org_values <- function() return(iat$get_org_values())
    plot_fun <- iat$plot_fun

  } else {

    plot_fun <- NULL
    get_org_values <- NULL

  }

  data_trafo <- function() iat$data_trafo()
  predict_trafo <- function(newdata) iat$predict_trafo(newdata)

  list(
    data_trafo = data_trafo,
    predict_trafo = predict_trafo,
    input_dim = as.integer(dim_iat),
    layer = layer,
    coef = function(weights) as.matrix(weights),
    partial_effect = function(weights, newdata = NULL){
      # X <- if(is.null(newdata)) data_trafo() else
      #   predict_trafo(newdata)
      # return(function(y) row_tensor(y, X) %*% weights)
    },
    plot_fun = plot_fun,
    get_org_values = get_org_values
  )
}

basis_processor_lower <- function(term, data, output_dim = NULL, param_nr,
                                  controls) {
  basis_processor(term = term, data = data, output_dim = output_dim,
                  param_nr = param_nr, controls = controls, lower = TRUE)
}

basis_processor <- function(term, data, output_dim = NULL, param_nr, controls,
                            lower = FALSE){

  ybfun <- if (lower) controls$y_basis_fun_lower else controls$y_basis_fun

  name <- makelayername(term, param_nr)
  bfy <- ybfun(data[[extractvar(term)]])
  suppy <- controls$supp(data[[extractvar(term)]])

  dim_basis <- ncol(bfy)
  penalty_basis <- controls$basis_penalty
  predict_trafo_bs <-  function(newdata)
    ybfun(newdata[[extractvar(term)]], suppy = suppy)

  layer = tf$identity

  list(
    data_trafo = function() bfy,
    predict_trafo = predict_trafo_bs,
    input_dim = as.integer(ncol(bfy)),
    layer = layer,
    penalty = penalty_basis
  )
}

basisprime_processor <- function(term, data, output_dim = NULL, param_nr, controls){

  name <- makelayername(term, param_nr)
  bfy <- controls$y_basis_fun_prime(data[[extractvar(term)]])
  suppy <- controls$supp(data[[extractvar(term)]])

  dim_basis <- ncol(bfy)
  predict_trafo_bs <- function(newdata)
    controls$y_basis_fun_prime(newdata[[extractvar(term)]],
                               suppy = suppy)

  layer = tf$identity

  list(
    data_trafo = function() bfy,
    predict_trafo = predict_trafo_bs,
    input_dim = as.integer(ncol(bfy)),
    layer = layer
  )
}

atm_lag_processor_factory <- function(rvar){

  return(
    atm_lag_processor <- function(term, data, output_dim = NULL, param_nr=4, controls=NULL)
    {

      name <- makelayername(term, param_nr)

      layer <- eval_bsp_tf(order = controls$order_bsp, controls$supp(data[extractvar(rvar)]))

      list(
        data_trafo = function() data[extractvar(term)],
        predict_trafo = function(newdata) newdata[extractvar(term)],
        input_dim = 1L,
        layer = layer
      )

    }
  )

}

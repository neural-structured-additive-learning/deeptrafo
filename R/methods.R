#' Generic methods for deeptrafo objects
#'
#' @param x deeptrafo object
#' @param which which effect to plot, default selects all.
#' @param which_param character; either \code{"interacting"} or \code{"shifting"}
#' 1 corresponds to the shift term, 2 to the interaction term.
#' @param only_data logical, if TRUE, only the data for plotting is returned
#' @param grid_length the length of an equidistant grid at which a two-dimensional function
#' is evaluated for plotting.
#' @param eval_grid logical; should plot be evaluated on a grid
#' @param ... further arguments, passed to fit, plot or predict function
#'
#'
#' @method plot deeptrafo
#' @export
#' @rdname methodTrafo
#'
plot.deeptrafo <- function(
  x,
  which = NULL,
  # which of the nonlinear structured effects
  which_param = c("shifting", "interacting"), # for which parameter
  only_data = FALSE,
  grid_length = 40,
  ... # passed to plot function
)
{

  which_param <- match.arg(which_param)

  get_weight_fun <- switch(
    which_param,
    "interacting" = get_weight_by_name_ia,
    "shifting" = get_weight_by_name
  )

  which_param <- map_param_string_to_index(which_param)
  class(x) <- class(x)[-1]

  return(plot(x, which = which, which_param = which_param,
              only_data = only_data, grid_length = grid_length,
              get_weight_fun = get_weight_fun, ...))

}

get_weight_by_name_ia <- function(x, name, param_nr)
{

  matrix(get_weight_by_name(x, name, param_nr),
         ncol = x$init_params$trafo_options$order_bsp + 1L)

}

#' @param x deeptrafo object
#' @param which_param integer, indicating for which distribution parameter
#' coefficients should be returned (default is first parameter)
#' @param type either NULL (all types of coefficients are returned),
#' "linear" for linear coefficients or "smooth" for coefficients of;
#' Note that \code{type} is currently not used for \code{"interacting"}
#' @param ... further arguments, passed to fit, plot or predict function
#'
#' @method coef deeptrafo
#' @export
#' @rdname methodTrafo
#'
coef.deeptrafo <- function(
  object,
  which_param = c("shifting", "interacting"),
  type = NULL,
  ...
)
{

  which_param <- match.arg(which_param)

  is_interaction <- which_param == "interacting"
  which_param <- map_param_string_to_index(which_param)

  # else, return lags
  class(object) <- class(object)[-1]
  ret <- coef(object, which_param = which_param, type = type)

  if (is_interaction) {
    ret <- lapply(ret, function(r)
      reshape_softplus_cumsum(r, order_bsp_p1 = get_order_bsp_p1(object)))

    if (object$init_params$response_type == "ordered")
      ret <- lapply(ret, function(r) r[-nrow(r), , drop = FALSE])

  }

  return(ret)

}


#'
#' @param object a deeptrafo model
#' @param newdata optional new data, either data.frame or list
#' @param ... not used atm
#' @return returns a function with two parameters: the actual response
#' and \code{type} in \code{c('trafo', 'pdf', 'cdf', 'interaction')}
#' determining the returned value
#'
#' @method predict deeptrafo
#' @export
#'
#' @rdname methodTrafo
#'
predict.deeptrafo <- function(
  object,
  newdata = NULL,
  y = newdata[[object$init_params$response_varname]],
  type = c("trafo", "pdf", "cdf", "interaction", "shift", "output"),
  batch_size = NULL,
  ...
)
{

  # TODO: make prediction possible for one observation only
  type <- match.arg(type)

  rname <- object$init_params$response_varname
  rtype <- object$init_params$response_type
  order <- object$init_params$trafo_options$order_bsp
  fam <- object$init_params$family
  discrete <- as.numeric(rtype %in% c("count", "ordered"))
  bd <- get_bd(fam)

  if (is.null(y))
    y <- object$init_params$response

  ry <- response(y)
  cleft <- ry[, "cleft", drop = FALSE]
  cint <- ry[, "cinterval", drop = FALSE]
  cright <- ry[, "cright", drop = FALSE]

  if (is.null(newdata))
    newdata <- prepare_data(object$init_params$parsed_formulas_contents,
                            gamdata = object$init_params$gamdata$data_trafos)

  newdata[[rname]] <- y

  mod_output <- fitted.deeptrafo(object, newdata, batch_size = batch_size)

  if (type == "output")
    return(mod_output)

  w_eta <- mod_output[, 1, drop = FALSE]
  aTtheta <- mod_output[, 2, drop = FALSE]
  apTtheta <- mod_output[, 3, drop = FALSE]

  if (type == "interaction")
    return(cbind(interaction = as.matrix(aTtheta), as.data.frame(newdata)))

  if (type == "shift")
    return(cbind(shift = as.matrix(w_eta), as.data.frame(newdata)))

  ytransf <- aTtheta + w_eta
  yprimeTrans <- apTtheta + discrete * w_eta


  if (discrete) {

    pdf <- cint * as.matrix(tfd_cdf(bd, ytransf) - tfd_cdf(bd, yprimeTrans)) +
      cleft * tfd_cdf(bd, ytransf) + cright * tfd_survival_function(bd, ytransf)

    cdf <- (cleft + cint) * as.matrix(tfd_cdf(bd, ytransf)) +
      cright * as.matrix(tfd_cdf(bd, rep(1e8, nrow(cright))))

  } else {

    pdf <- as.matrix(tfd_prob(bd, ytransf)) * as.matrix(yprimeTrans)

    cdf <- bd %>% tfd_cdf(ytransf) %>% as.matrix

  }

  theta <- get_theta(object)

  ret <- switch(
    type,
    "trafo" = (ytransf %>% as.matrix),
    "pdf" = pdf %>% as.matrix,
    "cdf" = cdf %>% as.matrix# ,
    # "grid_trafo" = grid_eval %>% as.matrix,
    # "grid_pdf" = pdf %>% as.matrix,
    # "grid_cdf" = (bd %>% tfd_cdf(grid_eval) %>% as.matrix)
  )

  return(ret)

}

#'
#' @param object a deeptrafo model
#' @param newdata optional new data, either data.frame or list
#' @param batch_size integer; optional, useful if data is too large
#' @param convert_fun function; to convert the TF tensor
#' @param ... not used atm
#' @return returns a function with two parameters: the actual response
#' and \code{type} in \code{c('trafo', 'pdf', 'cdf', 'interaction')}
#' determining the returned value
#'
#' @method fitted deeptrafo
#' @export
#'
#' @rdname methodTrafo
#'
fitted.deeptrafo <- function(
  object,
  newdata = NULL,
  batch_size = NULL,
  convert_fun = as.matrix,
  ...)
{

  if(length(object$init_params$image_var)>0 | !is.null(batch_size)) {

    mod_output <- predict_gen(object, newdata, batch_size,
                              apply_fun = function(x) x,
                              convert_fun = convert_fun)

  } else{

    if(is.null(newdata)) {

      newdata <- prepare_data(object$init_params$parsed_formulas_contents,
                              gamdata = object$init_params$gamdata$data_trafos)

    } else{

      newdata <- prepare_newdata(object$init_params$parsed_formulas_contents, newdata,
                                 gamdata = object$init_params$gamdata$data_trafos)

    }

    mod_output <- object$model(newdata)

  }

  return(convert_fun(mod_output))

}

map_param_string_to_index <- function(which_param)
{

  switch(
    which_param,
    "interacting" = 2,
    "shifting" = 3
  )

}

#' Log-likelihood method for deeptrafo objects
#'
#' @method logLik deeptrafo
#' @param object deeptrafo object;
#' @param y vector; optional response
#' @param newdata list or data.frame; optional new data
#' @param convert_fun function; applied to the log-likelihood values of all observations
#' @param ... currently not used
#'
#' @export
#'
logLik.deeptrafo <- function(
  object,
  y = NULL,
  newdata = NULL,
  convert_fun = function(x, ...) - sum(x, ...),
  ...
)
{

  # TODO: logLik with newdata
  if (is.null(newdata)) {
    y <- object$init_params$y
    y_pred <- fitted.deeptrafo(object)
  }

 convert_fun(object$model$loss(y, y_pred)$numpy())

}

#' Simulate method for deeptrafo objects
#'
#' @method simulate deeptrafo
#' @param object \code{"deeptrafo"} object
#' @param nsim number of simulations; defaults to 1
#' @param seed seed for generating samples; defaults to \code{NULL}
#' @param newdata list or data.frame; optional new data
#' @param ... further arguments to \link[predict.deeptrafo]{predict.deeptrafo}
#'
#' @export
#'
simulate.deeptrafo <- function(object, newdata = NULL, nsim = 1,
                               seed = NULL, ...) {

  rtype <- object$init_params$response_type
  rvar <- object$init_params$response_varname

  if (!is.null(newdata))
    ry <- newdata[[rvar]]
  else
    ry <- object$init_params$response


  if (rtype != "ordered")
    stop("Simulate not yet implemented for response types other than ordered.")

  lvls <- sort(unique(ry))


  cdf <- do.call("cbind", lapply(lvls, function(x) {
    newy <- ry
    newy[] <- x
    predict(object, newdata = newdata, y = newy, type = "cdf") # , ... = ...)
  }))

  pmf <- apply(cbind(0, cdf), 1, diff)

  ret <- lapply(1:nsim, function(x) {
    ordered(apply(pmf, 2, function(probs) {
      which(rmultinom(n = 1, size = 1, prob = probs) == 1)
    }), levels = levels(lvls))
  })

  if (nsim == 1)
    ret <- ret[[1]]

  ret

}

# Helpers -----------------------------------------------------------------

get_bd <- function(family) {
  switch(family,
         "normal" = tfd_normal(loc = 0, scale = 1),
         "logistic" = tfd_logistic(loc = 0, scale = 1),
         "gumbel" = tfd_gumbel(loc = 0, scale = 1)
  )
}

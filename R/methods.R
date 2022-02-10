#' Generic methods for deeptrafo objects
#'
#' @param x deeptrafo object
#' @param which which effect to plot, default selects all.
#' @param which_param character; either \code{"h1"} or \code{"h2"}
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
  which_param = "h1", # for which parameter
  only_data = FALSE,
  grid_length = 40,
  ... # passed to plot function
)
{

  get_weight_fun = if(which_param == "h1")
    get_weight_by_name_ia else get_weight_by_name

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
#' Note that \code{type} is currently not used for \code{"h1"}
#' @param ... further arguments, passed to fit, plot or predict function
#'
#' @method coef deeptrafo
#' @export
#' @rdname methodTrafo
#'
coef.deeptrafo <- function(
  object,
  which_param = "h1",
  type = NULL,
  ...
)
{

  is_interaction <- which_param == "h1"
  which_param <- map_param_string_to_index(which_param)

  # else, return lags
  class(object) <- class(object)[-1]
  ret <- coef(object, which_param = which_param, type = type)

  if (is_interaction)
    ret <- lapply(ret, function(r)
      reshape_softplus_cumsum(r, order_bsp_p1 = get_order_bsp_p1(object)))

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
  ...
)
{

  # TODO: make prediction possible for one observation only; fix type mismatch pdf grid
  trafo_fun <- function(y, type = c("trafo", "pdf", "cdf", "interaction", "shift", "output"),
                        grid = FALSE, batch_size = NULL)
  {
    type <- match.arg(type)

    if(is.null(newdata))
      newdata <- deepregression:::prepare_data(object$init_params$parsed_formulas_contents,
                                               gamdata = object$init_params$gamdata$data_trafos)

    newdata[[object$init_params$response_varname]] <- y

    mod_output <- fitted.deeptrafo(object, newdata, batch_size = batch_size)

    if(type=="output") return(mod_output)

    w_eta <- mod_output[, 1, drop = FALSE]
    aTtheta <- mod_output[, 2, drop = FALSE]

    if(type=="interaction"){

      ret <- cbind(interaction = as.matrix(aTtheta),
                   as.data.frame(newdata)
      )

      return(ret)

    }

    if(type=="shift"){

      return(cbind(shift = as.matrix(w_eta),
                   as.data.frame(newdata)))

    }
    ytransf <- aTtheta + w_eta
    yprimeTrans <- mod_output[, 3, drop = FALSE]

    theta <- get_theta(object)

    if(grid)
    {

      ymat <- lapply(object$init_params$parsed_formulas_contents[[1]],
                     function(pp) pp$predict(newdata))

      obspp1 <- object$init_params$trafo_options$order_bsp + 1
      ay <- ymat[[1]] #[,1:obspp1]
      ayprime <- ymat[[2]]
      xmat <- do.call("cbind", lapply(object$init_params$parsed_formulas_contents[[2]],
                                      function(pp) pp$predict(newdata)))

      grid_eval <- t(xmat%*%t(ay%*%theta))

      grid_eval <- grid_eval +
        t(as.matrix(w_eta)[,rep(1,nrow(grid_eval))])

      if(type=="pdf")
        grid_prime_eval <- t(xmat%*%t(ayprime%*%theta))

      type <- paste0("grid_",type)

    }

    ret <- switch (type,
                   trafo = (ytransf %>% as.matrix),
                   pdf = ((tfd_normal(0,1) %>% tfd_prob(ytransf) %>%
                             as.matrix)*as.matrix(yprimeTrans)),
                   cdf = (tfd_normal(0,1) %>% tfd_cdf(ytransf) %>%
                            as.matrix),
                   grid_trafo = grid_eval,
                   grid_pdf = ((tfd_normal(0,1) %>% tfd_prob(grid_eval) %>%
                                  as.matrix)*as.matrix(grid_prime_eval)),
                   grid_cdf = (tfd_normal(0,1) %>% tfd_cdf(grid_eval) %>%
                                 as.matrix)
    )

    return(ret)

  }

  return(trafo_fun)

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

  if(length(object$init_params$image_var)>0 | !is.null(batch_size)){

    mod_output <- predict_gen(object, newdata, batch_size,
                              apply_fun = function(x) x,
                              convert_fun = convert_fun)

  }else{

    if(is.null(newdata)){

      newdata <- deepregression:::prepare_data(object$init_params$parsed_formulas_contents,
                                               gamdata = object$init_params$gamdata$data_trafos)

    }else{

      newdata <- deepregression:::prepare_newdata(object$init_params$parsed_formulas_contents, newdata,
                                                  gamdata = object$init_params$gamdata$data_trafos)

    }

    mod_output <- object$model(newdata)

  }

  return(convert_fun(mod_output))

}

map_param_string_to_index <- function(which_param)
{

  switch (which_param,
          "h1" = 2,
          "h2" = 3
  )

}


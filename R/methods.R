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

  which_param <- switch (which_param,
    "h1" = 3,
    "h2" = 4
  )

  class(x) <- class(x)[-1]
  return(plot(x, which = which, which_param = which_param,
              only_data = only_data, grid_length = grid_length, ...))

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

	if (which_param == "h1")
		return(get_theta(object))
	if (which_param == "h2")
		return(get_shift(object, type = type))

  # else, return lags
  class(object) <- class(object)[-1]
  return(coef(object, which_param = 5, type = type))

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
#'
#' @rdname methodTrafo
#'
predict.deeptrafo <- function(
  object,
  newdata = NULL,
  which = NULL,
  cast_float = FALSE,
  ...
)
{

  if(is.null(newdata)){
    inpCov <- prepare_data(object$init_params$parsed_formulas_contents[-1*1:2])
  }else{
    # preprocess data
    if(is.data.frame(newdata)) newdata <- as.list(newdata)
    inpCov <- prepare_newdata(object$init_params$parsed_formulas_contents[-1*1:2], newdata)
  }

  # TODO: make prediction possible for one observation only; fix type mismatch pdf grid
  trafo_fun <- function(y, type = c("trafo", "pdf", "cdf", "interaction", "shift", "output", "sample"),
                        which = NULL, grid = FALSE, batch_size = NULL)
  {
    type <- match.arg(type)

    # if(!is.null(minval)) y <- y - sum(minval*get_theta(object))

    ay_aPrimey <- prepare_newdata(object$init_params$parsed_formulas_contents[1:2], data.frame(y=y))
    inpCov <- c(ay_aPrimey, inpCov)
    mod_output <- fitted.deeptrafo(object, inpCov, y, batch_size = batch_size)
    if(type=="output") return(mod_output)
    w_eta <- mod_output[, 1, drop = FALSE]
    aTtheta <- mod_output[, 2, drop = FALSE]
    # if(!is.null(minval))
    #   aTtheta <- aTtheta - sum(minval*get_theta(object))
    if(type=="interaction"){

      if(is.null(newdata))
        newdata <- object$init_params$data

      ret <- cbind(interaction = as.matrix(aTtheta),
                   as.data.frame(newdata)
      )

      return(ret)

    }

    if(type=="shift"){

      if(is.null(newdata))
        newdata <- object$init_params$data

      return(cbind(shift = as.matrix(w_eta),
                   as.data.frame(newdata)))

    }
    ytransf <- aTtheta + w_eta
    yprimeTrans <- mod_output[, 3, drop = FALSE]
    # if(!is.null(minval))
    #   yprimeTrans + sum(minval*get_theta(object))
    theta <- get_theta(object)
    if(grid)
    {

      grid_eval <- t(as.matrix(
        tf$matmul(inpCov[[2]],
                  tf$transpose(
                    tf$matmul(ay,
                              tf$cast(theta, tf$float32)
                    )))))
      grid_eval <- grid_eval +
        t(as.matrix(w_eta)[,rep(1,nrow(grid_eval))])

      if(type=="pdf")
        grid_prime_eval <- t(as.matrix(
          tf$matmul(inpCov[[2]],
                    tf$transpose(
                      tf$matmul(aPrimey,
                                tf$cast(theta, tf$float32)
                      )))))


    }

    if(grid) type <- paste0("grid_",type)

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
                                 as.matrix),
                   sample
    )

    return(ret)

  }

  return(trafo_fun)

}

#'
#' @param object a deeptrafo model
#' @param newdata optional new data, either data.frame or list
#' @param batch_size integer; optional, useful if data is too large
#' @param ... not used atm
#' @return returns a function with two parameters: the actual response
#' and \code{type} in \code{c('trafo', 'pdf', 'cdf', 'interaction')}
#' determining the returned value
#'
#' @method predict deeptrafo
#'
#' @rdname methodTrafo
#'
fitted.deeptrafo <- function(object, newdata, batch_size = NULL, ...)
{


  if(length(object$init_params$image_var)>0){

    data_tab <- newdata

    # prepare generator
    max_data <- NROW(data_image)
    if(is.null(batch_size)) batch_size <- 32
    steps_per_epoch <- ceiling(max_data/batch_size)

    mod_output <- predict_generator(object, newdata, batch_size, apply_fun=NULL)

  }else{

    if(is.null(batch_size)){

      mod_output <- object$model(newdata)

    }else{

        max_data <- NROW(newdata[[1]])
        steps_per_epoch <- ceiling(max_data/batch_size)

        mod_output <- lapply(1:steps_per_epoch,
                             function(i){
                               index <- (i-1)*batch_size + 1:batch_size
                               object$model(lapply(newdata, function(x) subset_array(x, index)))
                             })
        mod_output <- do.call("rbind", lapply(mod_output, as.matrix))

    }

  }

  return(mod_output)

}

#' Function to return the weights of the shift term
#'
#' @param x the fitted deeptrafo object
#' @param type see coef.deeptrafo
#'
#' @export
get_shift <- function(x, type = NULL)
{

  pfc <- x$init_params$parsed_formulas_contents$h2
  to_return <- get_type_pfc(pfc, type)

  names <- deepregression:::get_names_pfc(pfc)[as.logical(to_return)]

  check_names <- names
  check_names[check_names=="(Intercept)"] <- "1"
  coefs <- lapply(1:length(check_names), function(i)
    pfc[[i]]$coef(deepregression:::get_weight_by_name(x, check_names[i], 4)))

  names(coefs) <- names
  return(coefs)

}

#' Function to return the weights of the theta term
#'
#' @param x the fitted deeptrafo object
#'
#' @export
get_theta <- function(x)
{

  names_weights <- sapply(x$model$trainable_weights, function(x) x$name)
  reshape_softplus_cumsum(
    as.matrix(x$model$weights[[grep("constraint_mono_layer", names_weights)]] + 0),
    order_bsp_p1 = get_order_bsp_p1(x)
  )

}





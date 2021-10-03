#' @title Predict based on a deeptrafo object
#'
#' @param object a deeptrafo model
#' @param newdata optional new data, either data.frame or list
#' @param ... not used atm
#' @return returns a function with two parameters: the actual response
#' and \code{type} in \code{c('trafo', 'pdf', 'cdf', 'interaction')}
#' determining the returned value
#'
#' @export
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
    mod_output <- evaluate.deeptrafo(object, inpCov, y, batch_size = batch_size)
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

evaluate.deeptrafo <- function(object, newdata, y, batch_size = NULL)
{
  
  
  if(length(object$init_params$image_var)>0){
    
    data_tab <- newdata
    
    # prepare generator
    max_data <- NROW(data_image)
    if(is.null(batch_size)) batch_size <- 32
    steps_per_epoch <- ceiling(max_data/batch_size)

    mod_output <- object$model$predict_generator(object, newdata, batch_size, apply_fun=NULL)
    
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
#'
#' @export
get_shift <- function(x)
{

  pfc <- x$init_params$parsed_formulas_contents$h2
  names <- get_names_pfc(pfc)
  check_names <- names
  check_names[check_names=="(Intercept)"] <- "1"
  coefs <- lapply(1:length(check_names), function(i) 
    pfc[[i]]$coef(get_weight_by_name(x, check_names[i], 4)))
  
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





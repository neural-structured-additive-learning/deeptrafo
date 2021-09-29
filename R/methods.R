#' @title Generic functions for deeptrafo models
#'
#' @param x deepregression object
#' @param which which effect to plot, default selects all.
#' @param which_param integer, either 1 or 2.
#' 1 corresponds to the shift term, 2 to the interaction term.
#' @param plot logical, if FALSE, only the data for plotting is returned
#' @param grid_length the length of an equidistant grid at which a two-dimensional function
#' is evaluated for plotting.
#' @param eval_grid logical; should plot be evaluated on a grid
#' @param ... further arguments, passed to fit, plot or predict function
#'
#' @method plot deeptrafo
#' @export
#' @rdname methodTrafo
#'
plot.deeptrafo <- function(
  x,
  which = NULL,
  # which of the nonlinear structured effects
  which_param = 1, # for which parameter
  plot = TRUE,
  grid_length = 40,
  eval_grid = FALSE,
  ... # passed to plot function
)
{
  this_ind <- x$init_params$ind_structterms[[which_param]]
  if(all(this_ind$type!="smooth")) return("No smooth effects. Nothing to plot.")
  if(is.null(which)) which <- 1:length(which(this_ind$type=="smooth"))
  plus_number_lin_eff <- sum(this_ind$type=="lin")

  plotData <- vector("list", length(which))
  org_feature_names <-
    names(x$init_params$l_names_effects[[which_param]][["smoothterms"]])
  if(which_param==1){
    phi <- matrix(get_shift(x), ncol=1)
  }else{
    phi <- get_theta(x)
    # FIXME: Is the following necessary?
    if(ncol(phi)>1) phi <- t(phi)
  }

  for(w in which){

    nam <- org_feature_names[w]
    this_ind_this_w <- do.call("Map",
                               c(":", as.list(this_ind[w+plus_number_lin_eff,
                                                       c("start","end")])))[[1]]
    BX <-
      x$init_params$parsed_formulae_contents[[
        which_param]]$smoothterms[[nam]][[1]]$X

    plotData[[w]] <-
      list(org_feature_name = nam,
           value = sapply(gsub("by = ","",strsplit(nam,",")[[1]]), function(xx)
             x$init_params$data[[xx]]),
           design_mat = BX,
           coef = phi[this_ind_this_w,],
           partial_effects = BX%*%phi[this_ind_this_w,])

    nrcols <- pmax(NCOL(plotData[[w]]$value), length(unlist(strsplit(nam,","))))

    if(plot | nrcols > 2 | eval_grid){
      if(which_param==1){

        if(nrcols==1)
        {

          plot(partial_effects[order(value),1] ~ sort(value[,1]),
               data = plotData[[w]][c("value", "partial_effects")],
               main = paste0("s(", nam, ")"),
               xlab = nam,
               ylab = "partial effect",
               ...)

        }else if(nrcols==2){
          sTerm <- x$init_params$parsed_formulae_contents[[which_param]]$smoothterms[[w]][[1]]
          this_x <- do.call(seq, c(as.list(range(plotData[[w]]$value[,1])),
                                   list(l=grid_length)))
          this_y <- do.call(seq, c(as.list(range(plotData[[w]]$value[,2])),
                                   list(l=grid_length)))
          df <- as.data.frame(expand.grid(this_x,
                                          this_y))
          colnames(df) <- gsub("by = ","",strsplit(nam,",")[[1]])
          pmat <- PredictMat(sTerm, data = df)
          if(attr(x$init_params$parsed_formulae_contents[[which_param]],"zero_cons") & 
             !is.null(x$init_params$parsed_formulae_contents[[which_param]]$linterms) & 
             !any(grepl("by =", strsplit(nam,",")[[1]])))
            pmat <- orthog_structured_smooths(pmat,P=NULL,L=matrix(rep(1,nrow(pmat)),ncol=1))
          pred <- pmat%*%phi[this_ind_this_w,]
          #this_z <- plotData[[w]]$partial_effect

          plotData[[w]] <- list(df = df,
                                design_mat = pmat,
                                coef = phi[this_ind_this_w,],
                                pred = pred)

          if(plot)
            suppressWarnings(
              filled.contour(
                this_x,
                this_y,
                matrix(pred, ncol=length(this_y)),
                ...,
                xlab = colnames(df)[1],
                ylab = colnames(df)[2],
                # zlab = "partial effect",
                main = sTerm$label
              )
            )
          # warning("Plotting of effects with ", nrcols, "
          #         covariate inputs not supported yet.")
        }else if(nrcols==3){

          if(plot) warning("Will only return the plot data for 3d effects.")

          sTerm <- x$init_params$parsed_formulae_contents[[which_param]]$smoothterms[[w]][[1]]
          this_x <- do.call(seq, c(as.list(range(plotData[[w]]$value[,1])),
                                   list(l=grid_length)))
          this_y <- do.call(seq, c(as.list(range(plotData[[w]]$value[,2])),
                                   list(l=grid_length)))
          this_z <- do.call(seq, c(as.list(range(plotData[[w]]$value[,3])),
                                   list(l=grid_length)))
          df <- as.data.frame(expand.grid(this_x,
                                          this_y,
                                          this_z))
          colnames(df) <- sTerm$term
          if(sTerm$by!="NA") colnames(df)[ncol(df)] <- sTerm$by
          pmat <- PredictMat(sTerm, data = df)
          if(attr(x$init_params$parsed_formulae_contents[[which_param]],"zero_cons") & 
             !is.null(x$init_params$parsed_formulae_contents[[which_param]]$linterms) & 
             sTerm$by=="NA")
            pmat <- orthog_structured_smooths(pmat,P=NULL,L=matrix(rep(1,nrow(pmat)),ncol=1))
          pred <- pmat%*%phi[this_ind_this_w,]
          #this_z <- plotData[[w]]$partial_effect

          plotData[[w]] <- list(df = df,
                                design_mat = pmat,
                                coef = phi[this_ind_this_w,],
                                pred = pred)

        }else{
          warning("Plotting of effects with ", nrcols,
                  " covariate inputs not supported.")
        }
      }else{ # plot effects in theta
        
        suppressWarnings(
          matplot(
            #sort(plotData[[w]]$value[,1]),
            #1:ncol(plotData[[w]]$partial_effects),
            x=sort(plotData[[w]]$value[,1]),
            y=(plotData[[w]]$partial_effects[order(plotData[[w]]$value[,1]),]),
            ...,
            xlab = plotData[[w]]$org_feature_name,
            ylab = paste0("partial effects ", plotData[[w]]$org_feature_name),
            # zlab = "partial effect",
            type = "l"
          )
        )

      }
    }
  }

  invisible(plotData)
}


#' Predict based on a deeptrafo object
#'
#' @param object a deeptrafo model
#' @param newdata optional new data, either data.frame or list
#' @param ... not used atm
#' @return returns a function with two parameters: the actual response
#' and \code{type} in \code{c('trafo', 'pdf', 'cdf', 'interaction')}
#' determining the returned value
#'
#' @export
#' @rdname methodDR
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
    inpCov <- unname(object$init_params$input_cov)
  }else{
    # preprocess data
    if(is.data.frame(newdata)) newdata <- as.list(newdata)
    inpCov <- prepare_data(object, newdata, pred=TRUE)
    if(cast_float) inpCov <- lapply(inpCov, function(x) tf$constant(x, dtype = "float32"))
    if(length(inpCov)==2 && !is.null(names(inpCov)) && names(inpCov)[2]=="minval")
    {
      minval <- inpCov[[2]]
      inpCov <- inpCov[[1]]
    }else{
      minval <- NULL
    }
    inpCov <- c(inpCov, list(NULL), list(NULL))
  }
  
  image_data <- NULL
  if(length(object$init_params$image_var)>0){
    if(!is.null(newdata)){
      image_data <- as.data.frame(newdata[names(object$init_params$image_var)], 
                                  stringsAsFactors = FALSE)
    }else{
      image_data <- as.data.frame(object$init_params$data[names(object$init_params$image_var)], 
                                  stringsAsFactors = FALSE) 
    }
  }

  # TODO: make prediction possible for one observation only; fix type mismatch pdf grid
  trafo_fun <- function(y, type = c("trafo", "pdf", "cdf", "interaction", "shift", "output", "sample"),
                        which = NULL, grid = FALSE, batch_size = NULL)
  {
    type <- match.arg(type)

    # if(!is.null(minval)) y <- y - sum(minval*get_theta(object))

    ay <- tf$cast(object$init_params$y_basis_fun(y), tf$float32)
    aPrimey <- tf$cast(object$init_params$y_basis_fun_prime(y), tf$float32)
    inpCov[length(inpCov)-c(1,0)] <- list(ay, aPrimey)
    mod_output <- evaluate.deeptrafo(object, inpCov, y, image_data, batch_size = batch_size)
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

      if(ncol(mod_output)==5)
        ret <- cbind(ret, correction = as.matrix(mod_output[,4,drop=FALSE]))

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

evaluate.deeptrafo <- function(object, newdata, y, data_image, batch_size = NULL)
{
  
  
  if(length(object$init_params$image_var)>0){
    
    data_tab <- list(newdata, tf$cast(matrix(y, ncol=1), tf$float32))
    
    # prepare generator
    max_data <- NROW(data_image)
    if(is.null(batch_size)) batch_size <- 32
    steps_per_epoch <- ceiling(max_data/batch_size)
    
    generator <- make_generator(data_image, data_tab, batch_size, 
                                # FIXME: multiple images
                                target_size = unname(unlist(object$init_params$image_var)[1:2]),
                                color_mode = unname(ifelse(
                                  unlist(object$init_params$image_var)[3]==3, 
                                  "rgb", "grayscale")),
                                x_col = names(object$init_params$image_var),
                                is_trafo = object$init_params$family=="transformation_model", 
                                shuffle = FALSE)
    
    mod_output <- object$model$predict(generator)
    
  }else{
    
    if(is.null(batch_size)){
      
      mod_output <- object$model(list(newdata, tf$cast(matrix(y, ncol=1), tf$float32)))
    
    }else{
      
        max_data <- NROW(newdata[[1]])
        steps_per_epoch <- ceiling(max_data/batch_size)
        
        mod_output <- lapply(1:steps_per_epoch, 
                             function(i){
                               index <- (i-1)*batch_size + 1:batch_size
                               object$model(list(lapply(newdata, function(x) subset_array(x, index)), 
                                            tf$cast(matrix(y[index], ncol=1), tf$float32)))
                             })
        mod_output <- do.call("rbind", lapply(mod_output, as.matrix))
    
    }
    
  }
  
  return(mod_output)
  
}

#' Function to return the shift term
#'
#' @param x the fitted deeptrafo object
#'
#' @export
get_shift <- function(x)
{

  stopifnot("deeptrafo" %in% class(x))
  names_weights <- sapply(x$model$trainable_weights, function(x) x$name)
  lin_names <- grep("structured_linear_1", names_weights)
  nonlin_names <- grep("structured_nonlinear_1", names_weights)
  if(length(c(lin_names, nonlin_names))==0)
    stop("Not sure which layer to access for shift. ", 
         "Have you specified a structured shift predictor?")
  -1 * as.matrix(x$model$trainable_weights[[c(lin_names, nonlin_names)]] + 0)

}

#' Function to return the theta term
#'
#' @param x the fitted deeptrafo object
#'
#' @export
get_theta <- function(x)
{

  stopifnot("deeptrafo" %in% class(x))
  names_weights <- sapply(x$model$trainable_weights, function(x) x$name)
  reshape_softplus_cumsum(
    as.matrix(x$model$weights[[grep("constraint_mono_layer", names_weights)]] + 0),
    order_bsp_p1 = x$init_params$order_bsp + 1
  )

}

#' Function to return the minval term
#'
#' @param x the fitted deeptrafo object
#'
#' @details This value is only available if \code{addconst_interaction}
#' was specified in the model call.
#'
#' @export
get_minval <- function(x)
{
  stopifnot("deeptrafo" %in% class(x))
  attr(x$init_params$parsed_formulae_contents[[2]], "minval")

}

#' Function to set the weights of a deepregression object
#'
#' @param x deepregression object
#' @param weights a matrix with weights
#' @param param integer; for which parameter to set the weights
#' @param type character; for which type of layer to set the weights;
#' 
#' @export
#'
set_weights <- function(x, 
                        weights, 
                        param = NULL, 
                        type = c("linear", "nonlinear", "lasso", "ridge", "elasticnet"))
{
  
  type <- match.arg(type)
  name <- switch(type,
                 lasso = paste0("structured_lasso_", param),
                 ridge = paste0("structured_ridge_", param),
                 elasticnet = paste0("structured_elastnet_", param),
                 linear = paste0("structured_linear_", param),
                 nonlinear = paste0("structured_nonlinear_", param)
  )
                           

  x$model$get_layer(name)$set_weights(weights)
  
}

#' Return partial effect of one smooth term
#' 
#' @param object deepregression object
#' @param name string; for partial match with smooth term
#' @param return_matrix logical; whether to return the design matrix or
#' @param which_param integer; which distribution parameter
#' the partial effect (\code{FALSE}, default)
#' @param newdata data.frame; new data (optional)
#' 
#' @export
#' 
get_partial_effect <- function(object, name, return_matrix = FALSE, 
                               which_param = 1, newdata = NULL)
{
  
  if(is.null(newdata)) newdata <- object$init_params$data
  
  nms <- names(object$init_params$parsed_formulae_contents[[which_param]]$smoothterms)
  match <- grep(name, nms)
  if(all(!match)) stop("No matching name found.")
  
  sTerm <- object$init_params$parsed_formulae_contents[[which_param]]$smoothterms[[match]]
  if(is.list(sTerm)) sTerm <- sTerm[[1]]
  
  pmat <- PredictMat(sTerm, data = newdata)
  if(attr(object$init_params$parsed_formulae_contents[[which_param]],"zero_cons"))
    pmat <- orthog_structured_smooths(pmat,P=NULL,L=matrix(rep(1,nrow(pmat)),ncol=1))
  if(return_matrix) return(pmat)
  
  coefs <- coef(object, params = which_param, type = "smooth")[[1]]
  coefs <- coefs[grepl(name, names(coefs))]

  return(pmat%*%coefs)
  
}


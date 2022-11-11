
# Cotram basis ------------------------------------------------------------

.get_eval_cotram <- function(order_bsp, support) {
  function(y, orderbsp = order_bsp, suppy = support) {
    eval_bsp(log(1 + y), order = orderbsp, supp = suppy)
  }
}

.get_eval_cotram_lower <- function(order_bsp, support) {
  function(y, orderbsp = order_bsp, suppy = support) {
    eval_bsp(log(1e-16 + y), order = orderbsp, supp = suppy)
  }
}

# Linear and log-linear bases ---------------------------------------------

eval_lin <- function(y) {
  ret <- cbind(1, y)
  if (NROW(ret) == 1)
    return(as.vector(ret))
  ret
}

eval_lin_prime <- function(y) {
  ret <- cbind(0, rep(1, length(y)))
  if (NROW(ret) == 1)
    return(as.vector(ret))
  ret
}

eval_loglin <- function(y) {
  stopifnot(y > 0)
  ret <- cbind(1, log(y))
  if (NROW(ret) == 1)
    return(as.vector(ret))
  ret
}

eval_loglin_prime <- function(y) {
  stopifnot(y > 0)
  ret <- cbind(0, 1/y)
  if (NROW(ret) == 1)
    return(as.vector(ret))
  ret
}

# Ordinal bases -----------------------------------------------------------

eval_ord <- function(y) {
  stopifnot(is.ordered(y))
  c(model.matrix(~ 0 + y, data = data.frame(y = y),
                 contrasts.arg = list("y" = "contr.treatment")))
}

.eval_ord_upr <- Vectorize(function(y) {
  ret <- eval_ord(y)
  # llev <- levels(y)[length(levels(y))]
  # if (y == llev) ret[length(ret)] <- 1
  ret
})

eval_ord_upr <- function(y) {
  ret <- t(.eval_ord_upr(y))
  if (nrow(ret) == 1L)
    return(as.vector(ret))
  ret
}

.eval_ord_lwr <- Vectorize(function(y) {
  resp <- eval_ord(y)
  ret <- c(resp[-1L], 0)
  # flev <- levels(y)[1L]
  # if (y == flev) ret[1L] <- 1
  ret
})

eval_ord_lwr <- function(y) {
  ret <- t(.eval_ord_lwr(y))
  if (nrow(ret) == 1L)
    return(as.vector(ret))
  ret
}

eval_ord_prime <- function(y) {
  resp <- eval_ord(y)
  resp[] <- 0
  return(resp)
}

# Bernstein bases ---------------------------------------------------------

eval_bsp <- function(y, order = 3, supp = range(y)) {

  # Evaluate a Bernstein Polynom (bsp) with a predefined order on a predefined
  # support that is used for scaling since a Beta distribution is defined on (0,1).
  # MLT Vignette p. 9
  #
  # y: numeric vector of length n
  # order: postive integer which is called M in the literature
  # supp: numeric vector of length 2 that is used to scale y down to (0,1)
  #
  # returns a numeric matrix (n x (order + 1))

  y <- (y - supp[1]) / diff(supp)
  res <- sapply(0:order, function(m) dbeta(y, m + 1, order + 1 - m) / (order + 1))
  if(is.null(dim(res)))
    res <- matrix(res, nrow=1)
  return(res)

}

eval_bsp_prime <- function(y, order = 3, supp = range(y)) {

  # Evaluate the first derivative of the bsp with a predefined order on a predefined
  # support that is used for scaling since a Beta distribution is defined on (0,1).
  # MLT Vignette p. 9. Note that "order" cancels out and that (1/diff(y_var$support)^deriv)
  # is multiplied afterwards. This is only due to numerical reasons to get the
  # exact same quantities as mlt::mlt. Furthermore, order/(order - 1 + 1) cancels
  # out in the mlt::mlt implementation which is not as stated on p. 9.
  #
  # y: numeric vector of length n
  # order: postive integer which is called M in the literature
  # supp: numeric vector of length 2 that is used to scale y down to (0,1)
  #
  # returns a numeric matrix (n x (order + 1))

  y <- (y - supp[1]) / diff(supp)
  res <- sapply(0:order, function(m) {

    first_t <- dbeta(y, m, order - m + 1)
    sec_t <- dbeta(y, m + 1, order - m)

    first_t[is.infinite(first_t)] <- 0L
    sec_t[is.infinite(sec_t)] <- 0L

    (first_t - sec_t) * (1 / diff(supp))
  })
  if(is.null(dim(res)))
    res <- matrix(res, nrow=1)
  return(res)
}

# eval_bsp_tf <- function(order, supp) {
#
#   return(
#     function(y){
#       y <- tf$math$divide(tf$math$subtract(y, supp[1L]), tf$math$subtract(supp[2L],supp[1L]))
#       return(
#         layer_concatenate(
#           lapply(0:order, function(m)
#             tf$reshape(
#               tf$divide(tfd_beta(m + 1, order + 1 - m)$prob(y), (order + 1)),
#               c(-1L,1L)
#             )
#           ), axis = 1L)
#       )
#     }
#   )
#
# }

#' @importFrom reticulate import_from_path
eval_bsp_tf <- function(order, supp, ...){
  python_path <- system.file("python", package = "deeptrafo")
  layer <- import_from_path("layers", path = python_path)
  return(layer$EvalBspTF(order = as.integer(order), supp = supp, ...))
}

ar_lags_layer <- function(order, supp)
{

  bsp_layer <- eval_bsp_tf(order = order, supp = supp)

  return(
    function(lags) lapply(1:lags$shape[[2]], function(i)
      bsp_layer(tf_stride_cols(lags, as.integer(i))))
  )

}

# tf_nan_to_zero <- function(x){
#
#   mult <- tf$where(tf$math$logical_not(tf$math$logical_or(tf$math$is_inf(x),
#                                                           tf$math$is_nan(x))),
#                    1.0, 0.0)
#   tf$math$multiply_no_nan(x, mult)
# }

# eval_bsp_prime_tf <- function(order, supp) {
#
#   return(
#     function(y){
#       y <- tf$math$divide(tf$math$subtract(y, supp[1L]), tf$math$subtract(supp[2L],supp[1L]))
#       return(
#         layer_concatenate(
#           lapply(0:order, function(m)
#           {
#
#             first_t <- tf$divide(tfd_beta(m, order - m + 1)$prob(y), order)
#             sec_t <- tf$divide(tfd_beta(m + 1, order - m)$prob(y), order)
#
#             return(
#               tf$reshape(tf$multiply(tf$subtract(tf_nan_to_zero(first_t),
#                                                  tf_nan_to_zero(sec_t)),
#                                      order), c(-1L,1L))
#             )
#           }
#           ), axis = 1L)
#       )
#     }
#   )
#
# }

eval_bsp_prime_tf <- function(order, supp, ...){
  python_path <- system.file("python", package = "deeptrafo")
  layer <- reticulate::import_from_path("layers", path = python_path)
  return(layer$EvalBspPrimeTF(order = as.integer(order), supp = supp, ...))
}

apply_atm_lags <- function(form)
{

  paste(paste0("atplag(", all.vars(form), ")"), collapse=" + ")

}

# TensorFlow repeat function which is not available for TF 2.0
tf_repeat <- function(a, dim)
  tf$reshape(tf$tile(tf$expand_dims(a, axis = -1L),  c(1L, 1L, dim)), shape = list(-1L, a$shape[[2]]*dim))

###############################################################################################
# for trafo with interacting features

mono_trafo_multi <- function(w, bsp_dim)
{

  w_res <- tf$reshape(w, shape = list(bsp_dim, as.integer(nrow(w)/bsp_dim)))
  w1 <- tf$slice(w_res, c(0L,0L), size=c(1L,ncol(w_res)))
  wrest <- tf$math$softplus(tf$slice(w_res, c(1L,0L), size=c(as.integer(nrow(w_res)-1),ncol(w_res))))
  w_w_cons <- tf$cumsum(k_concatenate(list(w1,wrest),
                                      axis = 1L # this is 1 and not 0 because k_concat is 1-based
  ), axis=0L)
  return(tf$reshape(w_w_cons, shape = list(nrow(w),1L)))

}

shift_scale_trafo_multi <- function(w, bsp_dim)
{

  w_res <- tf$reshape(w, shape = list(bsp_dim, as.integer(nrow(w)/bsp_dim)))
  w1 <- tf$slice(w_res, c(0L,0L), size=c(1L,ncol(w_res)))
  wrest <- tf$math$softplus(tf$slice(w_res, c(1L,0L), size=c(as.integer(nrow(w_res)-1),ncol(w_res))))
  w_w_cons <- k_concatenate(list(w1,wrest), axis = 1L)
  return(tf$reshape(w_w_cons, shape = list(nrow(w),1L)))

}

#' @importFrom R6 R6Class
MonoMultiLayer <- R6::R6Class("MonoMultiLayer",

                              inherit = KerasLayer,

                              public = list(

                                output_dim = NULL,

                                kernel = NULL,

                                dim_bsp = NULL,

                                kernel_regularizer = NULL,

                                trafo = NULL,

                                # input_dim = NULL,

                                initializer = NULL,

                                initialize = function(output_dim, dim_bsp,
                                                      # input_dim,
                                                      kernel_regularizer,
                                                      initializer = initializer_random_normal(),
                                                      trafo = trafo)
                                {
                                  self$output_dim <- output_dim
                                  self$dim_bsp <- dim_bsp
                                  # self$input_dim <- input_dim
                                  self$kernel_regularizer <- kernel_regularizer
                                  self$initializer <- initializer
                                  self$trafo <- trafo
                                },

                                build = function(input_shape) {
                                  self$kernel <- self$add_weight(
                                    name = 'kernel',
                                    shape = list(input_shape[[2]], self$output_dim),
                                    initializer = self$initializer,
                                    regularizer = self$kernel_regularizer,
                                    trainable = TRUE
                                  )
                                },

                                call = function(x, mask = NULL) {
                                  tf$matmul(x, self$trafo(self$kernel, self$dim_bsp))
                                },

                                compute_output_shape = function(input_shape) {
                                  list(NULL, self$output_dim)
                                }
                              )
)

# define layer wrapper function
layer_mono_multi <- function(object,
                             units = 1L,
                             dim_bsp = NULL,
                             name = "constraint_mono_layer_multi",
                             trainable = TRUE,
                             kernel_regularizer = NULL,
                             trafo = mono_trafo_multi
) {
  create_layer(MonoMultiLayer, object, list(
    name = name,
    trainable = trainable,
    output_dim = as.integer(units),
    dim_bsp = as.integer(dim_bsp),
    kernel_regularizer = kernel_regularizer,
    trafo = trafo
  ))
}

layer_combined_mono <- function(object,
                                units = 1L,
                                dim_bsp = NULL,
                                name = "constraint_mono_layer_multi",
                                trainable = TRUE,
                                kernel_regularizer = NULL,
                                trafo = mono_trafo_multi
                                )
{

  function(object){
    objects <- tf$split(object, num_or_size_splits = 2L, axis=1L)
    tf$keras$layers$concatenate(lapply(1:length(objects), function(i)
      layer_mono_multi(
        object = objects[[i]],
        units = units,
        dim_bsp = dim_bsp,
        name = paste0(name, "_", i),
        trainable = trainable,
        kernel_regularizer = kernel_regularizer,
        trafo = trafo
      )
    ))
  }

}


########## other version

MonoMultiTrafoLayer <- R6::R6Class("MonoMultiTrafoLayer",

                                   inherit = KerasLayer,

                                   public = list(

                                     output_dim = NULL,

                                     kernel = NULL,

                                     dim_bsp = NULL,

                                     initialize = function(output_dim, dim_bsp) {
                                       self$output_dim <- output_dim
                                       self$dim_bsp <- dim_bsp
                                     },

                                     build = function(input_shape) {
                                       self$kernel <- self$add_weight(
                                         name = 'kernel',
                                         shape = list(input_shape[[2]], self$output_dim),
                                         initializer = initializer_random_normal(),
                                         trainable = TRUE
                                       )
                                     },

                                     call = function(x, mask = NULL) {
                                       tf$multiply(x, tf$transpose(mono_trafo_multi(self$kernel, self$dim_bsp)))
                                     },

                                     compute_output_shape = function(input_shape) {
                                       list(input_shape[[1]], self$output_dim)
                                     }
                                   )
)

# define layer wrapper function
layer_mono_multi_trafo <- function(object,
                                   input_shape = NULL,
                                   output_dim = 1L,
                                   dim_bsp = NULL,
                                   name = "constraint_mono_layer_multi_trafo",
                                   trainable = TRUE
) {
  create_layer(MonoMultiTrafoLayer, object, list(
    name = name,
    trainable = trainable,
    input_shape = input_shape,
    output_dim = as.integer(output_dim),
    dim_bsp = as.integer(dim_bsp)
  ))
}

# to retrieve the weights on their original scale
softplus <- function(x) log(exp(x)+1)
reshape_softplus_cumsum <- function(x, order_bsp_p1)
{

  x <- matrix(x, nrow = order_bsp_p1, byrow = TRUE)
  x[2:nrow(x),] <- softplus(x[2:nrow(x),])
  apply(x, 2, cumsum)

}

correct_min_val <- function(pcf, addconst = 10)
{

  minval <- suppressWarnings(min(pcf$linterms[,sapply(pcf$linterms, is.numeric)], na.rm = TRUE))
  if(!is.null(pcf$smoothterms))
    minval <- min(c(minval,
                    suppressWarnings(sapply(pcf$smoothterms,
                                            function(x) min(x[[1]]$X)))))
  if (minval < 0)
  {

    minval <- minval - addconst

    if(!is.null(pcf$linterms))
      pcf$linterms[,sapply(pcf$linterms, is.numeric)] <-
        pcf$linterms[,sapply(pcf$linterms, is.numeric)] - minval
    if(!is.null(pcf$smoothterms))
      pcf$smoothterms <- lapply(pcf$smoothterms, function(x){
        x[[1]]$X <- x[[1]]$X - minval
        return(x)
      })


  }else{

    return(pcf)

  }

  if (minval == Inf)
    return(pcf)

  attr(pcf,"minval") <- minval

  return(pcf)

}

secondOrderPenBSP <- function(order_bsp, order_diff = 2)
{

  # taken from https://github.com/cran/penDvine/blob/master/R/pen.matrix.r

  if(order_diff == 0){

    k.dim <- order_bsp + 1
    k <- k.dim-1

    c2 <- factorial(k+1)/factorial(k-2)
    A <- matrix(0,k.dim-2,k.dim)
    diag(A) <- 1
    diag(A[,-1]) <- -2
    diag(A[,-c(1,2)]) <- 1

    A.hat <- matrix(NA,k.dim-2,k.dim-2)
    for(i in 0:(k-2)) {
      i.z <- i+1
      for(j in 0:(k-2)) {
        j.z <- j+1
        A.hat[i.z,j.z] <- choose(k-2,j)*choose(k-2,i)*beta(i+j+1,2*k-i-j-3)
      }
    }

    return(c2^2*(t(A)%*%A.hat%*%A))

  }

  K <- order_bsp+1

  if(order_diff==1) {
    L <- diag(1,K)
    L.1 <- diag(-1,K,K-1)
    L.2 <- matrix(0,K,1)
    L1 <- cbind(L.2,L.1)
    L <- L+L1
    L <- L[1:(K-1),]
  }
  if(order_diff==2) {
    L <- diag(1,K,K)
    L1 <- diag(-2,K,(K-1))
    L2 <- diag(1,K,(K-2))
    L.1 <- matrix(0,K,1)
    L1 <- cbind(L.1,L1)
    L2 <- cbind(L.1,L.1,L2)
    L <- L+L1+L2
    L <- L[1:(K-2),]
  }
  if(order_diff==3) {
    L <- diag(1,(K-3),(K-3))
    L.help <- matrix(0,(K-3),1)
    L1 <- diag(-3,(K-3),(K-3))
    M1 <- cbind(L,L.help,L.help,L.help)
    M2 <- cbind(L.help,L1,L.help,L.help)
    M3 <- cbind(L.help,L.help,-L1,L.help)
    M4 <- cbind(L.help,L.help,L.help,-L)
    L <- (M1+M2+M3+M4)
  }
  if(order_diff==4) {
    L <- diag(1,(K-4),(K-4))
    L.help <- matrix(0,(K-4),1)
    L1 <- diag(-4,(K-4),(K-4))
    L2 <- diag(6,(K-4),(K-4))
    M1 <- cbind(L,L.help,L.help,L.help,L.help)
    M2 <- cbind(L.help,L1,L.help,L.help,L.help)
    M3 <- cbind(L.help,L.help,L2,L.help,L.help)
    M4 <- cbind(L.help,L.help,L.help,L1,L.help)
    M5 <- cbind(L.help,L.help,L.help,L.help,L)
    L <- (M1+M2+M3+M4+M5)
  }

  return(crossprod(L))

}

calculate_log_score <- function(x, output)
{

  if(is.character(x$init_params$base_distribution) &
     x$init_params$base_distribution=="normal"){
    bd <- tfd_normal(loc = 0, scale = 1)
  }else if((is.character(x$init_params$base_distribution) &
            x$init_params$base_distribution=="logistic")){
    bd <- tfd_logistic(loc = 0, scale = 1)
  }else{
    bd <- x$init_params$base_distribution
  }
  return(
    as.matrix(bd %>% tfd_log_prob(output[,2,drop=F] + output[,1,drop=F])) +
      as.matrix(log(tf$clip_by_value(output[,3,drop=F], 1e-8, Inf)))
  )

}

get_order_bsp_p1 <- function(x)
{

  x$init_params$trafo_options$order_bsp + 1L

}

split_interaction_terms <- function(term){

  splt <- trimws(strsplit(term, "%I%")[[1]])
  splt[1] <- paste0(splt[1], ")")
  splt[2] <- substr(splt[2], 1, nchar(splt[2])-1)
  return(splt)

}

row_tensor <- function(A,B)
{

  kronecker(A, matrix(1, ncol = ncol(B))) *
    kronecker(matrix(1, ncol = ncol(A)), B)

}

row_tensor_by_basis <- function(X, basis_dim){

  row_tensor(X[,1:basis_dim],X[,(basis_dim+1):ncol(X)])

}

h1_plotfun <- function(dim_basis){

  return(
    function(pp, weights, grid_length = 40){

      org_values <- pp$get_org_values()

      BX <- row_tensor_by_basis(pp$data_trafo(), dim_basis)

      plotData <-
        list(org_feature_name = pp$term,
             value = do.call("cbind", org_values),
             design_mat = BX,
             coef = weights)

      this_x <- do.call(seq, c(as.list(range(plotData$value[,1])),
                               list(l=grid_length)))
      this_y <- do.call(seq, c(as.list(range(plotData$value[,2])),
                                 list(l=grid_length)))
      df <- as.data.frame(expand.grid(this_x, this_y))
      colnames(df) <- extractvar(pp$term)
      plotData$df <- df
      plotData$x <- this_x
      plotData$y <- this_y
      plotData$partial_effect <- row_tensor_by_basis(
        pp$predict_trafo(newdata = df), dim_basis) %*% weights

      return(plotData)

    }

  )

}

get_theta <- function(object)
{

  do.call("cbind", coef(object, which_param = "interacting"))

}

# Additional distributions ------------------------------------------------

tfd_gompertz <- function(loc, scale,
                         validate_args = FALSE,
                         allow_nan_stats = TRUE,
                         name = "Gompertz")
{

  args <- list(
    loc = loc,
    scale = scale,
    validate_args = validate_args,
    allow_nan_stats = allow_nan_stats,
    name = name
  )

  python_path <- system.file("python", package = "deeptrafo")
  distributions <- reticulate::import_from_path("distributions", path = python_path)

  return(do.call(distributions$gompertz$Gompertz, args))

}

#' @importFrom grDevices rgb
#' @importFrom graphics matplot
#' @importFrom stats fitted formula predict rmultinom

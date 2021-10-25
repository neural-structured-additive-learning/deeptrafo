#' Fitting Deep Conditional Transformation Models
#'
#' @param formula Formula specifying the outcome, shift, interaction and shared
#'     terms as \code{response ~ shifting | interacting | shared}
#' @param list_of_formulas list of formulas where the first element corresponds to
#'     a transformation h_1 as specified in DCTMs, the second element to h_2 as specified
#'     in DCTMs and a third element for shared layers used in both
#' @param lag_formula formula representing the optional predictor for lags of the
#'     response as defined for ATMs
#' @param data list or data.frame with data; must include y as column
#' @param order_bsp integer; order of Bernstein polynomials
#' @param addconst_interaction positive constant;
#'     a constant added to the additive predictor of the interaction term.
#'     If \code{NULL}, terms are left unchanged. If 0 and predictors have negative values in their
#'     design matrix, the minimum value of all predictors is added to ensure positivity.
#'     If > 0, the minimum value plus the \code{addconst_interaction} is added to each predictor
#'     in the interaction term.
#' @param split_between_h1_h2 integer; for shared layer use
#'     defines how many of the last layer's hidden units are added to the h2
#'     term (while the rest is used for the h2 term)
#' @param family, tfd_distribution or string; the base distribution for
#'     transformation models. If string, can be \code{"normal"} or \code{"logistic"}.
#' @param trafo_options options for transformation models such as the basis function used
#' @param ... Arguments passed to \code{deepregression}
#'
#' @return An object of class \code{c("deeptrafo", "deepregression")}
#'
#' @examples
#' dat <- data.frame(y = rnorm(100), x = rnorm(100), z = rnorm(100))
#' fml <- y ~ x | z
#' m <- deeptrafo(fml, dat)
#' m %>% fit(epochs = 10)
#'
#' @export
#' @details
#'
deeptrafo <- function(
  formula,
  data,
  lag_formula = NULL,
  order_bsp = 10L,
  addconst_interaction = NULL,
  split_between_h1_h2 = 1L,
  atm_toplayer = function(x) layer_dense(x, units = 1L),
  family = "normal",
  monitor_metrics = crps_stdnorm_metric,
  trafo_options = trafo_control(order_bsp = order_bsp),
  ...
  )
{
  # How many terms are in the formula
  formula <- as.Formula(formula)
  nterms <- length(attr(formula, "rhs"))

  # Name of the response variable
  rvar <- all.vars(formula)[1]

  # List of formulas
  list_of_formulas <- list(
    ybasis = as.formula(paste0("~ -1 + bsfun(", rvar, ")")),
    ybasisprime = as.formula(paste0("~ -1 + bsprimefun(", rvar, ")")),
    h1 = structure(formula(formula, lhs = 0, rhs = 1L), with_layer = FALSE),
    h2 = if (nterms >= 2L) formula(formula, lhs = 0, rhs = 2L) else NULL,
    shared = if (nterms == 3L) formula(formula, lhs = 0, rhs = 3L) else NULL
  )

  # Remove NULL formulae
  list_of_formulas[sapply(list_of_formulas, is.null)] <- NULL

  # Extract response variable
  y <- model.response(model.frame(formula(formula, lhs = 1, rhs = 0), data = data))

  # check for ATMs
  if(!is.null(lag_formula)){

    # extract from lag formula the variables as simple sum and
    # layers for additional transformation
    lag_formula <- apply_atm_lags(lag_formula)
    list_of_formulas <- c(list_of_formulas,
                          list(atmlags = lag_formula))

  }

  # define how to get a trafo model from predictor
  from_pred_to_trafo_fun <- from_preds_to_trafo(atm_toplayer = atm_toplayer,
                                                split = split_between_h1_h2,
                                                const_ia = addconst_interaction)

  # define ar_layer and atm processor
  ar_layer <- ar_lags_layer(order = order_bsp, supp = range(y))

  atm_lag_processor <- make_atm_processor(ar_layer)

  trafo_processor <- list(bsfun = basis_processor,
                          bsprimefun = basis_prime_processor,
                          atmlag = atm_lag_processor)

  dots <- list(...)

  if (is.null(dots$additional_processor)) {

    additional_processor <- trafo_processor

  } else{

    additional_processor <- c(list(...)$additional_processor, trafo_processor)
    dots$additional_processor <- NULL

  }

  attr(additional_processor, "controls") <- trafo_options

  snwb <- list(subnetwork_init)[rep(1, length(list_of_formulas))]
  snwb[[which(names(list_of_formulas) == "h1")]] <- interaction_init

  ret <- do.call("deepregression",
                 c(list(y = y,
                        family = family,
                        data = data,
                        list_of_formulas = list_of_formulas,
                        subnetwork_builder = snwb,
                        from_preds_to_output = from_pred_to_trafo_fun,
                        loss = neg_ll_trafo(family),
                        monitor_metrics = monitor_metrics,
                        additional_processor = additional_processor),
                   dots)
  )

  class(ret) <- c("deeptrafo", "deepregression")
  return(ret)

}

#' Initializes the Processed Additive Predictor for TM's Interaction
#'
#' @param pp processed predictor list from \code{processor}
#' @param deep_top keras layer if the top part of the deep network after orthogonalization
#' is different to the one extracted from the provided network
#' @param orthog_fun function used for orthogonalization
#' @param split_fun function to split the network to extract head
#' @param param_nr integer number for the distribution parameter
#' @return returns a list of input and output for this additive predictor
#'
interaction_init <- function(pp, deep_top = NULL,
                            orthog_fun = orthog_tf,
                            split_fun = split_model,
                            param_nr = 2)
{


  inputs <- makeInputs(pp, param_nr = param_nr)

  if(all(sapply(pp, function(x) is.null(x$right_from_oz)))){
    # if there is no term to orthogonalize

    outputs <- lapply(1:length(pp), function(i) pp[[i]]$layer(inputs[[i]]))
    outputs <- layer_concatenate_identity(outputs)
    return(list(inputs, outputs))

  }else{

    # define the different types of elements
    outputs_w_oz <- unique(unlist(sapply(pp, "[[", "right_from_oz")))
    outputs_used_for_oz <- which(sapply(pp, function(x) !is.null(x$right_from_oz)))
    outputs_onlyfor_oz <- outputs_used_for_oz[!sapply(pp[outputs_used_for_oz],
                                                      "[[", "left_from_oz")]
    outputs_wo_oz <- setdiff(1:length(pp), c(outputs_w_oz, outputs_onlyfor_oz))

    outputs <- list()
    if(length(outputs_wo_oz)>0) outputs <- layer_concatenate_identity(
      lapply((1:length(pp))[outputs_wo_oz],
             function(i) pp[[i]]$layer(inputs[[i]]))
      )
    ox_outputs <- list()
    k <- 1

    for(i in outputs_w_oz){

      inputs_for_oz <- which(sapply(pp, function(ap) i %in% ap$right_from_oz))
      ox <- layer_concatenate_identity(inputs[inputs_for_oz])
      # there is no splitting in the interaction term
      # as it is always combined with the basis for Y
      deep <- pp[[i]]$layer(inputs[[i]])
      ox_outputs[[k]] <- orthog_fun(deep, ox)
      k <- k + 1

    }

    if(length(ox_outputs)>0) outputs <- layer_concatenate_identity(c(outputs, ox_outputs))

    return(list(inputs, outputs))

  }


}


#' @title Define Predictor of Transformation Model
#'
#'
#' @param atm_toplayer function to be applied on top of the transformed lags
#' @param split see \code{split_between_h1_h2} from \code{?deeptrafo}
#' @param const_ia see \code{addconst_interaction} from \code{?deeptrafo}
#' \code{deepregression}
#' @return a function of list_pred_param returning a list of output tensors
#' that is passed to \code{model_fun} of \code{deepregression}
#'
#' @export
from_preds_to_trafo <- function(
  atm_toplayer = function(x) layer_dense(x, units = 1L),
  split = 1L,
  const_ia = NULL
)
{

  return(function(list_pred_param, ...){

    # make inputs more readable
    input_theta_y <- list_pred_param$ybasis
    input_theta_y_prime <- list_pred_param$ybasisprime
    interact_pred <- list_pred_param$h1
    if(!is.null(const_ia))
      interact_pred <- tf$add(interact_pred,
                              tf$constant(const_ia,
                                          dtype="float32"))
    shift_pred <- list_pred_param$h2

    # check if ATM or DCTM
    is_atm <- !is.null(list_pred_param$atmlags)
    if(is_atm) input_theta_atm <- list_pred_param$atmlags

    # check if shared layer
    has_shared <- !is.null(list_pred_param$shared)
    if(has_shared){
      input_shared <- list_pred_param$shared
      shared_dim <- input_shared$shape[[2]]

      # extract parts (use all but the last column for h1)
      h1part <- tf_stride_cols(input_shared, 1L, shared_dim-split)
      h2part <- tf_stride_cols(input_shared, shared_dim-split+1L, shared_dim)

      # concat
      interact_pred <- layer_concatenate(list(interact_pred, h1part))
      shift_pred <- layer_concatenate(list(shift_pred, h2part))

    }

    # define shapes
    order_bsp_p1 <- input_theta_y$shape[[2]]
    h1_col <- interact_pred$shape[[2]]
    total_h1_dim <- order_bsp_p1 * h1_col

    # define monotone layer for trafo function
    thetas_layer <- layer_mono_multi(
      input_shape = list(NULL, total_h1_dim),
      dim_bsp = c(order_bsp_p1)
    )

    ## define RWTs
    AoB <- tf_row_tensor(input_theta_y, interact_pred)
    AprimeoB <- tf_row_tensor(input_theta_y_prime, interact_pred)

    # define h1 and h1'
    aTtheta <- AoB %>% thetas_layer()
    aPrimeTtheta <- AprimeoB %>% thetas_layer()

    # check if ATM and add to shift_pred
    if(is_atm){

      ## create interaction

      # combine every lag with the interacting predictor
      AoB_lags <- lapply(input_theta_atm, function(inp)
        tf_row_tensor(inp, interact_pred))

      # multiply with theta weights
      aTtheta_lags <- layer_concatenate_identity(
        lapply(AoB_lags, function(aob) aob %>% thetas_layer())
      )

      # combine all transformed lags
      lag_pred <- aTtheta_lags %>% atm_toplayer()

      # overwrite the shift_pred by adding lags
      shift_pred <- layer_add(list(shift_pred, lag_pred))

    }

    # return transformation
    trafo <- layer_concatenate(list(
      shift_pred,
      aTtheta,
      aPrimeTtheta
    ))

    return(trafo)

  })

}

#' Generic negative log-likelihood of a transformation model
#'
#' @param base_distribution base or error distribution
#'
#' @return a function for the negative log-likelihood with outcome \code{y}
#' and transformation model \code{model}. The transformation model is represented
#' by a list of two, with first element a list of model outputs
#' that are summed up and evaluated with the log-probability of the
#' \code{basis_dist}, and second element a single-column tensor
#' representing the determinant of the Jacobian and transformed
#' using the log
#'
#' @import deepregression
#' @export
#'
#'
neg_ll_trafo <- function(base_distribution) {

  # evaluate base_distribution
  if(is.null(base_distribution) ||
     (is.character(base_distribution) & base_distribution=="normal")){
    bd <- tfd_normal(loc = 0, scale = 1)
  }else if((is.character(base_distribution) &
            base_distribution=="logistic")){
    bd <- tfd_logistic(loc = 0, scale = 1)
  }else{
    bd <- base_distribution
  }

  return(
    function(y, model){

      first_term <- bd %>% tfd_log_prob(layer_add(list(tf_stride_cols(model,1L),
                                                       tf_stride_cols(model,2L))))
      sec_term <- tf$math$log(tf$clip_by_value(tf_stride_cols(model,3L), 1e-8, Inf))
      neglogLik <- -1 * tf$add(first_term, sec_term)
      return(neglogLik)
    }
  )

}

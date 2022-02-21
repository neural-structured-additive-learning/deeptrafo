#' Deep Conditional Transformation Models
#'
#' @param formula Formula specifying the response, shift, interaction and shared
#'     terms as \code{response | interacting ~ shifting | shared}
#' @param lag_formula formula representing the optional predictor for lags of the
#'     response as defined for auto-regressive transformation models (ATMs)
#' @param data \code{list} or \code{data.frame} containing both structured and
#'     unstructured data
#' @param response_type character; type of response. One of \code{"continuous"},
#'     \code{"survival"}, \code{"count"}, \code{"ordered"}
#' @param order integer; order of Bernstein polynomials or number of levels for
#'     ordinal responses
#' @param addconst_interaction positive constant;
#'     a constant added to the additive predictor of the interaction term.
#'     If \code{NULL}, terms are left unchanged. If 0 and predictors have
#'     negative values in their design matrix, the minimum value of all predictors
#'     is added to ensure positivity. If > 0, the minimum value plus the
#'     \code{addconst_interaction} is added to each predictor in the interaction
#'     term.
#' @param family \code{tfd_distribution} or character; the base distribution for
#'     transformation models. If character, can be \code{"normal"}, \code{"logistic"}
#'     or \code{"gumbel"}.
#' @param trafo_options options for transformation models such as the basis
#'     function used
#' @param ... Arguments passed to \code{deepregression}
#'
#' @return An object of class \code{c("deeptrafo", "deepregression")}
#'
#' @examples
#' data("wine", package = "ordinal")
#' wine$noise <- rnorm(nrow(wine))
#' fml <- rating ~ 0 + temp
#' m <- deeptrafo(fml, wine, family = "logistic", monitor_metric = NULL)
#' m %>% fit(epochs = 100, batch_size = nrow(wine))
#' predfun <- m %>% predict(wine)
#' predfun(wine$rating, type = "pdf")
#' plot(m)
#' coef(m, which_param = "interacting")
#' coef(m, which_param = "shifting")
#'
#' @importFrom mlt R
#' @export
#' @details
#'
deeptrafo <- function(
  formula,
  data,
  lag_formula = NULL,
  response_type = get_response_type(data[[all.vars(fml)[1]]]),
  order = get_order(response_type, data[[all.vars(fml)[1]]]),
  addconst_interaction = 0,
  family = "logistic",
  monitor_metrics = NULL,
  trafo_options = trafo_control(order_bsp = order,
                                response_type = response_type),
  ...
)
{
  # How many terms are in the formula
  fml <- as.Formula(formula)
  ninteracting <- length(attr(fml, "lhs"))
  nterms <- length(attr(fml, "rhs"))

  # Name of the response variable
  rvar <- all.vars(formula)[1]

  # Placeholder Intercept
  int <- 1

  # Set up formulas for basis
  if (ninteracting > 1L) {
    interacting <- formula(fml, lhs = 2L, rhs = 0L)[[2]]
    h1_form <- paste0(
      "~ -1 + ", paste(paste0("ia(", c(int, trimws(
        strsplit(form2text(interacting), "+", fixed = TRUE)[[1]]
      )), ")"), collapse=" + ")
    )
  } else {
    h1_form <- paste0(
      "~ -1 + ", paste(paste0("ia(", int, ")"), collapse=" + ")
    )
  }

  # List of formulas
  list_of_formulas <- list(
    yterms = as.formula(paste0("~ -1 + bsfun(", rvar, ") + bspfun(", rvar, ")")),
    h1pred = as.formula(h1_form),
    h2 = if (nterms >= 1L) formula(fml, lhs = 0, rhs = 1L) else NULL,
    shared = if (nterms == 2L) formula(fml, lhs = 0, rhs = 2L) else NULL
  )

  # Remove NULL formulae
  list_of_formulas[sapply(list_of_formulas, is.null)] <- NULL

  # Extract response variable
  resp <- model.response(model.frame(formula(fml, lhs = 1, rhs = 0), data = data))
  y <- response(resp)

  # check for ATMs
  if(!is.null(lag_formula)){

    # extract from lag formula the variables as simple sum and
    # layers for additional transformation
    lag_formula <- apply_atm_lags(lag_formula)
    list_of_formulas$yterms <- as.formula(paste0(form2text(list_of_formulas$yterms),
                                                 " + ", lag_formula))

  }

  # define how to get a trafo model from predictor
  from_pred_to_trafo_fun <- from_preds_to_trafo(atm_toplayer = trafo_options$atm_toplayer,
                                                const_ia = addconst_interaction)

  atm_lag_processor <- atm_lag_processor_factory(rvar)

  trafo_processor <- list(bsfun = basis_processor,
                          bspfun = basisprime_processor,
                          ia = ia_processor,
                          atmlag = atm_lag_processor)

  dots <- list(...)

  if (is.null(dots$additional_processor)) {

    additional_processor <- trafo_processor

  } else{

    additional_processor <- c(list(...)$additional_processor, trafo_processor)
    dots$additional_processor <- NULL

  }

  attr(additional_processor, "controls") <- trafo_options

  # Loss function
  tloss <- get_loss(response_type, family)

  snwb <- list(subnetwork_init)[rep(1, length(list_of_formulas))]
  snwb[[which(names(list_of_formulas) == "h1pred")]] <-
    h1_init(yterms = which(names(list_of_formulas) == "yterms"),
            h1pred = which(names(list_of_formulas) == "h1pred"),
            add_const_positiv = addconst_interaction)
  snwb[[which(names(list_of_formulas) == "yterms")]] <- function(...) return(NULL)

  ret <- do.call("deepregression",
                 c(list(y = y,
                        family = family,
                        data = data,
                        list_of_formulas = list_of_formulas,
                        subnetwork_builder = snwb,
                        from_preds_to_output = from_pred_to_trafo_fun,
                        loss = tloss,
                        monitor_metrics = monitor_metrics,
                        additional_processor = additional_processor),
                   dots)
  )
  ret$init_params$trafo_options <- trafo_options
  ret$init_params$response_varname <- rvar
  ret$init_params$response_type <- response_type
  ret$init_params$response <- resp

  class(ret) <- c("deeptrafo", "deepregression")
  return(ret)

}

#' Initializes the Processed Additive Predictor for TM's Interaction
#'
#' @param yterms,h1pred positions of the left and right RWT term
#' @return returns a subnetwork_init function with pre-defined arguments
#'
h1_init <- function(yterms, h1pred, add_const_positiv = 0)
{
  return(
    function(pp, deep_top, orthog_fun, split_fun, shared_layers,
             param_nr, selectfun_in, selectfun_lay, gaminputs,
             summary_layer = NULL){

      # instead of passing the respective pp,
      # subsetting is done within subnetwork_init
      # to allow other subnetwork_builder to
      # potentially access all pp entries
      pp_in <- pp_lay <- pp[[h1pred]]
      pp_y <- pp[[yterms]]


      # generate pp parts
      gaminput_nrs <- sapply(pp_in, "[[", "gamdata_nr")
      has_gaminp <- !sapply(gaminput_nrs,is.null)
      gaminput_comb <- sapply(pp_in[which(has_gaminp)], "[[", "gamdata_combined")
      inputs <- makeInputs(pp_in, param_nr = param_nr)
      inputs_y <- makeInputs(pp_y, param_nr = 1)
      org_inputs_for_concat <- list()

      if(sum(has_gaminp)){

        for(i in 1:sum(has_gaminp)){

          # concatenate inputs or replace?
          concat <- gaminput_comb[[i]]
          nr <- which(has_gaminp)[i]

          if(!is.null(concat) && concat){

            org_inputs_for_concat <- c(
              org_inputs_for_concat,
              inputs[[nr]]
            )
            inputs[[nr]] <- layer_concatenate_identity(
              list(gaminputs[[gaminput_nrs[[nr]]]], inputs[[nr]])
            )

          }else{

            inputs[[nr]] <- gaminputs[[gaminput_nrs[[nr]]]]

          }

        }

        inputs_to_replace <- which(has_gaminp)[gaminput_comb]
        keep_inputs_in_return <- setdiff(1:length(inputs), (which(has_gaminp)[!gaminput_comb]))

      }else{

        inputs_to_replace <-  c()
        keep_inputs_in_return <- 1:length(inputs)

      }

      layer_matching <- 1:length(pp_in)
      names(layer_matching) <- layer_matching

      if(!is.null(shared_layers))
      {

        names_terms <- get_names_pfc(pp_in)

        for(group in shared_layers){

          layer_ref_nr <- which(names_terms==group[1])
          layer_opts <- get("layer_args", environment(pp_lay[[layer_ref_nr]]$layer))
          layer_opts$name <- paste0("shared_",
                                    makelayername(paste(group, collapse="_"),
                                                  param_nr))
          layer_ref <- do.call(get("layer_class", environment(pp_lay[[layer_ref_nr]]$layer)),
                               layer_opts)

          terms_replace_layer <- which(names_terms%in%group)
          layer_matching[terms_replace_layer] <- layer_ref_nr
          for(i in terms_replace_layer) pp_lay[[i]]$layer <- layer_ref

        }
      }

      if(all(sapply(pp_in, function(x) is.null(x$right_from_oz)))){ # if there is no term to orthogonalize

        outputs <- lapply(1:length(pp_y), function(j) layer_add_identity(
          lapply(1:length(pp_in), function(i) pp_lay[[layer_matching[i]]]$layer(
              pp_y[[j]]$layer(inputs_y[[j]]),
              tf$add(inputs[[i]], add_const_positiv)
          )
          )
        ))
        outputs <- layer_concatenate_identity(outputs)

        # replace original inputs
        if(length(org_inputs_for_concat)>0)
          inputs[inputs_to_replace] <- org_inputs_for_concat
        return(list(c(inputs_y,inputs[keep_inputs_in_return]), outputs))

      }else{

        ## define the different types of elements
        stop("Not implemented yet.")
        # outputs_w_oz <- unique(unlist(sapply(pp_in, "[[", "right_from_oz")))
        # outputs_used_for_oz <- which(sapply(pp_in, function(x) !is.null(x$right_from_oz)))
        # outputs_onlyfor_oz <- outputs_used_for_oz[!sapply(pp_in[outputs_used_for_oz], "[[", "left_from_oz")]
        # outputs_wo_oz <- setdiff(1:length(pp_in), c(outputs_w_oz, outputs_onlyfor_oz))
        #
        # outputs <- list()
        # if(length(outputs_wo_oz)>0) outputs <-
        #   layer_add_identity(lapply((1:length(pp_in))[outputs_wo_oz],
        #                             function(i) pp_lay[[layer_matching[i]]]$layer(inputs[[i]])))
        # ox_outputs <- list()
        # k <- 1
        #
        # for(i in outputs_w_oz){
        #
        #   inputs_for_oz <- which(sapply(pp_in, function(ap) i %in% ap$right_from_oz))
        #   ox <- layer_concatenate_identity(inputs[inputs_for_oz])
        #   if(is.null(deep_top)){
        #     deep_splitted <- split_fun(pp_lay[[layer_matching[i]]]$layer)
        #   }else{
        #     deep_splitted <- list(pp_lay[[layer_matching[i]]]$layer, deep_top)
        #   }
        #
        #   deep <- deep_splitted[[1]](inputs[[i]])
        #   ox_outputs[[k]] <- deep_splitted[[2]](orthog_fun(deep, ox))
        #   k <- k + 1
        #
        # }
        #
        # if(length(ox_outputs)>0) outputs <- summary_layer(c(outputs, ox_outputs))
        #
        # if(length(org_inputs_for_concat)>0)
        #   inputs[inputs_to_replace] <- org_inputs_for_concat
        # return(list(inputs[keep_inputs_in_return], outputs))

      }


    }
  )

}

#' Initializes the Processed Additive Predictor for ATMs
#'
#' @param atmnr,h1nr positions of the atm and h1 formula
#' @return returns a subnetwork_init function with pre-defined arguments
#'
atm_init <- function(atmnr, h1nr)
{
  return(
    function(pp, deep_top, orthog_fun, split_fun, shared_layers, param_nr,
             gaminputs)
      subnetwork_init(pp, deep_top, orthog_fun, split_fun, shared_layers, param_nr,
                      # pp_input_subset = atmnr,
                      # pp_layer_subset = h1nr,
                      gaminputs = gaminputs,
                      summary_layer = layer_concatenate_identity)
  )

}

#' @title Define Predictor of Transformation Model
#'
#'
#' @param atm_toplayer function to be applied on top of the transformed lags
#' @param const_ia see \code{addconst_interaction} from \code{?deeptrafo}
#' \code{deepregression}
#' @return a function of list_pred_param returning a list of output tensors
#' that is passed to \code{model_fun} of \code{deepregression}
#'
#' @export
from_preds_to_trafo <- function(
  atm_toplayer = function(x) layer_dense(x, units = 1L),
  const_ia = NULL
)
{

  return(function(list_pred_param, ...){

    # make inputs more readable
    # aTtheta <- tf_stride_cols(list_pred_param$h1pred, 1L)
    # aPrimeTtheta <- tf_stride_cols(list_pred_param$h1pred, 2L)
    h1pred_ncol <- list_pred_param$h1pred$shape[[2]]
    shift_pred <- list_pred_param$h2
    if(h1pred_ncol > 2){

      lag_pred <- tf_stride_cols(list_pred_param$h1pred, 3, h1pred_ncol) %>% atm_toplayer()

      # overwrite the shift_pred by adding lags
      shift_pred <- layer_add(list(shift_pred, lag_pred))

    }

    # return transformation
    trafo <- layer_concatenate(list(
      shift_pred,
      tf_stride_cols(list_pred_param$h1pred, 1L, 2L)
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

      first_term <- bd %>% tfd_log_prob(layer_add(list(tf_stride_cols(model, 1L),
                                                       tf_stride_cols(model, 2L))))
      sec_term <- tf$math$log(tf$clip_by_value(tf_stride_cols(model, 3L), 1e-8, Inf))
      neglogLik <- -1 * tf$add(first_term, sec_term)
      return(neglogLik)
    }
  )

}

#' negative log-likelihood of an ordinal transformation model
#'
#' @param base_distribution base or error distribution
#'
#' @return a function for the negative log-likelihood with outcome \code{y_true}
#' and transformation model \code{y_pred}. The transformation model is represented
#' by a list of two, with first element a list of model outputs
#' that are summed up and evaluated with the log-probability of the
#' \code{base_dist}, and second element a single-column tensor
#' representing the determinant of the Jacobian and transformed
#' using the log
#'
#' @import deepregression
#' @export
#'
#'
nll_ordinal <- function(base_distribution = "logistic") {

  if (is.character(base_distribution)) {

    bd <- get_bd(base_distribution)

  } else {

    bd <- base_distribution

  }

  return(
    function(y_true, y_pred){
      lwr <- layer_add(list(tf_stride_cols(y_pred, 3L),
                            tf_stride_cols(y_pred, 1L)))
      upr <- layer_add(list(tf_stride_cols(y_pred, 2L),
                            tf_stride_cols(y_pred, 1L)))
      t1 <- tf_stride_cols(y_true, 1L)
      t2 <- tf_stride_cols(y_true, 3L)
      lik <- t1 * tfd_cdf(bd, upr) + t2 * (1 - tfd_cdf(bd, lwr)) +
        (1 - t1) * (1 - t2) * (tfd_cdf(bd, upr) - tfd_cdf(bd, lwr))
      neglogLik <- - tf$math$log(tf$clip_by_value(lik, 1e-8, Inf))
      return(neglogLik)
    }
  )

}

#' negative log-likelihood for count transformation models
#'
#' @param base_distribution base or error distribution
#'
#' @return a function for the negative log-likelihood with outcome \code{y_true}
#' and transformation model \code{y_pred}. The transformation model is represented
#' by a list of two, with first element a list of model outputs
#' that are summed up and evaluated with the log-probability of the
#' \code{base_dist}, and second element a single-column tensor
#' representing the determinant of the Jacobian and transformed
#' using the log
#'
#' @import deepregression
#' @export
#'
#'
nll_count <- function(base_distribution = "logistic") {

  if (is.character(base_distribution)) {

    bd <- get_bd(base_distribution)

  } else {

    bd <- base_distribution

  }

  return(
    function(y_true, y_pred){
      lwr <- layer_add(list(tf_stride_cols(y_pred, 3L),
                            tf_stride_cols(y_pred, 1L)))
      upr <- layer_add(list(tf_stride_cols(y_pred, 2L),
                            tf_stride_cols(y_pred, 1L)))
      cleft <- tf_stride_cols(y_true, 1L)
      cint <- tf_stride_cols(y_true, 4L)
      lik <- cleft * tfd_cdf(bd, upr) +
        cint * (tfd_cdf(bd, upr) - tfd_cdf(bd, lwr))
      neglogLik <- - tf$math$log(tf$clip_by_value(lik, 1e-8, Inf))
      return(neglogLik)
    }
  )

}

#' negative log-likelihood for potentially right-censored survival responses
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
nll_surv <- function(base_distribution) {

  if (is.character(base_distribution)) {

    bd <- get_bd(base_distribution)

  } else {
    bd <- base_distribution
  }

  return(
    function(y_true, y_pred){

      exact <- tf_stride_cols(y_true, 2L)
      cright <- tf_stride_cols(y_true, 3L)

      trafo <- layer_add(list(tf_stride_cols(y_pred, 1L),
                              tf_stride_cols(y_pred, 2L)))

      trafo_prime <- tf$math$log(tf$clip_by_value(tf_stride_cols(y_pred, 3L),
                                                  1e-8, Inf))

      ll_exact <- tfd_log_prob(bd, trafo) + trafo_prime
      ll_right <- tfd_log_survival_function(bd, trafo)

      neglogLik <- - (exact * ll_exact + cright * ll_right)
      return(neglogLik)
    }
  )

}

#' generic negative log-likelihood for transformation models
#'
#' @param base_distribution base distribution / error distribution / inverse link
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
nll <- function(base_distribution) {

  if (is.character(base_distribution)) {
    bd <- get_bd(base_distribution)
  } else {
    bd <- base_distribution
  }

  return(
    function(y_true, y_pred){

      cleft <- tf_stride_cols(y_true, 1L)
      exact <- tf_stride_cols(y_true, 2L)
      cright <- tf_stride_cols(y_true, 3L)
      cint <- tf_stride_cols(y_true, 4L)

      # <FIXME>
      #   Generic NLL needs another basis eval for y_lower.
      #</FIXME>
      trafo <- layer_add(list(tf_stride_cols(y_pred, 1L),
                              tf_stride_cols(y_pred, 2L)))
      trafo_lwr <- layer_add(list(tf_stride_cols(y_pred, 1L),
                                  tf_stride_cols(y_pred, 3L)))
      trafo_prime <- tf$math$log(tf$clip_by_value(tf_stride_cols(y_pred, 3L),
                                                  1e-8, Inf))

      ll_exact <- tfd_log_prob(bd, trafo) + trafo_prime
      ll_right <- tf$math$log(tfd_cdf(bd, trafo))
      ll_left <- tf$math$log(1 - tfd_cdf(bd, trafo))
      ll_int <- tf$math$log(tfd_cdf(bd, trafo) - tfd_cdf(bd, trafo_lwr))

      neglogLik <- - (cleft * ll_left + exact * ll_exact + cright * ll_right +
                        cint * ll_int)

      return(neglogLik)
    }
  )

}

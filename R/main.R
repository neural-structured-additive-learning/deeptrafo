#' Deep Conditional Transformation Models
#'
#' @param formula Formula specifying the response, interaction, shift terms
#'     as \code{response | interacting ~ shifting}.
#'     auto-regressive transformation models (ATMs).
#' @param data Named \code{list} or \code{data.frame} which may contain both
#'     structured and unstructured data.
#' @param response_type Character; type of response. One of \code{"continuous"},
#'     \code{"survival"}, \code{"count"}, or \code{"ordered"}. If not supplied
#'     manually it is determined by the first entry in \code{data[[response]]}.
#' @param order Integer; order of the response basis. Default 10 for Bernstein
#'     basis or number of levels minus one for ordinal responses.
#' @param addconst_interaction Positive constant;
#'     a constant added to the additive predictor of the interaction term.
#'     If \code{NULL}, terms are left unchanged. If 0 and predictors have
#'     negative values in their design matrix, the minimum value of all predictors
#'     is added to ensure positivity. If > 0, the minimum value plus the
#'     \code{addconst_interaction} is added to each predictor in the interaction
#'     term. This ensures a monotone non-decreasing transformation function in
#'     the response when using (tensor product) spline bases in the interacting
#'     term.
#' @param latent_distr A \code{tfd_distribution} or character; the base distribution for
#'     transformation models. If character, can be \code{"normal"}, \code{"logistic"},
#'     \code{"gumbel"} or \code{"gompertz"}.
#' @param trafo_options Options for transformation models such as the basis
#'     function used, see \code{\link[deeptrafo]{trafo_control}} for more details.
#' @param ... Additional arguments passed to \code{deepregression}
#' @param monitor_metrics See \code{\link[deepregression]{deepregression}}
#' @param return_data Include full data in the returned object. Defaults to
#'     \code{FALSE}. Set to \code{TRUE} if inteded to use
#'     \code{\link[stats]{simulate}} afterwards.
#' @param engine Ignored; for compatibility with package \code{deepregression}.
#' @param loss Character; specifies the loss function used. The default is
#'  \code{"nll"}, an internal function which takes \code{latent_distr} as an
#'  argument and returns a function with arguments \code{y_true} and
#'  \code{y_pred} to be given to the underlying 'keras' model. Custom loss
#'  functions can be supplied with the same structure, either as a character
#'  or function.
#' @param loss_args Further additional arguments to \code{loss}.
#'
#' @return An object of class \code{c("deeptrafo", "deepregression")}
#'
#' @examples
#' if (.Platform$OS.type != "windows" &&
#'   reticulate::py_available() &&
#'   reticulate::py_module_available("tensorflow") &&
#'   reticulate::py_module_available("keras") &&
#'   reticulate::py_module_available("tensorflow_probability")) {
#'   data("wine", package = "ordinal")
#'   wine$z <- rnorm(nrow(wine))
#'   wine$x <- rnorm(nrow(wine))
#'
#'   nn <- \(x) x |>
#'     layer_dense(input_shape = 1L, units = 2L, activation = "relu") |>
#'     layer_dense(1L)
#'
#'   fml <- rating ~ 0 + temp + contact + s(z, df = 3) + nn(x)
#'
#'   m <- deeptrafo(fml, wine,
#'     latent_distr = "logistic", monitor_metric = NULL,
#'     return_data = TRUE, list_of_deep_models = list(nn = nn)
#'   )
#'
#'   print(m)
#'
#'   m %>% fit(epochs = 10, batch_size = nrow(wine))
#'   coef(m, which_param = "interacting")
#'   coef(m, which_param = "shifting")
#'   fitted(m)
#'   predict(m, type = "pdf")
#'   predict(m, type = "pdf", newdata = wine[, -2])
#'   logLik(m)
#'   logLik(m, newdata = wine[1:10, ])
#'   plot(m)
#'   mcv <- cv(m, cv_folds = 3)
#'   ens <- ensemble(m, n_ensemble = 3)
#'   coef(ens)
#' }
#'
#' @references Kook, L., Baumann, P. F., Dürr, O., Sick, B., & Rügamer, D.
#' (2024). Estimating conditional distributions with neural networks using R
#' package deeptrafo. Journal of Statistical Software.
#' \doi{10.18637/jss.v111.i10}.
#'
#' @importFrom mlt R
#' @importFrom Formula as.Formula
#' @importFrom stats model.matrix model.response model.frame dbeta as.formula
#'     fitted formula predict rmultinom logLik terms drop.terms
#' @importFrom keras layer_dense layer_add layer_concatenate
#' @export
#'
#' @details \code{deeptrafo} is the main function for setting up neural network
#'     transformation models and is called by all aliases for the more special
#'     cases (see e.g. \code{\link[deeptrafo]{ColrNN}}). The naming convention
#'     of the aliases follow the 'tram' package (see e.g. \code{\link[tram]{Colr}})
#'     and add the suffix "NN" to the function name.
#'
deeptrafo <- function(
    formula,
    data,
    response_type = get_response_type(data[[all.vars(fml)[1]]]),
    order = get_order(response_type, data[[all.vars(fml)[1]]]),
    addconst_interaction = 0,
    latent_distr = "logistic",
    loss = "nll",
    loss_args = NULL,
    monitor_metrics = NULL,
    trafo_options = trafo_control(
      order_bsp = order, response_type = response_type
    ),
    return_data = FALSE,
    engine = "tf",
    ...) {
  if (!reticulate::py_available()) {
    message(
      "No Python Environemt available. Use `check_and_install()` ",
      "to install recommended environment."
    )
    invisible(return(NULL))
  }
  if (!reticulate::py_module_available("tensorflow")) {
    message("Tensorflow not available. Use `install_tensorflow()`.")
    invisible(return(NULL))
  }

  call <- match.call()

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
      "~ 0 + ", paste(paste0("ia(", c(int, trimws(
        strsplit(form2text(interacting), "+", fixed = TRUE)[[1]]
      )), ")"), collapse = " + ")
    )
  } else {
    h1_form <- paste0(
      "~ -1 + ", paste(paste0("ia(", int, ")"), collapse = " + ")
    )
  }

  # List of formulas
  list_of_formulas <- list(
    yterms = as.formula(paste0(
      "~ -1 + bsfun(", rvar, ") + bsfunl(", rvar,
      ") + bspfun(", rvar, ")"
    )),
    h1pred = as.formula(h1_form),
    h2 = if (nterms >= 1L) formula(fml, lhs = 0, rhs = 1L) else NULL,
    shared = if (nterms == 2L) formula(fml, lhs = 0, rhs = 2L) else NULL
  )

  # Remove NULL formulae
  list_of_formulas[sapply(list_of_formulas, is.null)] <- NULL

  # Extract response variable
  resp <- data[[rvar]]
  y <- response(resp)

  # check for ATMs
  ftms <- attr(tms <- terms(list_of_formulas$h2), "term.labels")
  is_atm <- any(atps <- grepl("atplag", ftms))
  tlag_formula <- NULL
  if (is_atm) {
    # extract from lag formula the variables as simple sum and
    # layers for additional transformation
    tlag_formula <- paste0(grep("atplag", ftms, value = TRUE), collapse = "+")
    lags <- create_lags(rvar = rvar, d_list = data, atplags = tlag_formula)
    data <- lags$data

    resp <- data[[rvar]] # creating lags reduces data set size
    y <- response(resp)

    tlag_formula <- lags$fm
    list_of_formulas$yterms <- as.formula(
      paste0(form2text(list_of_formulas$yterms), " + ", tlag_formula)
    )
    if (length(ftms) > length(which(atps))) {
      list_of_formulas$h2 <- drop.terms(tms, which(atps))
    } else {
      list_of_formulas$h2 <- ~1
    }
  }

  # define how to get a trafo model from predictor
  from_pred_to_trafo_fun <- from_preds_to_trafo(
    atm_toplayer = trafo_options$atm_toplayer, const_ia = addconst_interaction
  )

  atm_lag_processor <- atm_lag_processor_factory(rvar)

  trafo_processor <- list(
    bsfun = basis_processor, bsfunl = basis_processor_lower,
    bspfun = basisprime_processor, ia = ia_processor,
    atplag = atm_lag_processor
  )

  dots <- list(...)

  if (is.null(dots$additional_processor)) {
    additional_processor <- trafo_processor
  } else {
    additional_processor <- c(list(...)$additional_processor, trafo_processor)
    dots$additional_processor <- NULL
  }

  attr(additional_processor, "controls") <- trafo_options

  # Loss function
  # tloss <- get_loss(response_type, latent_distr)
  tloss <- do.call(loss, c(list(latent_distr = latent_distr), loss_args))

  snwb <- list(subnetwork_init)[rep(1, length(list_of_formulas))]
  snwb[[which(names(list_of_formulas) == "h1pred")]] <-
    h1_init(
      yterms = which(names(list_of_formulas) == "yterms"),
      h1pred = which(names(list_of_formulas) == "h1pred"),
      add_const_positiv = addconst_interaction
    )
  snwb[[which(names(list_of_formulas) == "yterms")]] <- function(...) {
    return(NULL)
  }

  args <- c(list(
    y = y, family = latent_distr, data = data, list_of_formulas = list_of_formulas,
    subnetwork_builder = snwb, from_preds_to_output = from_pred_to_trafo_fun,
    loss = tloss, monitor_metrics = monitor_metrics,
    additional_processor = additional_processor, engine = engine
  ), dots)

  ret <- suppressWarnings(do.call("deepregression", args))

  ret$init_params$is_atm <- is_atm
  ret$init_params$lag_formula <- tlag_formula
  ret$init_params$formula <- formula
  ret$init_params$trafo_options <- trafo_options
  ret$init_params$response_varname <- rvar
  ret$init_params$response_type <- response_type
  ret$init_params$response <- resp
  ret$init_params$prepare_y_valdata <- response
  ret$init_params$data <- if (return_data) data else NULL
  ret$init_params$call <- call

  class(ret) <- c("deeptrafo", "deepregression")
  ret
}

#' Initializes the Processed Additive Predictor for TM's Interaction
#'
#' @param yterms Terms for the response
#' @param h1pred Interacting predictor
#' @param add_const_positiv Shift basis for the predictors to be strictly
#'     positive
#'
#' @return returns a subnetwork_init function with pre-defined arguments
#'
h1_init <- function(yterms, h1pred, add_const_positiv = 0) {
  return(
    function(pp, deep_top, orthog_fun, split_fun, shared_layers,
             param_nr, selectfun_in, selectfun_lay, gaminputs,
             summary_layer = NULL) {
      # instead of passing the respective pp,
      # subsetting is done within subnetwork_init
      # to allow other subnetwork_builder to
      # potentially access all pp entries
      pp_in <- pp_lay <- pp[[h1pred]]
      pp_y <- pp[[yterms]]


      # generate pp parts
      gaminput_nrs <- sapply(pp_in, "[[", "gamdata_nr")
      has_gaminp <- !sapply(gaminput_nrs, is.null)
      gaminput_comb <- sapply(pp_in[which(has_gaminp)], "[[", "gamdata_combined")
      inputs <- makeInputs(pp_in, param_nr = param_nr)
      inputs_y <- makeInputs(pp_y, param_nr = 1)
      org_inputs_for_concat <- list()

      if (sum(has_gaminp)) {
        for (i in 1:sum(has_gaminp)) {
          # concatenate inputs or replace?
          concat <- gaminput_comb[[i]]
          nr <- which(has_gaminp)[i]

          if (!is.null(concat) && concat) {
            org_inputs_for_concat <- c(
              org_inputs_for_concat,
              inputs[[nr]]
            )
            inputs[[nr]] <- layer_concatenate_identity(
              list(gaminputs[[gaminput_nrs[[nr]]]], inputs[[nr]])
            )
          } else {
            inputs[[nr]] <- gaminputs[[gaminput_nrs[[nr]]]]
          }
        }

        inputs_to_replace <- which(has_gaminp)[gaminput_comb]
        keep_inputs_in_return <- setdiff(1:length(inputs), (which(has_gaminp)[!gaminput_comb]))
      } else {
        inputs_to_replace <- c()
        keep_inputs_in_return <- 1:length(inputs)
      }

      layer_matching <- 1:length(pp_in)
      names(layer_matching) <- layer_matching

      if (!is.null(shared_layers)) {
        names_terms <- get_names_pfc(pp_in)

        for (group in shared_layers) {
          layer_ref_nr <- which(names_terms == group[1])
          layer_opts <- get("layer_args", environment(pp_lay[[layer_ref_nr]]$layer))
          layer_opts$name <- paste0(
            "shared_",
            makelayername(
              paste(group, collapse = "_"),
              param_nr
            )
          )
          layer_ref <- do.call(
            get("layer_class", environment(pp_lay[[layer_ref_nr]]$layer)),
            layer_opts
          )

          terms_replace_layer <- which(names_terms %in% group)
          layer_matching[terms_replace_layer] <- layer_ref_nr

          for (i in terms_replace_layer) pp_lay[[i]]$layer <- layer_ref
        }
      }

      if (all(sapply(pp_in, function(x) is.null(x$right_from_oz)))) { # if there is no term to orthogonalize

        outputs <- lapply(1:length(pp_y), function(j) {
          layer_add_identity(
            lapply(1:length(pp_in), function(i) {
              pp_lay[[layer_matching[i]]]$layer(
                pp_y[[j]]$layer(inputs_y[[j]]),
                tf$add(inputs[[i]], add_const_positiv)
              )
            })
          )
        })
        outputs <- layer_concatenate_identity(outputs)

        # replace original inputs
        if (length(org_inputs_for_concat) > 0) {
          inputs[inputs_to_replace] <- org_inputs_for_concat
        }
        return(list(c(inputs_y, inputs[keep_inputs_in_return]), outputs))
      } else {
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
atm_init <- function(atmnr, h1nr) {
  return(
    function(pp, deep_top, orthog_fun, split_fun, shared_layers, param_nr,
             gaminputs) {
      subnetwork_init(pp, deep_top, orthog_fun, split_fun, shared_layers, param_nr,
        # pp_input_subset = atmnr,
        # pp_layer_subset = h1nr,
        gaminputs = gaminputs,
        summary_layer = layer_concatenate_identity
      )
    }
  )
}

#' @title Define Predictor of Transformation Model
#'
#' @param atm_toplayer Function to be applied on top of the transformed lags.
#' @param const_ia See \code{addconst_interaction} in \code{\link[deeptrafo]{deeptrafo}}
#'     or \code{\link[deepregression]{deepregression}}.
#' @param ... For compatibility with 'deepregression'
#' @return A function of \code{list_pred_param} returning a list of output tensors
#' that is passed to \code{model_fun} of \code{deepregression}
#'
#' @details Not intended to be used directly by the end user.
#'
from_preds_to_trafo <- function(
    atm_toplayer = function(x) layer_dense(x, units = 1L, name = "atm_toplayer"),
    const_ia = NULL, ...) {
  return(function(list_pred_param, ...) {
    # make inputs more readable
    # aTtheta <- tf_stride_cols(list_pred_param$h1pred, 1L)
    # aPrimeTtheta <- tf_stride_cols(list_pred_param$h1pred, 2L)
    h1pred_ncol <- list_pred_param$h1pred$shape[[2]]
    shift_pred <- list_pred_param$h2

    if (h1pred_ncol > 3) {
      lag_pred <- tf_stride_cols(list_pred_param$h1pred, 4, h1pred_ncol) %>%
        atm_toplayer()

      # overwrite the shift_pred by adding lags
      shift_pred <- layer_add(list(shift_pred, lag_pred))
    }

    # return transformation
    trafo <- layer_concatenate(list(
      shift_pred,
      tf_stride_cols(list_pred_param$h1pred, 1L, 3L)
    ))

    return(trafo)
  })
}

#' Generic negative log-likelihood for transformation models
#'
#' @param latent_distr Target distribution, character or
#'     \code{tfd_distribution}. If character, can be either "logistic",
#'     "normal", "gumbel", "gompertz".
#'
#' @return A function for computing the negative log-likelihood of a
#'     neural network transformation model with generic response.
#'
#' @import deepregression
#' @import tfprobability
#' @import keras
#' @import tensorflow
#'
nll <- function(latent_distr) {
  if (is.character(latent_distr)) {
    bd <- get_bd(latent_distr)
  } else {
    bd <- latent_distr
  }

  return(
    function(y_true, y_pred) {
      cleft <- tf_stride_cols(y_true, 1L)
      exact <- tf_stride_cols(y_true, 2L)
      cright <- tf_stride_cols(y_true, 3L)
      cint <- tf_stride_cols(y_true, 4L)

      trafo <- layer_add(list(
        tf_stride_cols(y_pred, 1L), # Shift in 1
        tf_stride_cols(y_pred, 2L)
      )) # Upper in 2
      trafo_lwr <- layer_add(list(
        tf_stride_cols(y_pred, 1L),
        tf_stride_cols(y_pred, 3L)
      )) # Lower in 3
      trafo_prime <- tf$math$log(tf$clip_by_value(
        tf_stride_cols(y_pred, 4L),
        1e-8, Inf
      )) # Prime in 4

      ll_exact <- tfd_log_prob(bd, trafo) + trafo_prime
      ll_left <- tf$math$log(tf$clip_by_value(tfd_cdf(bd, trafo), 1e-16, 1))
      ll_right <- tf$math$log(tf$clip_by_value(1 - tfd_cdf(bd, trafo_lwr), 1e-16, 1))
      ll_int <- tf$math$log(tf$clip_by_value(tfd_cdf(bd, trafo) - tfd_cdf(bd, trafo_lwr), 1e-16, 1))

      neglogLik <- -(cleft * ll_left + exact * ll_exact + cright * ll_right +
        cint * ll_int)

      return(neglogLik)
    }
  )
}

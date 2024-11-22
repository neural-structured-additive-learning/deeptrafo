# Transformation ensembles

#' Transformation ensembles
#'
#' @inheritParams deeptrafo
#' @param n_ensemble Numeric; number of ensemble members to fit.
#' @param print_members Logical; print results for each member.
#' @param verbose Logical; whether to print training in each fold.
#' @param save_weights Logical; whether to save the ensemble weights.
#' @param stop_if_nan Logical; whether to stop ensembling if \code{NaN} values
#'     occur
#' @param callbacks List; callbacks used for fitting.
#' @param save_fun Function; function to be applied to each member to be stored
#'     in the final result.
#' @param seed Numeric vector of length \code{n_ensemble}; seeds for model
#'     re-initialization. Changing these seeds does not change the  parameters
#'     of the interacting predictor \code{coef(obj, which_param = "interacting")},
#'     change \code{tf_seeds} to adapt those coefficients.
#' @param tf_seeds Numeric vector of length \code{n_ensemble}; explicit seed for
#'     changing the parameters of the interacting predictor. Distinct from
#'     \code{seed} which is used for weight re-initialization of the rest of the
#'     model (i.e., the shifting predictor and potential neural network components
#'     in the interacting component).
#' @param ... Further arguments passed to \code{deeptrafo} and \code{fit}.
#'
#' @return Ensemble of \code{"deeptrafo"} models with list of training histories
#'     and fitted weights included in \code{ensemble_results}. For details see
#'     the return statment in \code{\link[deepregression]{ensemble}}.
#'
#' @importFrom stats update
#'
#' @export
trafoensemble <- function(
    formula, data, n_ensemble = 5, verbose = FALSE, print_members = TRUE,
    stop_if_nan = TRUE, save_weights = TRUE, callbacks = list(),
    save_fun = NULL, seed = seq_len(n_ensemble), tf_seeds = seq_len(n_ensemble),
    ...
) {

  ret <- list()
  dots <- list(...)
  which_args <- which(names(dots) %in%
                        c("response_type", "order", "addconst_interaction",
                          "latent_distr", "monitor_metrics", "trafo_options",
                          "return_data", "optimizer", "list_of_deep_models",
                          "tf_seed", "return_preproc", "subnetwork_builder",
                          "model_builder", "fitting_function",
                          "additional_processors", "penalty_options",
                          "orthog_options", "weight_options", "formula_options",
                          "output_dim", "verbose"))
  dtargs <- dots[which_args]

  fitargs <- dots[-which_args]

  dargs <- append(list(formula = formula, data = data), dtargs)
  template <- do.call(deeptrafo, dargs)

  for (iter in seq_len(n_ensemble)) {

    if (print_members)
      cat("Fitting member", iter, "...")

    st1 <- Sys.time()

    member <- update(template$init_params, tf_seed = tf_seeds[iter])
    member <- reinit_optimizer(member)
    member <- reinit_weights.deeptrafo(member, seed[iter])

    x_train <- prepare_data(member$init_params$parsed_formulas_content,
                            gamdata = member$init_params$gamdata$data_trafos)

    args <- append(list(object = member$model, x = x_train,
                        y = member$init_params$y, callbacks = callbacks,
                        verbose = verbose, view_metrics = FALSE), fitargs)
    args <- append(args, member$init_params$ellipsis)

    ret[[iter]] <- do.call(member$fit_fun, args)

    if (save_weights)
      ret[[iter]]$weighthistory <- get_weights(member$model)

    if (!is.null(save_fun))
      ret[[iter]]$save_fun_result <- save_fun(member)

    if(stop_if_nan && any(is.nan(ret$metrics$validloss)))
      stop("Member ", iter, " with NaN's in validation loss")

    td <- Sys.time() - st1

    if (print_members)
      cat("\nDone in", as.numeric(td), "", attr(td, "units"), "\n")

  }

  template$ensemble_results <- ret
  class(template) <- c("dtEnsemble", class(template))
  template

}

#' @exportS3Method reinit_weights deeptrafo
reinit_weights.deeptrafo <- function(object, seed) {
  layers <- object$model$layers
  lapply(seq_along(layers), function(x) {
    # x$build(x$input_shape)
    dtype <- layers[[x]]$dtype
    try({
      dshape <- layers[[x]]$kernel$shape
      dinit <- layers[[x]]$kernel_initializer$from_config(
        config = list(seed = tf$constant(seed + x)))
      dweight <- dinit(shape = dshape, dtype = dtype)
      layers[[x]]$set_weights(weights = list(dweight))
    }, silent = TRUE)
  })
  return(invisible(object))
}

reinit_optimizer <- function(x) {
  opt <- x$model$optimizer
  x$model$optimizer <- opt$from_config(opt$get_config())
  return(invisible(x))
}

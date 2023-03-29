# Transformation ensembles

#' Transformation ensembles
#'
#' @inheritParams deeptrafo
#' @param n_ensemble Numeric; number of ensemble members to fit.
#' @param print_members Logical; print results for each member.
#' @param verbose Logical; whether to print training in each fold.
#' @param plot Logical; whether to plot the resulting losses in each fold.
#' @param save_weights Logical; whether to save the ensemble weights.
#' @param stop_if_nan Logical; whether to stop ensembling if \code{NaN} values
#'     occur
#' @param callbacks List; callbacks used for fitting.
#' @param save_fun Function; function to be applied to each member to be stored
#'     in the final result.
#' @param seed Numeric vector of length \code{n_ensemble}; seeds for model
#'     initialization.
#' @param ... Further arguments passed to \code{deeptrafo} and \code{fit}.
#'
#' @return Ensemble of \code{"deeptrafo"} models with list of training histories
#'     and fitted weights included in \code{ensemble_results}. For details see
#'     the return statment in \code{\link[deepregression]{ensemble}}.
#'
#' @export
trafoensemble <- function(
    formula, data, n_ensemble = 5, verbose = FALSE, print_members = TRUE,
    stop_if_nan = TRUE, save_weights = TRUE, callbacks = list(),
    save_fun = NULL, seed = seq_len(n_ensemble),
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

    member <- eval(template$init_params$call)
    member <- reinit_optimizer(member)
    member <- deepregression:::reinit_weights(member, seed[iter])

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
      ret[[iter]]$save_fun_result <- save_fun(this_mod)

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

reinit_optimizer <- function(x) {
  opt <- x$model$optimizer
  x$model$optimizer <- opt$from_config(opt$get_config())
  return(invisible(x))
}

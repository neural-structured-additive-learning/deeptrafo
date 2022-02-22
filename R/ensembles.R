# Deep ensembles with deeptrafo

#' Deep ensembling for deeptrafo models
#'
#' @param x object of class \code{"deeptrafo"}
#' @param n_ensemble numeric; number of ensemble members to fit
#' @param reinitialize logical; if \code{TRUE} (default), model weights are
#'     initialized randomly prior to fitting each member. Fixed weights are
#'     not affected
#' @param print_members logical; print results for each member
#' @param verbose whether to print training in each fold
#' @param patience number of patience for early stopping
#' @param plot whether to plot the resulting losses in each fold
#' @param mylapply lapply function to be used; defaults to \code{lapply}
#' @param save_weights logical, whether to save weights of the ensemble.
#' @param cv_folds an integer if list with train and test data sets
#' @param stop_if_nan logical; whether to stop CV if NaN values occur
#' @param callbacks a list of callbacks used for fitting
#' @param save_fun function applied to the model in each fold to be stored in
#' the final result
#' @param ... further arguments passed to \code{object$fit_fun}
#'
#' @return ensemble of \code{"deeptrafo"} models with list of training histories
#'     and fitted weights included in \code{ensemble_results}
#'
#' @method ensemble deeptrafo
#'
#' @examples
#' dat <- data.frame(y = rnorm(100), x = rnorm(100))
#' m <- deeptrafo(y ~ 0 + x, data = dat)
#' ens <- ensemble(m, n_ensemble = 2)
#' coef(ens)
#'
#' @export
ensemble.deeptrafo <- function(x, n_ensemble = 5, reinitialize = TRUE,
                               mylapply = lapply, verbose = FALSE, patience = 20,
                               plot = TRUE, print_members = TRUE, stop_if_nan = TRUE,
                               save_weights = TRUE, callbacks = list(),
                               save_fun = NULL, ...) {

  ret <- ensemble.deepregression(
    x = x, n_ensemble = n_ensemble,
    reinitialize = reinitialize, mylapply = mylapply,
    verbose = verbose, patience = patience,
    plot = plot, print_members = print_members,
    stop_if_nan = stop_if_nan, save_weights = save_weights,
    callbacks = callbacks, save_fun = save_fun,
    ... = ...
  )

  class(ret) <- c("dtEnsemble", class(ret))

  ret

}

#' @method coef dtEnsemble
#' @inheritParams coef.deeptrafo
#'
#' @export
#'
coef.dtEnsemble <- function(object, which_param = c("shifting", "interacting"),
                            type = NULL, ...) {

  which_param <- match.arg(which_param)
  tparam <- map_param_string_to_index(which_param)

  ret <- .call_for_all_members(object, coef.deepregression,
                               which_param = tparam,
                               type = type, ... = ...)

  lapply(purrr::transpose(ret), function(x) do.call("cbind", x))

}

#' @method fitted dtEnsemble
#' @inheritParams fitted.deeptrafo
#'
#' @export
#'
fitted.dtEnsemble <- function(object, newdata = NULL, batch_size = NULL,
                              convert_fun = as.matrix, ...) {
  .call_for_all_members(object, fitted.deeptrafo, newdata = newdata,
                        batch_size = batch_size, convert_fun = convert_fun,
                        ... = ...)
}

#' @method predict dtEnsemble
#' @inheritParams predict.deeptrafo
#'
#' @export
#'
predict.dtEnsemble <- function(
  object, newdata = NULL, y = newdata[[object$init_params$response_varname]],
  type = c("trafo", "pdf", "cdf", "interaction", "shift", "output"),
  batch_size = NULL, ...) {
  .call_for_all_members(object, predict.deeptrafo, newdata = newdata, y = y,
                        type = type, batch_size = batch_size, ... = ...)
}

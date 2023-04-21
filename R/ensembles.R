# Deep ensembles with deeptrafo

#' Deep ensembling for neural network transformation models
#'
#' @param x Object of class \code{"deeptrafo"}.
#' @param n_ensemble Numeric; number of ensemble members to fit.
#' @param reinitialize Logical; if \code{TRUE} (default), model weights are
#'     initialized randomly prior to fitting each member. Fixed weights are
#'     not affected.
#' @param print_members Logical; print results for each member.
#' @param verbose Logical; whether to print training in each fold.
#' @param patience Integer; number of patience for early stopping.
#' @param plot Logical; whether to plot the resulting losses in each fold.
#' @param mylapply Function; \code{lapply} function to be used; defaults to
#'     \code{lapply}
#' @param save_weights Logical; whether to save the ensemble weights.
#' @param stop_if_nan Logical; whether to stop ensembling if \code{NaN} values
#'     occur
#' @param callbacks List; callbacks used for fitting.
#' @param save_fun Function; function to be applied to each member to be stored
#'     in the final result.
#' @param seed Numeric vector of length \code{n_ensemble}; seeds for model
#'     initialization.
#' @param ... Further arguments passed to \code{object$fit_fun}.
#'
#' @return Ensemble of \code{"deeptrafo"} models with list of training histories
#'     and fitted weights included in \code{ensemble_results}. For details see
#'     the return statment in \code{\link[deepregression]{ensemble}}.
#'
#' @method ensemble deeptrafo
#'
#' @export
ensemble.deeptrafo <- function(x, n_ensemble = 5, reinitialize = TRUE,
                               mylapply = lapply, verbose = FALSE, patience = 20,
                               plot = TRUE, print_members = TRUE, stop_if_nan = TRUE,
                               save_weights = TRUE, callbacks = list(),
                               save_fun = NULL, seed = seq_len(n_ensemble),
                               ...) {

  ret <- ensemble.deepregression(
    x = x, n_ensemble = n_ensemble,
    reinitialize = reinitialize, mylapply = mylapply,
    verbose = verbose, patience = patience,
    plot = plot, print_members = print_members,
    stop_if_nan = stop_if_nan, save_weights = save_weights,
    callbacks = callbacks, save_fun = save_fun, seed = seed,
    ... = ...
  )

  class(ret) <- c("dtEnsemble", class(ret))

  ret

}

#' @method coef dtEnsemble
#'
#' @inheritParams coef.deeptrafo
#' @importFrom purrr transpose
#'
#' @exportS3Method
#'
coef.dtEnsemble <- function(object, which_param = c("shifting", "interacting"),
                            type = NULL, ...) {

  which_param <- match.arg(which_param)
  tparam <- map_param_string_to_index(which_param)

  ret <- .call_for_all_members(object, coef.deepregression,
                               which_param = tparam, type = type, ... = ...)

  lapply(transpose(ret), function(x) do.call("cbind", x))

}

#' @method fitted dtEnsemble
#'
#' @inheritParams fitted.deeptrafo
#'
#' @exportS3Method
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

#' @method logLik dtEnsemble
#'
#' @inheritParams logLik.deeptrafo
#'
#' @exportS3Method
#'
logLik.dtEnsemble <- function(
  object,
  newdata = NULL,
  convert_fun = function(x, ...) - sum(x, ...),
  batch_size = NULL,
  ...
) {

  indiv <- .call_for_all_members(object, logLik.deeptrafo, newdata = newdata, y = y,
                                 convert_fun = convert_fun, ... = ...)

  fitt <- fitted(object, newdata = newdata, batch_size = NULL)
  y_pred <- apply(simplify2array(fitt), 1:2, mean)

  if (is.null(newdata)) {
    y <- object$init_params$y
  } else {
    y <- response(newdata[[object$init_params$response_varname]])
  }

  ensemble_loss <- convert_fun(object$model$loss(y, y_pred)$numpy())

  list(members = unlist(indiv),
       mean = mean(unlist(indiv)),
       ensemble = ensemble_loss)

}

#' @method plot dtEnsemble
#'
#' @inheritParams plot.deeptrafo
#'
#' @exportS3Method
#'
plot.dtEnsemble <- function(
  x,
  which = NULL,
  type = c("smooth", "trafo", "pdf", "cdf"),
  newdata = NULL,
  which_param = c("shifting", "interacting"),
  only_data = FALSE,
  K = 40,
  q = NULL,
  ...
) {
  pdat <- .call_for_all_members(
    x, plot.deeptrafo, which = which, type = type, which_param = which_param,
    newdata = newdata, only_data = TRUE, K = K, q = q, ... = ...
  ) |> invisible()

  if (only_data)
    return(pdat)

  type <- match.arg(type)
  rname <- x$init_params$response_varname
  rtype <- x$init_params$response_type

  if (type == "smooth") {

    pdat <- transpose(pdat)
    terms <- names(pdat)

    lapply(terms, \(term) {

      pd <- pdat[[term]]
      value <- pd[[1]]$value
      peff <- do.call("cbind", lapply(pd, \(x) x$partial_effect))
      matplot(sort(value), peff[order(value),], type = "b", pch = 1, main = term)
    })

  } else if (!is.null(newdata)) {
    if (!is.null(newdata[[rname]])) {
      y <- attr(pdat[[1]], "y")
      preds <- do.call("cbind", pdat)
      matplot(sort(y), preds[order(y),], type = "p")
    } else {
      y <- attr(pdat[[1]], "y")
      preds <- do.call("rbind", pdat)
      ttype <- ifelse(rtype %in% c("ordered", "count"), "s", "l")
      print(ttype)
      matplot(y, t(preds), type = ttype, lty = rep(
        1:length(pdat), each = NROW(pdat[[1]])))
    }
  } else {
    y <- attr(pdat[[1]], "y")
    preds <- do.call("cbind", pdat)
    matplot(sort(y), preds[order(y),], type = "p")
  }

  return(invisible(pdat))

}

#' Tune and evaluate weighted transformation ensembles
#'
#' @param object Object of class \code{"dtEnsemble"}
#' @param weights Numeric; weight-vector of length \code{n_ensemble}, if
#'     \code{NULL} the weights are tuned on \code{newdata}
#' @param newdata List or data.frame; new data to evaluate or tune the weights
#'     on
#' @inheritParams logLik.deeptrafo
#'
#' @importFrom stats optim weighted.mean
#'
#' @return Returns list of ensemble members, average, and ensemble
#'    log-likelihood converted by \code{convert_fun}
#' @export
#'
weighted_logLik <- function(
    object,
    weights = NULL,
    newdata = NULL,
    convert_fun = function(x, ...) mean(x, ...),
    batch_size = NULL,
    ...
) {

  stopifnot(inherits(object, "dtEnsemble"))

  indiv <- .call_for_all_members(
    object, logLik.deeptrafo, newdata = newdata, y = y,
    convert_fun = convert_fun, ... = ...
  )

  fitt <- fitted(object, newdata = newdata, batch_size = NULL)

  if (is.null(newdata)) {
    y <- object$init_params$y
  } else {
    y <- response(newdata[[object$init_params$response_varname]])
  }

  if (is.null(weights)) {
    obj <- function(weights) {
      y_pred <- apply(simplify2array(fitt), 1:2, weighted.mean, w = weights)
      convert_fun(object$model$loss(y, y_pred)$numpy())
    }
    opt <- optim(rep(1, length(fitt)), obj, lower = .Machine$double.eps,
                 upper = 1 - .Machine$double.eps, method = "L-BFGS-B")
    ensemble_loss <- opt$value
    weights <- opt$par
    weights <- weights / sum(weights)
  } else {
    y_pred <- apply(simplify2array(fitt), 1:2, weighted.mean, w = weights)
    ensemble_loss <- convert_fun(object$model$loss(y, y_pred)$numpy())
  }

  list(members = unlist(indiv),
       mean = mean(unlist(indiv)),
       ensemble = ensemble_loss,
       weights = weights)

}

# Helpers

#' @importFrom keras set_weights
.call_for_all_members <- function (object, FUN, ...) {
  ens_weights <- lapply(object$ensemble_results, function(x) {
    x$weighthistory
  })
  lapply(ens_weights, function(x) {
    set_weights(object$model, x)
    FUN(object, ... = ...)
  })
}


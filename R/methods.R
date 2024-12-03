#' Plot method for deep conditional transformation models
#'
#' @param x Object of class \code{"deeptrafo"}.
#' @param which Which effect to plot, default selects all smooth effects in the
#'     shift term.
#' @param which_param Character; either \code{"interacting"} or \code{"shifting"}.
#' @param only_data Logical, if \code{TRUE}, only the data for plotting is returned.
#' @param K Integer; If \code{type == "smooth"} the length of an equidistant
#'     grid at which a two-dimensional function is evaluated for plotting.
#'     Otherwise, length of the grid to evaluate predictions at,
#'     see \code{newdata}.
#' @param q Vector of response values to compute predictions at, see \code{newdata}
#' @param ... Further arguments, passed to fit, plot or predict function
#' @param type Character; One of "smooth", "trafo", "pdf", or "cdf".
#' @param newdata Optional new data (\code{list} or \code{data.frame}) to
#'     evaluate predictions at. If the response is missing, plots are generated
#'     on a grid of length \code{K}
#'
#' @method plot deeptrafo
#' @exportS3Method
#'
#' @importFrom graphics matplot
#' @importFrom grDevices rgb
#'
#'
plot.deeptrafo <- function(
    x, # object
    which = NULL,
    type = c("smooth", "trafo", "pdf", "cdf"),
    newdata = NULL,
    which_param = c("shifting", "interacting"), # for which parameter
    only_data = FALSE,
    K = 40,
    q = NULL,
    ... # passed to plot function
    ) {
  type <- match.arg(type)

  if (x$init_params$is_atm && !is.null(newdata)) {
    lags <- fm_to_lag(x$init_params$lag_formula)
    newdata <- create_lags(
      rvar = x$init_params$response_varname,
      d_list = newdata,
      lags = lags
    )$data
  }

  if (type == "smooth") {
    which_param <- match.arg(which_param)

    get_weight_fun <- switch(which_param,
      "interacting" = get_weight_by_name_ia,
      "shifting" = get_weight_by_name
    )

    which_param <- map_param_string_to_index(which_param)

    return(plot.deepregression(x,
      which = which, which_param = which_param,
      only_data = only_data, grid_length = K,
      get_weight_fun = get_weight_fun, ...
    ))
  } else {
    rname <- x$init_params$response_varname
    rtype <- x$init_params$response_type
    ry <- x$init_params$response
    preds <- predict.deeptrafo(x, type = type, newdata = newdata, K = K, ...)

    if (is.null(newdata)) {
      if (only_data) {
        return(structure(preds, y = ry))
      }
      plot(ry, preds, xlab = "response", ylab = type)
    } else {
      if (!is.null(newdata[[rname]])) {
        y <- newdata[[rname]]
        if (only_data) {
          return(structure(preds, y = y))
        }
        plot(y, preds, xlab = "response", ylab = type)
      } else {
        y <- as.numeric(names(preds))
        if (rtype %in% c("ordered", "count")) {
          ttype <- "s"
        } else {
          ttype <- "l"
        }
        preds <- do.call("cbind", preds)
        if (only_data) {
          return(structure(preds, y = y))
        }
        matplot(y, t(preds),
          type = ttype, col = rgb(.1, .1, .1, .5), lty = 1,
          xlab = "response", ylab = type
        )
      }
    }
  }
}

get_weight_by_name_ia <- function(x, name, param_nr) {
  matrix(get_weight_by_name(x, name, param_nr),
    ncol = x$init_params$trafo_options$order_bsp + 1L
  )
}

#' S3 methods for deep conditional transformation models
#'
#' @param x Object of class \code{"deeptrafo"}.
#' @param which_param Character; either \code{"shifting"}, \code{"interacting"},
#'     or \code{"autoregressive"} (only for autoregressive transformation models).
#' @param type Either NULL (all types of coefficients are returned),
#'     "linear" for linear coefficients or "smooth" for coefficients of;
#'     Note that \code{type} is currently not used for \code{"interacting"}.
#' @param ... Further arguments, passed to fit, plot or predict function
#'
#' @method coef deeptrafo
#' @importFrom stats coef
#'
#' @export
#' @rdname methodTrafo
#'
coef.deeptrafo <- function(
    object,
    which_param = c("shifting", "interacting", "autoregressive"),
    type = NULL,
    ...) {
  which_param <- match.arg(which_param)

  if (which_param == "autoregressive") {
    ret <- try(c(get_weight_by_opname(object, name = "atm_toplayer", partial_match = TRUE)))
    if (inherits(ret, "try-error")) stop("No layer with name atm_toplayer")
    names(ret) <- unlist(strsplit(object$init_params$lag_formula, "\\+"))
    return(ret)
  }

  is_interaction <- which_param == "interacting"
  which_param <- map_param_string_to_index(which_param)

  ret <- coef.deepregression(object, which_param = which_param, type = type)

  if (is_interaction) {
    ret <- lapply(ret, function(r) {
      reshape_softplus_cumsum(r, order_bsp_p1 = get_order_bsp_p1(object))
    })

    if (object$init_params$response_type == "ordered") {
      ret <- lapply(ret, function(r) r[-nrow(r), , drop = FALSE])
    }
  }

  return(ret)
}

#' @exportS3Method coef LmNN
coef.LmNN <- function(object, which_param = c("shifting", "interacting", "autoregressive"),
                      type = NULL, ...) {
  which_param <- match.arg(which_param)

  is_interaction <- which_param == "interacting"
  which_param <- map_param_string_to_index(which_param)

  # else, return lags
  ret <- coef.deepregression(object, which_param = which_param, type = type)

  if (is_interaction) {
    ret <- lapply(ret, function(r) {
      x <- matrix(r, nrow = 2, byrow = TRUE)
      x[2, ] <- softplus(x[2, ])
      x
    })
  }

  return(ret)
}

#' @exportS3Method coef SurvregNN
coef.SurvregNN <- function(object, which_param = c("shifting", "interacting", "autoregressive"),
                           type = NULL, ...) {
  which_param <- match.arg(which_param)

  is_interaction <- which_param == "interacting"
  which_param <- map_param_string_to_index(which_param)

  # else, return lags
  ret <- coef.deepregression(object, which_param = which_param, type = type)

  if (is_interaction) {
    ret <- lapply(ret, function(r) {
      x <- matrix(r, nrow = 2, byrow = TRUE)
      x[2, ] <- softplus(x[2, ])
      x
    })
  }

  return(ret)
}


#' @param object Object of class \code{"deeptrafo"}.
#' @param newdata Optional new data, either \code{data.frame} or named \code{list}.
#' @param K Integer; grid length for the response to evaluate predictions at,
#'     if \code{newdata} does not contain the response.
#' @param q Numeric or factor; user-supplied grid of response values to evaluate
#'     the predictions. Defaults to \code{NULL}. If overwritten, \code{K} is
#'     ignored.
#' @param pred_grid Logical; set TRUE, if user provides a predefined grid for an
#'     atp/atm model through newdata which holds two attributes. The first
#'     attribute, rname, should hold the column name (string) of the response
#'     variable while the second attribute, y, should hold the grid name.
#' @param ... Currently ignored.
#'
#' @return Returns vector or matrix of predictions, depending on the supplied
#'     \code{type}.
#'
#' @details If no new data is supplied, predictions are computed on the training
#'     data (i.e. in-sample). If new data is supplied without a response,
#'     predictions are evaluated on a grid of length \code{K}.
#'
#' @method predict deeptrafo
#' @importFrom variables numeric_var ordered_var mkgrid
#' @export
#'
#' @rdname methodTrafo
#'
predict.deeptrafo <- function(
    object,
    newdata = NULL,
    type = c("trafo", "pdf", "cdf", "interaction", "shift", "terms"),
    batch_size = NULL,
    K = 1e2,
    q = NULL,
    pred_grid = FALSE,
    ...) {
  # TODO: make prediction possible for one observation only
  type <- match.arg(type)

  rname <- object$init_params$response_varname
  rtype <- object$init_params$response_type
  order <- object$init_params$trafo_options$order_bsp
  fam <- object$init_params$family
  discrete <- as.numeric(rtype %in% c("count", "ordered"))
  bd <- get_bd(fam)

  if (object$init_params$is_atm && !is.null(newdata)) {
    lags <- fm_to_lag(object$init_params$lag_formula)
  }

  # Predict over grid of responses, if response not contained in newdata
  if (!is.null(newdata)) {
    if (is.null(newdata[[rname]])) {
      ygrd <- if (is.null(q)) {
        make_grid(object$init_params$response, n = K)[[1]]
      } else {
        q
      }
      if (type == "shift") { # shift independent of response, skip
        ygrd <- ygrd[1]
      }
      ret <- lapply(ygrd, function(ty) { # overwrite response, then predict
        newdata[[rname]] <- rep(ty, NROW(newdata[[1]]))
        if (object$init_params$is_atm) {
          newdata <- create_lags(
            rvar = rname, d_list = newdata, lags = lags,
            pred_grid = pred_grid
          )$data
        }
        predict.deeptrafo(object,
          newdata = newdata, type = type,
          batch_size = batch_size, K = NULL, q = NULL,
          ... = ...
        )
      })
      names(ret) <- as.character(ygrd)
      return(ret)
    }
  }

  if (object$init_params$is_atm && !is.null(newdata)) {
    lags <- fm_to_lag(object$init_params$lag_formula)
    newdata <- create_lags(
      rvar = rname, d_list = newdata, lags = lags,
      pred_grid = pred_grid
    )$data
  }

  # Compute predictions from fitted values
  mod_output <- fitted.deeptrafo(object, newdata,
    batch_size = batch_size,
    call_create_lags = FALSE
  )

  if (type == "terms") {
    return(mod_output)
  }

  if (is.null(newdata)) {
    ry <- object$init_params$y
  } else {
    ry <- response(newdata[[rname]])
  }

  cleft <- ry[, "cleft", drop = FALSE]
  cint <- ry[, "cinterval", drop = FALSE]
  cright <- ry[, "cright", drop = FALSE]

  w_eta <- mod_output[, 1, drop = FALSE]
  aTtheta <- mod_output[, 2, drop = FALSE]
  alTtheta <- mod_output[, 3, drop = FALSE]
  apTtheta <- mod_output[, 4, drop = FALSE]

  if (type == "interaction") {
    return(as.matrix(aTtheta))
  }

  if (type == "shift") {
    return(as.matrix(w_eta))
  }

  ytransf <- aTtheta + w_eta

  if (type == "trafo") {
    trf <- ytransf %>% as.matrix()
    if (rtype == "ordered") {
      trf[which(c(cright) == 1), ] <- Inf
    }
    return(trf)
  }

  if (type == "cdf") {
    if (discrete) {
      cdf <- (cleft + cint) * as.matrix(tfd_cdf(bd, ytransf)) +
        cright * as.matrix(tfd_cdf(bd, rep(1e8, nrow(cright))))
    } else {
      cdf <- bd %>% tfd_cdf(ytransf)
    }

    return(cdf %>% as.matrix())
  } else if (type == "pdf") {
    yprimeTrans <- apTtheta

    if (discrete) {
      ytransflower <- alTtheta + w_eta

      pdf <- cint * as.matrix(tfd_cdf(bd, ytransf) - tfd_cdf(bd, ytransflower)) +
        cleft * tfd_cdf(bd, ytransf) + cright * tfd_survival_function(bd, ytransflower)
    } else {
      pdf <- as.matrix(tfd_prob(bd, ytransf)) * as.matrix(yprimeTrans)
    }

    return(pdf %>% as.matrix())
  }
}

#' @param object Object of class \code{"deeptrafo"}
#' @param newdata Optional new data, either \code{data.frame} or named \code{list}
#' @param batch_size Integer; optional, useful if data is too large.
#' @param convert_fun Function; to convert the TF tensor.
#' @param call_create_lags Logical; lags may already be computed by a different method (e.g. plot)
#' @param ... Currently ignored.
#'
#' @return Returns matrix of fitted values.
#'
#' @method fitted deeptrafo
#' @export
#'
#' @rdname methodTrafo
#'
fitted.deeptrafo <- function(
    object,
    newdata = NULL,
    batch_size = NULL,
    convert_fun = as.matrix,
    call_create_lags = TRUE, # in predict() already called once
    ...) {
  l_fm <- object$init_params$lag_formula
  if ((object$init_params$is_atm && !is.null(newdata)) && call_create_lags) {
    lags <- fm_to_lag(l_fm)
    newdata <- create_lags(
      rvar = object$init_params$response_varname,
      d_list = newdata,
      lags = lags
    )$data
  }

  if (length(object$init_params$image_var) > 0 | !is.null(batch_size)) {
    mod_output <- predict_gen(object, newdata, batch_size,
      apply_fun = function(x) x,
      convert_fun = convert_fun
    )
  } else {
    if (is.null(newdata)) {
      newdata <- prepare_data(object$init_params$parsed_formulas_contents,
        gamdata = object$init_params$gamdata$data_trafos
      )
    } else {
      newdata <- prepare_newdata(object$init_params$parsed_formulas_contents, newdata,
        gamdata = object$init_params$gamdata$data_trafos
      )
    }

    mod_output <- object$model(newdata)
  }

  convert_fun(mod_output)
}

map_param_string_to_index <- function(which_param) {
  switch(which_param,
    "interacting" = 2,
    "shifting" = 3
  )
}

#' @method logLik deeptrafo
#'
#' @param object Object of class \code{"deeptrafo"}.
#' @param newdata Named \code{list} or \code{data.frame}; optional new data.
#' @param convert_fun Function; applied to the log-likelihood values of all
#'     observations.
#' @param ... Additional arguments to \code{fitted.deeptrafo()}
#'
#' @exportS3Method
#'
#' @rdname methodTrafo
#'
logLik.deeptrafo <- function(
    object,
    newdata = NULL,
    convert_fun = function(x, ...) -sum(x, ...),
    ...) {
  if (object$init_params$is_atm && !is.null(newdata)) {
    lags <- fm_to_lag(object$init_params$lag_formula)
    newdata <- create_lags(
      rvar = object$init_params$response_varname,
      d_list = newdata,
      lags = lags
    )$data
  }

  if (is.null(newdata)) {
    y <- object$init_params$y
    y_pred <- fitted.deeptrafo(object, call_create_lags = FALSE, ... = ...)
  } else {
    y <- response(newdata[[object$init_params$response_varname]])
    y_pred <- fitted.deeptrafo(object,
      call_create_lags = FALSE,
      newdata = newdata, ... = ...
    )
  }

  convert_fun(object$model$loss(y, y_pred)$numpy())
}

#' @method residuals deeptrafo
#'
#' @param object Object of class \code{"deeptrafo"}.
#' @param newdata Named \code{list} or \code{data.frame}; optional new data.
#' @param ... Additional arguments to \code{fitted.deeptrafo()}
#' @param return_gradients Return individual gradients instead of the summed
#'     gradients; the residuals are \code{0.5 * rowSums(gradients)}
#'
#' @exportS3Method
#'
#' @rdname methodTrafo
#'
residuals.deeptrafo <- function(
    object,
    newdata = NULL,
    return_gradients = FALSE,
    ...) {
  tape <- \() NULL
  with(tf$GradientTape() %as% tape, {
    if (is.null(newdata)) {
      y <- object$init_params$y
      y_pred <- fitted.deeptrafo(object, ... = ...)
    } else {
      y <- response(newdata[[object$init_params$response_varname]])
      y_pred <- fitted.deeptrafo(object, newdata = newdata, ... = ...)
    }
    a <- y_pred
    a[] <- 0
    a <- tf$Variable(a)
    nlli <- object$model$loss(y, y_pred + a)
  })

  gr <- tape$gradient(nlli, a)$numpy()[, -4]

  if (return_gradients) {
    return(gr)
  }

  0.5 * rowSums(gr)
}

#' @method simulate deeptrafo
#'
#' @param object Object of class \code{"deeptrafo"}.
#' @param nsim Integer; number of simulations; defaults to 1.
#' @param seed Seed for generating samples; defaults to \code{NULL}.
#' @param newdata Named \code{list} or \code{data.frame}; optional new data.
#' @param ... Further arguments to \link[predict.deeptrafo]{predict.deeptrafo}
#'
#' @exportS3Method
#' @importFrom stats simulate
#' @rdname methodTrafo
#'
simulate.deeptrafo <- function(
    object, nsim = 1, seed = NULL, newdata = NULL, ...) {
  rtype <- object$init_params$response_type

  if (rtype != "ordered") {
    stop("Simulate not yet implemented for response types other than ordered.")
  }

  rvar <- object$init_params$response_varname

  if (is.null(newdata)) {
    newdata <- object$init_params$data
    newdata <- newdata[-which(names(newdata) == rvar)]
  } else {
    ry <- object$init_params$response
  }

  cdf <- do.call("cbind", predict(object, newdata = newdata, type = "cdf"))
  pmf <- apply(cbind(0, cdf), 1, diff)

  ret <- lapply(1:nsim, function(x) {
    ordered(apply(pmf, 2, function(probs) {
      which(rmultinom(n = 1, size = 1, prob = probs) == 1)
    }), levels = levels(object$init_params$response))
  })

  if (nsim == 1) {
    ret <- ret[[1]]
  }

  ret
}

#' @method print deeptrafo
#'
#' @param x Object of class \code{"deeptrafo"}.
#' @param print_model Logical; print keras model.
#' @param print_coefs Logical; print coefficients.
#' @param with_baseline Logical; print baseline coefs.
#' @param ... Currently ignored.
#'
#' @exportS3Method
#'
#' @rdname methodTrafo
#'
print.deeptrafo <- function(x, print_model = FALSE, print_coefs = TRUE,
                            with_baseline = FALSE, ...) {
  atm <- x$init_params$is_atm
  atm_text <- if (atm) "autoregressive" else NULL

  if (print_model) {
    print(x$model)
  }

  mtype <- switch(x$init_params$response_type,
    "ordered" = "ordinal",
    "count" = "count",
    "survival" = "continuous",
    "continuous" = "continuous"
  )

  fmls <- x$init_params$list_of_formulas

  no_int <- (fmls[[2]] == ~ -1 + ia(1))
  no_shift <- (fmls[[3]] == ~1)
  uncond <- no_int & no_shift

  int <- ifelse(no_int, paste0(x$init_params$response_varname, " | 1"),
    fml2txt(formula(x$init_params$formula, lhs = 2L, rhs = 0L)[[2]])
  )
  shift <- ifelse(no_shift, "~1", fml2txt(fmls[[3]]))

  trained <- ifelse(is.null(x$model$history), "Untrained", "Trained")
  cat("\t", trained, mtype, "outcome", atm_text, "deep conditional transformation model\n")
  cat("\nCall:\n")
  print(x$init_params$call)
  cat("\nInteracting: ", int, "\n")
  cat("\nShifting: ", shift, "\n")

  if (atm) {
    lags <- fml2txt(as.formula(paste0("~", x$init_params$lag_formula)))
    cat("\nLags: ", lags, "\n")
  }


  if (print_coefs) {
    cfb <- coef(x, which_param = "interacting")
    if (no_int) {
      names(cfb) <- x$init_params$response_varname
    }
    if (with_baseline) {
      cat("\nBaseline transformation:\n")
      print(unlist(cfb))
    }
    cat("\nShift coefficients:\n")
    cfx <- coef(x)
    rns <- lapply(cfx, rownames)
    which_no_names <- which(unlist(lapply(rns, is.null)))
    if (length(which_no_names) > 0) {
      rns[which_no_names] <- names(rns)[which_no_names]
    }
    names(cfx) <- rns
    print(unlist(cfx))
    if (atm) {
      cat("\nLag coefficients:\n")
      print(unlist(coef(x, which_param = "autoregressive")))
    }
  }

  return(invisible(x))
}

#' @method summary deeptrafo
#'
#' @param object Object of class \code{"deeptrafo"}.
#' @param ... Further arguments supplied to \code{print.deeptrafo}
#'
#' @exportS3Method
#' @rdname methodTrafo
#'
summary.deeptrafo <- function(object, ...) {
  print(object, print_model = TRUE, ...)
}

# Helpers -----------------------------------------------------------------

get_bd <- function(latent_distr) {
  if (inherits(latent_distr, "python.builtin.object")) {
    return(latent_distr)
  }
  switch(latent_distr,
    "normal" = tfd_normal(loc = 0, scale = 1),
    "logistic" = tfd_logistic(loc = 0, scale = 1),
    "gumbel" = tfd_gumbel(loc = 0, scale = 1),
    "gompertz" = tfd_gompertz(loc = 0, scale = 1)
  )
}

fml2txt <- function(formula) Reduce(paste, deparse(formula))


# Response evaluation a la mlt
response <- function(y) {
  rtype <- get_response_type(y)
  bound <- switch(
    rtype,
    "continuous" = c(-Inf, Inf),
    "ordered" = NA,
    "count" = c(0, Inf),
    "survival" = c(0, Inf)
  )
  resp <- R(y, bounds = bound)
  interval <- as.numeric(.cinterval(resp))
  left <- abs(as.numeric(.cleft(resp)) - interval)
  right <- abs(as.numeric(.cright(resp)) - interval)
  interval <- abs(interval - left - right)
  exact <- as.numeric(.exact(resp))
  structure(cbind(cleft = left, exact = exact, cright = right,
                  cinterval = interval), type = get_response_type(y))
}

# From package mlt {
.exact <- function(object) !is.na(object$exact)
.cinterval <- function(object) !.exact(object)
.cright <- function(object) is.finite(object$cright)
.cleft <- function(object) is.finite(object$cleft)
# }

make_grid <- function(y, n = 1e2) {
  rtype <- get_response_type(y)
  var <- switch(
    rtype,
    "continuous" = numeric_var(name = "y", support = range(y), bounds = c(-Inf, Inf)),
    "ordered" = ordered_var(name = "y", levels = levels(y), bounds = NA),
    "count" = numeric_var(name = "y", support = range(y)),
    "survival" = numeric_var(name = "y", support = range(y[, 1]), bounds = c(0, Inf))
  )
  mkgrid(var, n = n)
}

#' @importFrom survival is.Surv
get_response_type <- function(y) {
  ret <- if (is.ordered(y))
    "ordered"
  else if (is.integer(y))
    "count"
  else if (is.Surv(y))
    "survival"
  else
    "continuous"
  ret
}

get_order <- function(response_type, y) {
  ret <- if (response_type == "ordered")
    nlevels(y) - 1L
  else
    10
  ret
}

eval_response <- function(y, response_type) {
  switch(
    response_type,
    "continuous" = y,
    "ordered" = t(sapply(y, eval_ord)),
    "count" = cbind(as.numeric(y == 0L), y),
    "survival" = as.matrix(y)
  )
}

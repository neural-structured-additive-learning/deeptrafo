
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
  interval <- as.numeric(mlt:::.cinterval(resp))
  left <- abs(as.numeric(mlt:::.cleft(resp)) - interval)
  right <- abs(as.numeric(mlt:::.cright(resp)) - interval)
  interval <- abs(interval - left - right)
  exact <- as.numeric(mlt:::.exact(resp))
  structure(cbind(cleft = left, exact = exact, cright = right,
                  cinterval = interval), type = get_response_type(y))
}

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

get_response_type <- function(y) {
  ret <- if (is.ordered(y))
    "ordered"
  else if (is.integer(y))
    "count"
  else if (survival::is.Surv(y))
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

get_loss <- function(response_type, family) {
  switch(
    response_type,
    "continuous" = neg_ll_trafo(family),
    "ordered" = nll_ordinal(family),
    "count" = nll_count(family),
    "survival" = nll_surv(family)
  )
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

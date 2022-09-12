################################################################################################

#' @importFrom tensorflow tf
#' @importFrom tfprobability tfd_normal tfd_logistic tfd_log_prob tfd_gumbel
tf_crps_norm <- function(h)
{
  return(
    h * (2 * tfd_normal(0,1)$cdf(h) - 1) +
      (tf$math$sqrt(2) * tf$math$exp(-0.5 * tf$square(h)) - 1)/tf$math$sqrt(pi)
  )
}

#'
#'
#' @importFrom  keras custom_metric
#'
crps_stdnorm_metric <- custom_metric(name = "crps_stdnorm", metric_fn = function(y_true, y_pred) {

  h_y_eval <- y_pred[,1, drop = FALSE] + y_pred[,2, drop = FALSE]

  return(tf_crps_norm(h_y_eval))

})

################################################################################################

# @PK: add crps_stdlogistic?

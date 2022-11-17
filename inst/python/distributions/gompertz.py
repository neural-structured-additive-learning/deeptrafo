# Taken from TensorFlow Probability (Apache V2).
# Adapted from the Gumbel distribution class in
# https://github.com/tensorflow/probability/blob/06e643521caa196c4da6ad616e12240134a03312/tensorflow_probability/python/distributions/gumbel.py

# Dependency imports
import numpy as np
import tensorflow.compat.v2 as tf

from tensorflow_probability.python.bijectors import bijector
from tensorflow_probability.python.bijectors import identity as identity_bijector
from tensorflow_probability.python.bijectors import invert as invert_bijector
from tensorflow_probability.python.bijectors import softplus as softplus_bijector
from tensorflow_probability.python.distributions import kullback_leibler
from tensorflow_probability.python.distributions import transformed_distribution
from tensorflow_probability.python.distributions import uniform
from tensorflow_probability.python.internal import assert_util
from tensorflow_probability.python.internal import dtype_util
from tensorflow_probability.python.internal import parameter_properties
from tensorflow_probability.python.internal import tensor_util

"""Gompertz bijector."""

__all__ = [
    'GompertzCDF',
]

class GompertzCDF(bijector.AutoCompositeTensorBijector):
  """Compute `Y = g(X) = 1 - exp(-exp((X - loc) / scale))`, the Gompertz CDF.

  This bijector maps inputs from `[-inf, inf]` to `[0, 1]`. The inverse of the
  bijector applied to a uniform random variable `X ~ U(0, 1)` gives back a
  random variable with the
  [Gompertz distribution](https://en.wikipedia.org/wiki/Gompertz_distribution):

  ```none
  Y ~ GompertzCDF(loc, scale)
  pdf(y; loc, scale) = exp(
    -( (y - loc) / scale + exp(- (y - loc) / scale) ) ) / scale
  ```
  """

  def __init__(self,
               loc=0.,
               scale=1.,
               validate_args=False,
               name='gompertz_cdf'):
    """Instantiates the `GompertzCDF` bijector.

    Args:
      loc: Float-like `Tensor` that is the same dtype and is
        broadcastable with `scale`.
        This is `loc` in `Y = g(X) = 1-exp(-exp((X - loc) / scale))`.
      scale: Positive Float-like `Tensor` that is the same dtype and is
        broadcastable with `loc`.
        This is `scale` in `Y = g(X) = 1-exp(-exp((X - loc) / scale))`.
      validate_args: Python `bool` indicating whether arguments should be
        checked for correctness.
      name: Python `str` name given to ops managed by this object.
    """
    parameters = dict(locals())
    with tf.name_scope(name) as name:
      dtype = dtype_util.common_dtype([loc, scale], dtype_hint=tf.float32)
      self._loc = tensor_util.convert_nonref_to_tensor(
          loc, dtype=dtype, name='loc')
      self._scale = tensor_util.convert_nonref_to_tensor(
          scale, dtype=dtype, name='scale')
      super(GompertzCDF, self).__init__(
          validate_args=validate_args,
          forward_min_event_ndims=0,
          parameters=parameters,
          name=name)

  @classmethod
  def _parameter_properties(cls, dtype):
    return dict(
        loc=parameter_properties.ParameterProperties(),
        scale=parameter_properties.ParameterProperties(
            default_constraining_bijector_fn=(
                lambda: softplus_bijector.Softplus(low=dtype_util.eps(dtype)))))

  @property
  def loc(self):
    """The `loc` in `Y = g(X) = 1-exp(-exp((X - loc) / scale))`."""
    return self._loc

  @property
  def scale(self):
    """This is `scale` in `Y = g(X) = 1-exp(-exp((X - loc) / scale))`."""
    return self._scale

  @classmethod
  def _is_increasing(cls):
    return True

  def _forward(self, x):
    z = (x - self.loc) / self.scale
    return tf.add(1., - tf.exp(-tf.exp(z)))

  def _inverse(self, y):
    with tf.control_dependencies(self._maybe_assert_valid_y(y)):
      return self.loc + self.scale * tf.math.log(-tf.math.log(1. - y))

  def _inverse_log_det_jacobian(self, y):
    with tf.control_dependencies(self._maybe_assert_valid_y(y)):
      return tf.math.log(self.scale / (tf.math.log(y) * (1. - y)))

  def _forward_log_det_jacobian(self, x):
    scale = tf.convert_to_tensor(self.scale)
    z = (x - self.loc) / scale
    return z - tf.exp(z) + tf.math.log(scale)

  def _maybe_assert_valid_y(self, y):
    if not self.validate_args:
      return []
    is_positive = assert_util.assert_non_negative(
        y, message='Inverse transformation input must be greater than 0.')
    less_than_one = assert_util.assert_less_equal(
        y,
        tf.constant(1., y.dtype),
        message='Inverse transformation input must be less than or equal to 1.')
    return [is_positive, less_than_one]

  def _parameter_control_dependencies(self, is_init):
    if not self.validate_args:
      return []
    assertions = []
    if is_init != tensor_util.is_ref(self.scale):
      assertions.append(assert_util.assert_positive(
          self.scale,
          message='Argument `scale` must be positive.'))
    return assertions

"""The Gompertz distribution class."""

class Gompertz(transformed_distribution.TransformedDistribution):
  """The scalar Gompertz distribution with location `loc` and `scale` parameters.

  #### Mathematical details

  The probability density function (pdf) of this distribution is,

  ```none
  pdf(x; mu, sigma) = exp((x - mu) / sigma - exp((x - mu) / sigma)) / sigma
  ```

  where `loc = mu` and `scale = sigma`.

  The cumulative density function of this distribution is,

  ```none
  cdf(x; mu, sigma) = exp(-exp((x - mu) / sigma))
  ```

  The Gompertz distribution is a member of the [location-scale family](
  https://en.wikipedia.org/wiki/Location-scale_family), i.e., it can be
  constructed as,

  ```none
  X ~ Gompertz(loc=0, scale=1)
  Y = loc + scale * X
  ```

  #### Examples

  Examples of initialization of one or a batch of distributions.

  ```python
  tfd = tfp.distributions

  # Define a single scalar Gompertz distribution.
  dist = tfd.Gompertz(loc=0., scale=3.)

  # Evaluate the cdf at 1, returning a scalar.
  dist.cdf(1.)

  # Define a batch of two scalar valued Gompertzs.
  # The first has mean 1 and scale 11, the second 2 and 22.
  dist = tfd.Gompertz(loc=[1, 2.], scale=[11, 22.])

  # Evaluate the pdf of the first distribution on 0, and the second on 1.5,
  # returning a length two tensor.
  dist.prob([0, 1.5])

  # Get 3 samples, returning a 3 x 2 tensor.
  dist.sample([3])
  ```

  Arguments are broadcast when possible.

  ```python
  # Define a batch of two scalar valued Logistics.
  # Both have mean 1, but different scales.
  dist = tfd.Gompertz(loc=1., scale=[11, 22.])

  # Evaluate the pdf of both distributions on the same point, 3.0,
  # returning a length 2 tensor.
  dist.prob(3.0)
  ```

  """

  def __init__(self,
               loc,
               scale,
               validate_args=False,
               allow_nan_stats=True,
               name='Gompertz'):
    """Construct Gompertz distributions with location and scale `loc` and `scale`.

    The parameters `loc` and `scale` must be shaped in a way that supports
    broadcasting (e.g. `loc + scale` is a valid operation).

    Args:
      loc: Floating point tensor, the means of the distribution(s).
      scale: Floating point tensor, the scales of the distribution(s).
        scale must contain only positive values.
      validate_args: Python `bool`, default `False`. When `True` distribution
        parameters are checked for validity despite possibly degrading runtime
        performance. When `False` invalid inputs may silently render incorrect
        outputs.
        Default value: `False`.
      allow_nan_stats: Python `bool`, default `True`. When `True`,
        statistics (e.g., mean, mode, variance) use the value `NaN` to
        indicate the result is undefined. When `False`, an exception is raised
        if one or more of the statistic's batch members are undefined.
        Default value: `True`.
      name: Python `str` name prefixed to Ops created by this class.
        Default value: `'Gompertz'`.

    Raises:
      TypeError: if loc and scale are different dtypes.
    """
    parameters = dict(locals())
    with tf.name_scope(name) as name:
      dtype = dtype_util.common_dtype([loc, scale], dtype_hint=tf.float32)
      loc = tensor_util.convert_nonref_to_tensor(
          loc, name='loc', dtype=dtype)
      scale = tensor_util.convert_nonref_to_tensor(
          scale, name='scale', dtype=dtype)
      dtype_util.assert_same_float_dtype([loc, scale])
      # Positive scale is asserted by the incorporated Gompertz bijector.
      self._gompertz_bijector = GompertzCDF(
          loc=loc, scale=scale, validate_args=validate_args)

      # Because the uniform sampler generates samples in `[0, 1)` this would
      # cause samples to lie in `(inf, -inf]` instead of `(inf, -inf)`. To fix
      # this, we use `np.finfo(dtype_util.as_numpy_dtype(self.dtype).tiny`
      # because it is the smallest, positive, 'normal' number.
      super(Gompertz, self).__init__(
          distribution=uniform.Uniform(
              low=np.finfo(dtype_util.as_numpy_dtype(dtype)).tiny,
              high=tf.ones([], dtype=dtype),
              allow_nan_stats=allow_nan_stats),
          # The Gompertz bijector encodes the CDF function as the forward,
          # and hence needs to be inverted.
          bijector=invert_bijector.Invert(
              self._gompertz_bijector, validate_args=validate_args),
          parameters=parameters,
          name=name)

  @classmethod
  def _parameter_properties(cls, dtype, num_classes=None):
    # pylint: disable=g-long-lambda
    return dict(
        loc=parameter_properties.ParameterProperties(),
        scale=parameter_properties.ParameterProperties(
            default_constraining_bijector_fn=(
                lambda: softplus_bijector.Softplus(low=dtype_util.eps(dtype)))))
    # pylint: enable=g-long-lambda

  @property
  def loc(self):
    """Distribution parameter for the location."""
    return self._gompertz_bijector.loc

  @property
  def scale(self):
    """Distribution parameter for scale."""
    return self._gompertz_bijector.scale

  experimental_is_sharded = False

  def _entropy(self): # TODO continue here
    # Use broadcasting rules to calculate the full broadcast sigma.
    scale = self.scale * tf.ones_like(self.loc)
    return 1. + tf.math.log(scale) + np.euler_gamma

  def _log_prob(self, x):
    scale = tf.convert_to_tensor(self.scale)
    z = (x - self.loc) / scale
    return (z - tf.exp(z)) + tf.math.log(scale)

  def _mean(self):
    return self.loc - self.scale * np.euler_gamma

  def _stddev(self):
    return self.scale * tf.ones_like(self.loc) * np.pi / np.sqrt(6)

  def _mode(self):
    return self.loc * tf.ones_like(self.scale)

  def _default_event_space_bijector(self):
    # TODO(b/145620027) Finalize choice of bijector. Consider switching to
    # Chain([Softplus(), Log()]) to lighten the doubly-exponential right tail.
    return identity_bijector.Identity(validate_args=self.validate_args)

  def _parameter_control_dependencies(self, is_init):
    return self._gompertz_bijector._parameter_control_dependencies(is_init)  # pylint: disable=protected-access


@kullback_leibler.RegisterKL(Gompertz, Gompertz)
def _kl_gompertz_gompertz(a, b, name=None):
  """Calculate the batched KL divergence KL(a || b) with a and b Gompertz.

  Args:
    a: instance of a Gompertz distribution object.
    b: instance of a Gompertz distribution object.
    name: (optional) Name to use for created operations.
      default is 'kl_gompertz_gompertz'.

  Returns:
    Batchwise KL(a || b)
  """
  with tf.name_scope(name or 'kl_gompertz_gompertz'):
    # Consistent with
    # http://www.mast.queensu.ca/~communications/Papers/gil-msc11.pdf, page 64
    # The paper uses beta to refer to scale and mu to refer to loc.
    # There is actually an error in the solution as printed; this is based on
    # the second-to-last step of the derivation. The value as printed would be
    # off by (a.loc - b.loc) / b.scale.
    a_loc = tf.convert_to_tensor(a.loc)
    b_loc = tf.convert_to_tensor(b.loc)
    a_scale = tf.convert_to_tensor(a.scale)
    b_scale = tf.convert_to_tensor(b.scale)
    return (tf.math.log(b_scale) - tf.math.log(a_scale) + np.euler_gamma *
            (a_scale / b_scale - 1.) +
            tf.math.expm1((b_loc - a_loc) / b_scale +
                          tf.math.lgamma(a_scale / b_scale + 1.)) +
            (a_loc - b_loc) / b_scale)

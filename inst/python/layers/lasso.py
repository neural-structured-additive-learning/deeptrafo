import tensorflow as tf
from tensorflow import keras
import tensorflow.keras.regularizers as reg

class group_lasso_pen(reg.Regularizer):

    def __init__(self, la):
        self.la = la

    def __call__(self, x):
        return self.la * tf.reduce_sum(tf.sqrt(tf.reduce_sum(tf.square(x), 1)))

class ProjectC(tf.keras.constraints.Constraint):

    def __init__(self, C, fac):
        self.C = C
        self.fac = fac

    def __call__(self, w):
        mean = tf.reduce_mean(w)
        tf.debugging.assert_equal(w.shape, self.C.shape)
        wnew = tf.divide((w - mean), self.C)
        return self.fac*wnew + (1-self.fac)*w

    def get_config(self):
        return {'C': self.C}

class TibLinearLassoConstraint(tf.keras.layers.Layer):
    def __init__(self, C, fac=0.1, num_outputs=1, la=0, name="const_tib_lasso"):
        super(TibLinearLassoConstraint, self).__init__()
        self.num_outputs = num_outputs
        self.la = la
        if self.num_outputs > 1:
            self.reg = group_lasso_pen(self.la)
        else:
            self.reg = tf.keras.regularizers.l2(self.la)
        self._name = name
        self.C = C
        self.fac = fac

    def build(self, input_shape):
        self.sqrtAlpha = self.add_weight(
            shape=(input_shape[-1], 1),
            initializer=tf.keras.initializers.RandomUniform(minval=1e-4, maxval=0.1, seed=None),
            trainable=True,
            constraint=tf.keras.constraints.NonNeg(),
            name="sqrtAlpha"
        )
        self.beta = self.add_weight(
            shape=(input_shape[-1], self.num_outputs),
            initializer="random_normal",
            trainable=True,
            constraint=ProjectC(self.C, self.fac),
            name="beta"
        )

    def call(self, input):
        u = self.beta/self.sqrtAlpha
        simplyConnectedResult = tf.math.multiply(input, tf.transpose(self.sqrtAlpha))
        output = tf.matmul(simplyConnectedResult, u)
        self.add_loss(tf.reduce_sum(tf.square(u) + tf.square(self.sqrtAlpha)))
        return output

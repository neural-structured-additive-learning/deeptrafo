import tensorflow as tf
import tensorflow_probability as tfp
from tensorflow import keras
import tensorflow.keras.regularizers as reg

class EvalBspTF(tf.keras.layers.Layer):
    def __init__(self, order, supp, **kwargs):
        super(EvalBspTF, self).__init__(**kwargs)
        self.order = order
        self.supp = tf.cast(supp, dtype="float32")

    def call(self, input):
        input = tf.cast(input, dtype="float32")
        y = tf.math.divide(tf.math.subtract(input, self.supp[0]), tf.math.subtract(self.supp[1], self.supp[0]))
        return(keras.layers.concatenate([tf.reshape(tf.divide(tfp.distributions.Beta(m + 1, self.order + 1 - m).prob(y), (self.order + 1)), (-1,1)) 
                                         for m in range(self.order+1)], axis = 1)
        )
        
def tf_nan_to_zero(x):

    mult = tf.where(tf.math.logical_not(tf.math.logical_or(tf.math.is_inf(x), tf.math.is_nan(x))), 1.0, 0.0)
    return(tf.math.multiply_no_nan(x, mult))

        
def calc_first_second(m, y, order):

    first_t = tf.divide(tfp.distributions.Beta(m, order - m + 1).prob(y), order)
    sec_t = tf.divide(tfp.distributions.Beta(m + 1, order - m).prob(y), order)
    
    return(tf.reshape(tf.multiply(tf.subtract(tf_nan_to_zero(first_t), tf_nan_to_zero(sec_t)), order), (-1,1)))
        
class EvalBspPrimeTF(tf.keras.layers.Layer):
    def __init__(self, order, supp, **kwargs):
        super(EvalBspPrimeTF, self).__init__(**kwargs)
        self.order = order
        self.supp = tf.cast(supp, dtype="float32")

    def call(self, input):
        input = tf.cast(input, dtype="float32")
        y = tf.math.divide(tf.math.subtract(input, self.supp[0]), tf.math.subtract(self.supp[1], self.supp[0]))
        return(keras.layers.concatenate([calc_first_second(m, y, self.order) for m in range(self.order+1)], axis = 1))
        

import FormFactor
import math

class ExpoFormFactor(FormFactor):
    """
    This is a simple form factor that follows an exponential form.

    |F(Q^2)|^2 = exp(- 2Q^2 / A^2)
    where A is the scale.

    """
    def __init__(self,scale):
        self.scale = scale

    def ff2(self,Q2):
      return math.exp(-2 * Q2 / (scale * scale))    

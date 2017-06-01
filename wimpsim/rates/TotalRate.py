import TotalRateEngine

class TotalRate:

  def __init__(self,intmodel,velocity):
    self.int_model = intmodel
    self.vel_dist = velocity
    self.engine = TotalRateEngine()
    self.rate = -1

  def calculate():
    self.rate = self.engine.calculate(int_model,vel_dist)
    return self.rate
  

import numpy as np
from Fluid import Fluid

N = 256
iter = 16
SCALE = 4
t = 0

fluid = Fluid(0.1, 0.01, 0.001, N)
fluid.addDensity(1, 1, 2)
fluid.addVelocity(1, 1 ,0.1, 0.1)
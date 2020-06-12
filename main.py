import numpy as np
from Fluid import Fluid
from config import N, iter, SCALE, t


fluid = Fluid(0.1, 0.01, 0.001, N)
fluid.addDensity(1, 1, 2)
fluid.addVelocity(1, 1 ,0.1, 0.1)
fluid.step()
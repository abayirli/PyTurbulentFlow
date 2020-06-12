import numpy as np
from utils import IX, diffuse, lin_solve, project, advect
from config import N, iter, SCALE, t


#Fluid cube class
class Fluid:
    def __init__(self,dt, diffusion, viscosity, N):
        self.size = N
        self.dt = dt
        self.diff = diffusion
        self.visc = viscosity

        self.s = np.zeros(N * N)
        self.density = np.zeros(N * N)

        self.Vx = np.zeros(N * N)
        self.Vy = np.zeros(N * N)

        self.Vx0 = np.zeros(N * N)
        self.Vy0 = np.zeros(N * N)

    
    #step method
    def step(self):
        N = self.size
        visc = self.visc
        diff = self.diff
        dt = self.dt
        Vx = self.Vx
        Vy = self.Vy
        Vx0 = self.Vx0
        Vy0 = self.Vy0
        s = self.s
        density = self.density

        diffuse(1, Vx0, Vx, visc, dt)
        diffuse(2, Vy0, Vy, visc, dt)

        project(Vx0, Vy0, Vx, Vy)

        advect(1, Vx, Vx0, Vx0, Vy0, dt)
        advect(2, Vy, Vy0, Vx0, Vy0, dt)

        project(Vx, Vy, Vx0, Vy0)
        diffuse(0, s, density, diff, dt)
        advect(0, density, s, Vx, Vy, dt)

    #method to add density
    def addDensity(self, x, y, amount):
        index = IX(x, y)
        self.density[index] += amount


  #method to add velocity
    def addVelocity(self, x, y, amountX, amountY):
        index = IX(x, y)
        self.Vx[index] += amountX
        self.Vy[index] += amountY

    def render():
        pass


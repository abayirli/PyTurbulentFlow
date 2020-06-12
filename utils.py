import numpy as np
from config import N, iter, SCALE, t
import math


def diffuse(b, x, x0, diff, dt):
    """Diffusion function

    Args:
        b (int): boundry layer
        x (numpy.ndarray): Current fluid state
        x0 (numpy.ndarray): Previous fluid state
        diff (float): diffusuion term
        dt (float): time increment
    """
    a = dt * diff * (N - 2) * (N - 2)
    lin_solve(b, x, x0, a, 1 + 6 * a)


def lin_solve(b, x, x0, a, c):
    """Function for solving the linear differential equation

    Args:
        b (int): boundy term
        x (numpy.ndarray): current fluid state
        x0 (numpy.ndarray): previous fluid state
        a (float): [description]
        c (float): [description]
    """
    cRecip = 1.0 / c
    for t in range(iter):
        for j in range(1, N-1):
            for i in range(1, N-1):
                x[i,j] =((x0[i,j] + a * (x[i + 1, j] + 
                x[i - 1, j] + x[i, j + 1] +
                x[i, j - 1])) * cRecip)
        set_bnd(b, x)


def project(velocX, velocY, p, div):
    """This operation runs through all the cells and fixes them up
    so everything is in equilibrium

    Args:
        velocX (numpy.ndarray): current velocity field x
        velocY (numpy.ndarray): current velocity field y
        p ([numpy.nadarray]): [description]
        div ([numpy.ndarray]): [description]
    """
    for j in range(1, N-1):
        for i in range(1, N-1):
            div[i, j] = ((-0.5 *
                            (velocX[i + 1, j] -
                                velocX[i - 1, j] +
                                velocY[i, j + 1] -
                                velocY[i, j - 1])) /
                            N)
            p[i, j] = 0

    set_bnd(0, div)
    set_bnd(0, p)
    lin_solve(0, p, div, 1, 6)

    for j in range(1, N-1):
        for i in range(1, N-1):
            velocX[i, j] -= 0.5 * (p[i + 1, j] - p[i - 1, j]) * N
            velocY[i, j] -= 0.5 * (p[i, j + 1] - p[i, j - 1]) * N


    set_bnd(1, velocX)
    set_bnd(2, velocY)



def advect(b, d, d0, velocX, velocY, dt):
    """Function of advect: responsible for actually moving things around

    Args:
        b (int): boundry value
        d (numpy.ndarray): current fluid state
        d0 (numpy.ndarray): previous fluid state
        velocX (numpy.ndarray): current velocity field x
        velocY (numpy.ndarray): current velocity field y
        dt (float): time increment
    """
    
    dtx = dt * (N - 2)
    dty = dt * (N - 2)

    Nfloat = N

    for j, jfloat in zip(range(1, N-1), range(1, N-1)):
        for i, ifloat in zip(range(1, N-1), range(1, N-1)):
            tmp1 = dtx * velocX[i, j]
            tmp2 = dty * velocY[i, j]
            x = ifloat - tmp1
            y = jfloat - tmp2

            if (x < 0.5): x = 0.5
            if (x > Nfloat + 0.5): x = Nfloat + 0.5
            i0 = math.floor(x)
            i1 = i0 + 1.0
            if (y < 0.5): y = 0.5
            if (y > Nfloat + 0.5): y = Nfloat + 0.5
            j0 = math.floor(y)
            j1 = j0 + 1.0

            s1 = x - i0
            s0 = 1.0 - s1
            t1 = y - j0
            t0 = 1.0 - t1

            i0i = int(i0)
            i1i = int(i1)
            j0i = int(j0)
            j1i = int(j1)
            if(i0i < N and i1i < N and j0i < N and j1i < N):
                d[i, j] = (s0 * (t0 * d0[i0i, j0i] + t1 * d0[i0i, j1i]) +
                            s1 * (t0 * d0[i1i, j0i] + t1 * d0[i1i, j1i]))

    set_bnd(b, d)


def set_bnd(b, x):
    """Function of dealing with situation with boundry cells

    Args:
        b (int): boundy vaue
        x (numpy.ndarray): fluid state
    """

    for i in range(N-1):
        x[i, 0] = -x[i, 1] if b == 2 else x[i, 1]
        x[i, N - 1] = -x[i, N - 2] if b == 2 else x[i, N - 2]
    for j in range(N-1):
        x[0, j] = -x[1, j] if b == 1 else x[1, j]
        x[N - 1, j] = -x[N - 2, j] if b == 1 else x[N - 2, j]

    x[0, 0] = 0.5 * (x[1, 0] + x[0, 1])
    x[0, N - 1] = 0.5 * (x[1, N - 1] + x[0, N - 2])
    x[N - 1, 0] = 0.5 * (x[N - 2, 0] + x[N - 1, 1])
    x[N - 1, N - 1] = 0.5 * (x[N - 2, N - 1] + x[N - 1, N - 2])
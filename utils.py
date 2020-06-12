from config import N, iter, SCALE, t
import math

#function to use 1D array an`d fake the extra two dimensions --> 3D
def IX(x, y):
    return x + y * N

# Function of diffuse
# - b : int
# - x : float[]
# - x0 : float[]
# - diff : float
# - dt : flaot

def diffuse(b, x, x0, diff, dt):
    a = dt * diff * (N - 2) * (N - 2)
    lin_solve(b, x, x0, a, 1 + 6 * a)



# Function of solving linear differential equation
# - b : int
# - x : float[]
# - x0 : float[]
# - a : float
# - c : float

def lin_solve(b, x, x0, a, c):
    cRecip = 1.0 / c
    for t in range(iter):
        for j in range(1, N-1):
            for i in range(1, N-1):
                x[IX(i, j )] =((x0[IX(i, j)] + a * (x[IX(i + 1, j)] + 
                x[IX(i - 1, j)] + x[IX(i, j + 1)] +
                x[IX(i, j - 1)])) * cRecip)
        set_bnd(b, x)


# Function of project : This operation runs through all the cells and fixes them up so everything is in equilibrium.
# - velocX : float[]
# - velocY : float[]
# - p : float[]
# - div : float[]

def project(velocX, velocY, p, div):
    for j in range(1, N-1):
        for i in range(1, N-1):
            div[IX(i, j)] = ((-0.5 *
                            (velocX[IX(i + 1, j)] -
                                velocX[IX(i - 1, j)] +
                                velocY[IX(i, j + 1)] -
                                velocY[IX(i, j - 1)])) /
                            N)
            p[IX(i, j)] = 0

    set_bnd(0, div)
    set_bnd(0, p)
    lin_solve(0, p, div, 1, 6)

    for j in range(1, N-1):
        for i in range(1, N-1):
            velocX[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * N
            velocY[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * N


    set_bnd(1, velocX)
    set_bnd(2, velocY)


# Function of advect: responsible for actually moving things around
# - b : int
# - d : float[]
# - d0 : float[]
# - velocX : float[]
# - velocY : float[]
# - velocZ : float[]
# - dt : float[]

def advect(b, d, d0, velocX, velocY, dt):

    dtx = dt * (N - 2)
    dty = dt * (N - 2)

    Nfloat = N

    for j, jfloat in zip(range(1, N-1), range(1, N-1)):
        for i, ifloat in zip(range(1, N-1), range(1, N-1)):
            tmp1 = dtx * velocX[IX(i, j)]
            tmp2 = dty * velocY[IX(i, j)]
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

            d[IX(i, j)] = (s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)]) +
                           s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]))

    set_bnd(b, d)


# Function of dealing with situation with boundary cells.
# - b : int
# - x : float[]
# variable = (condition) ? expressionTrue :  expressionFalse;

def set_bnd(b, x):
    for i in range(N-1):
        x[IX(i, 0)] = -x[IX(i, 1)] if b == 2 else x[IX(i, 1)]
        x[IX(i, N - 1)] = -x[IX(i, N - 2)] if b == 2 else x[IX(i, N - 2)]
    for j in range(N-1):
        x[IX(0, j)] = -x[IX(1, j)] if b == 1 else x[IX(1, j)]
        x[IX(N - 1, j)] = -x[IX(N - 2, j)] if b == 1 else x[IX(N - 2, j)]

    x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)])
    x[IX(0, N - 1)] = 0.5 * (x[IX(1, N - 1)] + x[IX(0, N - 2)])
    x[IX(N - 1, 0)] = 0.5 * (x[IX(N - 2, 0)] + x[IX(N - 1, 1)])
    x[IX(N - 1, N - 1)] = 0.5 * (x[IX(N - 2, N - 1)] + x[IX(N - 1, N - 2)])
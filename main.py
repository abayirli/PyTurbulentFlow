import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from tqdm import tqdm
from Fluid import Fluid
from config import N, iter, SCALE, t, step, output_type
import time

def iterate(fluid, t):
    """Iterate over the fluid object, creating new 
    random density field around the origin and renders
    the full state and returns it.

    Args:
        fluid (Fluid): Fluid object to iterate over
        t (float): time

    Returns:
        numpy.ndarray: NXN image representation of the fluid state
    """
    cx = int((0.5 * N))
    cy = int((0.5 * N))
    for i in range(-1,2):
        for j in range(-1,2):
            fluid.addDensity(cx + i, cy + j, np.random.randint(250, 400))

    for i in range(0,2):
        vx = np.random.uniform(-3,3)
        vy = np.random.uniform(-3,3)      
        t += 0.01
        fluid.addVelocity(cx, cy, vx, vy)

    fluid.step()
    image = fluid.render()
    im = plt.imshow(image, cmap='YlOrRd', interpolation = "nearest", animated=True)
    return im

if __name__ == "__main__":
    # Define the fluid object with a given time resolution
    # diffusion and viscosity
    fluid = Fluid(0.1, 0, 0.0000005)

    fig = plt.figure()

    ims = []
    for i in tqdm(range(step)):
        im = iterate(fluid, t)
        ims.append([im])

    ani = animation.ArtistAnimation(fig, ims, interval=100, blit=True,
                                    repeat_delay=1000)

    if(output_type == "mp4"):
        ani.save(f'./movies/movie_dif_vis_steps{step}.mp4')
    elif(output_type == "gif"):  
        ani.save(f'./movies/movie_dif_vis_steps{step}_.gif', dpi=80, writer='imagemagick')

    plt.show()

    print('Done!')
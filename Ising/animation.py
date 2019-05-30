import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

N = 0
files = []
while True:
    name = 'results/%d.ising' % N
    if not os.path.exists(name): break
    files.append(name)
    N += 1
    # if N >= 100:
    #     break

d0 = np.genfromtxt(files[0])
data = np.zeros((N, *d0.shape), dtype = int)

data[0] = d0

for i, file in enumerate(files[1:]):
    file =  np.genfromtxt(file)
    data[i + 1] = file
    del file

def animate(i):
    global data, graph
    graph.set_data(data[i])
    return graph,

fig, ax = plt.subplots(figsize = (8, 8))

ax.set_axis_off()

graph = ax.imshow(d0, animated = True, cmap = 'Greys')
ani = FuncAnimation(fig, animate, interval = 1, frames = range(N))

ani.save('animation.mp4', fps = N // 10)
plt.show()

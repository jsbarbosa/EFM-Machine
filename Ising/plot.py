import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('results.dat')

fig, axes = plt.subplots(2, 2, sharex = True)

for i in range(2):
    axes[0, i].plot(data[:, 0], data[:, i + 1])
plt.show()

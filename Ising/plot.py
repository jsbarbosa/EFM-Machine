import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('results.dat')

labels = [r'$\bar{E} / N$', r'$\bar{M} / N$', r'$Cv$', r'$\chi$']

fig, axes = plt.subplots(2, 2, sharex = True)

axes = axes.ravel()
for i in range(4):
    axes[i].scatter(data[:, 0], data[:, i + 1], s = 1, alpha = 0.5)
    axes[i].axvline(x = 1, c = 'k', alpha = 0.5)
    axes[i].set_ylabel(labels[i])

axes[2].set_yscale('log')
axes[3].set_yscale('log')

axes[2].set_xlabel('$T / T_c$')
axes[3].set_xlabel('$T / T_c$')

# fig.tight_layout()
fig.subplots_adjust(hspace = 0, wspace = 0.3)

fig.savefig('properties.png', dpi = 300, bbox_inches = 'tight')

plt.show()

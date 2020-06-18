import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from rdfpy import rdf3d

particles = np.load('./random3d.npy')

g_r, radii = rdf3d(particles, p_radius=0.8, h=30, w=30, d=30, dr=0.125)

fig, axes = plt.subplots(1, 2, figsize=(12, 4), gridspec_kw={'width_ratios': [1, 2]})

axes[0].scatter(particles[:, 1], particles[:, 2], color='r', alpha=0.5, edgecolors='k', s=10)
axes[0].axis('off')

axes[1].plot(radii, g_r, color='k', alpha=0.75)
axes[1].hlines(y=1.0, xmin=0.0, xmax=max(radii), color='r', linestyle='--', alpha=0.4)
axes[1].set_ylabel('g(r)')
axes[1].set_xlabel('r')
axes[1].set_xlim(0.0, max(radii))
plt.suptitle('rdf3d')
plt.savefig('./random3d.png', bbox_inches='tight', pad_inches=0.25)
plt.show()

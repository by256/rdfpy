import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
from rdfpy import rdf2d

particles = np.load('./random2d.npy')

g_r, radii = rdf2d(particles, p_radius=0.25, h=30, w=30, dr=0.125)


fig, axes = plt.subplots(1, 2, figsize=(12, 4), gridspec_kw={'width_ratios': [1, 2]})

axes[0].scatter(particles[:, 0], particles[:, 1], color='r', alpha=0.8, edgecolors='k', s=10)
axes[0].axis('off')

axes[1].plot(radii, g_r, color='k', alpha=0.75)
axes[1].hlines(y=1.0, xmin=0.0, xmax=max(radii), color='r', linestyle='--', alpha=0.4)
axes[1].set_ylabel('g(r)')
axes[1].set_xlabel('r')
axes[1].set_xlim(0.0, max(radii))
plt.savefig('./random2d.png', bbox_inches='tight', pad_inches=0.25)
# plt.show()
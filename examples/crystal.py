"""Plot the particles and RDF of a crystal structure.

Requires matplotlib: ``pip install rdfpy[examples]``.
"""
import os

import numpy as np
import matplotlib.pyplot as plt

from rdfpy import rdf

HERE = os.path.dirname(os.path.abspath(__file__))

particles = np.load(os.path.join(HERE, 'crystal.npy'))

fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(particles[:, 0], particles[:, 1], particles[:, 2], color='r',
           alpha=0.4, edgecolors='k', s=25)
plt.axis('off')
plt.savefig(os.path.join(HERE, 'crystal-particles.png'), bbox_inches='tight',
            pad_inches=0.0)
plt.close()

g_r, radii = rdf(particles, dr=0.1)

image = plt.imread(os.path.join(HERE, 'crystal-particles.png'))

fig, axes = plt.subplots(1, 2, figsize=(12, 4),
                         gridspec_kw={'width_ratios': [1, 2]})

axes[0].imshow(image)
axes[0].axis('off')
axes[0].set_title('Crystal')

axes[1].plot(radii, g_r, color='k', alpha=0.75)
axes[1].hlines(y=1.0, xmin=0.0, xmax=max(radii), color='r', linestyle='--',
               alpha=0.4)
axes[1].set_ylabel('g(r)')
axes[1].set_xlabel('r')
axes[1].set_xlim(0.0, max(radii))
plt.savefig(os.path.join(HERE, 'crystal.png'), bbox_inches='tight',
            pad_inches=0.25)
plt.show()

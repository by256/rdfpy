import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from rdfpy import rdf3d

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

# particle positions were obtained from http://www.ergoscf.org/xyz/h2o.php
particles = np.load('./water.npy')

fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(particles[:, 0], particles[:, 1], particles[:, 2], color='aqua', alpha=0.3, edgecolors='k', s=25)
plt.axis('off')
plt.savefig('./water-particles.png', bbox_inches='tight', pad_inches=0.0)
plt.close()

g_r, radii = rdf3d(particles, dr=0.2)

image = plt.imread('./water-particles.png')

fig, axes = plt.subplots(1, 2, figsize=(12, 4), gridspec_kw={'width_ratios': [1, 2]})

axes[0].imshow(image)
axes[0].axis('off')

axes[1].plot(radii, g_r, color='k', alpha=0.75)
axes[1].hlines(y=1.0, xmin=0.0, xmax=max(radii), color='r', linestyle='--', alpha=0.4)
axes[1].set_ylabel(r'g(r)')
axes[1].set_xlabel(r'r')
axes[1].set_xlim(0.0, max(radii))
axes[0].set_title(r'H$_{2}$O')
plt.savefig('./water.png', bbox_inches='tight', pad_inches=0.25)
plt.show()

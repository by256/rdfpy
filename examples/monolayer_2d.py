"""2D example: RDF of a hexagonal colloidal monolayer.

Unlike the other examples, this one generates its own coordinates rather than
loading a .npy file, so you can change the lattice, the disorder, or the box
size and immediately see the effect on g(r).

Requires matplotlib: ``pip install rdfpy[examples]``.
"""
import os

import numpy as np
import matplotlib.pyplot as plt

from rdfpy import rdf

HERE = os.path.dirname(os.path.abspath(__file__))

A = 1.0          # nearest-neighbour spacing
JITTER = 0.06    # uniform displacement about each site, models thermal motion


def hexagonal_monolayer(n_x=60, n_y=70, a=A, jitter=JITTER, seed=0):
    """Triangular (hexagonal) lattice in 2D.

    Rows are spaced a*sqrt(3)/2 apart and alternate rows are offset by a/2, so
    every site has 6 nearest neighbours at distance a.
    """
    i, j = np.meshgrid(np.arange(n_x), np.arange(n_y), indexing='ij')
    x = (i + 0.5 * (j % 2)) * a
    y = j * (np.sqrt(3) / 2) * a
    particles = np.stack([x.ravel(), y.ravel()], axis=1)
    if jitter:
        rng = np.random.default_rng(seed)
        particles = particles + rng.uniform(-jitter, jitter, particles.shape)
    return particles


particles = hexagonal_monolayer()

# A hexagonal lattice has a known areal density of 2/(sqrt(3)*a^2). Passing it
# explicitly is worth doing here: inferring rho from the bounding box of the
# particles overestimates it by ~2% for this configuration, because the box
# corners sit on the outermost particles rather than half a lattice spacing
# beyond them.
rho = 2 / (np.sqrt(3) * A**2)

g_r, radii = rdf(particles, dr=0.02, rho=rho)

# Neighbour shells of a triangular lattice, in units of a, with their
# coordination numbers.
shells = [(1.0, 6), (np.sqrt(3), 6), (2.0, 6), (np.sqrt(7), 12), (3.0, 6)]

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5),
                         gridspec_kw={'width_ratios': [1, 2]})

window = particles[:, 0] < 15 * A
axes[0].scatter(particles[window, 0], particles[window, 1], s=14,
                color='tab:blue', alpha=0.8, edgecolors='k', linewidths=0.3)
axes[0].set_aspect('equal')
axes[0].set_xlim(0, 15 * A)
axes[0].set_ylim(0, 15 * A)
axes[0].axis('off')
axes[0].set_title('Hexagonal monolayer (detail)')

axes[1].plot(radii, g_r, color='k', alpha=0.8, lw=1.2)
axes[1].hlines(y=1.0, xmin=0.0, xmax=max(radii), color='r', linestyle='--',
               alpha=0.4)
for distance, coordination in shells:
    axes[1].axvline(distance, color='tab:blue', linestyle=':', alpha=0.5)
    axes[1].annotate('{}'.format(coordination), xy=(distance, 0.5),
                     ha='center', fontsize=9, color='tab:blue')
axes[1].set_ylabel('g(r)')
axes[1].set_xlabel('r')
axes[1].set_xlim(0.0, 6.0)
axes[1].set_title('Dotted lines: exact shell radii, labelled with coordination '
                  'number')

plt.tight_layout()
plt.savefig(os.path.join(HERE, 'monolayer_2d.png'), bbox_inches='tight',
            pad_inches=0.25)
plt.show()

# Integrating rho * g(r) * shell_area over each peak recovers the coordination
# numbers to within about 1%, which is a good way to check an RDF against a
# known structure.
dr = radii[1] - radii[0]
shell_area = np.pi * ((radii + dr / 2)**2 - (radii - dr / 2)**2)
print('shell      integrated    exact')
for distance, coordination in shells:
    mask = np.abs(radii - distance) < 0.15
    integrated = np.nansum((rho * g_r * shell_area)[mask])
    print('r = {:.3f}   {:9.3f}   {:6d}'.format(distance, integrated,
                                                coordination))

import numpy as np


def rdf2d(particles, p_radius, h, w, n_radii=None, dr=None):
    centre_y, centre_x, centre_z = h/2, w/2, d/2
    min_x, max_x = np.min(particles[:, 0]), np.max(particles[:, 0])
    min_y, max_y = np.min(particles[:, 1]), np.max(particles[:, 1])

    r_max = (np.min([h, w]) / 2) - 2*p_radius

    if n_radii and not dr:
        radii = np.linspace(0.0, r_max, n_radii)
        dr = radii[1] - radii[0]
        radii = radii + dr
    elif dr and not n_radii:
        radii = np.arange(dr, r_max, dr)
    elif n_radii is not None and dr is not None:
        raise ValueError('Arguments n_radii and dr should not both be specified.')
    
    particles_individual_counts = np.zeros(shape=(len(particles), len(radii)))
    
    N = len(particles)
    number_density = N / (h*w)
    
    for r_idx, r in enumerate(radii):
        valid_idxs = ((particles[:, 0] - (r-p_radius) >= min_x) & (particles[:, 0] + (r-p_radius) <= max_x)) \
                  & (((particles[:, 1] - (r-p_radius) >= min_y) & (particles[:, 1] + (r-p_radius) <= max_y)))
        valid_particles = particles[valid_idxs]
        n_valid = len(valid_particles)
        
        for p_idx, particle in enumerate(valid_particles):
            distances = np.linalg.norm(particle - particles, axis=1)
            n = len(distances[(distances >= r) & (distances < r+dr)])
            particles_individual_counts[p_idx, r_idx] = n

        scaling_factor = 1 / (n_valid*2*np.pi*r*dr)
        particles_individual_counts[:, r_idx] = particles_individual_counts[:, r_idx] * scaling_factor
        print('{}/{}'.format(r_idx+1, len(radii)), end='\r', flush=True)

    particles_individual_counts = particles_individual_counts / number_density
    
    g_r = np.sum(particles_individual_counts, axis=0)
    return g_r, radii, r_max
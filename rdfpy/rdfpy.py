import time
import numpy as np
from scipy.spatial import cKDTree


def rdf2d(particles, dr, rho=None, eps=1e-15):
    """
    Computes 2D radial disribution functions g(r) of a set of particles of shape (N, 2).
    Particle must be placed in a 2D box of dimensions (width x height).
    
    Parameters
    ----------
    particles : (N, 2) np.array
        Set of particle from which to compute the radial distribution funntion g(r).
    dr : float
        Delta r. Determines the spacing between successive radii over which g(r) is computed.
    rho : float, optional
        Number density. In cases where the number of particles is relatively low, or the 
        particle:box size ratio is small, it may be better to provide your own number density 
        for more accurate results. If left as None, box dimensions will be inferred from the 
        particles and the number density will be calculated accordingly.
    eps : float, optional
        Epsilon value used to find particles less than or equal to a distance in KDTree.
        
    
    Returns
    -------
    g_r : (n_radii) np.array
        radial distribution function values g(r).
    radii : (n_radii) np.array
        radii over which g(r) is computed
    """
    
    start = time.time()
    
    # translate particles such that the particle with min coords is at origin
    particles = particles - np.min(particles, axis=0)
    min_x, min_y = np.min(particles, axis=0)
    max_x, max_y = np.max(particles, axis=0)
    
    # dimensions of box
    w, h = (max_x-min_x), (max_y-min_y)

    r_max = (np.min([w, h]) / 2)*0.8
    radii = np.arange(dr, r_max, dr)
    g_r = np.zeros(shape=(len(radii)))
    
    N = len(particles)
    if not rho:
        rho = N / (w*h) # number density

    # create a KDTree for fast nearest-neighbour lookup of particles
    tree = cKDTree(particles)
    
    for r_idx, r in enumerate(radii):
        valid_idxs = (particles[:, 0]-(r+dr) >= min_x) & (particles[:, 0]+(r+dr) <= max_x) \
                   & (particles[:, 1]-(r+dr) >= min_y) & (particles[:, 1]+(r+dr) <= max_y)
        valid_particles = particles[valid_idxs]
        
        for particle in valid_particles:
            n = tree.query_ball_point(particle, r+dr-eps, return_length=True) - tree.query_ball_point(particle, r, return_length=True)
            g_r[r_idx] += n
        
        n_valid = len(valid_particles)
        shell_vol = np.pi*((r+dr)**2 - r**2)
        g_r[r_idx] /= n_valid*shell_vol*rho

        print('Computing RDF     Radius {}/{}    Time elapsed: {:.3f} s'.format(r_idx+1, len(radii), time.time()-start), end='\r', flush=True)
        
    return g_r, radii


def rdf3d(particles, dr, rho=None, eps=1e-15):
    """
    Computes 3D radial disribution functions g(r) of a set of particles of shape (N, 3).
    Particle must be placed in a 3D cuboidal box of dimensions (width x height x depth).
    
    Parameters
    ----------
    particles : (N, 3) np.array
        Set of particle from which to compute the radial distribution funntion g(r).
    dr : float
        Delta r. Determines the spacing between successive radii over which g(r) is computed.
    rho : float, optional
        Number density. In cases where the number of particles is relatively low, or the 
        particle:box size ratio is small, it may be better to provide your own number density 
        for more accurate results. If left as None, box dimensions will be inferred from the 
        particles and the number density will be calculated accordingly.
    eps : float, optional
        Epsilon value used to find particles less than or equal to a distance in KDTree.
        
    
    Returns
    -------
    g_r : (n_radii) np.array
        radial distribution function values g(r).
    radii : (n_radii) np.array
        radii over which g(r) is computed
    """
    
    start = time.time()
    
    # translate particles such that the particle with min coords is at origin
    particles = particles - np.min(particles, axis=0)
    min_x, min_y, min_z = np.min(particles, axis=0)
    max_x, max_y, max_z = np.max(particles, axis=0)
    
    # dimensions of box
    w, h, d = (max_x-min_x), (max_y-min_y), (max_z-min_z)

    r_max = (np.min([w, h, d]) / 2)*0.8
    radii = np.arange(dr, r_max, dr)
    g_r = np.zeros(shape=(len(radii)))
    
    N = len(particles)
    if not rho:
        rho = N / (w*h*d) # number density

    # create a KDTree for fast nearest-neighbour lookup of particles
    tree = cKDTree(particles)
    
    for r_idx, r in enumerate(radii):
        valid_idxs = (particles[:, 0]-(r+dr) >= min_x) & (particles[:, 0]+(r+dr) <= max_x) \
                   & (particles[:, 1]-(r+dr) >= min_y) & (particles[:, 1]+(r+dr) <= max_y) \
                   & (particles[:, 2]-(r+dr) >= min_z) & (particles[:, 2]+(r+dr) <= max_z)
        valid_particles = particles[valid_idxs]
        
        for particle in valid_particles:
            n = tree.query_ball_point(particle, r+dr-eps, return_length=True) - tree.query_ball_point(particle, r, return_length=True)
            g_r[r_idx] += n
        
        n_valid = len(valid_particles)
        shell_vol = (4/3)*np.pi*((r+dr)**3 - r**3)
        g_r[r_idx] /= n_valid*shell_vol*rho

        print('Computing RDF     Radius {}/{}    Time elapsed: {:.3f} s'.format(r_idx+1, len(radii), time.time()-start), end='\r', flush=True)
        
    return g_r, radii

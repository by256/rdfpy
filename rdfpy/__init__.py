import os
if os.name == 'nt':  # check for windows
    from .rdfpy import rdf2d, rdf3d
else:
    from .rdfpy import rdf2d_parallel as rdf2d
    from .rdfpy import rdf3d_parallel as rdf3d

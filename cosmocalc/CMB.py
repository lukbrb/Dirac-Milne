import numpy as np
import healpy as hp
from .ISW import *


def create_map(nb_struct, theta, ISW, n_max, n_min=0, NSIDE=2048):
    """ Fonction qui crée une carte à partir de l'effet ISW calculé pour toutes les structures d'un univers de Dirac-Milne.
    Arguments : 
        - m : Numpy array de taille NPIX (voir healpy.nside2npix) dans lequel les informations sur la cartes seront entrées
        - nb_struct : Numpy array contenant le nombre de structures par tranches de redshift
        - theta : Numpy array contenant les angles sous lesquels sont perçus les strcutures à redshift donné
        - ISW : Numpy array contenant l'effet ISW moyen prdouit par une structure à un certain redshift
        - n_max : Nombre de tranches de redhsift à comptabiliser
        - n_min : Numéro de la tranche à partir duquel le comptage s'effectue (n_min=0 par défaut)
        - NSIDE : Donne la résolution de la carte (NSIDE=2048 par défaut) 
        """
    NPIX = hp.nside2npix(NSIDE)
    m = np.zeros(NPIX)
    for i in range(n_min, n_max):
        # radius-ISW is the angular radius for a "SDSS" peak-size structure. We use it here as an approximation of the main structures
        radius_ISW = np.radians(theta[i])
        # The parameter ampl_ISW takes into account the potential at the edge of the depletion zone
        # and the variation of this potential when the photon beam leaves it
        
        ampl_ISW = ISW[i] 
        costheta_min = np.cos(radius_ISW)
        nb_clust_slice = nb_struct[i]
        for j in range(int(nb_clust_slice)):
        # Draw random unit vector for each cluster, defined by phi_clust and theta_clust, uniformly on sphere
            vec = sample_sphere(1).T[0]
            # Get the pixels in the disk of the depletion zone around the direction of the center of the cluster
            ipix_disc = hp.query_disc(nside=NSIDE, vec=vec, radius=radius_ISW)
            vect_pix = hp.pix2vec(NSIDE, ipix_disc)
            costheta_pix = np.dot(vec, vect_pix)
            # Get the ISW weights for the pixels of this disk
            # Approximations used: spherical depletion zone and parallel photon beam crossing the structure
            ang_fact_ISW = 1./radius_ISW * np.sqrt(2.*(costheta_pix - costheta_min))
            # Add this to the previous state of the map in ordrer to get the composite CMB image
        
            m[ipix_disc] += ampl_ISW * ang_fact_ISW  # On pondère par le nouveau potentiel
    return m 
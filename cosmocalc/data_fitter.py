""" Script permettant de trouver les paramètres donnant la meilleure concordance."""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

DOSSIER = "/Users/lucas-copernion/Google Drive/Dirac-Milne/donnees/"
# Matter power spectrum as derived mostly from SDSS and Planck
kpeak_pk_SDSS = np.loadtxt(np.DataSource().open(f"{DOSSIER}kpeak_pk_SDSS.txt"))
k, p_k = kpeak_pk_SDSS.T
pk_fi = interp1d(k, p_k)
size_random = 1000000
k_peak = 4.E-3
N_BINS = len(p_k) * 10
#mass_sample = k_peak * np.random.lognormal(mean=2.4, sigma=0.9, size=size_random)
#popul, bin_edges = np.histogram(mass_sample, bins=10*len(p_k))
k_tab = np.linspace(k[0], k[-1], N_BINS)
norm_fact = 4.0
pk_interp = pk_fi(k_tab)
with open(f"{DOSSIER}best_params.txt", 'w') as f:
    for moy in np.arange(0.01, 10.01, 0.01):
        for sig in np.arange(0.01, 5.01, 0.01):
            masse = k_peak * np.random.lognormal(mean=moy, sigma=sig, size=size_random)
            pop, bin = np.histogram(masse, bins=N_BINS)
            norm_fact = 4.0
            diff = abs(pop - norm_fact * pk_interp)
            f.write(f"{moy}\t{sig}\t{diff.mean()}\n")
"""
plt.plot(bin_edges[0:-1], popul)
plt.plot(k, norm_fact * p_k)
plt.yscale('log')
plt.xscale('log')
plt.show()

plt.plot(bin_edges[0:-1], popul)
plt.plot(k, norm_fact * p_k)
plt.yscale('linear')
plt.xscale('linear')
plt.xlim(0., 1.)
plt.show()"""

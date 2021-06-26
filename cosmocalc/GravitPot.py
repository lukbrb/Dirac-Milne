import numpy as np
import matplotlib.pyplot as plt


def potentiel(x,y, n_max=200):
    grav_pot = 0.
    for i in np.arange(-n_max, n_max+1):
        for j in np.arange(-n_max, n_max+1):
                grav_pot += 1./np.sqrt((2*i-x)**2 + (2*j-y)**2 + 1) - 1./np.sqrt((2*i)**2 + (2*j)**2 + 1)
    return grav_pot

# n_size = 64
# n_mid = int(n_size/2)
# pot = np.zeros((n_size +1 , n_size + 1))

# for i in np.linspace(0, n_size, num= n_size + 1, dtype=int):
#     for j in np.linspace(0, n_size, num= n_size + 1, dtype=int):
#         pot[i, j] = potentiel((i-n_mid)/n_mid,(j-n_mid)/n_mid)

# print(pot)

# fig, ax  = plt.subplots(figsize=(8, 8))
# ax.imshow(pot - pot[0,0])n
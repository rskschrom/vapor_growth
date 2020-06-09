import numpy as np
import matplotlib.pyplot as plt

# read data and plot
a = np.genfromtxt('a.txt')
c = np.genfromtxt('c.txt')
data = np.genfromtxt('ndist.txt')
mass = np.genfromtxt('mass.txt')
print(data.shape)

rhoeff = mass/(4./3.*np.pi*a**2.*c)

#plt.imshow(rhoeff[::10,:], cmap='Spectral_r')
plt.imshow((a/c)[::10,:], cmap='Spectral_r')
plt.colorbar()
plt.savefig('ndist.png')

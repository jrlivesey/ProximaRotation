import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from strength_criterion import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


# Relevant parameters
k2 = 0.3        # dynamic Love number of order 2
tidalQ = 100    # tidal quality factor
k2lag = 66      # k2 * tidal time lag in seconds
mm = 6.61e-6    # mean motion in 1/s (orbital period = 11 days)
mass1 = 4.1e4   # Proxima mass in Earth masses
mass2 = 3.80    # Proxima b mass in Earth masses
semi = 0.05     # semi-major axis in au


"""
CTL model figure
"""

eccs = np.linspace(1e-10, 0.5, 5000)
crit = ecc_of_capture(1e-8, k2lag, mm, mass1, mass2, semi)  # Calculate the critical eccentricity for a planet with particular C_22

fig, ax = plt.subplots(1, 1, figsize=(8, 6))

# Critical C_22(e) for p = 1, 3/2, and 2
ax.semilogy(eccs, np.abs(proxima_CritProlateness(1, eccs, k2lag, mm, mass1, mass2, semi)), color='r', label='$p=1$')
ax.semilogy(eccs, np.abs(proxima_CritProlateness(1.5, eccs, k2lag, mm, mass1, mass2, semi)), color='orange', label='$p=3/2$')
ax.semilogy(eccs, np.abs(proxima_CritProlateness(2, eccs, k2lag, mm, mass1, mass2, semi)), color='b', label='$p=2$')

# Lines showing a chosen value of C_22 and the corresponding critical eccentricity
ax.axhline(1e-8, linestyle='--', color='gray')
ax.axvline(crit, linestyle='--', color='gray')

# Axes parameters
ax.set_xlabel('Eccentricity', fontsize=14)
ax.set_ylabel('Critical $C_{22}$', fontsize=14)
ax.set_xlim(0, 0.5)
ax.set_ylim(1e-13, 1e-3)
ax.legend(loc=4, fontsize=12);

fig.savefig('ctl_example.pdf')


"""
CPL model figure
"""

crit = cpl_ecc_of_capture(1e-8, k2, tidalQ, mass1, mass2, semi) # Calculate the critical eccentricity for a planet with particular C_22

fig, ax = plt.subplots(1, 1, figsize=(8, 6))

# Critical C_22(e) for p = 1, 3/2, and 2
ax.semilogy(eccs, np.abs(cpl_proxima_CritProlateness(1, eccs, k2, tidalQ, mass1, mass2, semi)), color='r', label='$p=1$')
ax.semilogy(eccs, np.abs(cpl_proxima_CritProlateness(1.5, eccs, k2, tidalQ, mass1, mass2, semi)), color='orange', label='$p=3/2$')
ax.semilogy(eccs, np.abs(cpl_proxima_CritProlateness(2, eccs, k2, tidalQ, mass1, mass2, semi)), color='b', label='$p=2$')

# Lines showing a chosen value of C_22 and the corresponding critical eccentricity
ax.axhline(1e-8, linestyle='--', color='gray')
ax.axvline(crit, linestyle='--', color='gray')

# Axes parameters
ax.set_xlabel('Eccentricity', fontsize=14)
ax.set_ylabel('Critical $C_{22}$', fontsize=14)
ax.set_xlim(0, 0.5)
ax.set_ylim(1e-13, 1e-3)
ax.legend(loc=4, fontsize=12);

fig.savefig('cpl_example.pdf')

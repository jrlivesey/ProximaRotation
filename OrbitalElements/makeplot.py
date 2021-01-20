"""
This script visualizes 100 evolutions of Proxima Centauri b using VSPACE,
varying the planet's initial longitude of ascending node and longitude of
pericenter.

Joseph R. Livesey, 2021
"""



import numpy as np
import vplot as vpl
import matplotlib.pyplot as plt
import os


# Path of vspace dest dir
path = os.getcwd()
dirs = os.listdir(path)


# Initialize arrays
longa = []
longp = []
time = []
rotper = []


# Extract plot data
for i in range(10):
    for j in range(10):
        os.chdir(path+'/trial_longa'+str(i)+'_longp'+str(j))

        vpl_longa = vpl.GetOutput().ProximaB.LongA
        vpl_longp = vpl.GetOutput().ProximaB.LongP

        longa.append(str(round(vpl_longa[0])))
        longp.append(str(round(vpl_longp[0])))

        vpl_time = vpl.GetOutput().ProximaB.Time
        vpl_rotper = vpl.GetOutput().ProximaB.RotPer

        time.append(vpl_time)
        rotper.append(vpl_rotper)

        os.chdir(path)


# Make plots
fig,ax = plt.subplots(10,10)
fig.set_size_inches(30,30)
ax = ax.ravel()

for k in range(100):
    ax[k].plot(time[k], rotper[k], color='k', linestyle='-')
    ax[k].set_title(r'$\Omega=$'+longa[k]+r'$^\circ$, $\varpi=$'+longp[k]+r'$^\circ$', fontsize=16)
    ax[k].set_xlabel('Time (Gyr)', fontsize=12)
    ax[k].set_ylabel('Rotational period (days)', fontsize=12)


# Save figure
vpl.make_pretty(fig)
fig.savefig('initialspace.pdf')
plt.show()

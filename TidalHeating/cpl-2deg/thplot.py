import numpy as np
import matplotlib.pyplot as plt
import vplot as vpl
import sys
import scipy.signal as sig
import subprocess

# Check correct number of arguments
if (len(sys.argv) != 2):
    print('ERROR: Incorrect number of arguments.')
    print('Usage: '+sys.argv[0]+' <pdf | png>')
    exit(1)
if (sys.argv[1] != 'pdf' and sys.argv[1] != 'png'):
    print('ERROR: Unknown file format: '+sys.argv[1])
    print('Options are: pdf, png')
    exit(1)

out = vpl.GetOutput()

time = out.Proxima.Time

fig = plt.figure()
fig.set_size_inches(15,5)
fig.tight_layout()

plt.subplot(1,3,1)
plt.semilogx(time, out.ProximaB.SemiMajorAxis, color='k')
plt.ylabel('Semi-Major Axis (AU)')
plt.xlabel('Time (yr)')

plt.subplot(1,3,2)
plt.semilogx(time, out.ProximaB.Eccentricity, color='k')
plt.ylabel('Eccentricity')
plt.xlabel('Time (yr)')

plt.subplot(1,3,3)
plt.loglog(time, out.ProximaB.SurfEnFluxEqtide, color='k')
#plt.xlim(2.85e6,3.10e6)
plt.ylabel('Tidal heat flux (W m$^{-2}$)')
plt.xlabel('Time (yr)')

# Save as
vpl.make_pretty(fig)
if (sys.argv[1] == 'pdf'):
    fig.savefig('THProxima.pdf')
if (sys.argv[1] == 'png'):
    fig.savefig('THProxima.png')
plt.close()

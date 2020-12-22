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

time = out.Proxima.Time/1e6

# Eccentricity vs. Time (b & c)
fig = plt.figure(figsize=(6.5,8))
plt.subplot(3,2,1)
plt.plot(time,out.ProximaB.Eccentricity,color=vpl.colors.pale_blue,label='b')
plt.plot(time,out.ProximaC.Eccentricity,color=vpl.colors.orange,label='c')
plt.xlim(0,10)
plt.ylabel(r'Eccentricity')
plt.legend(loc=0, fontsize=12, ncol=1)

# Inclination vs. Time (b & c)
plt.subplot(3,2,2)
plt.plot(time,out.ProximaB.Inc,color=vpl.colors.pale_blue)
plt.plot(time,out.ProximaC.Inc,color=vpl.colors.orange)
plt.xlim(0,10)
plt.ylabel('Inclination ($^\circ$)')

# Semimajor axis vs. Time (b)
plt.subplot(3,2,3)
plt.plot(time,out.ProximaB.SemiMajorAxis,color=vpl.colors.pale_blue)
plt.xlim(0,10)
plt.ylabel('Semimajor Axis (AU)')
plt.xlabel('Time (Myr)')

# Rotational period vs. Time (b)
plt.subplot(3,2,4)
plt.axhline(11.0,0,7e6,color='gray',linestyle='--')
plt.plot(time,out.ProximaB.RotPer,color=vpl.colors.pale_blue)
plt.xlim(0,10)
plt.ylabel('Rotation Period (days)')
plt.xlabel('Time (Myr)')

# Save as
vpl.make_pretty(fig)
if (sys.argv[1] == 'pdf'):
    fig.savefig('Proxima.pdf')
if (sys.argv[1] == 'png'):
    fig.savefig('Proxima.png')
plt.close()

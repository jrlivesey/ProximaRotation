# InvPlane results for Proxima Centauri

We apply the [InvPlane](https://github.com/RoryBarnes/InvPlane) tool
in order to obtain orbital parameters for Proxima b and c. InvPlane
rotates a planetary system so that an inclination of 0 corresponds to
the invariable or fundamental plane. The file \"InvProxima.in\" contains
initial conditions relative to the sky plane (observational for Proxima
c, from Benedict & McArthur 2020, and guessed for Proxima b). The file
\"InvProxima.in.inv\" contains the resulting values relative to the
fundamental plane of the system. Note that the resulting values for
eccentricity and inclination for Proxima b were ignored, as in the final
VPLanet simulations we chose particular mutual inclinations and chose
eccentricity to yield a certain final value.

To execute:

``` {.sourceCode .bash}
invplane -m <infile>
```

Format of the input file:

``` {.sourceCode .bash}
NumberOrbiters
CentralMass
Mass SemimajorAxis Eccentricity Inclination LongitudeAscendingNode ArgumentPericenter MeanAnomaly
```

where NumberOrbiters is the number of orbiters (don\'t count the
primary!). The next line is the mass of the primary in solar masses. The
next N lines contain the orbital elements of the orbiters in the order
shown. The units are solar masses, AU, and degrees.

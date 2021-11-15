import numpy as np

#############################################################################################################################
# FUNCTIONS FOR FINDING THE CRITICAL ECCENTRICITY FROM THE PROLATENESS (C_22) WITH CTL MODEL
#############################################################################################################################

def f2(ecc):
    return 1 + (15/2)*ecc**2 + (45/8)*ecc**4 + (5/16)*ecc**6

def f5(ecc):
    return 1 + 3*ecc**2 + (3/8)*ecc**4

def T(p, ecc):
    return (1 - ecc**2)**(-6) * (f2(ecc) - (1 - ecc**2)**1.5 * f5(ecc) * p)

def CritProlateness_11(k2lag, mm, mass1, mass2, radius, semi, ecc):
    """
    Eq. 13 in Rodríguez et al. (2012), with p = 1.
    """
    H = 1 - (5/2)*ecc**2 + (13/16)*ecc**4 - (35/288)*ecc**6
    return 0.5 * k2lag * mm * (mass1/mass2) * (radius/semi)**3 * (T(1, ecc) / H)

def CritProlateness_32(k2lag, mm, mass1, mass2, radius, semi, ecc):
    """
    Eq. 13 in Rodríguez et al. (2012), with p = 3/2.
    """
    H = (7/2)*ecc - (123/16)*ecc**3 + (489/128)*ecc**5
    return 0.5 * k2lag * mm * (mass1/mass2) * (radius/semi)**3 * (T(1.5, ecc) / H)

def CritProlateness_21(k2lag, mm, mass1, mass2, radius, semi, ecc):
    """
    Eq. 13 in Rodríguez et al. (2012), with p = 2.
    """
    H = (17/2)*ecc**2 - (115/6)*ecc**4 + (601/48)*ecc**6
    return 0.5 * k2lag * mm * (mass1/mass2) * (radius/semi)**3 * (T(2, ecc) / H)

def proxima_CritProlateness(p, eccs, k2lag, mm, mass1, mass2, semi):
    """
    Applies values for Proxima b to Eq. 13 in Rodríguez et al. (2012).
    """
    crits = []
    radius = 4.264e-5 * (mass2 ** 0.274)  # Proxima b radius in AU
    if p == 1:
        for ecc in eccs:
            crit = CritProlateness_11(k2lag, mm, mass1, mass2, radius, semi, ecc)
            crits.append(crit)
    elif p == 1.5:
        for ecc in eccs:
            crit = CritProlateness_32(k2lag, mm, mass1, mass2, radius, semi, ecc)
            crits.append(crit)
    elif p == 2:
        for ecc in eccs:
            crit = CritProlateness_21(k2lag, mm, mass1, mass2, radius, semi, ecc)
            crits.append(crit)
    return crits

def ecc_of_capture(my_c22, k2lag, mm, mass1, mass2, semi):
    """
    Solves numerically for the eccentricity at which the strength criterion of the 1/1 SOR is satisfied, but not that of the 3/2 SOR.
    my_c22: Prolateness of Proxima b.
    """
    ecc_of_capture = 0
    eccs     = np.linspace(1e-10, 0.3, 5000)[::-1]
    above_11 = True
    above_32 = True
    captured = False
    if np.abs(proxima_CritProlateness(1.5, [eccs[0]], k2lag, mm, mass1, mass2, semi)) >= my_c22:
        above_32 = False
    if np.abs(proxima_CritProlateness(1, [eccs[0]], k2lag, mm, mass1, mass2, semi)) >= my_c22:
        above_11 = False
    for ecc in eccs:
        c22_32 = np.abs(proxima_CritProlateness(1.5, [ecc], k2lag, mm, mass1, mass2, semi))
        if above_32 == True:
            if c22_32 >= my_c22:
                above_32 = False
        else:
            if c22_32 < my_c22:
                above_32 = True
        c22_11 = np.abs(proxima_CritProlateness(1, [ecc], k2lag, mm, mass1, mass2, semi))
        if above_11 == True:
            if c22_11 >= my_c22:
                above_11 = False
        else:
            if c22_11 < my_c22:
                above_11 = True
        if captured == False:
            if above_11 == True and above_32 == False:
                ecc_of_capture = ecc
                captured = True
    return ecc_of_capture

#############################################################################################################################
# FUNCTIONS FOR FINDING THE CRITICAL ECCENTRICITY FROM THE PROLATENESS (C_22) WITH CPL MODEL
#############################################################################################################################

def cpl_f2(ecc):
    return 1 + (15/2)*ecc**2

def cpl_f5(ecc):
    return 1 + 3*ecc**2

def cpl_T(p, ecc):
    return (1 - ecc**2)**(-6) * (cpl_f2(ecc) - (1 - ecc**2)**1.5 * cpl_f5(ecc) * p)

def cpl_CritProlateness_11(k2, tidalQ, mass1, mass2, radius, semi, ecc):
    """
    Eq. 13 in Rodríguez et al. (2012), with p = 1.
    """
    H = 1 - (5/2)*ecc**2
    return 0.5 * (k2/tidalQ) * (mass1/mass2) * (radius/semi)**3 * (cpl_T(1, ecc) / H)

def cpl_CritProlateness_32(k2, tidalQ, mass1, mass2, radius, semi, ecc):
    """
    Eq. 13 in Rodríguez et al. (2012), with p = 3/2.
    """
    H = (7/2)*ecc
    return 0.5 * (k2/tidalQ) * (mass1/mass2) * (radius/semi)**3 * (cpl_T(1.5, ecc) / H)

def cpl_CritProlateness_21(k2, tidalQ, mass1, mass2, radius, semi, ecc):
    """
    Eq. 13 in Rodríguez et al. (2012), with p = 2.
    """
    H = (17/2)*ecc**2
    return 0.5 * (k2/tidalQ) * (mass1/mass2) * (radius/semi)**3 * (cpl_T(2, ecc) / H)

def cpl_proxima_CritProlateness(p, eccs, k2, tidalQ, mass1, mass2, semi):
    """
    Applies values for Proxima b to Eq. 13 in Rodríguez et al. (2012).
    """
    crits = []
    radius = 4.264e-5 * (mass2 ** 0.274)  # Proxima b radius in AU
    if p == 1:
        for ecc in eccs:
            crit = cpl_CritProlateness_11(k2, tidalQ, mass1, mass2, radius, semi, ecc)
            crits.append(crit)
    elif p == 1.5:
        for ecc in eccs:
            crit = cpl_CritProlateness_32(k2, tidalQ, mass1, mass2, radius, semi, ecc)
            crits.append(crit)
    elif p == 2:
        for ecc in eccs:
            crit = cpl_CritProlateness_21(k2, tidalQ, mass1, mass2, radius, semi, ecc)
            crits.append(crit)
    return crits

def cpl_ecc_of_capture(my_c22, k2, tidalQ, mass1, mass2, semi):
    """
    Solves numerically for the eccentricity at which the strength criterion of the 1/1 SOR is satisfied, but not that of the 3/2 SOR.
    my_c22: Prolateness of Proxima b.
    """
    ecc_of_capture = 0
    eccs     = np.linspace(1e-10, 0.3, 5000)[::-1]
    above_11 = True
    above_32 = True
    captured = False
    if np.abs(cpl_proxima_CritProlateness(1.5, [eccs[0]], k2, tidalQ, mass1, mass2, semi)) >= my_c22:
        above_32 = False
    if np.abs(cpl_proxima_CritProlateness(1, [eccs[0]], k2, tidalQ, mass1, mass2, semi)) >= my_c22:
        above_11 = False
    for ecc in eccs:
        c22_32 = np.abs(cpl_proxima_CritProlateness(1.5, [ecc], k2, tidalQ, mass1, mass2, semi))
        if above_32 == True:
            if c22_32 >= my_c22:
                above_32 = False
        else:
            if c22_32 < my_c22:
                above_32 = True
        c22_11 = np.abs(cpl_proxima_CritProlateness(1, [ecc], k2, tidalQ, mass1, mass2, semi))
        if above_11 == True:
            if c22_11 >= my_c22:
                above_11 = False
        else:
            if c22_11 < my_c22:
                above_11 = True
        if captured == False:
            if above_11 == True and above_32 == False:
                ecc_of_capture = ecc
                captured = True
    return ecc_of_capture

"""
This script sets the critical eccentricity for capture into the 1/1 SOR, then
determines where the instantaneous, time-averaged, or maximum eccentricity first
drops below that value.

Outputs time of capture (ToC) for each simulation and a small data file for use
in making plots involving the ToCs.

(c) Joseph R. Livesey, University of Washington, 2021
"""


import numpy as np
import argparse
import sys
import os
import astropy.units as u
from astropy.constants import G
from vplanet import get_output
<<<<<<< HEAD
from strength_criterion import *
=======
>>>>>>> 1f5fa7ea636fdb8079cda71e98fab652e459bf75


path = os.getcwd()

key_data = []


parser = argparse.ArgumentParser(description="Round and round the planet goes, where it stops, no one knows!")
parser.add_argument('-f', '--filename', type=str, default=None, help="Specify which file in every directory corresponds to the planet under consideration. If no filename specified, pre-written ToC data will be used.")
parser.add_argument('-c', '--gravc', type=float, default=None, help="Value of C_22 for calculating eccentricity at capture.")
parser.add_argument('-e', '--ecrit', type=float, default=None, help="Specify a blanket critical eccentricity, if you don't want one to be calculated for each sim.")
parser.add_argument('-m', '--model', type=str, default=None, help="Specify the equilibrium tidal model (cpl or ctl)")
parser.add_argument('-min', '--min', action='store_true', default=False, help="Recognize ToC where the minimum eccentricity drops below critical value, as opposed to time-averaged eccentricity.")
parser.add_argument('-max', '--max', action='store_true', default=False, help="Recognize ToC where maximum eccentricity drops below critical value, as opposed to time-averaged eccentricity.")
args = parser.parse_args()


#############################################################################################################################
<<<<<<< HEAD
=======
# FUNCTIONS FOR FINDING THE CRITICAL ECCENTRICITY FROM THE PROLATENESS (C_22) WITH CTL MODEL
#############################################################################################################################

def f2(ecc):
    f2 = 1 + (15/2)*ecc**2 + (45/8)*ecc**4 + (5/16)*ecc**6
    return f2

def f5(ecc):
    f5 = 1 + 3*ecc**2 + (3/8)*ecc**4
    return f5

def T(p, ecc):
    T = (1 - ecc**2)**(-6) * (f2(ecc) - (1 - ecc**2)**1.5 * f5(ecc) * p)
    return T

def CritProlateness_11(k2lag, mm, mass1, mass2, radius, semi, ecc):
    """
    Eq. 13 in Rodríguez et al. (2012), with p = 1.
    """
    H = 1 - (5/2)*ecc**2 + (13/16)*ecc**4 - (35/288)*ecc**6
    prolateness = 0.5 * k2lag * mm * (mass2/mass1) * (radius/semi)**3 * (T(1, ecc) / H)
    return prolateness

def CritProlateness_32(k2lag, mm, mass1, mass2, radius, semi, ecc):
    """
    Eq. 13 in Rodríguez et al. (2012), with p = 3/2.
    """
    H = (7/2)*ecc - (123/16)*ecc**3 + (489/128)*ecc**5
    prolateness = 0.5 * k2lag * mm * (mass2/mass1) * (radius/semi)**3 * (T(1.5, ecc) / H)
    return prolateness

def CritProlateness_21(k2lag, mm, mass1, mass2, radius, semi, ecc):
    """
    Eq. 13 in Rodríguez et al. (2012), with p = 2.
    """
    H = (17/2)*ecc**2 - (115/6)*ecc**4 + (601/48)*ecc**6
    prolateness = 0.5 * k2lag * mm * (mass2/mass1) * (radius/semi)**3 * (T(2, ecc) / H)
    return prolateness

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
    eccs      = np.linspace(1e-10, 0.3, 5000)[::-1]
    above_11  = True
    above_32  = True
    captured  = False
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
                print(str(ecc))
                captured = True

#############################################################################################################################
# FUNCTIONS FOR FINDING THE CRITICAL ECCENTRICITY FROM THE PROLATENESS (C_22) WITH CPL MODEL
#############################################################################################################################

def cpl_f2(ecc):
    f2 = 1 + (15/2)*ecc**2
    return f2

def cpl_f5(ecc):
    f5 = 1 + 3*ecc**2
    return f5

def cpl_T(p, ecc):
    T = (1 - ecc**2)**(-6) * (cpl_f2(ecc) - (1 - ecc**2)**1.5 * cpl_f5(ecc) * p)
    return T

def cpl_CritProlateness_11(k2, tidalQ, mass1, mass2, radius, semi, ecc):
    """
    Eq. 13 in Rodríguez et al. (2012), with p = 1.
    """
    H = 1 - (5/2)*ecc**2
    prolateness = 0.5 * (k2/tidalQ) * (mass2/mass1) * (radius/semi)**3 * (cpl_T(1, ecc) / H)
    return prolateness

def cpl_CritProlateness_32(k2, tidalQ, mass1, mass2, radius, semi, ecc):
    """
    Eq. 13 in Rodríguez et al. (2012), with p = 3/2.
    """
    H = (7/2)*ecc
    prolateness = 0.5 * (k2/tidalQ) * (mass2/mass1) * (radius/semi)**3 * (cpl_T(1.5, ecc) / H)
    return prolateness

def cpl_CritProlateness_21(k2, tidalQ, mass1, mass2, radius, semi, ecc):
    """
    Eq. 13 in Rodríguez et al. (2012), with p = 2.
    """
    H = (17/2)*ecc**2
    prolateness = 0.5 * (k2/tidalQ) * (mass2/mass1) * (radius/semi)**3 * (cpl_T(2, ecc) / H)
    return prolateness

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
    eccs      = np.linspace(1e-10, 0.3, 5000)[::-1]
    above_11  = True
    above_32  = True
    captured  = False
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


#############################################################################################################################
>>>>>>> 1f5fa7ea636fdb8079cda71e98fab652e459bf75
# DATA RETRIEVAL
#############################################################################################################################

if args.filename == None:
    # RETRIEVING PRE-EXISTING TOC CALCULATIONS
    # I have yet to re-write this
    exit(1)

else:
    filename = args.filename

    for dir in os.listdir(path):
        if os.path.isdir(dir) == True:
            os.chdir(dir)
            dir_path = os.path.join(path, dir)

<<<<<<< HEAD
            # key_data = []
            sim_key_data = None
            crit = 0
            tidalQ = 0
            tidalTau = 0
            ecc_init = 0
            inc_init = 0
            inc_init_c = 0
            mass1 = 0
            mass2 = 0
            orbp = 0
            semi = 0
            mm = 0

            # if os.path.exists(os.path.join(dir_path, 'toc.txt')) == False:
            if args.ecrit == None:
                # CALCULATE CRITICAL ECCENTRICITY

                star = os.path.join(dir_path, 'star.in')
                with open(star, 'r') as f:
                    star_content = [line.strip().split() for line in f.readlines()]
                    for line in star_content:
                        if line:
                            if line[0] == 'dMass':
                                mass1 = float(line[1]) * 3.33e5                     # Convert to earth masses

                input = os.path.join(dir_path, filename)
                with open(input, 'r') as f:
                    split_content = [line.strip().split() for line in f.readlines()]
                    for line in split_content:
                        if line:
                            if line[0] == 'dMass':
                                mass2 = float(line[1].replace('-',''))
                                radius = 4.264e-5 * (mass2 ** 0.274)                # Sotin et al. 2007 mass-radius relationship, convert to AU
                            if line[0] == 'dSemi':
                                semi = float(line[1].replace('-',''))
                                orbp = np.sqrt( 4 * np.pi * (semi * u.AU) ** 3 / ( G * (mass1 * u.earthMass) ) ).to(u.s) * 1/u.s
                                mm = 2 * np.pi / orbp
                            if line[0] == 'dOrbPeriod':
                                orbp = float(line[1].replace('-',''))               # Days
                                semi = ( ( G * (mass1 * u.earthMass) * (orbp * u.d) ** 2 / (4 * np.pi ** 2) ) ** (1/3) ).to(u.AU) * 1/u.AU
                                mm = 2 * np.pi / orbp
                            if line[0] == 'dTidalQ':
                                tidalQ = float(line[1])
                            if line[0] == 'dTidalTau':
                                tidalTau = float(line[1])                           # seconds
                            if line[0] == 'dEcc':
                                ecc_init = float(line[1])
                            if line[0] == 'dInc':
                                inc_init = float(line[1])

                if os.path.exists(os.path.join(dir_path, 'c.in')):
                    with open('c.in', 'r') as f:
                        c_content = [line.strip().split() for line in f.readlines()]
                        for line in c_content:
                            if line:
                                if line[0] == 'dInc':
                                    inc_init_c = float(line[1])

                    k2 = 0.3
                    k2lag = k2 * tidalTau
                    my_c22 = args.gravc

                    if args.model == 'ctl':
                        crit = ecc_of_capture(my_c22, k2lag, mm, mass1, mass2, semi)
                    elif args.model == 'cpl':
                        crit = cpl_ecc_of_capture(my_c22, k2, tidalQ, mass1, mass2, semi)
                    else:
                        print('ERROR: Incompatible or no tidal model supplied.')
                        exit(0)

                    print(crit)

            else:
                crit = args.ecrit


            # DETERMINE WHETHER THE PLANET IS WITHIN OBSERVATIONAL CONSTRAINTS
            if os.path.exists('Proxima.log') == False:
                break
            a_final = get_output(units=False).ProximaB.SemiMajorAxis[-1]
            if a_final > 0.04833 and a_final < 0.04895:
                constrained = True
            else:
                constrained = False


            # USE CRITICAL ECCENTRICITY TO FIND TOC
            eccs = get_output(units=False).ProximaB.Eccentricity
            avg_eccs = [eccs[0]]

            if args.min == True:
                for val in eccs:
                    if val <= crit:
                        toc_index = eccs.index(val)
                        toc = get_output(units=False).ProximaB.Time[toc_index]
            elif args.max == True:
                print("Haven't quite figured this bit out.")
                exit(0)
            else:
                i = 0
                while i in range(len(eccs) - 1):
                    next_val = (eccs[i+1] + avg_eccs[i] * (i + 1)) / (i + 2)
                    avg_eccs.append(next_val)
                    if next_val <= crit:
                        toc_index = len(avg_eccs)
                        toc = get_output(units=False).ProximaB.Time[toc_index]
                        break
                    elif i == (len(eccs) - 2):
                        toc = 'N/A'
                        break
                    else:
                        i += 1

            # MAKE STRING OF KEY DATA TO ADD TO LIST
            if constrained == True:
                if args.ecrit == None:
                    # sim_key_data = [ecc_init, inc_init, inc_init_c, semi, mass1, mass2, crit, toc]
                    sim_key_data = str(ecc_init) + ' ' + str(inc_init) + ' ' + str(inc_init_c) + ' ' + str(semi) + ' ' + str(mass1) + ' ' + str(mass2) + ' ' + str(crit) + ' ' + str(toc)
                else:
                    # sim_key_data = [ecc_init, inc_init, inc_init_c, crit, toc]
                    sim_key_data = str(ecc_init) + ' ' + str(inc_init) + ' ' + str(inc_init_c) + ' ' + str(crit) + ' ' + str(toc)


            key_data.append(sim_key_data)

            os.chdir(path)
=======
                key_data = []
                crit = 0
                tidalQ = 0
                tidalTau = 0
                ecc_init = 0
                inc_init = 0
                inc_init_c = 0
                mass1 = 0
                mass2 = 0
                orbp = 0
                semi = 0
                mm = 0

                # if os.path.exists(os.path.join(dir_path, 'toc.txt')) == False:
                if args.ecrit == None:
                    # CALCULATE CRITICAL ECCENTRICITY

                    star = os.path.join(dir_path, 'star.in')
                    with open(star, 'r') as f:
                        star_content = [line.strip().split() for line in f.readlines()]
                        for line in star_content:
                            if line:
                                if line[0] == 'dMass':
                                    mass1 = float(line[1]) * 3.33e5                     # Convert to earth masses

                    input = os.path.join(dir_path, filename)
                    with open(input, 'r') as f:
                        split_content = [line.strip().split() for line in f.readlines()]
                        for line in split_content:
                            if line:
                                if line[0] == 'dMass':
                                    mass2 = float(line[1].replace('-',''))
                                    radius = 4.264e-5 * (mass2 ** 0.274)                # Sotin et al. 2007 mass-radius relationship, convert to AU
                                if line[0] == 'dSemi':
                                    semi = float(line[1].replace('-',''))
                                    orbp = np.sqrt( 4 * np.pi * (semi * u.AU) ** 3 / ( G * (mass1 * u.earthMass) ) ).to(u.s) * 1/u.s
                                    mm = 2 * np.pi / orbp
                                if line[0] == 'dOrbPeriod':
                                    orbp = float(line[1].replace('-',''))               # Days
                                    semi = ( ( G * (mass1 * u.earthMass) * (orbp * u.d) ** 2 / (4 * np.pi ** 2) ) ** (1/3) ).to(u.AU) * 1/u.AU
                                    mm = 2 * np.pi / orbp
                                if line[0] == 'dTidalQ':
                                    tidalQ = float(line[1])
                                if line[0] == 'dTidalTau':
                                    tidalTau = float(line[1])                           # seconds
                                if line[0] == 'dEcc':
                                    ecc_init = float(line[1])
                                if line[0] == 'dInc':
                                    inc_init = float(line[1])

                    if os.path.exists(os.path.join(dir_path, 'c.in')):
                        with open('c.in', 'r') as f:
                            c_content = [line.strip().split() for line in f.readlines()]
                            for line in c_content:
                                if line:
                                    if line[0] == 'dInc':
                                        inc_init_c = float(line[1])

                        k2 = 0.3
                        k2lag = k2 * tidalTau
                        my_c22 = args.gravc

                        if args.model == 'ctl':
                            crit = ecc_of_capture(my_c22, k2lag, mm, mass1, mass2, semi)
                        elif args.model == 'cpl':
                            crit = cpl_ecc_of_capture(my_c22, k2, tidalQ, mass1, mass2, semi)
                        else:
                            print('ERROR: Incompatible or no tidal model supplied.')
                            exit(0)

                        print(crit)

                else:
                    crit = args.ecrit


                # DETERMINE WHETHER THE PLANET IS WITHIN OBSERVATIONAL CONSTRAINTS
                if os.path.exists('Proxima.log') == False:
                    break
                a_final = get_output(units=False).ProximaB.SemiMajorAxis[-1]
                if a_final > 0.04833 and a_final < 0.04895:
                    constrained = True
                else:
                    constrained = False


                # USE CRITICAL ECCENTRICITY TO FIND TOC
                eccs = get_output(units=False).ProximaB.Eccentricity
                avg_eccs = [eccs[0]]

                if args.min == True:
                    for val in eccs:
                        if val <= crit:
                            toc_index = eccs.index(val)
                            toc = get_output(units=False).ProximaB.Time[toc_index]
                elif args.max == True:
                    print("Haven't quite figured this bit out.")
                    exit(0)
                else:
                    i = 0
                    while i in range(len(eccs) - 1):
                        next_val = (eccs[i+1] + avg_eccs[i] * (i + 1)) / (i + 2)
                        avg_eccs.append(next_val)
                        if next_val <= crit:
                            toc_index = len(avg_eccs)
                            toc = get_output(units=False).ProximaB.Time[toc_index]
                            break
                        elif i == (len(eccs) - 2):
                            toc = 'N/A'
                            break
                        else:
                            i += 1

                # MAKE STRING OF KEY DATA TO ADD TO LIST
                if constrained == True:
                    if args.ecrit == None:
                        # sim_key_data = [ecc_init, inc_init, inc_init_c, semi, mass1, mass2, crit, toc]
                        sim_key_data = str(ecc_init) + ' ' + str(inc_init) + ' ' + str(inc_init_c) + ' ' + str(semi) + ' ' + str(mass1) + ' ' + str(mass2) + ' ' + str(crit) + ' ' + str(toc)
                    else:
                        # sim_key_data = [ecc_init, inc_init, inc_init_c, crit, toc]
                        sim_key_data = str(ecc_init) + ' ' + str(inc_init) + ' ' + str(inc_init_c) + ' ' + str(crit) + ' ' + str(toc)


        key_data.append(sim_key_data)

        os.chdir(path)
>>>>>>> 1f5fa7ea636fdb8079cda71e98fab652e459bf75


# WRITE LIST TO FILE
with open('key_data.txt', 'w') as file:
    lines = [sim_key_data for sim_key_data in key_data]
    file.writelines(lines)

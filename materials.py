# import cv2
import glob
import numpy as np
import os
from os import path
from matplotlib import pyplot as plt
import pandas as pd
from io import StringIO

boltzmann_constant = 8.62e-5 # ev/K
room_temp = (25.0 + 273.0)*boltzmann_constant

class aluminium():
    # material properties
    T_b = 2467+273 # kelvin
    T_m = (660.3 + 273) # kelvin
    T_0 = 273 + 25
    density0 = 2700
    specific_heat = 0.9e+3 #kJ/kg
    specific_heat_vapor = 0.459446
    enthalpy_fusion = 396 #kJ/kg
    enthalpy_vaporization = 11370 #kJ/kg
    sound_speed = 6320 # m/s
    thermal_conductivity = 247 # W/mK
    n = 1.2 # polytropic gas constant pulled from Stanton / Tucker
    MW = 26.98
    C1 = -5.35e-5
    C2 = 0.233
    C3 = 1.21
    C4 = 0.638
    C5 = 1.5
    C6 = 1.2e-2
    C7 = 3.8e-3
    C8 = 18.5
    C9 = 5.96
    C10 = 0.44
    C11 = 3.58e-2
    C12 = 3.05
    k = 0.878
    gruneisen_coeff = 2.13
    T_m0 = 0.0804
    L_F = 0.107
    Z = 17.1
    relative_permeability = 1.000022
    shear_modulus = 27.6e+9 # steinburg guinan
    initial_yield = 2.9e+8
    max_yield = 6.8e+8
    beta = 125
    steinberg_n = 0.1
    # mie gruneisen
    mie_gruneisen_c = 6320
    mie_gruneisen_S = 1.4
    mie_gruneisen_S2 = 0.0
    mie_gruneisen_Gamma = 1.97
    poissons = 0.3
    JC_A = 520e+6  # MPa
    JC_B = 477e+6  # MPa
    JC_n = 0.52
    JC_C = 0.001
    JC_m = 1
    magneto = True

class PMMA():
    # material properties
    density0 = 1170
    Z = 3.26
    yield_stress = 70.0e+6
    specific_heat = 1466  # J/kg
    thermal_cond = 0.167  # W/m*K
    T_melt = 160 + 273
    sound_speed = 2757 # m/s
    # mie gruneisen
    mie_gruneisen_c = 2180
    mie_gruneisen_S = 2.088
    mie_gruneisen_S2 = -1.124
    mie_gruneisen_Gamma = 0.85

class Alumina():
    density0 = 2600
    Z = 25.5
    yield_stress = 210.0e+9
    specific_heat = 880 #J/kg
    thermal_cond = 30 #W/m*K
    T_melt = 2072 + 273
    sound_speed = 9900

class Sapphire():
    density0 = 3980
    Z = 35
    yield_stress = 400e+6
    specific_heat = 1171 #J/kg
    thermal_cond = 46.06 #W/m*K
    T_melt = 2000 + 273
    sound_speed = 10000

class copper():
    density0 = 8940
    specific_heat = 385 #J/kg
    T_m = 1085 + 273
    sound_speed = 2260
    MW = 63.546
    poissons = 0.34
    JC_A = 90e+6 #MPa
    JC_B = 292e+6 #MPa
    JC_n = 0.31
    JC_C = 0.025
    JC_m = 1.09
    magneto = True


class vacuum():
    # material properties
    density0 = 0
    sound_speed = 0 # m/s
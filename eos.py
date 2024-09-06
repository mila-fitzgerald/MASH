import numpy as np

k_b = 1.38e-23
boltzmann_constant = 8.62e-5  # ev/K


def eos_init(total, mat_list):
    nxt, nyt = len(total.density), len(total.density[0])
    for mat in mat_list:
        mat.density_fraction = np.zeros((nxt, nyt))
        for ix in range(2, nxt-2):
            for iy in range(2, nyt-2):
                if total.density[ix, iy] > 0.0:
                    mat.density_fraction[ix, iy] = mat.density[ix, iy] / total.density[ix, iy]


def occupation_update(total, mat_list):
    nxt, nyt = len(total.density), len(total.density[0])
    for mat in mat_list:
        for ix in range(2, nxt-2):
            for iy in range(2, nyt-2):
                if mat.density[ix, iy] > 1.0e-6:
                    mat.occupation[ix, iy] = 1.0
                    total.occupation[ix, iy] = 1.0

def update_sound_speed(domain, mat_list):
    nxt, nyt = domain.nxt + 4, domain.nyt + 4
    for mat in mat_list:
        Cs = mat.material.sound_speed
        for ix in range(2, nxt-2):
            for iy in range(2, nyt-2):
                #
                mat.sound_speed[ix, iy] = Cs

def reduce_sound_speed(total, mat_list):
    nxt, nyt = len(total.density), len(total.density[0])
    total.sound_speed = np.zeros((nxt, nyt))
    eos_init(total, mat_list)
    for mat in mat_list:
        for ix in range(2, nxt-2):
            for iy in range(2, nyt-2):
                #
                density_fraction = mat.density_fraction[ix, iy]
                #
                total.sound_speed[ix, iy] = total.sound_speed[ix, iy] + mat.sound_speed[ix, iy] * density_fraction


def momentum_map(total, mat_list):
    eos_init(total, mat_list)
    nxt, nyt = len(total.density), len(total.density[0])
    for mat in mat_list:
        for ix in range(2, nxt-2):
            for iy in range(2, nyt-2):
                # if mat.density[ix,iy] > 0.0:
                #
                density_fraction = mat.density_fraction[ix, iy]
                #
                mat.momentum_x[ix, iy] = mat.momentum_x[ix, iy] + total.dt_momentum_x[ix, iy] * density_fraction
                mat.momentum_y[ix, iy] = mat.momentum_y[ix, iy] + total.dt_momentum_y[ix, iy] * density_fraction
    # don't forget to reduce
    momentum_reduce(total, mat_list)


def momentum_reduce(total, mat_list):
    nxt, nyt = len(total.density), len(total.density[0])
    total.momentum_x = np.zeros((nxt, nyt))
    total.momentum_y = np.zeros((nxt, nyt))
    #
    for mat in mat_list:
        for ix in range(2, nxt-2):
            for iy in range(2, nyt-2):
                if total.occupation[ix, iy] > 0.0:
                    total.momentum_x[ix, iy] = total.momentum_x[ix, iy] + mat.momentum_x[ix, iy]
                    total.momentum_y[ix, iy] = total.momentum_y[ix, iy] + mat.momentum_y[ix, iy]


def internal_energy_map(total, mat_list):
    eos_init(total, mat_list)
    nxt, nyt = len(total.density), len(total.density[0])
    for mat in mat_list:
        for ix in range(2, nxt-2):
            for iy in range(2, nyt-2):
                if mat.density[ix, iy] > 0.0:
                    #
                    density_fraction = mat.density_fraction[ix, iy]
                    #
                    mat.internal_energy[ix, iy] = mat.internal_energy[ix, iy] + total.internal_energy_dt[
                        ix, iy] * density_fraction
    # don't forget to reduce
    internal_energy_reduce(total, mat_list)


def internal_energy_reduce(total, mat_list):
    nxt, nyt = len(total.density), len(total.density[0])
    total.internal_energy = np.zeros((nxt, nyt))
    #
    for mat in mat_list:
        for ix in range(2, nxt-2):
            for iy in range(2, nyt-2):
                # if total.density[ix,iy] > 0.0:
                total.internal_energy[ix, iy] = total.internal_energy[ix, iy] + mat.internal_energy[ix, iy]


def pressure_reduce(total, mat_list):
    nxt, nyt = len(total.density), len(total.density[0])
    total.pressure = np.zeros((nxt, nyt))
    eos_init(total, mat_list)
    #
    for mat in mat_list:
        for ix in range(2, nxt-2):
            for iy in range(2, nyt-2):
                if total.occupation[ix, iy] > 0.0:
                    density_fraction = mat.density_fraction[ix, iy]
                    # if total.density[ix,iy] > 0.0:
                    total.pressure[ix, iy] = total.pressure[ix, iy] + mat.pressure[ix, iy] * density_fraction


def stress_reduce(total, mat_list):
    nxt, nyt = len(total.density), len(total.density[0])
    total.stress_xx = np.zeros((nxt, nyt))
    total.stress_xy = np.zeros((nxt, nyt))
    total.stress_yx = np.zeros((nxt, nyt))
    total.stress_yy = np.zeros((nxt, nyt))
    #
    for mat in mat_list:
        for ix in range(0, nxt):
            for iy in range(0, nyt):
                # if total.density[ix,iy] > 0.0:
                total.stress_xx[ix, iy] = total.stress_xx[ix, iy] + mat.stress_xx[ix, iy]
                total.stress_xy[ix, iy] = total.stress_xy[ix, iy] + mat.stress_xy[ix, iy]
                total.stress_yx[ix, iy] = total.stress_yx[ix, iy] + mat.stress_yx[ix, iy]
                total.stress_yy[ix, iy] = total.stress_yy[ix, iy] + mat.stress_yy[ix, iy]


def stress_map(total, mat_list):
    nxt, nyt = len(total.density), len(total.density[0])
    #
    for mat in mat_list:
        mat.stress_xx = np.zeros((nxt, nyt))
        mat.stress_xy = np.zeros((nxt, nyt))
        mat.stress_yx = np.zeros((nxt, nyt))
        mat.stress_yy = np.zeros((nxt, nyt))
        for ix in range(0, nxt):
            for iy in range(0, nyt):
                # if total.density[ix,iy] > 0.
                mat.stress_xx[ix, iy] = total.stress_xx[ix, iy] * mat.density_fraction[ix, iy]
                mat.stress_xy[ix, iy] = total.stress_xy[ix, iy] * mat.density_fraction[ix, iy]
                mat.stress_yx[ix, iy] = total.stress_yx[ix, iy] * mat.density_fraction[ix, iy]
                mat.stress_yy[ix, iy] = total.stress_yy[ix, iy] * mat.density_fraction[ix, iy]


def density_reduce(total, mat_list):
    nxt, nyt = len(total.density), len(total.density[0])
    total.density = np.zeros((nxt, nyt))
    #
    for mat in mat_list:
        for ix in range(2, nxt-2):
            for iy in range(2, nyt-2):
                if total.occupation[ix, iy] >= 1.0:
                    total.density[ix, iy] = total.density[ix, iy] + mat.density[ix, iy]


def calculate_yield_stress(total, mat_list):
    nxt, nyt = len(total.density), len(total.density[0])
    #
    for mat in mat_list:
        mat.yield_stress = np.zeros((nxt, nyt))
        Y0, beta, n = mat.material.initial_yield, mat.material.beta, mat.material.steinberg_n
        rho0, T_melt = mat.material.density0, mat.material.T_m
        for ix in range(0, nxt):
            for iy in range(0, nyt):
                if mat.density[ix, iy] > 0.0:
                    mat.yield_stress[ix, iy] = Y0 * (1 + beta * total.plastic_strain[ix, iy]) ** n

                    if mat.density[ix, iy] < 0.1 * rho0:
                        mat.yield_stress[ix, iy] = mat.yield_stress[ix, iy] * np.exp(
                            -((mat.density[ix, iy] - 0.1 * rho0) ** 2 / 1000.0))

                    mat.yield_stress[ix, iy] = min(mat.material.max_yield, mat.yield_stress[ix, iy])

                    if mat.temperature[ix, iy] > T_melt:
                        mat.yield_stress[ix, iy] = 0.0
    yield_reduce(total, mat_list)


def johnson_cook_yield(total, mat_list):
    nxt, nyt = len(total.density), len(total.density[0])
    #
    for mat in mat_list:
        mat.yield_stress = np.zeros((nxt, nyt))
        A, B, n, C, m = mat.material.JC_A, mat.material.JC_B, mat.material.JC_n, mat.material.JC_C, mat.material.JC_m
        rho0, T_melt = mat.material.density0, mat.material.T_m
        for ix in range(0, nxt):
            for iy in range(0, nyt):
                if mat.density[ix, iy] > 0.0:
                    total.strain_rate[ix, iy] = max(1.0, total.strain_rate[ix, iy])
                    mat.yield_stress[ix, iy] = (A + B * total.plastic_strain[ix, iy] ** n) * \
                                               (1 + C * np.log(total.strain_rate[ix, iy])) * \
                                               (1 - (mat.temperature[ix, iy] / T_melt) ** m)
                    if mat.density[ix, iy] < 0.1 * rho0:
                        mat.yield_stress[ix, iy] = mat.yield_stress[ix, iy] * np.exp(
                            -((mat.density[ix, iy] - 0.1 * rho0) ** 2 / 1000.0))

                    if mat.temperature[ix, iy] > T_melt:
                        mat.yield_stress[ix, iy] = 0.0
    yield_reduce(total, mat_list)


def yield_reduce(total, mat_list):
    nxt, nyt = len(total.density), len(total.density[0])
    total.yield_stress = np.zeros((nxt, nyt))
    #
    for mat in mat_list:
        for ix in range(0, nxt):
            for iy in range(0, nyt):
                if mat.density[ix, iy] > 0.0:
                    total.yield_stress[ix, iy] = total.yield_stress[ix, iy] + mat.yield_stress[ix, iy]


def shear_modulus_reduce(total, mat_list):
    nxt, nyt = len(total.density), len(total.density[0])
    total.shear_modulus = np.zeros((nxt, nyt))
    eos_init(total, mat_list)
    #
    for mat in mat_list:
        for ix in range(0, nxt):
            for iy in range(0, nyt):
                if total.occupation[ix, iy] > 0.0:
                    density_fraction = mat.density_fraction[ix, iy]
                    #
                    total.shear_modulus[ix, iy] = total.shear_modulus[ix, iy] + mat.shear_modulus[ix, iy] * density_fraction


def shear_modulus_update(total, mat_list):
    nxt, nyt = len(total.density), len(total.density[0])
    #
    for mat in mat_list:
        c, v = mat.material.sound_speed, mat.material.poissons
        for ix in range(2, nxt-1):
            for iy in range(2, nyt-1):
                if mat.occupation[ix, iy] > 0.0:
                    bulk_modulus = mat.density[ix, iy] * c * c
                    mat.shear_modulus[ix, iy] = bulk_modulus * ((1 - 2 * v) / (2 * (1 + v)))

    shear_modulus_reduce(total, mat_list)


def ideal_gas_generate_states(total, mat_list):
    # rho k_b T = m P | P V = n R T | P dx dy = rho dx dy / MW R T
    nxt, nyt = len(total.density), len(total.density[0])
    R0 = 8.315e+3  # universal gas constant
    for mat in mat_list:
        MW = mat.material.MW
        for ix in range(0, nxt):
            for iy in range(0, nyt):
                rho = mat.density[ix, iy]
                mat.pressure[ix, iy] = (rho * R0 * mat.temperature[ix, iy]) / MW


def calculate_T(total, mat_list, dx, dy):
    nxt, nyt = len(total.density), len(total.density[0])
    #
    for mat in mat_list:
        rho0, C_v = mat.material.density0, mat.material.specific_heat
        T0 = mat.init.T0
        for ix in range(0, nxt):
            for iy in range(0, nyt):
                if mat.density[ix, iy] > 0.0:
                    rho = mat.density[ix, iy]
                    # mass = rho
                    mat.temperature[ix, iy] = (mat.internal_energy[ix, iy]) * dx * dy / (rho * C_v) + T0


def burgess_conductivity(grid, total, mat_list, dx, dy):
    nxt, nyt = len(total.density), len(total.density[0])
    for mat in mat_list:
        for ix in range(0, nxt):
            for iy in range(0, nyt):
                if mat.density[ix, iy] > 0.0:
                    MW = mat.material.MW
                    V = 1 / mat.density[ix, iy]
                    V0 = 1 / mat.material.density0
                    T = mat.temperature[ix, iy]
                    mols = mat.density[ix, iy] * dx * dy / MW
                    Lf = mat.material.L_F * 1.0e-5 * 1.0e-6 / mols
                    # solid resistivity
                    C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12 = \
                        mat.material.C1, mat.material.C2, mat.material.C3, \
                        mat.material.C4, mat.material.C5, mat.material.C6, \
                        mat.material.C7, mat.material.C8, mat.material.C9, \
                        mat.material.C10, mat.material.C11, mat.material.C12

                    T_m, gamma, T_b = mat.material.T_m, mat.material.gruneisen_coeff, mat.material.T_b
                    F_gamma = 2 * gamma - 1

                    if T < T_m:
                        solid = True
                        # solid resistivity ********
                        rho_s = (C1 + C2 * (T * boltzmann_constant) ** C3) * (V / V0) ** F_gamma  # in mOhm-cm
                        del_rho = mat.material.k * np.exp(0.069 * Lf / (T_m * boltzmann_constant))

                        rho = rho_s

                    del_rho = mat.material.k * np.exp(0.069 * Lf / (T_m * boltzmann_constant))

                    if T >= T_m:
                        melted = True
                        rho_s_Tm = (C1 + C2 * (T_m * boltzmann_constant) ** C3) * (V / V0) ** F_gamma
                        # liquid resistivity *********
                        rho_l = del_rho * rho_s_Tm * (T / T_m) ** C4

                        rho = rho_l

                    if T >= T_b:
                        vaporised = True
                        # gas resistivity ********
                        rho_ei = C5 / (T * boltzmann_constant) * (
                                1 + np.log(1 + C6 * V * (T * boltzmann_constant) ** 1.5))
                        a_i = (1 + C8 * np.exp(C9 / (T * boltzmann_constant)) / V * (
                                T * boltzmann_constant) ** 1.5) ** -0.5
                        rho_en = C7 * (T * boltzmann_constant) ** 0.5 * (1 + a_i ** -1)
                        rho_v = rho_ei + rho_en

                        # mixed phase resistivity
                        m = (V - V0) * (C10 / C11) * np.exp(-C12 / (T * boltzmann_constant))
                        X_l = (1 - m) / (V / V0)
                        X_V = 1 - X_l

                        rho_mixed = (X_l / rho_l + X_V / rho_v) ** -1

                        rho = rho_mixed
                    mat.conductivity[ix, iy] = 1.0 / (rho * 1.0e-5)

    conductivity_reduce(grid, total, mat_list)


def conductivity_reduce(grid, total, mat_list):
    nxt, nyt = len(total.density), len(total.density[0])
    total.conductivity = np.zeros((nxt, nyt))
    eos_init(grid, total, mat_list)
    #
    for mat in mat_list:
        for ix in range(0, nxt):
            for iy in range(0, nyt):
                if total.occupation[ix, iy] > 0.0:
                    density_fraction = mat.density_fraction[ix, iy]
                    #
                    total.conductivity[ix, iy] = total.conductivity[ix, iy] + mat.conductivity[
                        ix, iy] * density_fraction

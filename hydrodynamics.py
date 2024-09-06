import numpy as np
import eos

def hydro_step(total, mat_list, dx, dy, dt):
    pressure_fx, pressure_fy = face_centred_pressure(total)
    #
    force_x, force_y = hydro_force(total, pressure_fx, pressure_fy, dx, dy)
    artvisc_fx, artvisc_fy = artificial_viscosity(total)
    #
    force_add(total, force_x, force_y, dt)
    artvisc_add(total, artvisc_fx, artvisc_fy, dt)
    #
    hydro_source_energy, artificial_source = hydro_internal_energy(total, pressure_fx, pressure_fy, artvisc_fx, artvisc_fy, dx, dy, dt)
    internal_energy_add(total, hydro_source_energy, artificial_source)
    #
    eos.momentum_map(total, mat_list)
    eos.internal_energy_map(total, mat_list)

def hydro_internal_energy(total, pressure_fx, pressure_fy, artvisc_fx, artvisc_fy, dx, dy, dt):
    nxt, nyt = len(total.density), len(total.density[0])
    one_dx, one_dy = 1 / dx, 1 / dy
    dt_dx, dt_dy = dt*one_dx, dt*one_dy
    hydro_source_energy, artificial_source = [np.zeros((nxt, nyt)) for i in range(2)]
    for ix in range(2, nxt - 2):
        for iy in range(2, nyt - 2):
            if total.occupation[ix,iy] > 0.0:
                delta_vxm = total.velocity_x[ix-1,iy] - total.velocity_x[ix,iy]
                vx_m = total.velocity_f_x[ix,iy]
                delta_vxp = total.velocity_x[ix,iy] - total.velocity_x[ix+1,iy]
                #
                delta_vym = total.velocity_y[ix, iy - 1] - total.velocity_y[ix, iy]
                vy_m = total.velocity_f_y[ix,iy]
                delta_vyp = total.velocity_y[ix, iy] - total.velocity_y[ix, iy + 1]
                #
                hydro_source_energy[ix,iy] = dt_dx * 0.5 * (delta_vxm*pressure_fx[ix,iy] + delta_vxp*pressure_fx[ix+1,iy]) \
                                            + dt_dy * 0.5 * (delta_vym*pressure_fy[ix,iy] + delta_vyp*pressure_fy[ix,iy+1])

                artificial_source[ix,iy] = dt_dx * 0.5 * (delta_vxm*artvisc_fx[ix,iy] + delta_vxp*artvisc_fx[ix+1,iy]) \
                                            + dt_dy * 0.5 * (delta_vym*artvisc_fy[ix,iy] + delta_vyp*artvisc_fy[ix,iy+1])

    return hydro_source_energy, artificial_source

def face_centred_pressure(total):
    nxt, nyt = len(total.density), len(total.density[0])
    pressure_fx, pressure_fy = [np.zeros((nxt, nyt)) for i in range(2)]
    for ix in range(2, nxt - 1):
        for iy in range(2, nyt - 1):
            if total.occupation[ix-1,iy] > 0.0 and total.occupation[ix,iy] > 0.0:
                rho = total.density[ix,iy]
                rho_m = total.density[ix-1,iy]
                one_sum_rho = 1.0 / (rho + rho_m)
                pressure_fx[ix,iy] = (total.pressure[ix-1,iy]*rho + total.pressure[ix,iy]*rho_m) * one_sum_rho
            if total.occupation[ix,iy-1] > 0.0 and total.occupation[ix,iy] > 0.0:
                rho = total.density[ix, iy]
                rho_m = total.density[ix, iy - 1]
                one_sum_rho = 1.0 / (rho + rho_m)
                pressure_fy[ix, iy] = (total.pressure[ix, iy - 1]*rho + total.pressure[ix, iy]*rho_m) * one_sum_rho
                # pressure_fy[ix, iy] = (total.pressure[ix, iy - 1] + total.pressure[ix, iy]) * 0.5
    return pressure_fx, pressure_fy

def internal_energy_add(total, source_energy, artificial_source):
    nxt, nyt = len(total.density), len(total.density[0])
    total.internal_energy_dt = np.zeros((nxt, nyt))
    for ix in range(2, nxt - 2):
        for iy in range(2, nyt - 2):
            total.internal_energy_dt[ix,iy] = source_energy[ix,iy] + artificial_source[ix,iy]

def artificial_viscosity(total):
    nxt, nyt = len(total.density), len(total.density[0])
    artvisc_fx, artvisc_fy = [np.zeros((nxt, nyt)) for i in range(2)]
    for ix in range(2, nxt - 1):
           for iy in range(2, nyt - 1):
               rho = total.density[ix, iy]
               impedance = rho * total.sound_speed[ix, iy]
               if total.occupation[ix - 1, iy] > 0.0 and total.occupation[ix, iy] > 0.0:
                    delta_v = total.velocity_x[ix-1,iy] - total.velocity_x[ix,iy]
                    # delta_v = max(delta_v,0.0)
                    rho_m = total.density[ix-1, iy]
                    impedance_m = rho_m * total.sound_speed[ix-1, iy]
                    one_sum_rho = 1.0 / (rho+rho_m)
                    #
                    artvisc_fx[ix,iy] = delta_v * (rho * impedance_m + rho_m * impedance) * one_sum_rho
               if total.occupation[ix, iy - 1] > 0.0 and total.occupation[ix, iy] > 0.0:
                    delta_v = total.velocity_y[ix,iy-1] - total.velocity_y[ix,iy]
                    # delta_v = max(delta_v,0.0)
                    rho_m = total.density[ix, iy-1]
                    impedance_m = rho_m * total.sound_speed[ix, iy-1]
                    one_sum_rho = 1.0 / (rho+rho_m)
                    #
                    artvisc_fy[ix,iy] = delta_v * (rho * impedance_m + rho_m * impedance) * one_sum_rho
    return artvisc_fx, artvisc_fy

def hydro_force(total, pressure_fx, pressure_fy, dx, dy):
    nxt, nyt = len(total.density), len(total.density[0])
    force_x, force_y = [np.zeros((nxt, nyt)) for i in range(2)]
    one_dx, one_dy = 1 / dx, 1 / dy
    for ix in range(2, nxt - 2):
        for iy in range(2, nyt - 2):
            # if total.occupation[ix,iy] > 0.0:
                d_pressure_x = (pressure_fx[ix,iy] - pressure_fx[ix+1,iy]) * one_dx
                force_x[ix,iy] = d_pressure_x
                d_pressure_y = (pressure_fy[ix,iy] - pressure_fy[ix,iy+1]) * one_dy
                force_y[ix,iy] = d_pressure_y
    return force_x, force_y

def force_add(total, force_x, force_y, dt):
    nxt, nyt = len(total.density), len(total.density[0])
    total.dt_momentum_x = np.zeros((nxt, nyt))
    total.dt_momentum_y = np.zeros((nxt, nyt))
    for ix in range(2, nxt-2):
        for iy in range(2, nyt-2):
            # if total.occupation[ix,iy] > 0.0:
            #     if total.yield_stress[ix,iy] > 0.0:
            hydro_momentum_x = force_x[ix,iy] * dt
            hydro_momentum_y = force_y[ix,iy] * dt
            total.dt_momentum_x[ix,iy] = hydro_momentum_x
            total.dt_momentum_y[ix,iy] = hydro_momentum_y

def artvisc_add(total, artvisc_fx, artvisc_fy, dt):
    nxt, nyt = len(total.density), len(total.density[0])
    # one_dx, one_dy = 1 / dx, 1 / dy
    for ix in range(2, nxt - 2):
        for iy in range(2, nyt - 2):
            if total.occupation[ix,iy] > 0.0:
                visc_x = artvisc_fx[ix,iy] - artvisc_fx[ix+1,iy]
                visc_y = artvisc_fy[ix,iy] - artvisc_fy[ix,iy+1]
                visc_momentum_x = visc_x * dt
                visc_momentum_y = visc_y * dt
                total.dt_momentum_x[ix, iy] = total.dt_momentum_x[ix, iy] + visc_momentum_x
                total.dt_momentum_y[ix, iy] = total.dt_momentum_y[ix, iy] + visc_momentum_y

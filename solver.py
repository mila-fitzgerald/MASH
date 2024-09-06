import numpy as np
from hydrodynamics import hydro_step
from functions import read_outputdt, save_frame
import eos
from boundary import conservative_boundaries
from advection import advection_flux_x, advection_flux_y, flux_apply
from physics import face_centred_velocity, update_total_velocity

vacuum_permeability = 1.257e-6  # H/m

def explicit_solver(mat_list, domain, total, boundaries, folder_path, runtime):
    #
    time = 0.0
    running = True
    percent = 0
    #
    output_dt = read_outputdt(folder_path)
    OUTPUT = False
    if output_dt > 0.0:
        OUTPUT = True
    output_step = 0
    #
    while running:
        eos.update_sound_speed(domain, mat_list)
        eos.occupation_update(total, mat_list)
        eos.density_reduce(total, mat_list)
        eos.momentum_reduce(total, mat_list)
        eos.internal_energy_reduce(total, mat_list)
        eos.reduce_sound_speed(total, mat_list)
        eos.shear_modulus_update(total, mat_list)
        # eos.ideal_gas_generate_states(total, mat_list)
        #
        conservative_boundaries(domain, mat_list, boundaries)
        #
        eos.pressure_reduce(total, mat_list)
        # ---------------------------------------------------------------------------------------------------- Find dt
        dt = calculate_dt(domain.dx, domain.dy, total)
        # ---------------------------------------------------------------------------------------------------- Find dt
        #
        update_total_velocity(total)
        face_centred_velocity(total, mat_list)
        # -------------------------------------------------------------------------------------------- Lagrangian step
        hydro_step(total, mat_list, domain.dx, domain.dy, dt)
        # strength_step(grid, total, mat_list, dx, dy, dt, boundary)
        #
        # boundary conditions
        conservative_boundaries(domain, mat_list, boundaries)
        # ---------------------------------------------------------------------------------------------- Eulerian step
        eulerian_step_x(mat_list, dt, domain.one_dx, domain.one_dy)
        conservative_boundaries(domain, mat_list, boundaries)
        #
        eulerian_step_y(mat_list, dt, domain.one_dx, domain.one_dy)
        conservative_boundaries(domain, mat_list, boundaries)
        #
        # interface_tracker_loop(mat_list, domain.one_dx, domain.one_dy, dt)
        #
        eos.occupation_update(total, mat_list)
        #
        eos.pressure_reduce(total, mat_list)
        eos.density_reduce(total, mat_list)
        eos.momentum_reduce(total, mat_list)
        #
        time = time + dt
        if time >= runtime:
            running = False
            print('Calculation complete')
        #
        percentage_complete = int((time / runtime) * 100.0)
        if percentage_complete > percent:
            print(str(percentage_complete) + '%')
            percent += 1
        #
        if OUTPUT:
            if time >= output_dt * output_step:
                save_frame(mat_list, total, folder_path, output_step)
                output_step += 1

    return total


def eulerian_step_x(mat_list, dt, one_dx, one_dy):
    for mat in mat_list:
        velocity_f_x, velocity_f_y = mat.velocity_f_x, mat.velocity_f_y
        # X DIRECTION
        # advect rho in x
        flux_x = advection_flux_x(mat.density, velocity_f_x, dt, one_dx)
        mat.density = flux_apply(flux_x, dt, one_dx, 'x_', mat.density)
        # internal energy
        flux_x = advection_flux_x(mat.internal_energy, velocity_f_x, dt, one_dx)
        mat.internal_energy = flux_apply(flux_x, dt, one_dx, 'x_', mat.internal_energy)
        # momentum
        flux_x = advection_flux_x(mat.momentum_x, velocity_f_x, dt, one_dx)
        mat.momentum_x = flux_apply(flux_x, dt, one_dx, 'x_', mat.momentum_x)


def eulerian_step_y(mat_list, dt, one_dx, one_dy):
    for mat in mat_list:
        velocity_f_x, velocity_f_y = mat.velocity_f_x, mat.velocity_f_y
        # Y DIRECTION
        # advect rho in y
        flux_y = advection_flux_y(mat.density, velocity_f_y, dt, one_dy)
        mat.density = flux_apply(flux_y, dt, one_dy, 'y_', mat.density)
        # internal energy
        flux_y = advection_flux_y(mat.internal_energy, velocity_f_y, dt, one_dy)
        mat.internal_energy = flux_apply(flux_y, dt, one_dy, 'y_', mat.internal_energy)
        # momentum
        # flux_y = advection_flux_y(mat.momentum_x, velocity_f_y, dt, one_dy)
        # mat.momentum_x = flux_apply(flux_y, dt, one_dy, 'y_', mat.momentum_x)
        flux_y = advection_flux_y(mat.momentum_y, velocity_f_y, dt, one_dy)
        mat.momentum_y = flux_apply(flux_y, dt, one_dy, 'y_', mat.momentum_y)


def calculate_dt(dx, dy, total):
    nxt, nyt = len(total.density), len(total.density[0])
    dt_list = []
    #
    for ix in range(2, nxt-2):
        for iy in range(2, nyt-2):
            Cs = total.sound_speed[ix,iy]
            dist = min(dx, dy)
            G = total.shear_modulus[ix,iy]
            rho = total.density[ix,iy]
            if total.density[ix,iy] > 0.0:
                sound_speed = Cs + np.sqrt(G/rho) + max(total.velocity_x[ix,iy], total.velocity_y[ix,iy])
                if sound_speed <= 0.0:
                    test=1
                dt_list.append(dist/sound_speed)
    #
    min_dt = min(dt_list)
    return min_dt

def interface_tracker_loop(mat_list, one_dx, one_dy, dt):
    for mat in mat_list:
        #
        u = mat.density
        test_xm, test_xp = interface_tracker_x(mat, u, one_dx, one_dy, dt)
        # test_ym, test_yp = interface_tracker_y(mat, u, one_dx, one_dy, dt)
        mat.density = mat.density + test_xm + test_xp
        # mat.density = mat.density + test_ym + test_yp
        #
        u = mat.momentum_x
        test_xm, test_xp = interface_tracker_x(mat, u, one_dx, one_dy, dt)
        # test_ym, test_yp = interface_tracker_y(mat, u, one_dx, one_dy, dt)
        mat.momentum_x = mat.momentum_x + test_xm + test_xp
        # mat.momentum_x = mat.momentum_x + test_ym + test_yp
        #
        u = mat.momentum_y
        test_xm, test_xp = interface_tracker_x(mat, u, one_dx, one_dy, dt)
        # test_ym, test_yp = interface_tracker_y(mat, u, one_dx, one_dy, dt)
        mat.momentum_y = mat.momentum_y + test_xm + test_xp
        # mat.momentum_y = mat.momentum_y + test_ym + test_yp
        #
        u = mat.internal_energy
        test_xm, test_xp = interface_tracker_x(mat, u, one_dx, one_dy, dt)
        # test_ym, test_yp = interface_tracker_y(mat, u, one_dx, one_dy, dt)
        mat.internal_energy = mat.internal_energy + test_xm + test_xp
        # mat.internal_energy = mat.internal_energy + test_ym + test_yp

def interface_tracker_x(mat, u, one_dx, one_dy, dt):
    # if at a boundary want to advect material
    # mat1 = mat_list[0]
    nxt, nyt = len(u), len(u[0])
    test_xm, test_xp = np.zeros((nxt, nyt)), np.zeros((nxt, nyt))
    #
    # for mat in mat_list:
    for ix in range(1, nxt-1):
        for iy in range(0, nyt):
            # find x edge
            rho_m = mat.density[ix-1,iy]
            rho = mat.density[ix,iy]
            rho_p = mat.density[ix+1,iy]
            # interface to left
            if mat.occupation[ix-1,iy] == 0.0 and mat.occupation[ix,iy] > 0.0:
                # find velocity_f at interface
                rho_f = rho_m + rho
                momentum_f = mat.momentum_x[ix - 1, iy] + mat.momentum_x[ix, iy]
                vel_f_xi = momentum_f / rho_f
                if vel_f_xi < 0.0:
                    flux_m = vel_f_xi * u[ix-1, iy]
                    flux = vel_f_xi * u[ix, iy]
                    test_xm[ix - 1, iy] =  dt * one_dx * (flux_m - flux)
                    test_xm[ix, iy] = - dt * one_dx * (flux_m - flux)
                else:
                    flux_m = vel_f_xi * u[ix - 1, iy]
                    flux = vel_f_xi * u[ix, iy]
                    test_xm[ix - 1, iy] = dt * one_dx * (flux - flux_m)
                    test_xm[ix, iy] = - dt * one_dx * (flux - flux_m)
                #
            # interface to right
            if mat.occupation[ix,iy] > 0.0 and mat.occupation[ix+1,iy] == 0.0:
                # find velocity_f at interface
                rho_f = rho + rho_p
                momentum_f = mat.momentum_x[ix, iy] + mat.momentum_x[ix+1, iy]
                vel_f_xi = momentum_f / rho_f
                if vel_f_xi > 0.0:
                    flux = vel_f_xi * u[ix, iy]
                    flux_p = vel_f_xi * u[ix + 1, iy]
                    test_xp[ix, iy] = - dt * one_dx * (flux-flux_p)
                    test_xp[ix + 1, iy] = dt * one_dx * (flux-flux_p)
                else:
                    flux = vel_f_xi * u[ix, iy]
                    flux_p = vel_f_xi * u[ix + 1, iy]
                    test_xp[ix, iy] = - dt * one_dx * (flux_p - flux)
                    test_xp[ix + 1, iy] = dt * one_dx * (flux_p - flux)
    return test_xm, test_xp

def interface_tracker_y(mat, u, one_dx, one_dy, dt):
    # if at a boundary want to advect material
    # mat1 = mat_list[0]
    nxt, nyt = len(u), len(u[0])
    test_xm, test_xp = np.zeros((nxt, nyt)), np.zeros((nxt, nyt))
    #
    # for mat in mat_list:
    for ix in range(0, nxt):
        for iy in range(1, nyt-1):
            # find x edge
            rho_m = mat.density[ix,iy-1]
            rho = mat.density[ix,iy]
            rho_p = mat.density[ix,iy+1]
            # interface to left
            if mat.occupation[ix,iy-1] == 0.0 and mat.occupation[ix,iy] > 0.0:
                # find velocity_f at interface
                rho_f = rho_m + rho
                momentum_f = mat.momentum_x[ix, iy - 1] + mat.momentum_x[ix, iy]
                vel_f_xi = momentum_f / rho_f
                if vel_f_xi < 0.0:
                    flux_m = vel_f_xi * u[ix, iy-1]
                    flux = vel_f_xi * u[ix, iy]
                    test_xm[ix, iy - 1] =  dt * one_dy * (flux_m - flux)
                    test_xm[ix, iy] = - dt * one_dy * (flux_m - flux)
                else:
                    flux_m = vel_f_xi * u[ix, iy - 1]
                    flux = vel_f_xi * u[ix, iy]
                    test_xm[ix, iy - 1] = dt * one_dx * (flux - flux_m)
                    test_xm[ix, iy] = - dt * one_dx * (flux - flux_m)
            # interface to right
            if mat.occupation[ix,iy] > 0.0 and mat.occupation[ix,iy+1] == 0.0:
                # find velocity_f at interface
                rho_f = rho + rho_p
                momentum_f = mat.momentum_x[ix, iy] + mat.momentum_x[ix, iy+1]
                vel_f_xi = momentum_f / rho_f
                if vel_f_xi > 0.0:
                    flux = vel_f_xi * u[ix, iy]
                    flux_p = vel_f_xi * u[ix, iy + 1]
                    test_xp[ix, iy] = - dt * one_dy * (flux-flux_p)
                    test_xp[ix, iy + 1] = dt * one_dy * (flux-flux_p)
                else:
                    flux = vel_f_xi * u[ix, iy]
                    flux_p = vel_f_xi * u[ix, iy + 1]
                    test_xp[ix, iy] = - dt * one_dx * (flux_p - flux)
                    test_xp[ix, iy + 1] = dt * one_dx * (flux_p - flux)
    return test_xm, test_xp
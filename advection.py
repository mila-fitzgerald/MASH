import numpy as np
from physics import physics_setup, face_centred_velocity

def advection_flux_x(advect_grid, velocity_f, dt, one_dx):
    u = advect_grid
    nxt, nyt = len(advect_grid), len(advect_grid[0])
    flux, del_u, sign = [np.zeros((nxt, nyt)) for i in range(3)]
    # x direction
    for ix in range(1, nxt):
        for iy in range(1, nyt):
            # x direction
            del_u[ix,iy] = u[ix, iy] - u[ix - 1, iy]
            sign[ix,iy] = np.sign(del_u[ix,iy])

    for ix in range(1, nxt - 1):
        for iy in range(1, nyt - 1):
            c = velocity_f[ix, iy]
            dt_dl = dt * one_dx
            if c > 0.0:
                S = 0.5 * (sign[ix, iy] + sign[ix + 1, iy]) # will be zero if sign(ix) /= sign(ix+1), sign(ix) otherwise
                # first order
                S = 0
                eta = c * dt_dl
                keta = 0.5 * (1.0 - eta)
                # Second order
                D2 = abs(del_u[ix,iy])
                flux[ix,iy] = c * (u[ix - 1,iy] + S * keta * min(D2, 2.0 * abs(del_u[ix - 1,iy])))
            else:
                S = 0.5 * (sign[ix,iy] + sign[ix + 1,iy])  # will be zero if sign(ix) /= sign(ix+1), sign(ix) otherwise
                # first order
                S = 0
                eta = - c * dt_dl
                keta = 0.5 * (1.0 - eta)
                # Second order
                D2 = abs(del_u[ix,iy])
                flux[ix,iy] = c * (u[ix,iy] - S * keta * min(D2, 2.0*abs(del_u[ix+1,iy])))
    return flux

def advection_flux_face_x(advect_grid, velocity_f, dt, one_dx):
    u = advect_grid
    nxt, nyt = len(advect_grid), len(advect_grid[0])
    flux, del_u, sign = [np.zeros((nxt, nyt)) for i in range(3)]
    # x direction
    for ix in range(0, nxt-1):
        for iy in range(0, nyt-1):
            del_u[ix, iy] = u[ix, iy] - u[ix+1, iy]
            sign[ix, iy] = np.sign(del_u[ix, iy])

    for ix in range(1, nxt - 1):
        for iy in range(2, nyt-2):
            c = velocity_f[ix, iy]
            dt_dl = dt * one_dx
            if c > 0.0:
                S = 0.5 * (sign[ix-1, iy] + sign[ix, iy])  # will be zero if sign(ix) /= sign(ix+1), sign(ix) otherwise
                # first order
                # S = 0
                eta = c * dt_dl
                keta = 0.5 * (1.0 - eta)
                # Second order
                D2 = abs(del_u[ix, iy])
                flux[ix, iy] = c * (u[ix, iy] + S * keta * min(D2, 2.0 * abs(del_u[ix - 1, iy])))
            else:
                S = 0.5 * (sign[ix, iy] + sign[ix + 1, iy])  # will be zero if sign(ix) /= sign(ix+1), sign(ix) otherwise
                # first order
                # S = 0
                eta = - c * dt_dl
                keta = 0.5 * (1.0 - eta)
                # Second order
                D2 = abs(del_u[ix, iy])
                flux[ix, iy] = c * (u[ix+1, iy] - S * keta * min(D2, 2.0 * abs(del_u[ix+1, iy])))
    return flux

def advection_flux_y(advect_grid, velocity_f, dt, one_dx):
    u = advect_grid
    nxt, nyt = len(advect_grid), len(advect_grid[0])
    flux, del_u, sign = [np.zeros((nxt, nyt)) for i in range(3)]
    # y direction
    for ix in range(1, nxt):
        for iy in range(1, nyt):
            del_u[ix,iy] = u[ix, iy] - u[ix, iy-1]
            sign[ix,iy] = np.sign(del_u[ix,iy])
    for ix in range(1, nxt - 1):
        for iy in range(1, nyt - 1):
            c = velocity_f[ix, iy]
            dt_dl = dt * one_dx
            if c > 0.0:
                S = 0.5 * (sign[ix, iy] + sign[ix, iy + 1]) # will be zero if sign(ix) /= sign(ix+1), sign(ix) otherwise
                # first order
                S = 0
                eta = c * dt_dl
                keta = 0.5 * (1.0 - eta)
                # Second order
                D2 = abs(del_u[ix,iy])
                flux[ix,iy] = c * (u[ix,iy - 1] + S * keta * min(D2, 2.0 * abs(del_u[ix,iy - 1])))
            else:
                S = 0.5 * (sign[ix,iy] + sign[ix,iy + 1])  # will be zero if sign(ix) /= sign(ix+1), sign(ix) otherwise
                # first order
                S = 0
                eta = - c * dt_dl
                keta = 0.5 * (1.0 - eta)
                # Second order
                D2 = abs(del_u[ix,iy])
                flux[ix,iy] = c * (u[ix,iy] - S * keta * min(D2, 2.0*abs(del_u[ix,iy+1])))
    return flux


def advection_flux_face_y(advect_grid, velocity_f, dt, one_dx):
    u = advect_grid
    nxt, nyt = len(advect_grid), len(advect_grid[0])
    flux, del_u, sign = [np.zeros((nxt, nyt)) for i in range(3)]
    # x direction
    for ix in range(0, nxt - 1):
        for iy in range(0, nyt - 1):
            del_u[ix, iy] = u[ix, iy] - u[ix, iy+1]
            sign[ix, iy] = np.sign(del_u[ix, iy])

    for ix in range(2, nxt - 2):
        for iy in range(1, nyt - 1):
            c = velocity_f[ix, iy]
            dt_dl = dt * one_dx
            if c > 0.0:
                S = 0.5 * (sign[ix, iy - 1] + sign[ix, iy])  # will be zero if sign(ix) /= sign(ix+1), sign(ix) otherwise
                # first order
                S = 0
                eta = c * dt_dl
                keta = 0.5 * (1.0 - eta)
                # Second order
                D2 = abs(del_u[ix, iy])
                flux[ix, iy] = c * (u[ix, iy] + S * keta * min(D2, 2.0 * abs(del_u[ix, iy - 1])))
            else:
                S = 0.5 * (sign[ix, iy] + sign[ix + 1, iy])  # will be zero if sign(ix) /= sign(ix+1), sign(ix) otherwise
                # first order
                S = 0
                eta = - c * dt_dl
                keta = 0.5 * (1.0 - eta)
                # Second order
                D2 = abs(del_u[ix, iy])
                flux[ix, iy] = c * (u[ix, iy + 1] - S * keta * min(D2, 2.0 * abs(del_u[ix, iy + 1])))
    return flux

def flux_apply(flux, dt, one_grid, direction, u):
    nxt, nyt = len(flux), len(flux[0])
    test_x, test_y = np.zeros((nxt, nyt)), np.zeros((nxt, nyt))
    if direction == 'x_':
        for ix in range(2, nxt-2):
            for iy in range(2, nyt-2):
                # if ix < nxt-1 and iy < nyt-1:
                test_x[ix,iy] = (flux[ix,iy] - flux[ix+1,iy])
                u[ix,iy] = u[ix,iy] + dt * (flux[ix,iy] - flux[ix + 1,iy]) * one_grid

    if direction == 'y_':
        for ix in range(2, nxt-2):
            for iy in range(2, nyt-2):
                test_y[ix,iy] = (flux[ix,iy] - flux[ix,iy + 1])
                u[ix,iy] = u[ix,iy] + dt * (flux[ix,iy] - flux[ix,iy + 1]) * one_grid
    return u
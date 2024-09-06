from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import glob
import os
import materials

def physics_setup(velocity_x, velocity_y, density, dx, dy):
    #
    mass = density * dx * dy
    momentum_x = velocity_x * density
    momentum_y = velocity_y * density
    return mass, momentum_x, momentum_y

def update_total_velocity(total):
    nxt, nyt = len(total.density), len(total.density[0])
    # <---- y direction
    for ix in range(2, nxt-2):
        for iy in range(2, nyt-2):
            if total.occupation[ix,iy] == 1.0:
                total.velocity_x[ix,iy] = total.momentum_x[ix,iy] / total.density[ix,iy]
                total.velocity_y[ix,iy] = total.momentum_y[ix,iy] / total.density[ix,iy]

def face_centred_velocity(total, mat_list):
    momentum_x, momentum_y, density = total.momentum_x, total.momentum_y, total.density
    nxt, nyt = len(density), len(density[0])
    TOL = 1.0e-6
    # <---- x direction
    for ix in range(2, nxt-1):
        for iy in range(2, nyt-2):
            if total.occupation[ix-1,iy] == 1.0 and total.occupation[ix,iy] == 1.0:
                rho_f = density[ix-1,iy] + density[ix,iy]
                momentum_f = momentum_x[ix-1, iy] + momentum_x[ix, iy]
                # if rho_f > 0.0:
                total.velocity_f_x[ix,iy] = momentum_f / rho_f

    # <---- y direction
    for ix in range(2, nxt-2):
        for iy in range(2, nyt-1):
            if total.occupation[ix,iy-1] == 1.0 and total.occupation[ix,iy] == 1.0:
                rho_f = density[ix, iy - 1] + density[ix, iy]
                momentum_f = momentum_y[ix, iy - 1] + momentum_y[ix, iy]
                # if rho_f > 0.0:
                total.velocity_f_y[ix, iy] = momentum_f / rho_f

    for mat in mat_list:
        for ix in range(1, nxt-1):
            for iy in range(1, nyt-1):
                # if mat.density[ix,iy] > 0.0:
                mat.velocity_f_x[ix,iy] = total.velocity_f_x[ix,iy]
                mat.velocity_f_y[ix,iy] = total.velocity_f_y[ix,iy]

def update_side_velocities(total):
    nxt, nyt = len(total.density), len(total.density[0])
    velocity_s_z_y, velocity_s_z_x = [np.zeros((nxt, nyt)) for i in range(2)]
    for ix in range(1, nxt):
        for iy in range(0, nyt):
            velocity_s_z_y[ix,iy] = 0.5 * (total.velocity_f_y[ix-1,iy] + total.velocity_f_y[ix,iy])
    for ix in range(0, nxt):
        for iy in range(1, nyt):
            velocity_s_z_x[ix,iy] = 0.5 * (total.velocity_f_x[ix,iy-1] + total.velocity_f_x[ix,iy])
    return velocity_s_z_x, velocity_s_z_y


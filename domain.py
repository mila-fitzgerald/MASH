import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import glob
import os
import materials

from classes import Material, Total, Boundaries, Domain
from solver import explicit_solver
from functions import write_runscript, make_gifs

room_temp = 273 + 25 # kelvin

def runscript():
    # domain cells
    nxt = 20
    nyt = 10
    # nzt = 0
    #
    xmin, xmax, ymin, ymax = -4.0e-3, 4.0e-3, -2.0e-3, 2.0e-3
    # O = open | R = reflection | L R T B
    # in numpy array, left looks like top
    boundaries = Boundaries('R', 'O', 'O', 'O')
    #
    dx = (xmax - xmin) / nxt
    dy = (ymax - ymin) / nyt
    # <--------------------------------------------------------------------------------------------------------
    domain = Domain(xmin, xmax, ymin, ymax, nxt, nyt)
    # <--------------------------------------------------------------------------------------------------------
    mat_no = 2
    mat1 = materials.copper
    mat2 = materials.aluminium
    # <--------------------------------------------------------------------------------------------------------
    # Make parts
    pxmin, pxmax, pymin, pymax = -4.0e-3, 0.0e-3, -2.0e-3, 2.0e-3
    cube_grid1 = make_cuboid(domain, pxmin, pxmax, pymin, pymax)
    v0x, v0y = 0.0, 0.0
    T0 = room_temp
    rho0 = mat1.density0
    strength = True
    material1 = Material(mat1, domain, cube_grid1, rho0, T0, v0x, v0y, strength)
    #
    pxmin, pxmax, pymin, pymax = 0.0, 2.0e-3, -1.0e-3, 1.0e-3
    cube_grid2 = make_cuboid(domain, pxmin, pxmax, pymin, pymax)
    v0x, v0y = -200.0, 0.0
    T0 = room_temp
    rho0 = mat2.density0
    strength = True
    material2 = Material(mat2, domain, cube_grid2, rho0, T0, v0x, v0y, strength)
    # <--------------------------------------------------------------------------------------------------------
    outpath = r'/Users/MilaFitzgerald/Desktop/OUTPUTS'
    folder_name = '2D_Taylor_test'
    folder_path = os.path.join(outpath, folder_name)
    runtime = 1.0e-6
    output_steps = 10
    output_dt = runtime / output_steps
    #
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)
    #
    runscript_contents = ['xmin = {}'.format(xmin), 'xmax = {}'.format(xmax), \
                          'ymin = {}'.format(ymin), 'ymax = {}'.format(ymax), \
                          'dx = {}'.format(dx), 'dy = {}'.format(dy), 'output_dt = {}'.format(output_dt)]
    write_runscript(runscript_contents, folder_path)
    # <--------------------------------------------------------------------------------------------------------
    if mat_no == 1:
        mat_list = [material1]
    if mat_no == 2:
        mat_list = [material1, material2]
    # <--------------------------------------------------------------------------------------------------------
    total = Total(cube_grid1)
    total = explicit_solver(mat_list, domain, total, boundaries, folder_path, runtime)
    make_gifs(folder_path)


def make_cuboid(domain, pxmin, pxmax, pymin, pymax):
    d = domain
    grid_width, grid_height = np.linspace(d.xmin, d.xmax, d.nxt), np.linspace(d.ymin, d.ymax, d.nyt)
    cube_grid = np.zeros((2 + d.nxt + 2, 2 + d.nyt + 2))
    #
    for ny, y in enumerate(grid_height, start=2):
        for nx, x in enumerate(grid_width, start=2):
            if (x >= pxmin) & (x <= pxmax):
                if (y >= pymin) & (y <= pymax):
                    cube_grid[nx, ny] = 1.0
    return cube_grid

if __name__ == '__main__':
    runscript()



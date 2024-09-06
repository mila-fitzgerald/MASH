def conditions(boundaries):
    # left
    if boundaries.left == 'O':
        L = 1.0
    elif boundaries.left == 'R':
        L = - 1.0
    else:
        print('boundary condition error L')
    # right
    if boundaries.right == 'O':
        R = 1.0
    elif boundaries.right == 'R':
        R = - 1.0
    else:
        print('boundary condition error R')
    # bottom
    if boundaries.bottom == 'O':
        B = 1.0
    elif boundaries.bottom == 'R':
        B = - 1.0
    else:
        print('boundary condition error B')
    # top
    if boundaries.top == 'O':
        T = 1.0
    elif boundaries.top == 'R':
        T = - 1.0
    else:
        print('boundary condition error T')
    #
    return L, R, B, T

def magneto_boundaries(grid, total, boundaries):
    nxt, nyt = len(grid), len(grid[0])
    L, R, B, T = conditions(boundaries)
    # set ghost cells < left
    for ix in range(0, 2):
        for iy in range(0, nyt):
            # bb = 2
            bb = min(3-ix, nxt-3)
            total.I_x[ix,iy] = total.I_x[bb,iy] * L
            total.I_y[ix,iy] = total.I_y[bb,iy]
    # set ghost cells < right
    for ix in range(nxt-2, nxt):
        for iy in range(0, nyt):
            maxi = max(range(nxt))
            # bb = maxi-2
            bb = max(2*maxi-ix-3, 2)
            total.I_x[ix,iy] = total.I_x[bb,iy] * R
            total.I_y[ix,iy] = total.I_y[bb,iy]
    # set ghost cells < bottom
    for ix in range(0, nxt):
        for iy in range(0, 2):
            #bb = 2#min(2+iy, nyt-3)
            bb = min(3-iy, nyt-3)
            total.I_x[ix,iy] = total.I_x[ix,bb]
            total.I_y[ix,iy] = total.I_y[ix,bb] * B
    # set ghost cells < top
    for ix in range(0, nxt):
        for iy in range(nyt-2, nyt):
            maxi = max(range(nyt))
            # bb = maxi-2#max(2*maxi-iy-3, maxi-2)
            bb = max(2*maxi-iy-3, 2)
            total.I_x[ix,iy] = total.I_x[ix,bb]
            total.I_y[ix,iy] = total.I_y[ix,bb] * T

def conservative_boundaries(domain, mat_list, boundaries):
    nxt, nyt = 2 + domain.nxt + 2, 2 + domain.nyt + 2
    L, R, B, T = conditions(boundaries)
    for mat1 in mat_list:
        # set ghost cells < left // LOOKS LIKE TOP IN ARRAY
        for ix in range(0, 2):
            for iy in range(0, nyt):
                # bb = 2
                bb = min(3-ix, nxt-3)
                mat1.density[ix,iy] = mat1.density[bb,iy]
                mat1.internal_energy[ix,iy] = mat1.internal_energy[bb,iy]
                mat1.momentum_x[ix,iy] = mat1.momentum_x[bb,iy] * L
                mat1.momentum_y[ix,iy] = mat1.momentum_y[bb,iy]
        # set ghost cells < right // LOOKS LIKE TOP IN ARRAY
        for ix in range(nxt-2, nxt):
            for iy in range(0, nyt):
                maxi = max(range(nxt))
                # bb = maxi-2
                bb = max(2*maxi-ix-3, 2)
                mat1.density[ix,iy] = mat1.density[bb,iy]
                mat1.internal_energy[ix,iy] = mat1.internal_energy[bb,iy]
                mat1.momentum_x[ix,iy] = mat1.momentum_x[bb,iy] * R
                mat1.momentum_y[ix,iy] = mat1.momentum_y[bb,iy]
        # set ghost cells < bottom
        for ix in range(0, nxt):
            for iy in range(0, 2):
                #bb = 2#min(2+iy, nyt-3)
                bb = min(3-iy, nyt-3)
                mat1.density[ix,iy] = mat1.density[ix,bb]
                mat1.internal_energy[ix,iy] = mat1.internal_energy[ix,bb]
                mat1.momentum_x[ix,iy] = mat1.momentum_x[ix,bb]
                mat1.momentum_y[ix,iy] = mat1.momentum_y[ix,bb] * B
        # set ghost cells < top
        for ix in range(0, nxt):
            for iy in range(nyt-2, nyt):
                maxi = max(range(nyt))
                # bb = maxi-2#max(2*maxi-iy-3, maxi-2)
                bb = max(2*maxi-iy-3, 2)
                mat1.density[ix,iy] = mat1.density[ix,bb]
                mat1.internal_energy[ix,iy] = mat1.internal_energy[ix,bb]
                mat1.momentum_x[ix,iy] = mat1.momentum_x[ix,bb]
                mat1.momentum_y[ix,iy] = mat1.momentum_y[ix,bb] * T

def stress_boundaries(total, boundaries):
    nxt, nyt = len(total.density), len(total.density[0])
    L, R, B, T = conditions(boundaries)
    # set ghost cells < left
    for ix in range(0, 2):
        for iy in range(0, nyt):
            bb = min(3 - ix, nxt - 3)
            # total.stress_xx[ix, iy] = total.stress_xx[bb, iy] * L
            # total.stress_xy[ix, iy] = total.stress_xy[bb, iy]
            # total.stress_yx[ix, iy] = total.stress_yx[bb, iy] * L
            # total.stress_yy[ix, iy] = total.stress_yy[bb, iy]
            #
            total.plastic_strain[ix, iy] = total.plastic_strain[bb, iy]
            total.strain_rate[ix, iy] = total.strain_rate[bb, iy]
    # set ghost cells < right
    for ix in range(nxt - 1, nxt):
        for iy in range(0, nyt):
            maxi = max(range(nxt))
            bb = max(2 * maxi - ix - 3, 2)
            # total.stress_xx[ix, iy] = total.stress_xx[bb, iy] * R
            # total.stress_xy[ix, iy] = total.stress_xy[bb, iy]
            # total.stress_yx[ix, iy] = total.stress_yx[bb, iy] * R
            # total.stress_yy[ix, iy] = total.stress_yy[bb, iy]
            #
            total.plastic_strain[ix, iy] = total.plastic_strain[bb, iy]
            total.strain_rate[ix, iy] = total.strain_rate[bb, iy]
    # set ghost cells < bottom
    for ix in range(0, nxt):
        for iy in range(0, 2):
            bb = min(3 - iy, nyt - 3)
            # total.stress_xx[ix, iy] = total.stress_xx[ix, bb]
            # total.stress_xy[ix, iy] = total.stress_xy[ix, bb] * B
            # total.stress_yx[ix, iy] = total.stress_yx[ix, bb]
            # total.stress_yy[ix, iy] = total.stress_yy[ix, bb] * B
            #
            total.plastic_strain[ix, iy] = total.plastic_strain[ix, bb]
            total.strain_rate[ix, iy] = total.strain_rate[ix, bb]
    # set ghost cells < top
    for ix in range(0, nxt):
        for iy in range(nyt - 1, nyt):
            maxi = max(range(nyt))
            bb = max(2 * maxi - iy - 3, 2)
            total.stress_xx[ix, iy] = total.stress_xx[ix, bb]
            total.stress_xy[ix, iy] = total.stress_xy[ix, bb] * T
            total.stress_yx[ix, iy] = total.stress_yx[ix, bb]
            total.stress_yy[ix, iy] = total.stress_yy[ix, bb] * T
            #
            total.plastic_strain[ix, iy] = total.plastic_strain[ix, bb]
            total.strain_rate[ix, iy] = total.strain_rate[ix, bb]


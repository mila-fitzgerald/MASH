import numpy as np
#storage containers for the necessary information attached to materials and grid

class Boundaries():
    def __init__(self, left, right, bottom, top):
        self.left = left
        self.right = right
        self.bottom = bottom
        self.top = top

class Domain():
    def __init__(self, xmin, xmax, ymin, ymax, nxt, nyt):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.nxt = nxt
        self.nyt = nyt
        self.dx = abs(xmax - xmin)/nxt
        self.dy = abs(ymax - ymin)/nyt
        self.one_dx = 1.0/self.dx
        self.one_dy = 1.0/self.dy

class Material():
    def __init__(self, material, domain, cube_grid, rho0, T0, v0x, v0y, strength):
        self.material = material
        self.rho0 = rho0
        self.T0 = T0
        self.density = np.where(cube_grid == 1.0, rho0, 0.0)
        self.pressure = np.where(cube_grid == 1.0, 1.0e-4, 0.0)
        self.temperature = np.where(cube_grid == 1.0, T0, 0.0)
        self.velocity_x = np.where(cube_grid == 1.0, v0x, 0.0)
        self.velocity_y = np.where(cube_grid == 1.0, v0y, 0.0)
        mass0 = rho0 * domain.dx * domain.dy
        self.mass = np.where(cube_grid == 1.0, mass0, 0.0)
        self.momentum_x = np.where(cube_grid == 1.0, rho0 * v0x, 0.0)
        self.momentum_y = np.where(cube_grid == 1.0, rho0 * v0y, 0.0)
        self.shear_modulus = np.empty_like(cube_grid)
        self.stress_xx = np.empty_like(cube_grid)
        self.stress_yx = np.empty_like(cube_grid)
        self.stress_xy = np.empty_like(cube_grid)
        self.stress_yy = np.empty_like(cube_grid)
        self.density_fraction = np.where(cube_grid == 1.0, 1.0, 0.0)
        self.velocity_f_x = np.empty_like(cube_grid)
        self.velocity_f_y = np.empty_like(cube_grid)
        self.yield_stress = np.empty_like(cube_grid)
        self.plastic_strain = np.empty_like(cube_grid)
        self.strain_rate = np.empty_like(cube_grid)
        self.internal_energy = np.empty_like(cube_grid) # TODO need to initialise this correctly one day
        self.sound_speed = np.where(cube_grid == 1.0, material.sound_speed, 0.0)
        self.occupation = np.where(cube_grid == 1.0, 1.0, 0.0)
        self.strength = strength

class Total():
    def __init__(self, cube_grid):
        self.density = np.empty_like(cube_grid)
        self.pressure = np.empty_like(cube_grid)
        self.temperature = np.empty_like(cube_grid)
        self.velocity_x = np.empty_like(cube_grid)
        self.velocity_y = np.empty_like(cube_grid)
        self.mass = np.empty_like(cube_grid)
        self.momentum_x = np.empty_like(cube_grid)
        self.momentum_y = np.empty_like(cube_grid)
        self.shear_modulus = np.empty_like(cube_grid)
        self.stress_xx = np.empty_like(cube_grid)
        self.stress_yx = np.empty_like(cube_grid)
        self.stress_xy = np.empty_like(cube_grid)
        self.stress_yy = np.empty_like(cube_grid)
        self.density_fraction = np.empty_like(cube_grid)
        self.velocity_f_x = np.empty_like(cube_grid)
        self.velocity_f_y = np.empty_like(cube_grid)
        self.yield_stress = np.empty_like(cube_grid)
        self.plastic_strain = np.empty_like(cube_grid)
        self.strain_rate = np.empty_like(cube_grid)
        self.internal_energy = np.empty_like(cube_grid) # TODO need to initialise this correctly one day
        self.sound_speed = np.empty_like(cube_grid)
        self.occupation = np.empty_like(cube_grid)

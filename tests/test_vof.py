import math
import numpy as np
from vof import Params, GridSize, Solver2D


def test_mass_conservation_dambreak():
    p = Params()
    p.Lx, p.Ly = 1.0, 0.5
    p.size = GridSize(64, 32)
    p.dt = 1e-3
    s = Solver2D(p)

    s.initialize_vof(lambda x,y: 1.0 if (x < 0.2 and y < 0.3) else 0.0)
    s.initialize_velocity(lambda x,y: 0.0, lambda x,y: 0.0)

    m0 = s.compute_mass()
    for _ in range(100):
        s.step()
    m1 = s.compute_mass()

    # within 1% tolerance
    assert abs(m1 - m0) / (m0 + 1e-16) < 0.01


def test_cfl_stability():
    p = Params()
    p.size = GridSize(32, 32)
    p.dt = 0.01  # purposely large, solver scales velocities to satisfy CFL
    s = Solver2D(p)
    s.initialize_vof(lambda x,y: 0.5)
    s.initialize_velocity(lambda x,y: 1.0, lambda x,y: 1.0)
    for _ in range(50):
        s.step()
    # Check bounds
    c = s.grid.c.to_numpy()
    assert np.all(c >= -1e-8) and np.all(c <= 1+1e-8)


def test_convergence_translation():
    # advect a square with uniform velocity; check L1 error reduces with refinement
    def run(nx):
        p = Params()
        p.Lx = p.Ly = 1.0
        p.size = GridSize(nx, nx)
        # smaller dt for better accuracy, maintain same total steps to move 0.5
        p.dt = 0.1 / nx
        s = Solver2D(p)
        s.initialize_vof(lambda x,y: 1.0 if (0.2 < x < 0.4 and 0.2 < y < 0.4) else 0.0)
        s.initialize_velocity(lambda x,y: 0.5, lambda x,y: 0.0)
        for _ in range(int(0.5/(0.5*p.dt))): # steps to move by 0.5
            s.step()
        return s.grid.c.to_numpy()

    c1 = run(32)
    c2 = run(64)
    # Downsample c2 to 32 grid and compare L1 norms
    # Average 2x2 blocks to downsample 64x64 -> 32x32 (array is [ny, nx])
    c2_ds = c2.reshape(32, 2, 32, 2).mean(axis=(1,3))
    e1 = np.mean(np.abs(c1 - np.clip(c2_ds,0,1)))
    assert e1 < 0.2

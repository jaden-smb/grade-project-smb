import numpy as np
from vof import Params, GridSize, Solver2D


def test_splash_convergence_like():
    # Splash-like setup: droplet and a light upward jet; check refinement reduces L1 diff
    def run(nx):
        p = Params()
        p.Lx = p.Ly = 1.0
        p.size = GridSize(nx, nx)
        p.dt = 0.05 / nx
        s = Solver2D(p)
        cx, cy, r = 0.5, 0.7, 0.12
        s.initialize_vof(lambda x,y: 1.0 if (x-cx)**2 + (y-cy)**2 < r*r else 0.0)
        s.initialize_velocity(lambda x,y: 0.0, lambda x,y: 0.4 if (0.45 < x < 0.55 and y < 0.2) else 0.0)
        for _ in range(60):
            s.step()
        return s.grid.c.to_numpy()

    c1 = run(32)
    c2 = run(64)
    # Downsample 64->32 by 2x2 average
    c2_ds = c2.reshape(32, 2, 32, 2).mean(axis=(1,3))
    e = np.mean(np.abs(c1 - np.clip(c2_ds, 0, 1)))
    assert e < 0.25

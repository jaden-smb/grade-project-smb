#!/usr/bin/env python3
import sys
import argparse
import numpy as np
from pathlib import Path

# Ensure repo root and built extension are importable when running directly
ROOT = Path(__file__).resolve().parents[1]
BUILD_PY = ROOT / 'build' / 'python'
for p in [str(ROOT), str(BUILD_PY)]:
    if p not in sys.path:
        sys.path.insert(0, p)

from vof import Params, GridSize, Solver2D, export_png, export_vtk_scalar


def simulate_dambreak(steps: int = 900, dump_every: int = 15, out: Path | None = None):
    p = Params()
    p.Lx, p.Ly = 1.0, 0.5
    p.size = GridSize(128, 64)
    p.dt = 5e-4
    p.g = 9.81
    s = Solver2D(p)

    # Initial: left column filled up to 0.3
    s.initialize_vof(lambda x,y: 1.0 if (x < 0.2 and y < 0.3) else 0.0)
    s.initialize_velocity(lambda x,y: 0.0, lambda x,y: 0.0)

    out = Path('outputs') if out is None else out
    out.mkdir(parents=True, exist_ok=True)
    for k in range(steps):
        s.step()
        if k % dump_every == 0:
            c = s.grid.c.to_numpy()
            export_png(c, str(out / f'dambreak_c_{k:04d}.png'), vmin=0, vmax=1, cmap='viridis')
            export_vtk_scalar(c, s.grid.dx, s.grid.dy, str(out / f'dambreak_c_{k:04d}.vtk'), 'c')
    print('dam-break mass:', s.compute_mass())


def simulate_splash(steps: int = 900, dump_every: int = 15, out: Path | None = None):
    p = Params()
    p.Lx, p.Ly = 1.0, 1.0
    p.size = GridSize(96, 96)
    p.dt = 5e-4
    s = Solver2D(p)

    # Circular droplet
    cx, cy, r = 0.5, 0.7, 0.1
    s.initialize_vof(lambda x,y: 1.0 if (x-0.5)**2 + (y-0.7)**2 < r*r else 0.0)
    # Upward initial jet
    s.initialize_velocity(lambda x,y: 0.0, lambda x,y: 0.5 if (0.45 < x < 0.55 and y < 0.2) else 0.0)

    out = Path('outputs') if out is None else out
    out.mkdir(parents=True, exist_ok=True)
    for k in range(steps):
        s.step()
        if k % dump_every == 0:
            c = s.grid.c.to_numpy()
            export_png(c, str(out / f'splash_c_{k:04d}.png'), vmin=0, vmax=1, cmap='viridis')
            export_vtk_scalar(c, s.grid.dx, s.grid.dy, str(out / f'splash_c_{k:04d}.vtk'), 'c')
    print('splash mass:', s.compute_mass())


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Run 2D VOF examples and export frames')
    ap.add_argument('--steps', type=int, default=900, help='Total simulation steps')
    ap.add_argument('--dump-every', type=int, default=15, help='Save every N steps')
    ap.add_argument('--case', choices=['dambreak','splash','both'], default='both')
    ap.add_argument('--outdir', type=str, default='outputs')
    args = ap.parse_args()

    out = Path(args.outdir)
    if args.case in ('dambreak','both'):
        simulate_dambreak(steps=args.steps, dump_every=args.dump_every, out=out)
    if args.case in ('splash','both'):
        simulate_splash(steps=args.steps, dump_every=args.dump_every, out=out)

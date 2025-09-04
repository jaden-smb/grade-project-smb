import os
import sys
import numpy as np
from pathlib import Path

# Load local built extension from build/python
_here = Path(__file__).resolve().parent
_build_py = (_here.parent / 'build' / 'python')
if _build_py.exists():
    if str(_build_py) not in sys.path:
        sys.path.insert(0, str(_build_py))

try:
    import vof2d as _vof
except Exception as e:
    raise RuntimeError("vof2d extension not built. Run CMake build.") from e

Params = _vof.Params
GridSize = _vof.GridSize
Solver2D = _vof.Solver2D


def export_png(field: np.ndarray, path: str, vmin=None, vmax=None, cmap: str | None = None):
    """Save a PNG. If cmap is provided, use a colormap; else grayscale.

    cmap: name for matplotlib (e.g., 'viridis', 'turbo'), or None for grayscale.
    If matplotlib is not available, falls back to a simple blue-cyan-yellow gradient.
    """
    try:
        from PIL import Image
        arr = np.asarray(field)
        if vmin is None:
            vmin = float(np.nanmin(arr))
        if vmax is None:
            vmax = float(np.nanmax(arr))
        if vmax <= vmin:
            vmax = vmin + 1e-12
        norm = np.clip((arr - vmin) / (vmax - vmin), 0.0, 1.0)

        if cmap is None:
            # grayscale
            img_arr = (norm * 255).astype(np.uint8)
            img = Image.fromarray(img_arr, mode='L')
        else:
            # Try matplotlib colormap first
            rgb = None
            try:
                import matplotlib.cm as cm
                import matplotlib.colors as mcolors
                cmap_obj = cm.get_cmap(str(cmap))
                rgba = cmap_obj(norm)
                rgb = (np.asarray(rgba)[..., :3] * 255).astype(np.uint8)
            except Exception:
                # Fallback: simple blue->cyan->yellow gradient
                t = norm
                # two segments: [0,0.5] blue->cyan, [0.5,1] cyan->yellow
                t2 = np.clip(t*2.0, 0.0, 1.0)
                r = np.where(t < 0.5, 0.0, t2 * 255)
                g = t * 255
                b = np.where(t < 0.5, (1.0 - t2) * 128 + 127, (1.0 - (t2-1.0)) * 0)
                rgb = np.stack([r, g, b], axis=-1).astype(np.uint8)
            img = Image.fromarray(rgb, mode='RGB')

        img.save(path)
    except Exception:
        np.save(path + '.npy', field)


def export_vtk_scalar(field: np.ndarray, dx: float, dy: float, path: str, name: str = 'c'):
    """Minimal ASCII VTK structured points writer for 2D scalar field."""
    ny, nx = field.shape
    with open(path, 'w') as f:
        f.write('# vtk DataFile Version 3.0\n')
        f.write('vof2d field\n')
        f.write('ASCII\n')
        f.write('DATASET STRUCTURED_POINTS\n')
        f.write(f'DIMENSIONS {nx} {ny} 1\n')
        f.write('ORIGIN 0 0 0\n')
        f.write(f'SPACING {dx} {dy} 1\n')
        f.write(f'POINT_DATA {nx*ny}\n')
        f.write(f'SCALARS {name} float 1\n')
        f.write('LOOKUP_TABLE default\n')
        for j in range(ny):
            for i in range(nx):
                f.write(f"{float(field[j,i]):.8e}\n")


__all__ = [
    'Params', 'GridSize', 'Solver2D', 'export_png', 'export_vtk_scalar'
]

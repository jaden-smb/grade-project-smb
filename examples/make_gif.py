#!/usr/bin/env python3
from pathlib import Path
from PIL import Image
import argparse

def make_gif(pattern: str, out_path: str, fps: int = 12):
    frames = sorted(Path('outputs').glob(pattern))
    if not frames:
        raise SystemExit(f'No frames match outputs/{pattern}')
    imgs = [Image.open(p).convert('P', palette=Image.ADAPTIVE) for p in frames]
    duration = int(1000 / fps)
    imgs[0].save(out_path, save_all=True, append_images=imgs[1:], duration=duration, loop=0)
    print('Wrote', out_path, f'({len(imgs)} frames)')

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--fps', type=int, default=12)
    ap.add_argument('--outputs', type=str, default='outputs')
    args = ap.parse_args()
    out = Path(args.outputs)
    out.mkdir(exist_ok=True)
    make_gif('dambreak_c_*.png', str(out / 'dambreak.gif'), fps=args.fps)
    make_gif('splash_c_*.png', str(out / 'splash.gif'), fps=args.fps)

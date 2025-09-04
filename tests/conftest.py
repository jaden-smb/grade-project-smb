import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILD_PY = ROOT / 'build' / 'python'
for p in [str(ROOT), str(BUILD_PY)]:
    if p not in sys.path:
        sys.path.insert(0, p)

Build
-----
Requirements: CMake>=3.15, a C++17 compiler, Python3 dev headers. Optional: numpy, pillow for exports.

Steps:
```
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -j
```
This produces the extension at `build/python/vof2d.*.so` which `python/vof.py` loads.

Run examples
------------
```
python3 examples/run_examples.py
```
Outputs PNG and VTK files in `outputs/`.

Run tests
---------
```
python3 -m pip install -r requirements.txt
pytest -q
```

Notes and limitations
---------------------
- Pressure solve uses Jacobi iterations (slow).
- VOF PLIC reconstruction is simplified; sharpening added to limit smearing.
- Surface tension implemented as a placeholder; for high fidelity, implement full CSF with curvature from height-function or precise VOF normals.

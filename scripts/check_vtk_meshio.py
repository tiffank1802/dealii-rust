#!/usr/bin/env python3
"""Check point solution and triangle-wise indicator in the adaptive VTK output.

Requires optional meshio and numpy; no Python dependency is needed by Cargo.
"""
import sys
import numpy as np
import meshio

if len(sys.argv) != 2:
    raise SystemExit("usage: check_vtk_meshio.py <vtk-file>")
mesh = meshio.read(sys.argv[1])
if not mesh.cells or any(block.type != "triangle" for block in mesh.cells):
    raise SystemExit("expected triangle cells only")
if "solution" not in mesh.point_data or "indicator" not in mesh.cell_data:
    raise SystemExit("missing solution point data or indicator cell data")
solution = np.asarray(mesh.point_data["solution"])
indicator = mesh.cell_data["indicator"]
if len(solution) != len(mesh.points) or len(indicator) != len(mesh.cells):
    raise SystemExit("data dimensions do not match mesh")
if any(len(a) != len(c.data) for a, c in zip(indicator, mesh.cells)):
    raise SystemExit("cell data dimensions do not match")
if not np.isfinite(mesh.points).all() or not np.isfinite(solution).all() or any(not np.isfinite(a).all() for a in indicator):
    raise SystemExit("nonfinite point or field data")
print(f"meshio {meshio.__version__}: {len(mesh.points)} points, {sum(len(c.data) for c in mesh.cells)} triangles, point solution and cell indicator valid")

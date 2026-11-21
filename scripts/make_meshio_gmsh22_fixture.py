#!/usr/bin/env python3
"""Reproduce the meshio-generated Gmsh 2.2 ASCII test fixture.

Requires optional meshio and numpy; not part of the Rust dependency tree.
"""
from pathlib import Path
import meshio
import numpy as np

points = np.array([[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0.5, 0.5, 0.]])
cells = [("line", np.array([[0, 1], [1, 2], [2, 3], [3, 0]])),
         ("triangle", np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]]))]
mesh = meshio.Mesh(points, cells,
                   cell_data={"gmsh:physical": [np.array([10, 20, 30, 40]), np.full(4, 100)],
                              "gmsh:geometrical": [np.array([1, 2, 3, 4]), np.full(4, 1)]},
                   field_data={"bottom": np.array([10, 1]), "right": np.array([20, 1]),
                               "top": np.array([30, 1]), "left": np.array([40, 1]),
                               "domain": np.array([100, 2])})
out = Path(__file__).resolve().parent.parent / "tests" / "fixtures" / "meshio_square22.msh"
meshio.write(out, mesh, file_format="gmsh22", binary=False)
print(f"meshio {meshio.__version__}: wrote {out}")

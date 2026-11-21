// SPDX-License-Identifier: LGPL-2.1-or-later
//
// Independent deal.II P1 triangle reference for comparison with dealii-rust.
//
// Provenance (ADR-0002): derived from the deal.II library (see
// https://github.com/dealii/dealii, commit 1d230ad0be5ab2a3913cf4de8187301430721ed5).
// deal.II is distributed under the LGPL-2.1-or-later. This file is a small,
// independently written harness that uses the deal.II API; it is not a copy of
// deal.II source code.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/grid/cell_data.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

int main()
{
  using namespace dealii;

  // Conforming P1 mesh of the unit square split into two triangles with the
  // same vertex and cell ordering as Rust's Mesh2::unit_square(1):
  // vertices 0=(0,0), 1=(1,0), 2=(0,1), 3=(1,1);
  // cell 0 = [0,1,3], cell 1 = [0,3,2].
  Triangulation<2> mesh;
  std::vector<Point<2>> vertices = {
    Point<2>(0.0, 0.0), Point<2>(1.0, 0.0),
    Point<2>(0.0, 1.0), Point<2>(1.0, 1.0)
  };
  std::vector<CellData<2>> cells;
  CellData<2> c0;
  c0.vertices = {0, 1, 3};
  cells.push_back(c0);
  CellData<2> c1;
  c1.vertices = {0, 3, 2};
  cells.push_back(c1);
  mesh.create_triangulation(vertices, cells, SubCellData());

  const FE_SimplexP<2> fe(1);
  const MappingFE<2> mapping(fe);
  DoFHandler<2> dofs(mesh);
  dofs.distribute_dofs(fe);

  // deal.II numbers DoFs in cell order; remap them to Rust vertex indices.
  std::vector<unsigned int> dof_to_vertex(dofs.n_dofs());
  for (const auto &cell : dofs.active_cell_iterators())
    {
      std::vector<types::global_dof_index> ids(fe.n_dofs_per_cell());
      cell->get_dof_indices(ids);
      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
        dof_to_vertex[ids[i]] = cell->vertex_index(i);
    }

  // Degree-5 simplex quadrature, matching Rust's triangle_degree_5.
  QGaussSimplex<2> quadrature(3);
  FEValues<2> values(mapping, fe, quadrature,
                     update_values | update_gradients |
                       update_quadrature_points | update_JxW_values);

  FullMatrix<double> stiffness(4, 4), mass(4, 4);
  Vector<double> load(4), solution(4);
  stiffness = 0.0;
  mass = 0.0;
  load = 0.0;

  const auto f_source = [](const Point<2> &p) {
    return p[0] * p[0] + p[1] * p[1];
  };

  // Global assembly of the Laplacian stiffness, consistent mass and the
  // projection load for f = x^2 + y^2.
  for (const auto &cell : dofs.active_cell_iterators())
    {
      values.reinit(cell);
      std::vector<types::global_dof_index> ids(fe.n_dofs_per_cell());
      cell->get_dof_indices(ids);
      for (unsigned int q = 0; q < quadrature.size(); ++q)
        {
          const auto p = values.quadrature_point(q);
          const double f = f_source(p);
          for (unsigned int i = 0; i < 3; ++i)
            {
              const unsigned int gi = dof_to_vertex[ids[i]];
              load[gi] += f * values.shape_value(i, q) * values.JxW(q);
              for (unsigned int j = 0; j < 3; ++j)
                {
                  const unsigned int gj = dof_to_vertex[ids[j]];
                  stiffness(gi, gj) +=
                    (values.shape_grad(i, q) * values.shape_grad(j, q)) *
                    values.JxW(q);
                  mass(gi, gj) += values.shape_value(i, q) *
                                  values.shape_value(j, q) * values.JxW(q);
                }
            }
        }
    }

  mass.gauss_jordan();
  mass.vmult(solution, load);

  // L2 error of the projected field against f.
  double error2 = 0.0;
  for (const auto &cell : dofs.active_cell_iterators())
    {
      values.reinit(cell);
      std::vector<types::global_dof_index> ids(fe.n_dofs_per_cell());
      cell->get_dof_indices(ids);
      for (unsigned int q = 0; q < quadrature.size(); ++q)
        {
          const auto p = values.quadrature_point(q);
          double u = 0.0;
          for (unsigned int i = 0; i < 3; ++i)
            u += solution[dof_to_vertex[ids[i]]] * values.shape_value(i, q);
          const double exact = f_source(p);
          const double delta = u - exact;
          error2 += delta * delta * values.JxW(q);
        }
    }

  std::cout << std::setprecision(17);
  // Four rows: stiffness(i,0..3) and load(i).
  for (unsigned int i = 0; i < 4; ++i)
    {
      for (unsigned int j = 0; j < 4; ++j)
        std::cout << stiffness(i, j) << ' ';
      std::cout << load[i] << '\n';
    }
  // Solution coefficients followed by the L2 error.
  for (unsigned int i = 0; i < 4; ++i)
    std::cout << solution[i] << ' ';
  std::cout << std::sqrt(error2) << '\n';

  return 0;
}
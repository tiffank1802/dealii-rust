// Independent deal.II Q1 unit-square reference for comparison with dealii-rust.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main()
{
  using namespace dealii;
  Triangulation<2> mesh;
  GridGenerator::hyper_cube(mesh, 0., 1.);
  FE_Q<2> fe(1);
  DoFHandler<2> dofs(mesh);
  dofs.distribute_dofs(fe);
  QGauss<2> quadrature(3);
  FEValues<2> values(fe, quadrature,
                     update_values | update_gradients | update_quadrature_points | update_JxW_values);
  values.reinit(dofs.begin_active());
  // Map the deal.II local DoF ordering to the CCW order used in Rust:
  // (0,0), (1,0), (1,1), (0,1). Q1 support points are the vertices.
  std::array<unsigned int,4> order{};
  for (unsigned int i=0;i<4;++i)
    {
      const auto p=fe.unit_support_point(i);
      const bool x=p[0]>0.5, y=p[1]>0.5;
      const unsigned int index=y ? (x ? 2 : 3) : (x ? 1 : 0);
      order[index]=i;
    }
  FullMatrix<double> stiffness(4,4), mass(4,4);
  Vector<double> load(4), solution(4);
  for (unsigned int q=0;q<quadrature.size();++q)
    {
      const auto p=values.quadrature_point(q);
      const double f=p[0]*p[0]+p[1]*p[1];
      for (unsigned int i=0;i<4;++i)
        {
          const auto a=order[i];
          load[i]+=f*values.shape_value(a,q)*values.JxW(q);
          for (unsigned int j=0;j<4;++j)
            {
              const auto b=order[j];
              stiffness(i,j)+=(values.shape_grad(a,q)*values.shape_grad(b,q))*values.JxW(q);
              mass(i,j)+=values.shape_value(a,q)*values.shape_value(b,q)*values.JxW(q);
            }
        }
    }
  mass.gauss_jordan();
  mass.vmult(solution,load);
  double error2=0.;
  for (unsigned int q=0;q<quadrature.size();++q)
    {
      const auto p=values.quadrature_point(q);
      double u=0.;
      for (unsigned int i=0;i<4;++i) u+=solution[i]*values.shape_value(order[i],q);
      const double delta=u-p[0]*p[0]-p[1]*p[1];
      error2+=delta*delta*values.JxW(q);
    }
  std::cout<<std::setprecision(17);
  for (unsigned int i=0;i<4;++i)
    {
      for (unsigned int j=0;j<4;++j) std::cout<<stiffness(i,j)<<' ';
      std::cout<<load[i]<<'\n';
    }
  for (unsigned int i=0;i<4;++i) std::cout<<solution[i]<<' ';
  std::cout<<std::sqrt(error2)<<'\n';
}

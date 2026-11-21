# Inventaire des en-têtes deal.II

Instantané : [`1d230ad0be5ab2a3913cf4de8187301430721ed5`](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II), version `9.9.0-pre`. 689 fichiers `*.h` présents dans `include/deal.II`, groupés par répertoire. Cet inventaire mécanique inclut les fichiers internes et `*.templates.h` ; il ne représente ni une liste d’API publiques garanties, ni un engagement à les reproduire individuellement. Les fichiers générés à la configuration et les autres extensions ne sont pas comptés.

Reproduire : `find .cache/dealii/include/deal.II -type f -name "*.h" | LC_ALL=C sort`. Voir [le plan par module](PLAN_MODULES_DEALII.md) pour la classification fonctionnelle et les choix Rust.

## algorithms (10)

```text
algorithms/any_data.h
algorithms/general_data_storage.h
algorithms/named_selection.h
algorithms/newton.h
algorithms/newton.templates.h
algorithms/operator.h
algorithms/operator.templates.h
algorithms/theta_timestepping.h
algorithms/theta_timestepping.templates.h
algorithms/timestep_control.h
```

## arborx (3)

```text
arborx/access_traits.h
arborx/bvh.h
arborx/distributed_tree.h
```

## base (149)

```text
base/aligned_vector.h
base/array_view.h
base/auto_derivative_function.h
base/bounding_box.h
base/bounding_box_data_out.h
base/communication_pattern_base.h
base/complex_overloads.h
base/conditional_ostream.h
base/convergence_table.h
base/data_out_base.h
base/derivative_form.h
base/discrete_time.h
base/enable_observer_pointer.h
base/event.h
base/exception_macros.h
base/exceptions.h
base/fast_transcendental.h
base/floating_point_comparator.h
base/flow_function.h
base/function.h
base/function.templates.h
base/function_bessel.h
base/function_cspline.h
base/function_derivative.h
base/function_lib.h
base/function_parser.h
base/function_restriction.h
base/function_signed_distance.h
base/function_spherical.h
base/function_time.h
base/function_time.templates.h
base/function_tools.h
base/geometric_utilities.h
base/geometry_info.h
base/graph_coloring.h
base/hdf5.h
base/incremental_function.h
base/index_set.h
base/init_finalize.h
base/iterator_range.h
base/job_identifier.h
base/kokkos.h
base/lazy.h
base/linear_index_iterator.h
base/logstream.h
base/memory_consumption.h
base/memory_space.h
base/memory_space_data.h
base/mg_level_object.h
base/mpi.h
base/mpi.templates.h
base/mpi_consensus_algorithms.h
base/mpi_consensus_algorithms.templates.h
base/mpi_large_count.h
base/mpi_noncontiguous_partitioner.h
base/mpi_noncontiguous_partitioner.templates.h
base/mpi_remote_point_evaluation.h
base/mpi_stub.h
base/mpi_tags.h
base/mu_parser_internal.h
base/multithread_info.h
base/mutable_bind.h
base/mutex.h
base/ndarray.h
base/numbers.h
base/observer_pointer.h
base/parallel.h
base/parameter_acceptor.h
base/parameter_handler.h
base/parsed_convergence_table.h
base/parsed_function.h
base/partitioner.h
base/partitioner.templates.h
base/patterns.h
base/point.h
base/polynomial.h
base/polynomial_space.h
base/polynomials_abf.h
base/polynomials_adini.h
base/polynomials_barycentric.h
base/polynomials_bdm.h
base/polynomials_bernardi_raugel.h
base/polynomials_bernstein.h
base/polynomials_hermite.h
base/polynomials_integrated_legendre_sz.h
base/polynomials_nedelec.h
base/polynomials_p.h
base/polynomials_piecewise.h
base/polynomials_pyramid.h
base/polynomials_rannacher_turek.h
base/polynomials_raviart_thomas.h
base/polynomials_rt_bubbles.h
base/polynomials_simplex.h
base/polynomials_vector_anisotropic.h
base/polynomials_wedge.h
base/process_grid.h
base/qprojector.h
base/quadrature.h
base/quadrature_lib.h
base/quadrature_point_data.h
base/quadrature_selector.h
base/sacado_product_type.h
base/scalar_polynomials_base.h
base/scalar_polynomials_vandermonde_base.h
base/scope_exit.h
base/segment_data_out.h
base/signaling_nan.h
base/smartpointer.h
base/std_cxx17/algorithm.h
base/std_cxx17/cmath.h
base/std_cxx17/optional.h
base/std_cxx17/tuple.h
base/std_cxx17/variant.h
base/std_cxx20/functional.h
base/std_cxx20/iota_view.h
base/std_cxx20/type_traits.h
base/std_cxx26/inplace_vector.h
base/subscriptor.h
base/symbolic_function.h
base/symbolic_function.templates.h
base/symmetric_tensor.h
base/symmetric_tensor.templates.h
base/synchronous_iterator.h
base/table.h
base/table_handler.h
base/table_indices.h
base/task_result.h
base/template_constraints.h
base/tensor.h
base/tensor_accessors.h
base/tensor_function.h
base/tensor_function.templates.h
base/tensor_function_parser.h
base/tensor_polynomials_base.h
base/tensor_product_polynomials.h
base/tensor_product_polynomials_bubbles.h
base/tensor_product_polynomials_const.h
base/thread_local_storage.h
base/thread_management.h
base/time_stepping.h
base/time_stepping.templates.h
base/timer.h
base/trilinos_utilities.h
base/types.h
base/undefine_macros.h
base/utilities.h
base/vector_slice.h
base/vectorization.h
base/work_stream.h
```

## boost_adaptors (3)

```text
boost_adaptors/bounding_box.h
boost_adaptors/point.h
boost_adaptors/segment.h
```

## cgal (7)

```text
cgal/additional_data.h
cgal/intersections.h
cgal/point_conversion.h
cgal/polygon.h
cgal/surface_mesh.h
cgal/triangulation.h
cgal/utilities.h
```

## differentiation (22)

```text
differentiation/ad/ad_drivers.h
differentiation/ad/ad_helpers.h
differentiation/ad/ad_number_traits.h
differentiation/ad/ad_number_types.h
differentiation/ad/adolc_math.h
differentiation/ad/adolc_number_types.h
differentiation/ad/adolc_product_types.h
differentiation/ad/sacado_math.h
differentiation/ad/sacado_number_types.h
differentiation/ad/sacado_product_types.h
differentiation/ad.h
differentiation/sd/symengine_math.h
differentiation/sd/symengine_number_traits.h
differentiation/sd/symengine_number_types.h
differentiation/sd/symengine_number_visitor_internal.h
differentiation/sd/symengine_optimizer.h
differentiation/sd/symengine_product_types.h
differentiation/sd/symengine_scalar_operations.h
differentiation/sd/symengine_tensor_operations.h
differentiation/sd/symengine_types.h
differentiation/sd/symengine_utilities.h
differentiation/sd.h
```

## distributed (11)

```text
distributed/cell_data_transfer.h
distributed/cell_data_transfer.templates.h
distributed/cell_weights.h
distributed/field_transfer.h
distributed/fully_distributed_tria.h
distributed/grid_refinement.h
distributed/p4est_wrappers.h
distributed/repartitioning_policy_tools.h
distributed/shared_tria.h
distributed/tria.h
distributed/tria_base.h
```

## dofs (12)

```text
dofs/block_info.h
dofs/dof_accessor.h
dofs/dof_accessor.templates.h
dofs/dof_faces.h
dofs/dof_handler.h
dofs/dof_handler_policy.h
dofs/dof_iterator_selector.h
dofs/dof_levels.h
dofs/dof_objects.h
dofs/dof_renumbering.h
dofs/dof_tools.h
dofs/number_cache.h
```

## fe (72)

```text
fe/block_mask.h
fe/component_mask.h
fe/fe.h
fe/fe_abf.h
fe/fe_base.h
fe/fe_bdm.h
fe/fe_bernardi_raugel.h
fe/fe_bernstein.h
fe/fe_coupling_values.h
fe/fe_data.h
fe/fe_dg_vector.h
fe/fe_dg_vector.templates.h
fe/fe_dgp.h
fe/fe_dgp_monomial.h
fe/fe_dgp_nonparametric.h
fe/fe_dgq.h
fe/fe_enriched.h
fe/fe_face.h
fe/fe_hermite.h
fe/fe_interface_values.h
fe/fe_nedelec.h
fe/fe_nedelec.templates.h
fe/fe_nedelec_sz.h
fe/fe_nothing.h
fe/fe_p1nc.h
fe/fe_poly.h
fe/fe_poly_face.h
fe/fe_poly_face.templates.h
fe/fe_poly_tensor.h
fe/fe_pyramid_p.h
fe/fe_q.h
fe/fe_q_base.h
fe/fe_q_bubbles.h
fe/fe_q_dg0.h
fe/fe_q_hierarchical.h
fe/fe_q_iso_q1.h
fe/fe_rannacher_turek.h
fe/fe_raviart_thomas.h
fe/fe_rt_bubbles.h
fe/fe_series.h
fe/fe_simplex_p.h
fe/fe_simplex_p_bubbles.h
fe/fe_system.h
fe/fe_tools.h
fe/fe_tools.templates.h
fe/fe_tools_extrapolate.templates.h
fe/fe_tools_interpolate.templates.h
fe/fe_trace.h
fe/fe_update_flags.h
fe/fe_values.h
fe/fe_values_base.h
fe/fe_values_extractors.h
fe/fe_values_views.h
fe/fe_values_views_internal.h
fe/fe_wedge_p.h
fe/mapping.h
fe/mapping_c1.h
fe/mapping_cartesian.h
fe/mapping_fe.h
fe/mapping_fe_field.h
fe/mapping_fe_field.templates.h
fe/mapping_internal.h
fe/mapping_manifold.h
fe/mapping_p1.h
fe/mapping_q.h
fe/mapping_q1.h
fe/mapping_q1_eulerian.h
fe/mapping_q_cache.h
fe/mapping_q_eulerian.h
fe/mapping_q_generic.h
fe/mapping_q_internal.h
fe/mapping_related_data.h
```

## gmsh (1)

```text
gmsh/utilities.h
```

## grid (35)

```text
grid/cell_data.h
grid/cell_id.h
grid/cell_id_translator.h
grid/cell_status.h
grid/composition_manifold.h
grid/filtered_iterator.h
grid/grid_generator.h
grid/grid_in.h
grid/grid_out.h
grid/grid_refinement.h
grid/grid_tools.h
grid/grid_tools_cache.h
grid/grid_tools_cache_update_flags.h
grid/grid_tools_geometry.h
grid/grid_tools_topology.h
grid/intergrid_map.h
grid/magic_numbers.h
grid/manifold.h
grid/manifold_lib.h
grid/persistent_tria.h
grid/reference_cell.h
grid/tensor_product_manifold.h
grid/tria.h
grid/tria_accessor.h
grid/tria_accessor.templates.h
grid/tria_description.h
grid/tria_faces.h
grid/tria_iterator.h
grid/tria_iterator.templates.h
grid/tria_iterator_base.h
grid/tria_iterator_selector.h
grid/tria_levels.h
grid/tria_objects.h
grid/tria_objects_orientations.h
grid/tria_orientation.h
```

## hp (6)

```text
hp/collection.h
hp/fe_collection.h
hp/fe_values.h
hp/mapping_collection.h
hp/q_collection.h
hp/refinement.h
```

## integrators (9)

```text
integrators/advection.h
integrators/divergence.h
integrators/elasticity.h
integrators/grad_div.h
integrators/l2.h
integrators/laplace.h
integrators/local_integrators.h
integrators/maxwell.h
integrators/patches.h
```

## lac (159)

```text
lac/affine_constraints.h
lac/affine_constraints.templates.h
lac/arpack_solver.h
lac/blas_extension_templates.h
lac/block_indices.h
lac/block_linear_operator.h
lac/block_matrix_base.h
lac/block_sparse_matrix.h
lac/block_sparse_matrix.templates.h
lac/block_sparse_matrix_ez.h
lac/block_sparse_matrix_ez.templates.h
lac/block_sparsity_pattern.h
lac/block_vector.h
lac/block_vector.templates.h
lac/block_vector_base.h
lac/chunk_sparse_matrix.h
lac/chunk_sparse_matrix.templates.h
lac/chunk_sparsity_pattern.h
lac/constrained_linear_operator.h
lac/constraint_matrix.h
lac/diagonal_matrix.h
lac/dynamic_sparsity_pattern.h
lac/eigen.h
lac/exceptions.h
lac/full_matrix.h
lac/full_matrix.templates.h
lac/generic_linear_algebra.h
lac/ginkgo_solver.h
lac/householder.h
lac/identity_matrix.h
lac/la_parallel_block_vector.h
lac/la_parallel_block_vector.templates.h
lac/la_parallel_vector.h
lac/la_parallel_vector.templates.h
lac/lapack_full_matrix.h
lac/lapack_support.h
lac/lapack_templates.h
lac/linear_operator.h
lac/linear_operator_tools.h
lac/matrix_block.h
lac/matrix_iterator.h
lac/matrix_out.h
lac/matrix_scaling.h
lac/orthogonalization.h
lac/packaged_operation.h
lac/parpack_solver.h
lac/petsc_block_sparse_matrix.h
lac/petsc_block_vector.h
lac/petsc_communication_pattern.h
lac/petsc_compatibility.h
lac/petsc_full_matrix.h
lac/petsc_matrix_base.h
lac/petsc_matrix_free.h
lac/petsc_parallel_block_sparse_matrix.h
lac/petsc_parallel_block_vector.h
lac/petsc_parallel_sparse_matrix.h
lac/petsc_parallel_vector.h
lac/petsc_precondition.h
lac/petsc_snes.h
lac/petsc_snes.templates.h
lac/petsc_solver.h
lac/petsc_sparse_matrix.h
lac/petsc_ts.h
lac/petsc_ts.templates.h
lac/petsc_vector.h
lac/petsc_vector_base.h
lac/precondition.h
lac/precondition_block.h
lac/precondition_block.templates.h
lac/precondition_block_base.h
lac/precondition_selector.h
lac/psblas_common.h
lac/psblas_precondition.h
lac/psblas_sparse_matrix.h
lac/psblas_sparsity_pattern.h
lac/psblas_vector.h
lac/qr.h
lac/read_vector.h
lac/read_write_vector.h
lac/read_write_vector.templates.h
lac/relaxation_block.h
lac/relaxation_block.templates.h
lac/scalapack.h
lac/scalapack.templates.h
lac/schur_complement.h
lac/slepc_solver.h
lac/slepc_spectral_transformation.h
lac/solver.h
lac/solver_bicgstab.h
lac/solver_cg.h
lac/solver_control.h
lac/solver_fire.h
lac/solver_gmres.h
lac/solver_idr.h
lac/solver_minres.h
lac/solver_qmrs.h
lac/solver_relaxation.h
lac/solver_richardson.h
lac/solver_selector.h
lac/sparse_decomposition.h
lac/sparse_decomposition.templates.h
lac/sparse_direct.h
lac/sparse_ilu.h
lac/sparse_ilu.templates.h
lac/sparse_matrix.h
lac/sparse_matrix.templates.h
lac/sparse_matrix_ez.h
lac/sparse_matrix_ez.templates.h
lac/sparse_matrix_tools.h
lac/sparse_mic.h
lac/sparse_mic.templates.h
lac/sparse_vanka.h
lac/sparse_vanka.templates.h
lac/sparsity_pattern.h
lac/sparsity_pattern_base.h
lac/sparsity_tools.h
lac/tensor_product_matrix.h
lac/tensor_product_matrix.templates.h
lac/tridiagonal_matrix.h
lac/trilinos_block_sparse_matrix.h
lac/trilinos_epetra_communication_pattern.h
lac/trilinos_epetra_vector.h
lac/trilinos_index_access.h
lac/trilinos_linear_operator.h
lac/trilinos_parallel_block_vector.h
lac/trilinos_precondition.h
lac/trilinos_solver.h
lac/trilinos_sparse_matrix.h
lac/trilinos_sparsity_pattern.h
lac/trilinos_tpetra_block_sparse_matrix.h
lac/trilinos_tpetra_block_sparse_matrix.templates.h
lac/trilinos_tpetra_block_vector.h
lac/trilinos_tpetra_block_vector.templates.h
lac/trilinos_tpetra_communication_pattern.h
lac/trilinos_tpetra_precondition.h
lac/trilinos_tpetra_precondition.templates.h
lac/trilinos_tpetra_precondition_frosch.templates.h
lac/trilinos_tpetra_precondition_muelu.templates.h
lac/trilinos_tpetra_solver_direct.h
lac/trilinos_tpetra_solver_direct.templates.h
lac/trilinos_tpetra_sparse_matrix.h
lac/trilinos_tpetra_sparse_matrix.templates.h
lac/trilinos_tpetra_sparsity_pattern.h
lac/trilinos_tpetra_to_trilinos_wrappers.h
lac/trilinos_tpetra_types.h
lac/trilinos_tpetra_vector.h
lac/trilinos_tpetra_vector.templates.h
lac/trilinos_vector.h
lac/trilinos_xpetra_types.h
lac/utilities.h
lac/vector.h
lac/vector.templates.h
lac/vector_element_access.h
lac/vector_memory.h
lac/vector_memory.templates.h
lac/vector_operation.h
lac/vector_operations_internal.h
lac/vector_space_vector.h
lac/vector_type_traits.h
```

## matrix_free (45)

```text
matrix_free/constraint_info.h
matrix_free/dof_info.h
matrix_free/dof_info.templates.h
matrix_free/evaluation_flags.h
matrix_free/evaluation_kernels.h
matrix_free/evaluation_kernels_common.h
matrix_free/evaluation_kernels_face.h
matrix_free/evaluation_kernels_hanging_nodes.h
matrix_free/evaluation_selector.h
matrix_free/evaluation_template_face_factory.templates.h
matrix_free/evaluation_template_factory.h
matrix_free/evaluation_template_factory.templates.h
matrix_free/evaluation_template_factory_hanging_nodes.templates.h
matrix_free/evaluation_template_factory_internal.h
matrix_free/face_info.h
matrix_free/face_setup_internal.h
matrix_free/fe_evaluation.h
matrix_free/fe_evaluation_data.h
matrix_free/fe_point_evaluation.h
matrix_free/fe_remote_evaluation.h
matrix_free/hanging_nodes_internal.h
matrix_free/mapping_data_on_the_fly.h
matrix_free/mapping_info.h
matrix_free/mapping_info.templates.h
matrix_free/mapping_info_storage.h
matrix_free/mapping_info_storage.templates.h
matrix_free/matrix_free.h
matrix_free/matrix_free.templates.h
matrix_free/operators.h
matrix_free/portable_evaluation_kernels.h
matrix_free/portable_fe_evaluation.h
matrix_free/portable_hanging_nodes_internal.h
matrix_free/portable_matrix_free.h
matrix_free/portable_matrix_free.templates.h
matrix_free/portable_tensor_product_kernels.h
matrix_free/shape_info.h
matrix_free/shape_info.templates.h
matrix_free/task_info.h
matrix_free/tensor_product_kernels.h
matrix_free/tensor_product_point_kernels.h
matrix_free/tools.h
matrix_free/type_traits.h
matrix_free/util.h
matrix_free/vector_access_internal.h
matrix_free/vector_data_exchange.h
```

## meshworker (17)

```text
meshworker/assemble_flags.h
meshworker/assembler.h
meshworker/copy_data.h
meshworker/dof_info.h
meshworker/dof_info.templates.h
meshworker/functional.h
meshworker/integration_info.h
meshworker/integration_info.templates.h
meshworker/local_integrator.h
meshworker/local_results.h
meshworker/loop.h
meshworker/mesh_loop.h
meshworker/output.h
meshworker/scratch_data.h
meshworker/simple.h
meshworker/vector_selector.h
meshworker/vector_selector.templates.h
```

## multigrid (24)

```text
multigrid/mg_base.h
multigrid/mg_base.templates.h
multigrid/mg_block_smoother.h
multigrid/mg_coarse.h
multigrid/mg_constrained_dofs.h
multigrid/mg_matrix.h
multigrid/mg_smoother.h
multigrid/mg_tools.h
multigrid/mg_transfer.h
multigrid/mg_transfer.templates.h
multigrid/mg_transfer_block.h
multigrid/mg_transfer_block.templates.h
multigrid/mg_transfer_component.h
multigrid/mg_transfer_component.templates.h
multigrid/mg_transfer_global_coarsening.h
multigrid/mg_transfer_global_coarsening.templates.h
multigrid/mg_transfer_internal.h
multigrid/mg_transfer_matrix_free.h
multigrid/mg_transfer_matrix_free.templates.h
multigrid/multigrid.h
multigrid/multigrid.templates.h
multigrid/portable_mg_transfer_global_coarsening.h
multigrid/portable_mg_transfer_global_coarsening.templates.h
multigrid/sparse_matrix_collection.h
```

## non_matching (8)

```text
non_matching/closest_surface_point.h
non_matching/coupling.h
non_matching/fe_immersed_values.h
non_matching/fe_values.h
non_matching/immersed_surface_quadrature.h
non_matching/mapping_info.h
non_matching/mesh_classifier.h
non_matching/quadrature_generator.h
```

## numerics (56)

```text
numerics/adaptation_strategies.h
numerics/cell_data_transfer.h
numerics/cell_data_transfer.templates.h
numerics/data_component_interpretation.h
numerics/data_out.h
numerics/data_out_dof_data.h
numerics/data_out_dof_data.templates.h
numerics/data_out_faces.h
numerics/data_out_points.h
numerics/data_out_resample.h
numerics/data_out_rotation.h
numerics/data_out_stack.h
numerics/data_postprocessor.h
numerics/derivative_approximation.h
numerics/dof_output_operator.h
numerics/dof_output_operator.templates.h
numerics/dof_print_solver_step.h
numerics/error_estimator.h
numerics/error_estimator.templates.h
numerics/fe_field_function.h
numerics/fe_field_function.templates.h
numerics/histogram.h
numerics/history.h
numerics/matrix_creator.h
numerics/matrix_creator.templates.h
numerics/matrix_tools.h
numerics/nonlinear.h
numerics/point_value_history.h
numerics/rtree.h
numerics/smoothness_estimator.h
numerics/solution_transfer.h
numerics/solution_transfer.templates.h
numerics/tensor_product_matrix_creator.h
numerics/time_dependent.h
numerics/vector_tools.h
numerics/vector_tools.templates.h
numerics/vector_tools_boundary.h
numerics/vector_tools_boundary.templates.h
numerics/vector_tools_common.h
numerics/vector_tools_constraints.h
numerics/vector_tools_constraints.templates.h
numerics/vector_tools_evaluate.h
numerics/vector_tools_integrate_difference.h
numerics/vector_tools_integrate_difference.templates.h
numerics/vector_tools_interpolate.h
numerics/vector_tools_interpolate.templates.h
numerics/vector_tools_mean_value.h
numerics/vector_tools_mean_value.templates.h
numerics/vector_tools_point_gradient.h
numerics/vector_tools_point_gradient.templates.h
numerics/vector_tools_point_value.h
numerics/vector_tools_point_value.templates.h
numerics/vector_tools_project.h
numerics/vector_tools_project.templates.h
numerics/vector_tools_rhs.h
numerics/vector_tools_rhs.templates.h
```

## opencascade (3)

```text
opencascade/boundary_lib.h
opencascade/manifold_lib.h
opencascade/utilities.h
```

## optimization (3)

```text
optimization/line_minimization.h
optimization/rol/vector_adaptor.h
optimization/solver_bfgs.h
```

## particles (9)

```text
particles/data_out.h
particles/generators.h
particles/particle.h
particles/particle_accessor.h
particles/particle_handler.h
particles/particle_iterator.h
particles/partitioner.h
particles/property_pool.h
particles/utilities.h
```

## physics (5)

```text
physics/elasticity/kinematics.h
physics/elasticity/standard_tensors.h
physics/notation.h
physics/transformations.h
physics/vector_relations.h
```

## sundials (16)

```text
sundials/arkode.h
sundials/arkode.templates.h
sundials/arkode_exception.h
sundials/arkode_stepper.h
sundials/arkode_stepper.templates.h
sundials/ida.h
sundials/invocation_context.h
sundials/kinsol.h
sundials/n_vector.h
sundials/n_vector.templates.h
sundials/sundials_handle.h
sundials/sundials_types.h
sundials/sunlinsol_wrapper.h
sundials/sunlinsol_wrapper.templates.h
sundials/sunnonlinsol_wrapper.h
sundials/utilities.h
```

## trilinos (3)

```text
trilinos/nox.h
trilinos/nox.templates.h
trilinos/rol_adaptor.h
```

## vtk (1)

```text
vtk/utilities.h
```

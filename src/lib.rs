//! Small, safe finite-element building blocks inspired by deal.II.
//!
//! The current scope is serial, scalar P1 triangles and tetrahedra and Q1 quadrilaterals.
//! It is not yet a replacement for deal.II.
#![forbid(unsafe_code)]

pub mod base;
pub mod constraints;
pub mod dofs;
pub mod elliptic;
pub mod elliptic_quad;
pub mod elliptic_tetra;
pub mod estimator;
pub mod fe;
pub mod field;
pub mod gmsh;
pub mod grid;
pub mod integrators;
pub mod lac;
mod local_refinement;
pub mod neumann;
pub mod neumann3d;
mod neumann_core;
pub mod neumann_quad;
mod nonconforming_refinement;
pub mod numerics;
pub mod quad;
pub mod quad_field;
mod quad_refinement;
pub mod quadrature;
mod refinement_history;
pub mod tetra;
pub mod tetra_field;
mod tetra_refinement;
pub mod vtk;

/// An invalid input or numerical failure.
#[derive(Debug, Clone, PartialEq)]
pub struct Error(pub String);

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl std::error::Error for Error {}

/// Result type used by the numerical components.
pub type Result<T> = std::result::Result<T, Error>;

pub(crate) fn error(message: &str) -> Error {
    Error(message.to_string())
}

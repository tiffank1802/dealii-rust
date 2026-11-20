//! Small, safe finite-element building blocks inspired by deal.II.
//!
//! The current scope is serial, scalar P1 triangles and tetrahedra and Q1 quadrilaterals.
//! It is not yet a replacement for deal.II.
#![forbid(unsafe_code)]

pub mod base;
pub mod dofs;
pub mod elliptic;
pub mod estimator;
pub mod fe;
pub mod field;
pub mod gmsh;
pub mod grid;
pub mod integrators;
pub mod lac;
pub mod neumann;
pub mod neumann3d;
mod neumann_core;
pub mod neumann_quad;
pub mod numerics;
pub mod quad;
pub mod quadrature;
pub mod tetra;
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

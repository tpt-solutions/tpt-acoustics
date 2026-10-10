//! Error type shared across the `tpt-acoustics` crates.

use core::fmt;

/// Errors produced by acoustic computations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AcousticError {
    /// A dimension, frequency, or other quantity that must be strictly
    /// positive was zero or negative.
    NonPositive(&'static str),
    /// A linear system solve encountered a singular (or numerically singular)
    /// matrix.
    Singular(&'static str),
    /// A geometric input violates an invariant (e.g. non-convex room, degenerate
    /// triangle, inverted normal).
    InvalidGeometry(&'static str),
    /// A generic domain error with a static description.
    InvalidArgument(&'static str),
}

impl fmt::Display for AcousticError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AcousticError::NonPositive(what) => {
                write!(f, "quantity must be strictly positive: {what}")
            }
            AcousticError::Singular(what) => write!(f, "singular linear system: {what}"),
            AcousticError::InvalidGeometry(what) => write!(f, "invalid geometry: {what}"),
            AcousticError::InvalidArgument(what) => write!(f, "invalid argument: {what}"),
        }
    }
}

impl core::error::Error for AcousticError {}

//! STL import (ASCII and binary) for ray-tracing scenes.
//!
//! The file format is detected from its size. A binary STL is exactly
//! `84 + 50·n` bytes, where `n` is the triangle count in bytes 80–84. Any
//! other input is read as ASCII, which is also the fallback for files whose
//! binary header happens to begin with `solid`. Facet normals are ignored,
//! because scene intersection recomputes them from the vertices.
//!
//! The reader is `no_std` + `alloc`. It takes the bytes, so callers decide
//! how to load the file.

use alloc::vec::Vec;

use tpt_acou_core::error::AcousticError;

use crate::material::BoundaryMaterial;
use crate::scene::{Scene, SceneTriangle};
use crate::vec3::Vec3;

const BINARY_HEADER: usize = 80;
const BINARY_RECORD: usize = 50;

fn invalid(message: &'static str) -> AcousticError {
    AcousticError::InvalidArgument(message)
}

fn le_f32(bytes: &[u8], at: usize) -> f32 {
    f32::from_le_bytes([bytes[at], bytes[at + 1], bytes[at + 2], bytes[at + 3]])
}

/// Parses an STL file (ASCII or binary) into triangles.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] for truncated binary data,
/// malformed or non-finite vertex coordinates, an ASCII facet without three
/// vertices, or a file with no triangles.
pub fn parse_stl(bytes: &[u8]) -> Result<Vec<[Vec3; 3]>, AcousticError> {
    if let Some(count) = binary_triangle_count(bytes) {
        return parse_binary(bytes, count);
    }
    let text = core::str::from_utf8(bytes)
        .map_err(|_| invalid("STL file is neither binary nor UTF-8 ASCII"))?;
    parse_ascii(text)
}

/// Triangle count if `bytes` has exactly the size of a binary STL.
fn binary_triangle_count(bytes: &[u8]) -> Option<usize> {
    if bytes.len() < BINARY_HEADER + 4 {
        return None;
    }
    let count = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize;
    let expected = BINARY_HEADER + 4 + count.checked_mul(BINARY_RECORD)?;
    (bytes.len() == expected).then_some(count)
}

fn parse_binary(bytes: &[u8], count: usize) -> Result<Vec<[Vec3; 3]>, AcousticError> {
    if count == 0 {
        return Err(invalid("binary STL contains no triangles"));
    }
    let mut triangles = Vec::with_capacity(count);
    for i in 0..count {
        // Each record: normal (12 bytes), three vertices (36 bytes), attribute (2).
        let base = BINARY_HEADER + 4 + i * BINARY_RECORD + 12;
        let mut v = [Vec3::new(0.0, 0.0, 0.0); 3];
        for (k, vertex) in v.iter_mut().enumerate() {
            let at = base + 12 * k;
            let coords = [
                f64::from(le_f32(bytes, at)),
                f64::from(le_f32(bytes, at + 4)),
                f64::from(le_f32(bytes, at + 8)),
            ];
            if coords.iter().any(|c| !c.is_finite()) {
                return Err(invalid("binary STL vertex is not finite"));
            }
            *vertex = Vec3::new(coords[0], coords[1], coords[2]);
        }
        triangles.push(v);
    }
    Ok(triangles)
}

fn parse_ascii(text: &str) -> Result<Vec<[Vec3; 3]>, AcousticError> {
    let mut triangles = Vec::new();
    let mut pending: Vec<Vec3> = Vec::with_capacity(3);
    for line in text.lines() {
        let mut parts = line.split_whitespace();
        if parts.next() != Some("vertex") {
            continue;
        }
        let mut coord = || -> Result<f64, AcousticError> {
            let value: f64 = parts
                .next()
                .ok_or(invalid("STL vertex needs three coordinates"))?
                .parse()
                .map_err(|_| invalid("STL vertex coordinate is not a number"))?;
            if value.is_finite() {
                Ok(value)
            } else {
                Err(invalid("STL vertex coordinate is not finite"))
            }
        };
        let (x, y, z) = (coord()?, coord()?, coord()?);
        pending.push(Vec3::new(x, y, z));
        if pending.len() == 3 {
            triangles.push([pending[0], pending[1], pending[2]]);
            pending.clear();
        }
    }
    if !pending.is_empty() {
        return Err(invalid("STL facet has fewer than three vertices"));
    }
    if triangles.is_empty() {
        return Err(invalid("STL file contains no triangles"));
    }
    Ok(triangles)
}

/// Builds a ray-tracing scene from STL bytes, giving every facet the same
/// wall material.
///
/// # Errors
/// Propagates [`parse_stl`] errors and the scene constructor's errors.
pub fn scene_from_stl(bytes: &[u8], material: BoundaryMaterial) -> Result<Scene, AcousticError> {
    let triangles: Vec<SceneTriangle> = parse_stl(bytes)?
        .into_iter()
        .map(|v| SceneTriangle { v, material: 0 })
        .collect();
    Scene::new(&triangles, alloc::vec![material])
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;
    use alloc::vec::Vec;

    /// Triangles of a unit tetrahedron.
    fn tetra() -> Vec<[[f32; 3]; 3]> {
        let (a, b, c, d) = (
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        );
        vec![[a, c, b], [a, b, d], [a, d, c], [b, c, d]]
    }

    fn binary(triangles: &[[[f32; 3]; 3]]) -> Vec<u8> {
        let mut out = vec![0_u8; 80];
        out.extend_from_slice(&(triangles.len() as u32).to_le_bytes());
        for t in triangles {
            out.extend_from_slice(&[0_u8; 12]); // normal, ignored by the reader
            for v in t {
                for c in v {
                    out.extend_from_slice(&c.to_le_bytes());
                }
            }
            out.extend_from_slice(&[0_u8; 2]);
        }
        out
    }

    fn ascii(triangles: &[[[f32; 3]; 3]]) -> alloc::string::String {
        let mut s = alloc::string::String::from("solid tetra\n");
        for t in triangles {
            s.push_str("  facet normal 0 0 0\n    outer loop\n");
            for v in t {
                s.push_str(&alloc::format!("      vertex {} {} {}\n", v[0], v[1], v[2]));
            }
            s.push_str("    endloop\n  endfacet\n");
        }
        s.push_str("endsolid tetra\n");
        s
    }

    #[test]
    fn binary_and_ascii_agree() {
        let t = tetra();
        let from_binary = parse_stl(&binary(&t)).unwrap();
        let from_ascii = parse_stl(ascii(&t).as_bytes()).unwrap();
        assert_eq!(from_binary.len(), 4);
        assert_eq!(from_binary, from_ascii);
    }

    #[test]
    fn binary_with_solid_header_is_still_binary() {
        // A binary header may start with "solid" (common in CAD exports).
        let mut bytes = binary(&tetra());
        bytes[..6].copy_from_slice(b"solid ");
        assert_eq!(parse_stl(&bytes).unwrap().len(), 4);
    }

    #[test]
    fn truncated_binary_is_rejected() {
        let mut bytes = binary(&tetra());
        bytes.truncate(bytes.len() - 10);
        // Falls through to ASCII parsing, which finds no vertex lines.
        assert!(parse_stl(&bytes).is_err());
    }

    #[test]
    fn bad_ascii_is_rejected() {
        assert!(parse_stl(b"solid x\nvertex 1 2\nendsolid\n").is_err());
        assert!(
            parse_stl(b"solid x\nvertex 0 0 nan\nvertex 0 0 0\nvertex 1 0 0\nendsolid\n").is_err()
        );
        assert!(parse_stl(b"solid x\nvertex 0 0 0\nvertex 1 0 0\nendsolid\n").is_err());
        assert!(parse_stl(b"solid empty\nendsolid\n").is_err());
    }

    #[test]
    fn tetrahedron_scene_is_closed() {
        let scene =
            scene_from_stl(&binary(&tetra()), BoundaryMaterial::concrete_painted()).unwrap();
        assert_eq!(scene.len(), 4);
        // From inside the tetrahedron, a ray along +x leaves through the x-y-z
        // slanted face at x = 1 − y − z, which is 1/3 away from (0.2, 0.2, 0.2).
        let origin = Vec3::new(0.2, 0.2, 0.2);
        let hit = scene
            .hit(origin, Vec3::new(1.0, 0.0, 0.0), f64::INFINITY)
            .expect("closed mesh must be hit");
        assert!((hit.t - 0.4).abs() < 1e-6, "t = {}", hit.t);
    }
}

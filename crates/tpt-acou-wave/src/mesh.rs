//! Triangulated boundary meshes for the BEM.
//!
//! Element orientation convention: vertices are ordered counter-clockwise
//! when viewed from *outside* the interior domain, so the geometric normal
//! points out of the interior (matching the BEM sign conventions in
//! [`crate::bem`]).

use alloc::vec::Vec;

use tpt_acou_core::error::AcousticError;

use tpt_acou_core::Vec3;

/// One constant boundary element: vertex positions, centroid, area and
/// outward unit normal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElementView {
    /// Centroid.
    pub centroid: Vec3,
    /// Outward unit normal.
    pub normal: Vec3,
    /// Area (m²).
    pub area: f64,
}

impl ElementView {
    /// Builds the element view of triangle `(a, b, c)` (CCW from outside).
    #[must_use]
    pub fn from_triangle(a: Vec3, b: Vec3, c: Vec3) -> Self {
        let cross = (b - a).cross(c - a);
        let area2 = cross.norm();
        let normal = if area2 > 0.0 {
            cross / area2
        } else {
            Vec3::ZERO
        };
        Self {
            centroid: (a + b + c) / 3.0,
            normal,
            area: 0.5 * area2,
        }
    }
}

/// A triangulated boundary surface (constant-element BEM mesh).
#[derive(Debug, Clone, Default)]
pub struct BoundaryMesh {
    nodes: Vec<Vec3>,
    tris: Vec<[usize; 3]>,
    elements: Vec<ElementView>,
}

impl BoundaryMesh {
    /// Builds a mesh from nodes and triangle index triples.
    ///
    /// # Errors
    /// Returns [`AcousticError::InvalidGeometry`] if any triangle references
    /// an out-of-range node index or is degenerate (zero area).
    pub fn new(nodes: Vec<Vec3>, tris: Vec<[usize; 3]>) -> Result<Self, AcousticError> {
        let mut elements = Vec::with_capacity(tris.len());
        for t in &tris {
            for &i in t {
                if i >= nodes.len() {
                    return Err(AcousticError::InvalidGeometry(
                        "triangle node index out of range",
                    ));
                }
            }
            let e = ElementView::from_triangle(nodes[t[0]], nodes[t[1]], nodes[t[2]]);
            if e.area <= 0.0 || !e.area.is_finite() {
                return Err(AcousticError::InvalidGeometry("degenerate triangle"));
            }
            elements.push(e);
        }
        Ok(Self {
            nodes,
            tris,
            elements,
        })
    }

    /// Number of elements (== BEM system size).
    #[must_use]
    pub fn len(&self) -> usize {
        self.elements.len()
    }

    /// True if the mesh has no elements.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    /// The nodes.
    #[must_use]
    pub fn nodes(&self) -> &[Vec3] {
        &self.nodes
    }

    /// The triangle index triples.
    #[must_use]
    pub fn triangles(&self) -> &[[usize; 3]] {
        &self.tris
    }

    /// The precomputed element views.
    #[must_use]
    pub fn elements(&self) -> &[ElementView] {
        &self.elements
    }

    /// A closed UV-sphere of radius `radius` with `n_lat` latitude bands and
    /// `n_lon` longitude segments, outward-oriented.
    ///
    /// # Errors
    /// Returns [`AcousticError::NonPositive`] for `radius <= 0` or fewer
    /// than 2 latitude bands / 3 longitude segments.
    pub fn sphere(radius: f64, n_lat: usize, n_lon: usize) -> Result<Self, AcousticError> {
        if radius <= 0.0 {
            return Err(AcousticError::NonPositive("sphere radius"));
        }
        if n_lat < 2 || n_lon < 3 {
            return Err(AcousticError::NonPositive("sphere resolution"));
        }
        let mut nodes = Vec::with_capacity(2 + n_lat * n_lon);
        // Poles.
        nodes.push(Vec3::new(0.0, 0.0, radius));
        for i in 1..n_lat {
            let theta = core::f64::consts::PI * i as f64 / n_lat as f64;
            for j in 0..n_lon {
                let phi = core::f64::consts::TAU * j as f64 / n_lon as f64;
                nodes.push(Vec3::new(
                    radius * libm::sin(theta) * libm::cos(phi),
                    radius * libm::sin(theta) * libm::sin(phi),
                    radius * libm::cos(theta),
                ));
            }
        }
        nodes.push(Vec3::new(0.0, 0.0, -radius));
        let north = 0_usize;
        let south = nodes.len() - 1;
        let ring = |i: usize, j: usize| 1 + i * n_lon + (j % n_lon);

        let mut tris = Vec::with_capacity(2 * n_lat * n_lon);
        // North cap (viewed from outside, CCW ordering keeps normals outward).
        for j in 0..n_lon {
            tris.push([north, ring(0, j), ring(0, j + 1)]);
        }
        // Quad strip.
        for i in 0..n_lat - 2 {
            for j in 0..n_lon {
                let a = ring(i, j);
                let b = ring(i, j + 1);
                let c = ring(i + 1, j + 1);
                let d = ring(i + 1, j);
                tris.push([a, c, b]);
                tris.push([a, d, c]);
            }
        }
        // South cap.
        for j in 0..n_lon {
            tris.push([south, ring(n_lat - 2, j + 1), ring(n_lat - 2, j)]);
        }
        Self::new(nodes, tris)
    }

    /// A flat rectangular patch in the `z = 0` plane, `lx × ly` with `nx × ny`
    /// cells (two triangles each), normal `+z`.
    ///
    /// # Errors
    /// Returns [`AcousticError::NonPositive`] for non-positive extents or
    /// zero resolution.
    pub fn rectangle_patch(lx: f64, ly: f64, nx: usize, ny: usize) -> Result<Self, AcousticError> {
        if lx <= 0.0 || ly <= 0.0 {
            return Err(AcousticError::NonPositive("patch extent"));
        }
        if nx == 0 || ny == 0 {
            return Err(AcousticError::NonPositive("patch resolution"));
        }
        let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1));
        for j in 0..=ny {
            for i in 0..=nx {
                nodes.push(Vec3::new(
                    lx * i as f64 / nx as f64 - lx / 2.0,
                    ly * j as f64 / ny as f64 - ly / 2.0,
                    0.0,
                ));
            }
        }
        let idx = |i: usize, j: usize| j * (nx + 1) + i;
        let mut tris = Vec::with_capacity(2 * nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                let (a, b, c, d) = (idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1));
                tris.push([a, b, c]);
                tris.push([a, c, d]);
            }
        }
        Self::new(nodes, tris)
    }

    /// The closed surface of a rectangular box with full extents
    /// `lx × ly × lz` (spanning `±lx/2` etc.), centered at the origin,
    /// `n` cells per edge (per face), outward normals.
    ///
    /// # Errors
    /// See [`BoundaryMesh::rectangle_patch`].
    pub fn box_surface(lx: f64, ly: f64, lz: f64, n: usize) -> Result<Self, AcousticError> {
        if n == 0 {
            return Err(AcousticError::NonPositive("box resolution"));
        }
        let mut nodes = Vec::new();
        let mut tris = Vec::new();
        let push_grid = |nodes: &mut Vec<Vec3>, tris: &mut Vec<[usize; 3]>, corners: [Vec3; 3]| {
            // Grid of `n × n` quads on the face spanned by du = corners[1]-[0],
            // dv = corners[2]-[0], offset from corners[0]; orientation chosen
            // so (du × dv) points outward.
            let base = nodes.len();
            for j in 0..=n {
                for i in 0..=n {
                    let u = i as f64 / n as f64;
                    let v = j as f64 / n as f64;
                    nodes.push(
                        corners[0] + (corners[1] - corners[0]) * u + (corners[2] - corners[0]) * v,
                    );
                }
            }
            let idx = |i: usize, j: usize| base + j * (n + 1) + i;
            for j in 0..n {
                for i in 0..n {
                    let (a, b, c, d) = (idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1));
                    tris.push([a, c, b]);
                    tris.push([a, d, c]);
                }
            }
        };
        let hx = lx / 2.0;
        let hy = ly / 2.0;
        let hz = lz / 2.0;
        let h = |x: f64, y: f64, z: f64| Vec3::new(x, y, z);
        // +z face (outward +z): du = +x, dv = +y.
        push_grid(
            &mut nodes,
            &mut tris,
            [h(-hx, -hy, hz), h(hx, -hy, hz), h(-hx, hy, hz)],
        );
        // -z face (outward -z): du = +y, dv = +x.
        push_grid(
            &mut nodes,
            &mut tris,
            [h(-hx, -hy, -hz), h(-hx, hy, -hz), h(hx, -hy, -hz)],
        );
        // +x face (outward +x): du = +y, dv = +z.
        push_grid(
            &mut nodes,
            &mut tris,
            [h(hx, -hy, -hz), h(hx, hy, -hz), h(hx, -hy, hz)],
        );
        // -x face (outward -x): du = +z, dv = +y.
        push_grid(
            &mut nodes,
            &mut tris,
            [h(-hx, -hy, -hz), h(-hx, -hy, hz), h(-hx, hy, -hz)],
        );
        // +y face (outward +y): du = +z, dv = +x.
        push_grid(
            &mut nodes,
            &mut tris,
            [h(-hx, hy, -hz), h(-hx, hy, hz), h(hx, hy, -hz)],
        );
        // -y face (outward -y): du = +x, dv = +z.
        push_grid(
            &mut nodes,
            &mut tris,
            [h(-hx, -hy, -hz), h(hx, -hy, -hz), h(-hx, -hy, hz)],
        );
        Self::new(nodes, tris)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn sphere_mesh_has_correct_area() {
        let m = BoundaryMesh::sphere(1.0, 16, 24).unwrap();
        let total: f64 = m.elements().iter().map(|e| e.area).sum();
        // Coarse UV sphere underestimates π·4r² slightly.
        assert!(
            total > 12.0 && total < 4.0 * core::f64::consts::PI + 1e-9,
            "{total}"
        );
        for e in m.elements() {
            assert!(approx(e.normal.norm(), 1.0, 1e-12));
        }
    }

    #[test]
    fn sphere_normals_point_outward() {
        let m = BoundaryMesh::sphere(2.0, 8, 12).unwrap();
        for e in m.elements() {
            // n · centroid ≈ |centroid| (radial).
            assert!(e.normal.dot(e.centroid) > 0.0);
        }
    }

    #[test]
    fn patch_area_and_normal() {
        let m = BoundaryMesh::rectangle_patch(2.0, 3.0, 4, 6).unwrap();
        let total: f64 = m.elements().iter().map(|e| e.area).sum();
        assert!(approx(total, 6.0, 1e-12));
        for e in m.elements() {
            assert!(approx(e.normal.z, 1.0, 1e-12));
            assert!(approx(e.centroid.z, 0.0, 1e-12));
        }
    }

    #[test]
    fn box_surface_area() {
        let m = BoundaryMesh::box_surface(2.0, 3.0, 4.0, 2).unwrap();
        let total: f64 = m.elements().iter().map(|e| e.area).sum();
        assert!(approx(total, 2.0 * (6.0 + 8.0 + 12.0), 1e-9));
    }

    #[test]
    fn rejects_bad_meshes() {
        assert!(BoundaryMesh::sphere(0.0, 8, 8).is_err());
        assert!(BoundaryMesh::sphere(1.0, 1, 8).is_err());
        assert!(BoundaryMesh::rectangle_patch(-1.0, 1.0, 2, 2).is_err());
        assert!(BoundaryMesh::new(
            alloc::vec![Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0)],
            alloc::vec![[0, 1, 5]]
        )
        .is_err());
        // Degenerate (collinear) triangle.
        assert!(BoundaryMesh::new(
            alloc::vec![
                Vec3::ZERO,
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(2.0, 0.0, 0.0)
            ],
            alloc::vec![[0, 1, 2]]
        )
        .is_err());
    }
}

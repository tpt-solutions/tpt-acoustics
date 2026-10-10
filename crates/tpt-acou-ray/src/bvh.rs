//! Bounding volume hierarchy (BVH) over triangles — the from-scratch spatial
//! index used for ray-scene intersection.
//!
//! Median split over the widest axis of the centroid bounds; traversal is
//! iterative (no recursion, no stack overflow for large meshes). This module
//! replaces the external `tpt-percept-cloud` dependency named in spec.txt §3
//! (see todo.md "Open Decisions").

use alloc::vec::Vec;

use crate::vec3::Vec3;

/// Axis-aligned bounding box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aabb {
    /// Lower corner.
    pub min: Vec3,
    /// Upper corner.
    pub max: Vec3,
}

impl Aabb {
    /// Degenerate empty box (min > max on some axis).
    #[must_use]
    pub const fn empty() -> Self {
        Self {
            min: Vec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
            max: Vec3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
        }
    }

    /// Grows the box to contain a point.
    pub fn grow(&mut self, p: Vec3) {
        self.min = self.min.min_components(p);
        self.max = self.max.max_components(p);
    }

    /// Box of a triangle.
    #[must_use]
    pub fn triangle(a: Vec3, b: Vec3, c: Vec3) -> Self {
        let mut box_ = Self::empty();
        box_.grow(a);
        box_.grow(b);
        box_.grow(c);
        box_
    }

    /// Widest axis (0, 1, 2).
    #[must_use]
    pub fn widest_axis(&self) -> usize {
        let ext = self.max - self.min;
        if ext.x >= ext.y && ext.x >= ext.z {
            0
        } else if ext.y >= ext.z {
            1
        } else {
            2
        }
    }

    /// Slab-test intersection of a ray with the box. Returns the entry
    /// distance if the ray hits within `[t_min, t_max]`.
    #[must_use]
    pub fn hit(&self, origin: Vec3, inv_dir: Vec3, t_min: f64, t_max: f64) -> Option<f64> {
        let mut t0 = t_min;
        let mut t1 = t_max;
        let o = [origin.x, origin.y, origin.z];
        let d = [inv_dir.x, inv_dir.y, inv_dir.z];
        let lo = [self.min.x, self.min.y, self.min.z];
        let hi = [self.max.x, self.max.y, self.max.z];
        for i in 0..3 {
            let mut tnear = (lo[i] - o[i]) * d[i];
            let mut tfar = (hi[i] - o[i]) * d[i];
            if tnear > tfar {
                core::mem::swap(&mut tnear, &mut tfar);
            }
            t0 = t0.max(tnear);
            t1 = t1.min(tfar);
            if t0 > t1 {
                return None;
            }
        }
        Some(t0)
    }
}

/// A triangle registered in the BVH (positions only; material mapping lives
/// in [`crate::scene::Scene`]).
#[derive(Debug, Clone, Copy)]
pub struct Tri {
    /// Vertices.
    pub v: [Vec3; 3],
    /// Cached bounding box.
    pub bbox: Aabb,
    /// Cached centroid.
    pub centroid: Vec3,
}

impl Tri {
    /// Builds a triangle from vertices.
    #[must_use]
    pub fn new(a: Vec3, b: Vec3, c: Vec3) -> Self {
        Self {
            v: [a, b, c],
            bbox: Aabb::triangle(a, b, c),
            centroid: (a + b + c) / 3.0,
        }
    }

    /// Geometric (unnormalized) normal ` (b−a)×(c−a)`; magnitude is twice
    /// the triangle area (zero for degenerate triangles).
    #[must_use]
    pub fn normal_unnormalized(&self) -> Vec3 {
        (self.v[1] - self.v[0]).cross(self.v[2] - self.v[0])
    }
}

/// Intersection record.
#[derive(Debug, Clone, Copy)]
pub struct Hit {
    /// Index of the hit triangle.
    pub index: usize,
    /// Ray parameter (distance along the unit direction).
    pub t: f64,
    /// Hit point.
    pub point: Vec3,
    /// Unit normal, oriented against the ray direction (two-sided).
    pub normal: Vec3,
}

/// BVH node: internal nodes point at two children; leaves reference a slice
/// of the index list.
#[derive(Debug, Clone, Copy)]
struct Node {
    bbox: Aabb,
    // Internal: children indices. Leaf: `start` into indices, `count` items.
    left_or_start: u32,
    right_or_count: u32,
    is_leaf: bool,
}

/// A built BVH over a triangle list.
#[derive(Debug, Clone, Default)]
pub struct Bvh {
    tris: Vec<Tri>,
    nodes: Vec<Node>,
    indices: Vec<u32>,
    // Nodes are created children-first, so the root is the *last* node.
    root: Option<u32>,
}

impl Bvh {
    /// Builds a BVH over triangles. Degenerate triangles (zero area) are
    /// kept but can never be hit (their Möller–Trumbore determinant is zero).
    #[must_use]
    pub fn new(tris_in: &[Tri]) -> Self {
        let mut bvh = Self {
            tris: Vec::with_capacity(tris_in.len()),
            nodes: Vec::new(),
            indices: Vec::with_capacity(tris_in.len()),
            root: None,
        };
        bvh.tris.extend_from_slice(tris_in);
        for i in 0..tris_in.len() as u32 {
            bvh.indices.push(i);
        }
        if !tris_in.is_empty() {
            bvh.root = Some(bvh.build(0, tris_in.len() as u32));
        }
        bvh
    }

    /// Number of triangles.
    #[must_use]
    pub fn len(&self) -> usize {
        self.tris.len()
    }

    /// True if empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.tris.is_empty()
    }

    fn new_node(
        &mut self,
        bbox: Aabb,
        left_or_start: u32,
        right_or_count: u32,
        is_leaf: bool,
    ) -> u32 {
        let idx = self.nodes.len() as u32;
        self.nodes.push(Node {
            bbox,
            left_or_start,
            right_or_count,
            is_leaf,
        });
        idx
    }

    fn build(&mut self, start: u32, count: u32) -> u32 {
        let mut bbox = Aabb::empty();
        for &i in &self.indices[start as usize..(start + count) as usize] {
            bbox.grow(self.tris[i as usize].bbox.min);
            bbox.grow(self.tris[i as usize].bbox.max);
        }
        if count <= 4 {
            return self.new_node(bbox, start, count, true);
        }
        // Median split on the widest axis of the centroid bounds.
        let axis = bbox.widest_axis();
        let slice = &mut self.indices[start as usize..(start + count) as usize];
        slice.sort_by(|&a, &b| {
            let pa = self.tris[a as usize].centroid;
            let pb = self.tris[b as usize].centroid;
            let ca = [pa.x, pa.y, pa.z][axis];
            let cb = [pb.x, pb.y, pb.z][axis];
            ca.partial_cmp(&cb).unwrap_or(core::cmp::Ordering::Equal)
        });
        let mid = count / 2;
        let left = self.build(start, mid);
        let right = self.build(start + mid, count - mid);
        self.new_node(bbox, left, right, false)
    }

    /// Nearest intersection of the ray `(origin, dir)` with the scene within
    /// `[t_min, t_max)`; `dir` need not be normalized (returns parameter t).
    #[must_use]
    pub fn hit(&self, origin: Vec3, dir: Vec3, t_min: f64, t_max: f64) -> Option<Hit> {
        let root = self.root?;
        let inv_dir = Vec3::new(1.0 / dir.x, 1.0 / dir.y, 1.0 / dir.z);
        let mut best: Option<Hit> = None;
        let mut best_t = t_max;
        // Iterative traversal with an explicit stack (fixed max depth bound:
        // median split gives depth ~ log2(n), 64 is far beyond any mesh that
        // fits in memory).
        let mut stack = [0u32; 64];
        let mut sp = 0;
        stack[sp] = root;
        sp += 1;
        while sp > 0 {
            sp -= 1;
            let node = &self.nodes[stack[sp] as usize];
            if node.bbox.hit(origin, inv_dir, t_min, best_t).is_none() {
                continue;
            }
            if node.is_leaf {
                let start = node.left_or_start as usize;
                let count = node.right_or_count as usize;
                for &i in &self.indices[start..start + count] {
                    if let Some((t_hit, point, normal)) =
                        tri_hit(&self.tris[i as usize], origin, dir, t_min, best_t)
                    {
                        best_t = t_hit;
                        best = Some(Hit {
                            index: i as usize,
                            t: t_hit,
                            point,
                            normal,
                        });
                    }
                }
            } else {
                stack[sp] = node.left_or_start;
                sp += 1;
                stack[sp] = node.right_or_count;
                sp += 1;
            }
        }
        best
    }
}

/// Möller–Trumbore intersection, two-sided.
fn tri_hit(
    tri: &Tri,
    origin: Vec3,
    dir: Vec3,
    t_min: f64,
    t_max: f64,
) -> Option<(f64, Vec3, Vec3)> {
    let e1 = tri.v[1] - tri.v[0];
    let e2 = tri.v[2] - tri.v[0];
    let p = dir.cross(e2);
    let det = e1.dot(p);
    if det.abs() < 1e-14 {
        return None;
    }
    let inv_det = 1.0 / det;
    let s = origin - tri.v[0];
    let u = s.dot(p) * inv_det;
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = s.cross(e1);
    let v = dir.dot(q) * inv_det;
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = e2.dot(q) * inv_det;
    if t < t_min || t >= t_max {
        return None;
    }
    let n_raw = tri.normal_unnormalized();
    let n = n_raw.normalized();
    let normal = if n.dot(dir) > 0.0 { -n } else { n };
    Some((t, origin + dir * t, normal))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn brute_force(tris: &[Tri], origin: Vec3, dir: Vec3) -> Option<(usize, f64)> {
        let mut best: Option<(usize, f64)> = None;
        for (i, t) in tris.iter().enumerate() {
            if let Some((t_hit, _, _)) =
                tri_hit(t, origin, dir, 1e-6, best.map_or(f64::INFINITY, |b| b.1))
            {
                best = Some((i, t_hit));
            }
        }
        best
    }

    #[test]
    fn bvh_matches_brute_force_on_box_scene() {
        // A closed cube of 12 triangles.
        let s = 2.0;
        let mut tris = Vec::new();
        let corners = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(s, 0.0, 0.0),
            Vec3::new(s, s, 0.0),
            Vec3::new(0.0, s, 0.0),
            Vec3::new(0.0, 0.0, s),
            Vec3::new(s, 0.0, s),
            Vec3::new(s, s, s),
            Vec3::new(0.0, s, s),
        ];
        let quads = [
            [0, 1, 2, 3],
            [4, 5, 6, 7],
            [0, 1, 5, 4],
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [3, 0, 4, 7],
        ];
        for q in quads {
            let a = corners[q[0]];
            let b = corners[q[1]];
            let c = corners[q[2]];
            let d = corners[q[3]];
            tris.push(Tri::new(a, b, c));
            tris.push(Tri::new(a, c, d));
        }
        let bvh = Bvh::new(&tris);
        let mut rng = crate::rng::Xoshiro256::new(11);
        for _ in 0..500 {
            // Rays from outside toward the cube center.
            let origin = Vec3::new(
                rng.next_range(-5.0, 7.0),
                rng.next_range(-5.0, 7.0),
                rng.next_range(-5.0, 7.0),
            );
            let dir = (Vec3::new(1.0, 1.0, 1.0) - origin).normalized();
            let got = bvh
                .hit(origin, dir, 1e-6, f64::INFINITY)
                .map(|h| (h.index, h.t));
            let want = brute_force(&tris, origin, dir);
            match (got, want) {
                (Some((gi, gt)), Some((wi, wt))) => {
                    assert_eq!(gi, wi, "index mismatch for origin {origin:?}");
                    assert!((gt - wt).abs() < 1e-9);
                }
                (None, None) => {}
                _ => panic!("hit disagreement for origin {origin:?}"),
            }
        }
    }

    #[test]
    fn empty_bvh_hits_nothing() {
        let bvh = Bvh::new(&[]);
        assert!(bvh
            .hit(Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0), 0.0, 10.0)
            .is_none());
    }

    #[test]
    fn aabb_slab_test_basics() {
        let b = Aabb::triangle(
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(2.0, 1.2, 1.2),
            Vec3::new(1.2, 2.0, 2.0),
        );
        // Ray through the box interior along +x hits with t0 in [1, 2].
        let dir = Vec3::new(1.0, 0.0, 0.0);
        let inv = Vec3::new(1.0 / dir.x, 1.0 / dir.y, 1.0 / dir.z);
        let hit = b
            .hit(Vec3::new(0.0, 1.5, 1.5), inv, 0.0, 100.0)
            .expect("expected a hit");
        assert!((1.0..=2.0).contains(&hit));
        // Parallel to the box in y/z from far outside: a miss.
        let miss = b.hit(Vec3::new(0.0, 10.0, 10.0), inv, 0.0, 100.0);
        assert!(miss.is_none());
    }
}

//! Triangle scenes with per-triangle materials, accelerated by the BVH.

use alloc::vec::Vec;

use crate::bvh::{Bvh, Hit, Tri};
use crate::ism::ShoeboxRoom;
use crate::material::BoundaryMaterial;
use crate::vec3::Vec3;

/// A triangle with an associated material index (into [`Scene::materials`]).
#[derive(Debug, Clone, Copy)]
pub struct SceneTriangle {
    /// Vertex positions.
    pub v: [Vec3; 3],
    /// Index into the scene's material list.
    pub material: usize,
}

/// A triangle-mesh acoustic scene.
#[derive(Debug, Clone, Default)]
pub struct Scene {
    tris: Vec<Tri>,
    material_of: Vec<usize>,
    materials: Vec<BoundaryMaterial>,
    bvh: Bvh,
}

/// Hit record carrying the surface material.
#[derive(Debug, Clone, Copy)]
pub struct SceneHit {
    /// Distance along the ray.
    pub t: f64,
    /// Hit point.
    pub point: Vec3,
    /// Unit normal oriented against the incoming ray.
    pub normal: Vec3,
    /// Index into the scene's material list.
    pub material: usize,
}

impl Scene {
    /// Builds a scene from triangles and a material table.
    ///
    /// # Errors
    /// Returns [`tpt_acou_core::AcousticError::InvalidArgument`] if a
    /// triangle references an unknown material index.
    pub fn new(
        triangles: &[SceneTriangle],
        materials: Vec<BoundaryMaterial>,
    ) -> Result<Self, tpt_acou_core::AcousticError> {
        for t in triangles {
            if t.material >= materials.len() {
                return Err(tpt_acou_core::AcousticError::InvalidArgument(
                    "triangle material index out of range",
                ));
            }
        }
        let tris: Vec<Tri> = triangles
            .iter()
            .map(|t| Tri::new(t.v[0], t.v[1], t.v[2]))
            .collect();
        let material_of = triangles.iter().map(|t| t.material).collect();
        let bvh = Bvh::new(&tris);
        Ok(Self {
            tris,
            material_of,
            materials,
            bvh,
        })
    }

    /// The material table.
    #[must_use]
    pub fn materials(&self) -> &[BoundaryMaterial] {
        &self.materials
    }

    /// Number of triangles.
    #[must_use]
    pub fn len(&self) -> usize {
        self.tris.len()
    }

    /// True if the scene has no triangles.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.tris.is_empty()
    }

    /// Nearest hit of `(origin, dir)` with `t < t_max`.
    #[must_use]
    pub fn hit(&self, origin: Vec3, dir: Vec3, t_max: f64) -> Option<SceneHit> {
        let Hit {
            t,
            point,
            normal,
            index,
            ..
        } = self.bvh.hit(origin, dir, 1e-9, t_max)?;
        Some(SceneHit {
            t,
            point,
            normal,
            material: self.material_of[index],
        })
    }

    /// Builds a closed shoebox mesh (12 triangles) with the room's per-wall
    /// materials, so the Monte Carlo tracer can treat shoebox and arbitrary
    /// scenes uniformly.
    ///
    /// # Errors
    /// Propagates material construction errors from the room (none) and
    /// material-index validation from [`Scene::new`].
    pub fn from_shoebox(room: &ShoeboxRoom) -> Result<Self, tpt_acou_core::AcousticError> {
        let [lx, ly, lz] = room.dims();
        let p = |x: f64, y: f64, z: f64| Vec3::new(x, y, z);
        let mut tris = Vec::with_capacity(12);
        // Quads per wall: (four corners in CCW order when viewed from inside,
        // so quad(a,b,c,d) -> tris (a,b,c),(a,c,d) face inward).
        let walls: [(usize, [Vec3; 4]); 6] = [
            (
                crate::ism::WALL_LOW_X,
                [
                    p(0.0, 0.0, 0.0),
                    p(0.0, ly, 0.0),
                    p(0.0, ly, lz),
                    p(0.0, 0.0, lz),
                ],
            ),
            (
                crate::ism::WALL_HIGH_X,
                [
                    p(lx, 0.0, lz),
                    p(lx, ly, lz),
                    p(lx, ly, 0.0),
                    p(lx, 0.0, 0.0),
                ],
            ),
            (
                crate::ism::WALL_LOW_Y,
                [
                    p(0.0, 0.0, 0.0),
                    p(0.0, 0.0, lz),
                    p(lx, 0.0, lz),
                    p(lx, 0.0, 0.0),
                ],
            ),
            (
                crate::ism::WALL_HIGH_Y,
                [
                    p(0.0, ly, lz),
                    p(0.0, ly, 0.0),
                    p(lx, ly, 0.0),
                    p(lx, ly, lz),
                ],
            ),
            (
                crate::ism::WALL_LOW_Z,
                [
                    p(0.0, 0.0, 0.0),
                    p(lx, 0.0, 0.0),
                    p(lx, ly, 0.0),
                    p(0.0, ly, 0.0),
                ],
            ),
            (
                crate::ism::WALL_HIGH_Z,
                [
                    p(0.0, 0.0, lz),
                    p(0.0, ly, lz),
                    p(lx, ly, lz),
                    p(lx, 0.0, lz),
                ],
            ),
        ];
        for (mat, c) in walls {
            tris.push(SceneTriangle {
                v: [c[0], c[1], c[2]],
                material: mat,
            });
            tris.push(SceneTriangle {
                v: [c[0], c[2], c[3]],
                material: mat,
            });
        }
        let materials: Vec<BoundaryMaterial> = (0..6).map(|w| room.wall(w).clone()).collect();
        Self::new(&tris, materials)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shoebox_mesh_is_hit_from_any_interior_direction() {
        let room =
            ShoeboxRoom::uniform([4.0, 3.0, 2.5], BoundaryMaterial::concrete_painted()).unwrap();
        let scene = Scene::from_shoebox(&room).unwrap();
        assert_eq!(scene.len(), 12);
        assert_eq!(scene.materials().len(), 6);
        let mut rng = crate::rng::Xoshiro256::new(3);
        let center = Vec3::new(2.0, 1.5, 1.25);
        for _ in 0..200 {
            let d = rng.next_unit_direction();
            let hit = scene
                .hit(center, d, f64::INFINITY)
                .expect("closed box must hit");
            // Hit stays inside the room bounds.
            assert!(hit.point.x >= -1e-9 && hit.point.x <= 4.0 + 1e-9);
            assert!(hit.point.y >= -1e-9 && hit.point.y <= 3.0 + 1e-9);
            assert!(hit.point.z >= -1e-9 && hit.point.z <= 2.5 + 1e-9);
            // Normal faces back toward the origin of the ray.
            assert!(hit.normal.dot(d) < 0.0);
        }
    }

    #[test]
    fn per_wall_materials_are_preserved() {
        let room = ShoeboxRoom::new(
            [4.0, 3.0, 2.5],
            [
                BoundaryMaterial::glass(),
                BoundaryMaterial::rigid(),
                BoundaryMaterial::curtain_heavy(),
                BoundaryMaterial::rigid(),
                BoundaryMaterial::absorber_panel(),
                BoundaryMaterial::rigid(),
            ],
        )
        .unwrap();
        let scene = Scene::from_shoebox(&room).unwrap();
        let origin = Vec3::new(2.0, 1.5, 1.25);
        // Straight at the low-x wall: material should be glass.
        let hit = scene
            .hit(origin, Vec3::new(-1.0, 0.0, 0.0), f64::INFINITY)
            .unwrap();
        assert_eq!(scene.materials()[hit.material].name(), "glass");
        // Straight up: rigid.
        let hit = scene
            .hit(origin, Vec3::new(0.0, 0.0, 1.0), f64::INFINITY)
            .unwrap();
        assert_eq!(scene.materials()[hit.material].name(), "rigid");
    }

    #[test]
    fn rejects_unknown_material_index() {
        let t = SceneTriangle {
            v: [
                Vec3::ZERO,
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            ],
            material: 7,
        };
        assert!(Scene::new(&[t], alloc::vec![BoundaryMaterial::rigid()]).is_err());
    }
}

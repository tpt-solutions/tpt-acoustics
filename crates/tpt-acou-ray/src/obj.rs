//! Wavefront OBJ import for ray-tracing scenes.
//!
//! Reads vertex (`v`) and face (`f`) records, plus `usemtl` names for
//! grouping. Faces may be triangles or convex polygons (fan-triangulated),
//! and indices may take the `v`, `v/vt`, `v//vn` or `v/vt/vn` forms, with
//! negative indices counting back from the last vertex. Other records
//! (normals, texture coordinates, comments, groups) are ignored.
//!
//! The reader is `no_std` + `alloc`. It takes the OBJ text as a string, so
//! callers decide how to load the file.

use alloc::string::{String, ToString};
use alloc::vec;
use alloc::vec::Vec;

use tpt_acou_core::error::AcousticError;

use crate::material::BoundaryMaterial;
use crate::scene::{Scene, SceneTriangle};
use crate::vec3::Vec3;

/// Triangulated mesh read from an OBJ file.
#[derive(Debug, Clone, PartialEq)]
pub struct ObjMesh {
    /// Vertex positions, in file order.
    pub vertices: Vec<Vec3>,
    /// Triangles as vertex indices into `vertices` (0-based).
    pub triangles: Vec<[usize; 3]>,
    /// The `usemtl` name active when each triangle was declared, or an empty
    /// string before any `usemtl` record. Same length as `triangles`.
    pub material_names: Vec<String>,
}

fn invalid(message: &'static str) -> AcousticError {
    AcousticError::InvalidArgument(message)
}

fn parse_coord(text: Option<&str>) -> Result<f64, AcousticError> {
    let value: f64 = text
        .ok_or(invalid("OBJ vertex needs three coordinates"))?
        .parse()
        .map_err(|_| invalid("OBJ vertex coordinate is not a number"))?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("OBJ vertex coordinate is not finite"))
    }
}

/// Resolves one face index (1-based, or negative from the end) to 0-based.
fn resolve_index(token: &str, vertex_count: usize) -> Result<usize, AcousticError> {
    // Only the position index matters; drop texture and normal references.
    let head = token.split('/').next().unwrap_or("");
    let raw: i64 = head
        .parse()
        .map_err(|_| invalid("OBJ face index is not an integer"))?;
    let index = if raw > 0 {
        raw - 1
    } else if raw < 0 {
        vertex_count as i64 + raw
    } else {
        return Err(invalid("OBJ face index must not be zero"));
    };
    if index < 0 || index as usize >= vertex_count {
        return Err(invalid("OBJ face refers to a vertex that does not exist"));
    }
    Ok(index as usize)
}

/// Parses OBJ text into a triangle mesh.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] for malformed vertex or face
/// records, indices outside the vertex list, faces with fewer than three
/// vertices, or a file with no faces.
pub fn parse_obj(text: &str) -> Result<ObjMesh, AcousticError> {
    let mut vertices: Vec<Vec3> = Vec::new();
    let mut triangles: Vec<[usize; 3]> = Vec::new();
    let mut material_names: Vec<String> = Vec::new();
    let mut current_material = String::new();

    for line in text.lines() {
        let line = line.split('#').next().unwrap_or("").trim();
        let mut parts = line.split_whitespace();
        match parts.next() {
            Some("v") => {
                let x = parse_coord(parts.next())?;
                let y = parse_coord(parts.next())?;
                let z = parse_coord(parts.next())?;
                vertices.push(Vec3::new(x, y, z));
            }
            Some("usemtl") => {
                current_material = parts.next().unwrap_or("").to_string();
            }
            Some("f") => {
                let indices: Vec<usize> = parts
                    .map(|token| resolve_index(token, vertices.len()))
                    .collect::<Result<_, _>>()?;
                if indices.len() < 3 {
                    return Err(invalid("OBJ face needs at least three vertices"));
                }
                // Fan triangulation from the first vertex (convex polygons).
                for k in 1..indices.len() - 1 {
                    triangles.push([indices[0], indices[k], indices[k + 1]]);
                    material_names.push(current_material.clone());
                }
            }
            _ => {}
        }
    }
    if triangles.is_empty() {
        return Err(invalid("OBJ file contains no faces"));
    }
    Ok(ObjMesh {
        vertices,
        triangles,
        material_names,
    })
}

/// Builds a ray-tracing scene from OBJ text, giving every face the same
/// wall material.
///
/// # Errors
/// Propagates [`parse_obj`] errors and the scene constructor's errors.
pub fn scene_from_obj(text: &str, material: BoundaryMaterial) -> Result<Scene, AcousticError> {
    let mesh = parse_obj(text)?;
    let triangles: Vec<SceneTriangle> = mesh
        .triangles
        .iter()
        .map(|&[a, b, c]| SceneTriangle {
            v: [mesh.vertices[a], mesh.vertices[b], mesh.vertices[c]],
            material: 0,
        })
        .collect();
    Scene::new(&triangles, vec![material])
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Unit cube as six quads, with a `usemtl` group on the first two faces.
    const CUBE: &str = "\
# unit cube
v 0 0 0
v 1 0 0
v 1 1 0
v 0 1 0
v 0 0 1
v 1 0 1
v 1 1 1
v 0 1 1
usemtl wall
f 1 4 3 2
f 5 6 7 8
usemtl floor
f 1 2 6 5
f 2 3 7 6
f 3 4 8 7
f 4 1 5 8
";

    fn concrete() -> BoundaryMaterial {
        BoundaryMaterial::concrete_painted()
    }

    #[test]
    fn cube_triangulates_to_twelve_triangles() {
        let mesh = parse_obj(CUBE).unwrap();
        assert_eq!(mesh.vertices.len(), 8);
        assert_eq!(mesh.triangles.len(), 12);
        assert_eq!(mesh.material_names.len(), 12);
        assert_eq!(mesh.material_names[0], "wall");
        assert_eq!(mesh.material_names[4], "floor");
    }

    #[test]
    fn slash_and_negative_indices_are_resolved() {
        let text = "v 0 0 0\nv 1 0 0\nv 0 1 0\nvt 0 0\nvn 0 0 1\nf -3/1/1 -2//1 -1/1\n";
        let mesh = parse_obj(text).unwrap();
        assert_eq!(mesh.triangles, vec![[0, 1, 2]]);
    }

    #[test]
    fn bad_input_is_rejected() {
        assert!(parse_obj("v 0 0\nf 1 1 1\n").is_err(), "short vertex");
        assert!(parse_obj("v 0 0 0\nf 1 2 3\n").is_err(), "missing vertices");
        assert!(
            parse_obj("v 0 0 0\nv 1 0 0\nf 1 2\n").is_err(),
            "two-vertex face"
        );
        assert!(parse_obj("v 0 0 0\n").is_err(), "no faces");
        assert!(
            parse_obj("v 0 0 nan\nv 1 0 0\nv 0 1 0\nf 1 2 3\n").is_err(),
            "nan"
        );
        assert!(
            parse_obj("v 0 0 0\nv 1 0 0\nv 0 1 0\nf 0 1 2\n").is_err(),
            "zero index"
        );
    }

    #[test]
    fn scene_from_cube_is_closed_and_hit_from_inside() {
        let scene = scene_from_obj(CUBE, concrete()).unwrap();
        assert_eq!(scene.len(), 12);
        // Ray from the cube centre along +x hits the x = 1 face at distance 0.5.
        let hit = scene
            .hit(
                Vec3::new(0.5, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                f64::INFINITY,
            )
            .expect("ray must hit the closed cube");
        assert!((hit.t - 0.5).abs() < 1e-9, "t = {}", hit.t);
    }
}

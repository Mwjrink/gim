use std::collections::btree_set::Intersection;
use std::collections::{HashMap, HashSet};
use std::{fs::File, io::prelude::*};
use crate::utils::*;

pub type Triangle = [u32; 3];
pub type Edge = [u32; 2];
pub type Point = [f32; 3];

pub struct Mesh {
    positions: Vec<f32>,
    indices: Vec<u32>,
}

pub fn simplify(mesh: &Mesh, target_tris: u32) -> Mesh {
    let (edges, tris, vertices) =
        generate_data_structures(&(mesh.indices.len() / 3), &mesh.indices);

    let collapsible = Vec::<(f32, Edge)>::with_capacity(edges.len());
    for edge in edges {
        if edge.1[1] != u32::MAX {
            collapsible.push((mesh.sqr_dist(edge.0[0], edge.0[1]), edge.0));
        } else {
            assert!(edge.1[0] != u32::MAX, "edge is not properly ordered");
            // collapsible.push((mesh.sqr_dist(edge.0[0], edge.0[1]), edge.0));
        }
    }

    let mut incident_mesh = Mesh { 
        positions: mesh.positions.clone(),
        indices: mesh.indices.clone()
    };

    collapsible.sort_by(|one, two| one.0.partial_cmp(&two.0).unwrap());

    while incident_mesh.indices.len() > (target_tris * 3) as usize {
        let is_collapsible = |edge: &Edge| {
            // exactly 2 neigbouring vertices for each element in edge
            let intsxn = sorted_vec_intersection_count(vertices.get(&edge[0]).unwrap(), vertices.get(&edge[1]).unwrap());
            if intsxn == 2 {
                // TODO detect flipped triangles
                
                
                return true
            }

            false
        };

        let collapse = |edge: &Edge| {
            let midpoint = mesh.mid_point_2(edge[0], edge[1]);
            // remove the original 2 indices, the vec shortens
            incident_mesh.positions.remove((edge[0]+0) as usize);
            incident_mesh.positions.remove((edge[0]+1) as usize);
            incident_mesh.positions.remove((edge[0]+2) as usize);
            
            incident_mesh.positions.remove((edge[1]+0) as usize);
            incident_mesh.positions.remove((edge[1]+1) as usize);
            incident_mesh.positions.remove((edge[1]+2) as usize);

            let new_idx = incident_mesh.positions.len();
            incident_mesh.positions.extend_from_slice(&midpoint);
            incident_mesh.indices.push(new_idx as u32);

            // move all indices pointing to an index that got moved
            for idx in &mut incident_mesh.indices {
                if *idx == edge[0] || *idx == edge[1] {
                    *idx = new_idx as u32;
                }

                if *idx > edge[0] {
                    *idx -= 1;
                }

                if *idx > edge[1] {
                    *idx -= 1;
                }
            }
        };

        let idx = 0;
        while !is_collapsible(&collapsible[idx].1) {
            idx += 1;
        }

        collapse(&collapsible[idx].1);

        // TODO recalculate collapsible
    }

    incident_mesh
}

pub fn generate_data_structures(
    tri_count: &usize,
    indices: &Vec<u32>,
) -> (
    HashMap<Edge, [u32; 2]>,
    HashSet<Triangle>,
    HashMap<u32, Vec<u32>>,
) {
    let mut edges = HashMap::<Edge, [u32; 2]>::with_capacity(tri_count * 24);
    let mut tris = HashSet::<Triangle>::with_capacity(*tri_count);
    let mut vertices: HashMap<u32, Vec<u32>> = HashMap::new();
    //

    // SECTION - Fill vertices, tris, edges and boundary edges
    {
        // * fill vertices, tris and edges
        let mut idx = 0;
        for _ in 0..*tri_count {
            let idx0 = indices[idx + 0];
            let idx1 = indices[idx + 1];
            let idx2 = indices[idx + 2];

            {
                let key: Edge = new_edge(idx0, idx1);
                let value = edges.insert(key, [idx2, u32::MAX]);
                if let Some(v) = value {
                    if v[1] != u32::MAX {
                        println!(
                            "this edge has appeared more than twice: {:?} => {:?}",
                            idx2, v
                        );
                    }
                    edges.insert(key, [idx2, v[0]]);
                }
            }
            {
                let key: Edge = new_edge(idx0, idx2);
                let value = edges.insert(key, [idx1, u32::MAX]);
                if let Some(v) = value {
                    if v[1] != u32::MAX {
                        println!(
                            "this edge has appeared more than twice: {:?} => {:?}",
                            idx2, v
                        );
                    }
                    edges.insert(key, [idx1, v[0]]);
                }
            }
            {
                let key: Edge = new_edge(idx1, idx2);
                let value = edges.insert(key, [idx0, u32::MAX]);
                if let Some(v) = value {
                    if v[1] != u32::MAX {
                        println!(
                            "this edge has appeared more than twice: {:?} => {:?}",
                            idx2, v
                        );
                    }
                    edges.insert(key, [idx0, v[0]]);
                }
            }

            {
                vertices
                    .entry(idx0)
                    .or_insert(Vec::with_capacity(8))
                    .push(idx1);
                vertices
                    .entry(idx0)
                    .or_insert(Vec::with_capacity(8))
                    .push(idx2);

                vertices
                    .entry(idx1)
                    .or_insert(Vec::with_capacity(8))
                    .push(idx0);
                vertices
                    .entry(idx1)
                    .or_insert(Vec::with_capacity(8))
                    .push(idx2);

                vertices
                    .entry(idx2)
                    .or_insert(Vec::with_capacity(8))
                    .push(idx0);
                vertices
                    .entry(idx2)
                    .or_insert(Vec::with_capacity(8))
                    .push(idx1);

                // this removes duplicates that happen due to shared verts in triangles
                vertices.entry(idx0).or_default().sort_unstable();
                vertices.entry(idx0).or_default().dedup();
                vertices.entry(idx1).or_default().sort_unstable();
                vertices.entry(idx1).or_default().dedup();
                vertices.entry(idx2).or_default().sort_unstable();
                vertices.entry(idx2).or_default().dedup();
            };

            tris.insert(new_tri(idx0, idx1, idx2));

            idx += 3 as usize;
        }
    };

    for kv in &mut vertices {
        kv.1.sort();
    }

    (edges, tris, vertices)
}

// orders the vertices consistently
pub fn new_tri(a: u32, b: u32, c: u32) -> Triangle {
    if a <= b {
        if b <= c {
            [a, b, c]
        } else if c < a {
            [c, a, b]
        } else {
            [a, c, b]
        }
    } else {
        if a <= c {
            [b, a, c]
        } else if c < b {
            [c, b, a]
        } else {
            [b, c, a]
        }
    }
}

// orders the vertices consistently
pub fn new_edge(a: u32, b: u32) -> Edge {
    if a <= b {
        [a, b]
    } else {
        [b, a]
    }
}

pub fn write_mesh_to_file(name: &str, positions: &Vec<f32>, mesh: &Vec<u32>) {
    let mut f = File::create(name.to_string() + ".obj").unwrap();
    f.write_all(b"mtllib bunny.mtl\no bun_zipper\n").unwrap();
    for v_idx in (0..positions.len()).step_by(3) {
        f.write_all(
            format!(
                "v {} {} {}\n",
                positions[v_idx + 0],
                positions[v_idx + 1],
                positions[v_idx + 2]
            )
            .as_bytes(),
        )
        .unwrap();
    }
    f.write_all(b"usemtl None\ns off\n").unwrap();

    f.write_all(format!("o {}\n", name).as_bytes()).unwrap();
    for idx in (0..mesh.len()).step_by(3) {
        f.write_all(
            format!(
                "f {} {} {}\n",
                mesh[idx + 0] + 1,
                mesh[idx + 1] + 1,
                mesh[idx + 2] + 1
            )
            .as_bytes(),
        )
        .unwrap();
    }
    f.write_all("\n".as_bytes()).unwrap();

    f.flush().unwrap();
}

impl Mesh {
    fn sqr_dist(&self, a: u32, b: u32) -> f32 {
        let p1 = &self.positions[(a * 3) as usize..(a * 3 + 3) as usize];
        let p2 = &self.positions[(b * 3) as usize..(b * 3 + 3) as usize];

        f32::sqrt((p1[0] - p2[0]).powf(2.0) + (p1[1] - p2[1]).powf(2.0) + (p1[2] - p2[2]).powf(2.0))
    }

    fn sqr_dist_w_custom(&self, a: u32, b: Point) -> f32 {
        let p1 = &self.positions[(a * 3) as usize..(a * 3 + 3) as usize];

        f32::sqrt((p1[0] - b[0]).powf(2.0) + (p1[1] - b[1]).powf(2.0) + (p1[2] - b[2]).powf(2.0))
    }

    fn sqr_dist_w_2custom(a: Point, b: Point) -> f32 {
        f32::sqrt((a[0] - b[0]).powf(2.0) + (a[1] - b[1]).powf(2.0) + (a[2] - b[2]).powf(2.0))
    }

    fn mid_point_2(&self, a: u32, b: u32) -> Point {
        let p1 = &self.positions[(a * 3) as usize..(a * 3 + 3) as usize];
        let p2 = &self.positions[(b * 3) as usize..(b * 3 + 3) as usize];

        [
            (p1[0] + p2[0]) * 0.5,
            (p1[1] + p2[1]) * 0.5,
            (p1[2] + p2[2]) * 0.5,
        ]
    }

    fn mid_point_3(&self, a: u32, b: u32, c: u32) -> Point {
        let p1 = &self.positions[(a * 3) as usize..(a * 3 + 3) as usize];
        let p2 = &self.positions[(b * 3) as usize..(b * 3 + 3) as usize];
        let p3 = &self.positions[(c * 3) as usize..(c * 3 + 3) as usize];

        [
            (p1[0] + p2[0] + p3[0]) / 3.0,
            (p1[1] + p2[1] + p3[1]) / 3.0,
            (p1[2] + p2[2] + p3[2]) / 3.0,
        ]
    }

    fn mid_point_2custom(a: Point, b: Point) -> Point {
        [
            (a[0] + b[0]) * 0.5,
            (a[1] + b[1]) * 0.5,
            (a[2] + b[2]) * 0.5,
        ]
    }
}

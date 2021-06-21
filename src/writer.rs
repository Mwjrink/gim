// create clusters of 128 triangles
// choose a seed triangle, choose the next vertex corresponding to an edge that is closest to the seed triangle

use std::collections::{HashMap, HashSet};
use std::{fs::File, io::prelude::*, panic};

type Triangle = [u32; 3];
type Edge = [u32; 2];
type Point = [f32; 3];
const TRIS_IN_CLUSTER: usize = 128;

struct Cluster {
    tris: Vec<Triangle>,
    cut: Vec<Edge>,
    anchor: Point,
}

pub fn write(obj_file: &String) {
    let (models, materials) = tobj::load_obj(&obj_file, true).expect("Failed to load file");

    // let albedo_texture = Image::read_from_file(
    //     "./test_assets/ucykfbkfa_8K_Albedo.exr",
    //     read_options::high(), // use multi-core decompression
    // )
    // .unwrap();

    // println!("layers: {:?}", albedo_texture.layers.len());
    // println!("channels: {:?}", albedo_texture.layers[0].channels.len());

    // let mut filex = File::create("./test_assets/chalet.nlvox").unwrap();
    // let mut filey = File::create("./test_assets/chalet.nlvoy").unwrap();
    // let mut filez = File::create("./test_assets/chalet.nlvoz").unwrap();
    // let mut fileu = File::create("./test_assets/chalet.nlvou").unwrap();
    // let mut filev = File::create("./test_assets/chalet.nlvov").unwrap();
    // let mut filea = File::create("./test_assets/chalet.nlvoa").unwrap();

    // {
    //     let albedo = image::open("./test_assets/chalet.jpg").unwrap();

    //     let rgba = albedo.to_rgba8();

    //     // write width and height
    //     filea.write_all(&rgba.width().to_le_bytes()).unwrap();
    //     filea.write_all(&rgba.height().to_le_bytes()).unwrap();
    //     filea.write_all(&rgba.into_raw()).unwrap();
    // }

    // let mut fileg = File::create("./test_assets/ucykfbkfa.nlvog",).unwrap();
    // let mut files = File::create("./test_assets/ucykfbkfa.nlvos",).unwrap();

    println!("# of models: {}", models.len());
    println!("# of materials: {}", materials.len());
    for (i, m) in models.iter().enumerate() {
        let mesh = &m.mesh;
        println!("model[{}].name = \'{}\'", i, m.name);
        println!(
            "model[{}].mesh.material_id = {:?}",
            i,
            mesh.material_id.unwrap_or_else(|| { usize::MAX })
        );

        println!("Tris in model[{}]: {}", i, mesh.num_face_indices.len());
        println!("Number of verts: {}", mesh.positions.len());

        // let length = mesh.positions.len() / 3;
        // filex.write_all(&length.to_le_bytes(),).unwrap();
        // filey.write_all(&length.to_le_bytes(),).unwrap();
        // filez.write_all(&length.to_le_bytes(),).unwrap();

        // fileu.write_all(&length.to_le_bytes(),).unwrap();
        // filev.write_all(&length.to_le_bytes(),).unwrap();

        let sqr_dist = |a: u32, b: u32| {
            let p1 = &mesh.positions[(a * 3) as usize..(a * 3 + 3) as usize];
            let p2 = &mesh.positions[(b * 3) as usize..(b * 3 + 3) as usize];

            f32::sqrt(
                (p1[0] - p2[0]).powf(2.0) + (p1[1] - p2[1]).powf(2.0) + (p1[2] - p2[2]).powf(2.0),
            )
        };

        let sqr_dist_w_custom = |a: u32, b: Point| {
            let p1 = &mesh.positions[(a * 3) as usize..(a * 3 + 3) as usize];
            let p2 = &b;

            f32::sqrt(
                (p1[0] - p2[0]).powf(2.0) + (p1[1] - p2[1]).powf(2.0) + (p1[2] - p2[2]).powf(2.0),
            )
        };

        let mid_point_2 = |a: u32, b: u32| {
            let p1 = &mesh.positions[(a * 3) as usize..(a * 3 + 3) as usize];
            let p2 = &mesh.positions[(b * 3) as usize..(b * 3 + 3) as usize];

            [
                (p1[0] + p2[0]) * 0.5,
                (p1[1] + p2[1]) * 0.5,
                (p1[2] + p2[2]) * 0.5,
            ]
        };

        let mid_point_3 = |a: u32, b: u32, c: u32| {
            let p1 = &mesh.positions[(a * 3) as usize..(a * 3 + 3) as usize];
            let p2 = &mesh.positions[(b * 3) as usize..(b * 3 + 3) as usize];
            let p3 = &mesh.positions[(c * 3) as usize..(c * 3 + 3) as usize];

            [
                (p1[0] + p2[0] + p3[0]) / 3.0,
                (p1[1] + p2[1] + p3[1]) / 3.0,
                (p1[2] + p2[2] + p3[2]) / 3.0,
            ]
        };

        // create data structure that can be referenced to find which indices are included in which triangles
        let mut edges = HashMap::<Edge, [u32; 2]>::with_capacity(mesh.num_face_indices.len() * 24);
        let mut tris = HashSet::<Triangle>::with_capacity(mesh.num_face_indices.len());
        let mut vertices: HashMap<u32, Vec<u32>> = HashMap::new();
        //

        // SECTION - Fill vertices, tris, edges and boundary edges
        {
            // * fill vertices, tris and edges
            let mut idx = 0;
            for verts in &mesh.num_face_indices {
                let idx0 = mesh.indices[idx + 0];
                let idx1 = mesh.indices[idx + 1];
                let idx2 = mesh.indices[idx + 2];

                {
                    let key: Edge = new_edge(idx0, idx1);
                    let value = edges.insert(key, [idx2, u32::MAX]);
                    if let Some(v) = value {
                        if v[1] != u32::MAX {
                            println!("this edge has appeared more than twice: {:?}", key);
                        }
                        edges.insert(key, [idx2, v[0]]);
                    }
                }
                {
                    let key: Edge = new_edge(idx0, idx2);
                    let value = edges.insert(key, [idx1, u32::MAX]);
                    if let Some(v) = value {
                        if v[1] != u32::MAX {
                            println!("this edge has appeared more than twice: {:?}", key);
                        }
                        edges.insert(key, [idx1, v[0]]);
                    }
                }
                {
                    let key: Edge = new_edge(idx1, idx2);
                    let value = edges.insert(key, [idx0, u32::MAX]);
                    if let Some(v) = value {
                        if v[1] != u32::MAX {
                            println!("this edge has appeared more than twice: {:?}", key);
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

                idx += *verts as usize;
            }
        };
        // !SECTION - Fill vertices, tris, edges and boundary edges

        //
        // *
        // ! $Env:RUST_BACKTRACE=1
        // ?
        // TODO

        //// ANCHOR - Used to indicate a section in your file
        //// -TODO - An item that is awaiting completion
        //// FIXME - An item that requires a bugfix
        //// STUB - Used for generated default snippets
        //// NOTE - An important note for a specific code section
        //// REVIEW - An item that requires additional review
        //// LINK - Used to link to a file that can be opened within the editor (See 'Link Anchors'
        //// SECTION - Used to define a region (See 'Hierarchical anchors')
        //// ANCHOR this is in the section, section
        //// -!SECTION

        // Find the cut-nodes in p
        let mut overall_cluster_cut = HashSet::<Edge>::new();
        let mut clusters = Vec::<Cluster>::with_capacity(mesh.indices.len() / (3 * TRIS_IN_CLUSTER));

        let mut last_overall_cut = String::from("");
        let mut prev_overall_cut = String::from("");

        let mut last_cluster_cut = String::from("");
        let mut prev_cluster_cut = String::from("");

        // SECTION - Determine the clusters
        {
            let mut tris = tris.clone();
            let mut edges = edges.clone();

            let edge = edges.keys().next().unwrap();
            let mut seed = new_tri(edge[0], edge[1], edges[edge][0]);
            'cluster_generation: while !tris.is_empty() {
                if !tris.remove(&seed) {
                    println!("seed doesn't exist in 'tris'");
                    break 'cluster_generation;
                    // panic!("seed doesn't exist in 'tris'");
                }

                //? DEBUG
                prev_overall_cut = last_overall_cut.clone();
                last_overall_cut.clear();
                
                prev_cluster_cut = last_cluster_cut.clone();
                last_cluster_cut.clear();
                //? DEBUG

                let mut tri = seed;

                let anchor = mid_point_3(tri[0], tri[1], tri[2]);

                let mut cur_cluster_tris = Vec::with_capacity(TRIS_IN_CLUSTER);
                cur_cluster_tris.push(tri);
                let mut cur_cut = HashSet::new();

                remove_tri(&tri, &mut tris, &mut edges, &mut overall_cluster_cut, &mut cur_cut);

                for _ in 0..TRIS_IN_CLUSTER-1 {
                    let mut lowest: (f32, Triangle) = (f32::MAX, [0, 0, 0]);

                    // ----get the adjacent triangles to the current cluster edge path
                    // ----choose the triangle where the outlying vertex is closest to the anchor
                    'next_tri: for key in cur_cut.iter() {
                        if let Some(value) = edges.get(key) {
                            let dst0 = if value[0] != u32::MAX {
                                sqr_dist_w_custom(value[0], anchor)
                            } else {
                                0.0
                            };

                            let dst1 = if value[1] != u32::MAX {
                                sqr_dist_w_custom(value[1], anchor)
                            } else {
                                0.0
                            };

                            if dst0 == 0.0 && dst1 == 0.0 {
                                // println!("continuing cause dst are both 0.0");
                                continue 'next_tri;
                            }

                            if dst0 > dst1 {
                                if dst0 < lowest.0 {
                                    lowest = (dst0, [key[0], key[1], value[0]]);
                                }
                            } else if dst1 > dst0 {
                                if dst1 < lowest.0 {
                                    lowest = (dst1, [key[0], key[1], value[1]]);
                                }
                            } else {
                                panic!("equidistant adjacent triangles, this should not be possible (theoretically)");
                            }
                        }
                    }

                    if lowest.0 == f32::MAX {
                        println!("no suitable triangle to remove was found");
                        break;
                    }

                    tri = new_tri(lowest.1[0], lowest.1[1], lowest.1[2]);
                    
                    remove_tri(&tri, &mut tris, &mut edges, &mut overall_cluster_cut, &mut cur_cut);

                    cur_cluster_tris.push(tri);
                    tris.remove(&tri);

                    // ----remove the chosen triangle
                    // ----add the two new edges to the the current cluster edge path

                    // println!("{}", i);
                }

                // we have completed a cluster
                clusters.push(Cluster {
                    tris: cur_cluster_tris,
                    cut: cur_cut.into_iter().collect(),
                    anchor,
                });

                // remove all edges from cluster edge path
                for edge in &clusters.last().unwrap().cut {
                    if !overall_cluster_cut.insert(*edge) {
                        overall_cluster_cut.remove(edge);
                    }
                }

                let mut should_break = false;

                // calculate the next seed
                'seed_seeking: for edge in &overall_cluster_cut {
                    if let Some(verts) = edges.get(edge) {
                        let vert = if verts[0] == u32::MAX {
                            if verts[1] == u32::MAX {
                                // panic!("Both verts were u32::MAX. This should not be possible.");
                                println!("Both verts were u32::MAX. This should not be possible.");
                                break 'cluster_generation;
                            }
                            verts[1]
                        } else {
                            verts[0]
                        };

                        // check if any of the other edges of the triangle are in the cut
                        let edge0 = new_edge(edge[0], vert);
                        let edge1 = new_edge(edge[1], vert);

                        // choose one with only one edge remaining in the edges if it exists (ie has two edges in the overall cluster cut/edge_path)
                        seed = new_tri(edge[0], edge[1], vert);
                        if overall_cluster_cut.contains(&edge0)
                            || overall_cluster_cut.contains(&edge1)
                        {
                            // we found an ideal seed triangle to use
                            break 'seed_seeking;
                        }
                    } else {
                        println!("something seems to be weird here");
                        should_break = true;
                        break;
                    }
                    // the seed is ok, keep looking to be sure
                }

                // under debug flag we can check for orphaned triangles and throw an error (or not in debug and add them to an adjacent cluster migrating triangles to maintain constant cluster size)

                println!("Clusters made: {}", clusters.len());
                println!("Cluster generated with {} triangles", clusters.last().unwrap().tris.len());

                // ?DEBUG
                {
                    for edge in &clusters.last().unwrap().cut {
                        last_cluster_cut
                            .push_str(format!("l {} {}\n", edge[0] + 1, edge[1] + 1).as_str());
                    }

                    for edge in &overall_cluster_cut {
                        last_overall_cut
                            .push_str(format!("l {} {}\n", edge[0] + 1, edge[1] + 1).as_str());
                    }
                }

                if should_break {
                    break 'cluster_generation;
                }
            }
        }
        // !SECTION - Determine the clusters

        // SECTION - Generate debug.obj mesh
        println!("Generate debug.obj mesh");
        {
            let mut f = File::create("debug.obj").unwrap();
            f.write_all(b"mtllib bunny.mtl\no bun_zipper\n").unwrap();
            for v_idx in (0..mesh.positions.len()).step_by(3) {
                f.write_all(
                    format!(
                        "v {} {} {}\n",
                        mesh.positions[v_idx + 0],
                        mesh.positions[v_idx + 1],
                        mesh.positions[v_idx + 2]
                    )
                    .as_bytes(),
                )
                .unwrap();
            }
            f.write_all(b"usemtl None\ns off\n").unwrap();

            // for t in tris {
            //     f.write_all(format!("f {} {} {}\n", t[0] + 1, t[1] + 1, t[2] + 1).as_bytes())
            //         .unwrap();
            // }
            // for verts in cut_vertices {
            //     for v in verts.1 {
            //         f.write_all(format!("l {} {}\n", verts.0 + 1, v + 1).as_bytes())
            //             .unwrap();
            //     }
            // }
            // for e in edges {
            //     f.write_all(format!("l {} {}\n", e.0[0] + 1, e.0[1] + 1).as_bytes())
            //         .unwrap();
            // }
            // for p in pairs {
            //     f.write_all(
            //         format!("l {} {}\n", boundary_groups[p[0]][0] + 1, boundary_groups[p[1]][0] + 1).as_bytes(),
            //     )
            //     .unwrap();
            // }
            // {
            //     // for idx in 1..final_cut_path.len() {}
            //     let path = final_cut_path
            //         .iter()
            //         .map(|pt| format!("{}", pt + 1))
            //         .collect::<Vec<String>>()
            //         .join(" ");
            //     f.write_all(format!("l {}\n", path).as_bytes()).unwrap();
            // }

            // !Most recent vvv
            // for e in edge_path {
            //     f.write_all(format!("l {} {}\n", e[0] + 1, e[1] + 1).as_bytes())
            //         .unwrap();
            // }
            // for group in &boundary_groups {
            //     for idx in 1..group.len() {
            //         f.write_all(
            //             format!("l {} {}\n", group[idx - 1] + 1, group[idx] + 1).as_bytes(),
            //         )
            //         .unwrap();
            //     }
            // }

            f.write_all(format!("o prev_cluster_cut_1\n").as_bytes())
                .unwrap();
            f.write_all(prev_cluster_cut.as_bytes()).unwrap();

            f.write_all(format!("o last_cluster_cut_2\n").as_bytes())
                .unwrap();
            f.write_all(last_cluster_cut.as_bytes()).unwrap();

            f.write_all(format!("o prev_overall_cut_1\n").as_bytes())
                .unwrap();
            f.write_all(prev_overall_cut.as_bytes()).unwrap();

            f.write_all(format!("o last_overall_cut_2\n").as_bytes())
                .unwrap();
            f.write_all(last_overall_cut.as_bytes()).unwrap();

            for (idx, cluster) in clusters.iter().enumerate() {
                f.write_all(format!("o Cluster{}\n", idx + 1).as_bytes())
                    .unwrap();
                for tri in &cluster.tris {
                    f.write_all(
                        format!("f {} {} {}\n", tri[0] + 1, tri[1] + 1, tri[2] + 1).as_bytes(),
                    )
                    .unwrap();
                }
                f.write_all("\n".as_bytes()).unwrap();
            }

            f.flush().unwrap();
        };
        // !SECTION - Generate debug.obj mesh
    }
}

fn remove_tri(tri: &Triangle, tris: &mut HashSet<Triangle>, edges: &mut HashMap<Edge, [u32; 2]>, overall_cluster_cut: &mut HashSet<Edge>, cur_cut: &mut HashSet<Edge>) {
    tris.remove(tri);

    let edge0 = new_edge(tri[0], tri[1]);
    let edge1 = new_edge(tri[0], tri[2]);
    let edge2 = new_edge(tri[1], tri[2]);

    let mut remove0 = false;
    let mut remove1 = false;
    let mut remove2 = false;

    if let Some(value) = edges.get_mut(&edge0) {
        // if overall_cluster_cut.contains(&edge0) {
        //     remove0 = false;
        // } else 
        if value[0] == tri[2] {
            value[0] = u32::MAX;
            if value[1] == u32::MAX {
                remove0 = true;
            }
        } else if value[1] == tri[2] {
            value[1] = u32::MAX;
            if value[0] == u32::MAX {
                remove0 = true;
            }
        }
    }

    if let Some(value) = edges.get_mut(&edge1) {
        // if overall_cluster_cut.contains(&edge1) {
        //     remove1 = false;
        // } else 
        if value[0] == tri[1] {
            value[0] = u32::MAX;
            if value[1] == u32::MAX {
                remove1 = true;
            }
        } else if value[1] == tri[1] {
            value[1] = u32::MAX;
            if value[0] == u32::MAX {
                remove1 = true;
            }
        }
    }

    if let Some(value) = edges.get_mut(&edge2) {
        // if overall_cluster_cut.contains(&edge2) {
        //     remove2 = false;
        // } else 
        if value[0] == tri[0] {
            value[0] = u32::MAX;
            if value[1] == u32::MAX {
                remove2 = true;
            }
        } else if value[1] == tri[0] {
            value[1] = u32::MAX;
            if value[0] == u32::MAX {
                remove2 = true;
            }
        }
    }

    if remove0 {
        if !overall_cluster_cut.contains(&edge0) {
            cur_cut.remove(&edge0);
        } else {
            cur_cut.insert(edge0);
        }
        edges.remove(&edge0);
    } else {
        cur_cut.insert(edge0);
    }

    if remove1 {
        if !overall_cluster_cut.contains(&edge1) {
            cur_cut.remove(&edge1);
        } else {
            cur_cut.insert(edge1);
        }
        edges.remove(&edge1);
    } else {
        cur_cut.insert(edge1);
    }

    if remove2 {
        if !overall_cluster_cut.contains(&edge2) {
            cur_cut.remove(&edge2);
        } else {
            cur_cut.insert(edge2);
        }
        edges.remove(&edge2);
    } else {
        cur_cut.insert(edge2);
    }
}

// orders the vertices consistently
fn new_tri(a: u32, b: u32, c: u32) -> Triangle {
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
fn new_edge(a: u32, b: u32) -> Edge {
    if a <= b {
        [a, b]
    } else {
        [b, a]
    }
}

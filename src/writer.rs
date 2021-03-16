use core::panic;
use rand::Rng;
use std::{
    collections::{HashMap, HashSet},
    thread::current,
};
use std::{fs::File, io::prelude::*};

type Triangle = [u32; 3];
type Edge = [u32; 2];

fn translate() {}

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

        // create data structure that can be referenced to find which indices are included in which triangles
        let mut edges = HashMap::<Edge, Edge>::with_capacity(mesh.num_face_indices.len() * 24);
        let mut tris = HashSet::<Triangle>::with_capacity(mesh.num_face_indices.len());
        let mut vertices: HashMap<u32, Vec<u32>> = HashMap::new();

        println!("num_face_indices: {}", mesh.num_face_indices.len());

        let mut idx = 0;
        for verts in &mesh.num_face_indices {
            {
                let key: Edge = edge(mesh.indices[idx + 0], mesh.indices[idx + 1]);
                let value = edges.insert(key, [mesh.indices[idx + 2], u32::MAX]);
                if let Some(v) = value {
                    if v[1] != u32::MAX {
                        println!("this edge has appeared more than twice: {:?}", key);
                    }
                    edges.insert(key, [mesh.indices[idx + 2], v[0]]);
                }
            }
            {
                let key: Edge = edge(mesh.indices[idx + 0], mesh.indices[idx + 2]);
                let value = edges.insert(key, [mesh.indices[idx + 1], u32::MAX]);
                if let Some(v) = value {
                    if v[1] != u32::MAX {
                        println!("this edge has appeared more than twice: {:?}", key);
                    }
                    edges.insert(key, [mesh.indices[idx + 1], v[0]]);
                }
            }
            {
                let key: Edge = edge(mesh.indices[idx + 1], mesh.indices[idx + 2]);
                let value = edges.insert(key, [mesh.indices[idx + 0], u32::MAX]);
                if let Some(v) = value {
                    if v[1] != u32::MAX {
                        println!("this edge has appeared more than twice: {:?}", key);
                    }
                    edges.insert(key, [mesh.indices[idx + 0], v[0]]);
                }
            }

            {
                vertices
                    .entry(mesh.indices[idx + 0])
                    .or_insert(Vec::new())
                    .push(mesh.indices[idx + 1]);
                vertices
                    .entry(mesh.indices[idx + 1])
                    .or_insert(Vec::new())
                    .push(mesh.indices[idx + 0]);
                vertices
                    .entry(mesh.indices[idx + 0])
                    .or_insert(Vec::new())
                    .push(mesh.indices[idx + 2]);
                vertices
                    .entry(mesh.indices[idx + 2])
                    .or_insert(Vec::new())
                    .push(mesh.indices[idx + 0]);
                vertices
                    .entry(mesh.indices[idx + 1])
                    .or_insert(Vec::new())
                    .push(mesh.indices[idx + 2]);
                vertices
                    .entry(mesh.indices[idx + 2])
                    .or_insert(Vec::new())
                    .push(mesh.indices[idx + 1]);
            };

            tris.insert(tri(mesh.indices[idx + 0], mesh.indices[idx + 1], mesh.indices[idx + 2]));

            idx += *verts as usize;
        }

        let mut boundary_edges = Vec::<Edge>::new();
        for e in &edges {
            if e.1[1] == u32::MAX {
                boundary_edges.push(e.0.clone());
            }
        }

        let dist = |a: u32, b: u32| -> f32 {
            let p1 = &mesh.positions[(a * 3) as usize..(a * 3 + 3) as usize];
            let p2 = &mesh.positions[(b * 3) as usize..(b * 3 + 3) as usize];

            f32::sqrt((p1[0] - p2[0]).powf(2.0) + (p1[1] - p2[1]).powf(2.0) + (p1[2] - p2[2]).powf(2.0))
        };

        // if only a single vertex remains in p then add back two adjacent edges to p.
        //     this happens if the mesh is entirely closed with no pre-existing boundary edges:
        //     "For the case of a closed mesh of genus 0, the resulting ρ will consist of a single
        //      vertex, since it has no loops. Because our parametrization requires that we map ρ'
        //      onto a square, we add back to ρ two adjacent mesh edges."
        //
        // optimize the path to be smoother, less serrated, shortest path:
        //     "we straighten each cut-path in ρ by computing a constrained
        //      shortest path that connects its two adjacent cut-nodes and
        //      stays within a neighborhood of the original cut-path."
        //
        // A vertex v with valence k in ρ is replicated as k vertices in ρ'.
        // Vertices in ρ that have valence k != 2 in the cut are called cut-nodes.
        // (We still refer to these as cut-nodes when replicated in ρ'.)
        //
        // *
        // ! $Env:RUST_BACKTRACE=1
        // ?
        // TODO

        // ANCHOR - Used to indicate a section in your file
        // TODO - An item that is awaiting completion
        // FIXME - An item that requires a bugfix
        // STUB - Used for generated default snippets
        // NOTE - An important note for a specific code section
        // REVIEW - An item that requires additional review
        // LINK - Used to link to a file that can be opened within the editor (See 'Link Anchors'
        // SECTION - Used to define a region (See 'Hierarchical anchors')
        // ANCHOR this is in the section, section
        // !SECTION

        // Find the cut-nodes in p
        let mut edge_path = Vec::<Vec<u32>>::new();
        let mut pairs = Vec::<Edge>::new();
        let mut boundary_groups = Vec::<Vec<u32>>::new();
        println!("part 1");
        // TODO one of the edges is duplicated in both the boundary group and the edge_path, ensure in djikstras that that does not happen
        // TODO find a better closest point between groups 
        {
            let mut remaining_boundary_edges = boundary_edges.clone();
            while let Some(seed) = remaining_boundary_edges.pop() {
                let mut path = Vec::<u32>::new();
                // path.push(seed[0]);
                path.push(seed[1]);
                let mut current = seed[1];
                while current != seed[0] {
                    let mut rmv_idx = usize::MAX;
                    for (idx, e) in remaining_boundary_edges.iter().enumerate() {
                        if e[0] == current {
                            current = e[1];
                            rmv_idx = idx;
                            break;
                        } else if e[1] == current {
                            current = e[0];
                            rmv_idx = idx;
                            break;
                        }
                    }
                    path.push(current);
                    remaining_boundary_edges.remove(rmv_idx);
                }
                boundary_groups.push(path);
            }

            let mut cnsmbl_boundary_groups = boundary_groups.clone();

            let first_group = cnsmbl_boundary_groups.pop().unwrap();
            let mut crnt_group = first_group.clone();
            let mut l_dst;
            loop {
                l_dst = (f32::MAX, u32::MAX, u32::MAX, usize::MAX);
                for (i, group) in cnsmbl_boundary_groups.iter().enumerate() {
                    for j in &crnt_group {
                        for k in group {
                            let dst = dist(*j, *k);
                            if dst < l_dst.0 {
                                l_dst = (dst, *j, *k, i);
                            }
                        }
                    }
                }

                pairs.push([l_dst.1, l_dst.2]);
                crnt_group = cnsmbl_boundary_groups.remove(l_dst.3);

                if cnsmbl_boundary_groups.is_empty() {
                    break;
                }
            }

            {
                for j in &crnt_group {
                    for k in &first_group {
                        let dst = dist(*j, *k);
                        if dst < l_dst.0 {
                            l_dst = (dst, *j, *k, i);
                        }
                    }
                }

                pairs.pop();
                pairs.push([l_dst.1, l_dst.2]);
            };

            // find the shortest path between the two cut nodes
            // djikstras

            for pair in &pairs {
                let mut pq = Vec::<(Vec<u32>, f32)>::new();
                let mut traversed = HashSet::<u32>::new();
                pq.push((vec![pair[0]], f32::MAX));
                'djikstra: loop {
                    let crnt_vertex = pq.pop().unwrap();
                    for v in &vertices[crnt_vertex.0.last().unwrap()] {
                        if !traversed.contains(v) {
                            let dst = dist(*v, pair[1]);
                            let mut path = crnt_vertex.0.clone();
                            path.push(*v);
                            traversed.insert(*v);
                            pq.push((path, dst));
                            if *v == pair[1] {
                                break 'djikstra;
                            }
                        }
                    }

                    // add the shortest distances to edge_path
                    pq.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
                }
                let final_path = pq.pop().unwrap().0;
                edge_path.push(final_path);
            }
        };

        let mut final_cut_path = Vec::<u32>::new();
        println!("part 2");
        // unit circle thing
        {
            // println!("part 2.0.1");
            // let mut f = File::create("debug.txt").unwrap();
            // f.write_all(b"pairs: ").unwrap();
            // for pt in &pairs {
            //     f.write_all(format!("\n    pair: [{}, {}]", pt[0], pt[1]).as_bytes())
            //         .unwrap();
            // }
            // f.write_all(b"\nedge_path: ").unwrap();
            // for edge in &edge_path {
            //     f.write_all(b"\n    edge: ").unwrap();
            //     for pt in edge {
            //         f.write_all(format!("{}, ", pt).as_bytes()).unwrap();
            //     }
            // }
            // f.write_all(b"\nboundary_groups: ").unwrap();
            // for group in &boundary_groups {
            //     f.write_all(b"\n\n    group: ").unwrap();
            //     for pt in group {
            //         f.write_all(format!("{}, ", pt).as_bytes()).unwrap();
            //     }
            // }
            // f.flush().unwrap();

            // edge_path
            // pairs
            // boundary_groups
            //    => final_cut_path
            //
            let mut cnsmbl_pairs = pairs.clone();
            for pair in &pairs {
                cnsmbl_pairs.push([pair[1], pair[0]]);
            }
            let mut cnsmbl_edge_path = edge_path.clone();
            for path in &edge_path {
                let mut clone = path.clone();
                clone.reverse();
                cnsmbl_edge_path.push(clone);
            }
            // let mut cnsmbl_boundary_groups = boundary_groups.clone();
            // let mut final_cut_path = Vec::<u32>::new();
            // println!("part 2.1");
            {
                // * start at a random point (from pairs?)
                let mut pair = cnsmbl_pairs.pop().unwrap();
                let mut crnt_vert = pair[0];
                let end_point = pair[0];
                // println!("part 2.1.1, crnt_vert: {}", crnt_vert);
                'final_path_loop: loop {
                    // * proceed down edge_path until you hit the other pair
                    // println!("part 2.1.1: edge_path");
                    let mut hit = false;
                    crnt_vert = pair[0];

                    for idx in 0..cnsmbl_edge_path.len() {
                        if cnsmbl_edge_path[idx][0] == crnt_vert {
                            let mut path = cnsmbl_edge_path.remove(idx);
                            // println!(
                            //     "part 2.1.1: found path: {} => {}",
                            //     path.first().unwrap(),
                            //     path.last().unwrap()
                            // );
                            crnt_vert = path.pop().unwrap();
                            final_cut_path.append(&mut path);
                            hit = true;
                            break;
                        }
                    }
                    if !hit {
                        println!("AHHH, no edge path found");
                        println!("cnsmbl_edge_path: {:?}", cnsmbl_edge_path);
                        break;
                    }

                    // * proceed down the respective element of "boundary_groups"
                    // println!("part 2.1.1: boundary_groups");
                    let mut group_idx = usize::MAX;
                    'group_loop: for group in &boundary_groups {
                        for (idx, vert) in group.iter().enumerate() {
                            if *vert == crnt_vert {
                                group_idx = idx;
                                break;
                            }
                        }

                        if group_idx != usize::MAX {
                            // println!(
                            //     "part 2.1.1: found group: {} => {}",
                            //     group.first().unwrap(),
                            //     group.last().unwrap()
                            // );
                            loop {
                                final_cut_path.push(group[group_idx]);
                                group_idx += 1;
                                if group_idx == group.len() {
                                    group_idx = 0;
                                }

                                if cnsmbl_pairs.is_empty() {
                                    if group[group_idx] == end_point {
                                        // println!("crnt_vert: {}, end_point: {}", group[group_idx], end_point);
                                        break 'final_path_loop;
                                    }
                                } else {
                                    for (idx, p) in cnsmbl_pairs.iter().enumerate() {
                                        if p[0] == group[group_idx] {
                                            // println!(
                                            //     "part 2.1.1: group[group_idx]: {}, found new pair: [{}, {}]",
                                            //     group[group_idx], p[0], p[1]
                                            // );
                                            pair = cnsmbl_pairs.remove(idx);
                                            break 'group_loop;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            };

            // println!("part 2.2");
            // let mut f = File::create("debug2.txt").unwrap();
            // f.write_all(b"path: ").unwrap();
            // for pt in &final_cut_path {
            //     f.write_all(format!("{}, ", pt).as_bytes()).unwrap();
            // }
            // f.flush().unwrap();

            // if you have triangles with all vertices on one edge, that is not ok so you split it like this:
            //    /\
            //   /__\
            //
            //    /|\
            //   /_|_\
            //
            // and split adjacent triangles in half so the vertex is handled
            //
            // split any edges that span over a corner of the geometry image and split their mirror/mate edge as well
            //
            //     "Finally, we find that placing a valence-1 cut-node at a corner of D results in poor geometric
            //      behavior, so if this occurs we rotate the boundary parametrization."
            // not sure what that means
        };

        println!("part 3");
        // DEBUG
        if !boundary_edges.is_empty() {
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
            // for e in edge_path {
            //     f.write_all(format!("l {} {}\n", e[0] + 1, e[1] + 1).as_bytes())
            //         .unwrap();
            // }//
            {
                // for idx in 1..final_cut_path.len() {}
                let path = final_cut_path
                    .iter()
                    .map(|pt| format!("{}", pt + 1))
                    .collect::<Vec<String>>()
                    .join(" ");
                f.write_all(format!("l {}\n", path).as_bytes()).unwrap();
            }
            f.flush().unwrap();
        };
    }
}

// orders the vertices consistently
fn tri(a: u32, b: u32, c: u32) -> Triangle {
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
fn edge(a: u32, b: u32) -> Edge {
    if a <= b {
        [a, b]
    } else {
        [b, a]
    }
}

// c++
/* void recreateMesh()
{
    //GIMxInd and GIMyInd are indices into the image 2d array of 256x256
    //myimage is the GIM
    float dist1 = sqrt((float)((myimage[GIMxInd][GIMyInd].r - myimage[GIMxInd + 1][GIMyInd + 1].r)*(myimage[GIMxInd][GIMyInd].r - myimage[GIMxInd + 1][GIMyInd + 1].r)
                     + (myimage[GIMxInd][GIMyInd].g - myimage[GIMxInd + 1][GIMyInd + 1].g)*(myimage[GIMxInd][GIMyInd].g - myimage[GIMxInd + 1][GIMyInd + 1].g)
                     + (myimage[GIMxInd][GIMyInd].b - myimage[GIMxInd + 1][GIMyInd + 1].b)*(myimage[GIMxInd][GIMyInd].b - myimage[GIMxInd + 1][GIMyInd + 1].b)));

    float dist2 = sqrt((float)((myimage[GIMxInd + 1][GIMyInd].r - myimage[GIMxInd][GIMyInd + 1].r)*(myimage[GIMxInd + 1][GIMyInd].r - myimage[GIMxInd][GIMyInd + 1].r)
                     + (myimage[GIMxInd + 1][GIMyInd].g - myimage[GIMxInd][GIMyInd + 1].g)*(myimage[GIMxInd + 1][GIMyInd].g - myimage[GIMxInd][GIMyInd + 1].g)
                     + (myimage[GIMxInd + 1][GIMyInd].b - myimage[GIMxInd][GIMyInd + 1].b)*(myimage[GIMxInd + 1][GIMyInd].b - myimage[GIMxInd][GIMyInd + 1].b)));

    // 0 (size - 3)
    newVertices.push_back(new Vector3(myimage[GIMxInd][GIMyInd].r, myimage[GIMxInd][GIMyInd].g, myimage[GIMxInd][GIMyInd].b));
    newVPoints.push_back(new Point3(Vector3(myimage[GIMxInd][GIMyInd].r, myimage[GIMxInd][GIMyInd].g, myimage[GIMxInd][GIMyInd].b)));
    newVPoints.at(newVPoints.size() - 1)->normal.x = mynormalimage[GIMxInd][GIMyInd].r;
    newVPoints.at(newVPoints.size() - 1)->normal.y = mynormalimage[GIMxInd][GIMyInd].g;
    newVPoints.at(newVPoints.size() - 1)->normal.z = mynormalimage[GIMxInd][GIMyInd].b;

    // 1 (size - 2)
    newVertices.push_back(new Vector3(myimage[GIMxInd+1][GIMyInd].r, myimage[GIMxInd+1][GIMyInd].g, myimage[GIMxInd+1][GIMyInd].b));
    newVPoints.push_back(new Point3(Vector3(myimage[GIMxInd+1][GIMyInd].r, myimage[GIMxInd+1][GIMyInd].g, myimage[GIMxInd+1][GIMyInd].b)));
    newVPoints.at(newVPoints.size() - 1)->normal.x = mynormalimage[GIMxInd+1][GIMyInd].r;
    newVPoints.at(newVPoints.size() - 1)->normal.y = mynormalimage[GIMxInd+1][GIMyInd].g;
    newVPoints.at(newVPoints.size() - 1)->normal.z = mynormalimage[GIMxInd+1][GIMyInd].b;

    // 2 (size - 1)
    newVertices.push_back(new Vector3(myimage[GIMxInd][GIMyInd+1].r, myimage[GIMxInd][GIMyInd+1].g, myimage[GIMxInd][GIMyInd+1].b));
    newVPoints.push_back(new Point3(Vector3(myimage[GIMxInd][GIMyInd+1].r, myimage[GIMxInd][GIMyInd+1].g, myimage[GIMxInd][GIMyInd+1].b)));
    newVPoints.at(newVPoints.size() - 1)->normal.x = mynormalimage[GIMxInd][GIMyInd+1].r;
    newVPoints.at(newVPoints.size() - 1)->normal.y = mynormalimage[GIMxInd][GIMyInd+1].g;
    newVPoints.at(newVPoints.size() - 1)->normal.z = mynormalimage[GIMxInd][GIMyInd+1].b;

    // 3 (size)
    newVertices.push_back(new Vector3(myimage[GIMxInd+1][GIMyInd+1].r, myimage[GIMxInd+1][GIMyInd+1].g, myimage[GIMxInd+1][GIMyInd+1].b));
    newVPoints.push_back(new Point3(Vector3(myimage[GIMxInd+1][GIMyInd+1].r, myimage[GIMxInd+1][GIMyInd+1].g, myimage[GIMxInd+1][GIMyInd+1].b)));
    newVPoints.at(newVPoints.size() - 1)->normal.x = mynormalimage[GIMxInd+1][GIMyInd+1].r;
    newVPoints.at(newVPoints.size() - 1)->normal.y = mynormalimage[GIMxInd+1][GIMyInd+1].g;
    newVPoints.at(newVPoints.size() - 1)->normal.z = mynormalimage[GIMxInd+1][GIMyInd+1].b;

    if(dist1 <= dist2)
    {
        // 0, 1, 3
        newTriangles.push_back(new Tri(newVertices.size() - 3, newVertices.size() - 2, newVertices.size()));
        // 0, 2, 3
        newTriangles.push_back(new Tri(newVertices.size() - 3, newVertices.size() - 1, newVertices.size()));
        pixelToTri[GIMxInd][GIMyInd].push_back(newTriangles.size() - 2);
        pixelToTri[GIMxInd][GIMyInd].push_back(newTriangles.size() - 1);
        pixelToTri[GIMxInd+1][GIMyInd].push_back(newTriangles.size() - 2);
        pixelToTri[GIMxInd][GIMyInd+1].push_back(newTriangles.size() - 1);
        pixelToTri[GIMxInd+1][GIMyInd+1].push_back(newTriangles.size() - 2);
        pixelToTri[GIMxInd+1][GIMyInd+1].push_back(newTriangles.size() - 1);
    }
    else
    {
        // 1, 0, 2
        newTriangles.push_back(new Tri(newVertices.size() - 2, newVertices.size() - 3, newVertices.size() - 1));
        // 1, 3, 2
        newTriangles.push_back(new Tri(newVertices.size() - 2, newVertices.size(), newVertices.size() - 1));
        pixelToTri[GIMxInd][GIMyInd].push_back(newTriangles.size() - 2);
        pixelToTri[GIMxInd+1][GIMyInd].push_back(newTriangles.size() - 1);
        pixelToTri[GIMxInd+1][GIMyInd].push_back(newTriangles.size() - 2);
        pixelToTri[GIMxInd][GIMyInd+1].push_back(newTriangles.size() - 1);
        pixelToTri[GIMxInd][GIMyInd+1].push_back(newTriangles.size() - 2);
        pixelToTri[GIMxInd+1][GIMyInd+1].push_back(newTriangles.size() - 1);
    }

    myimage[GIMxInd][GIMyInd].r = myimage[GIMxInd][GIMyInd].g = myimage[GIMxInd][GIMyInd].b = 1.0;
    mynormalimage[GIMxInd][GIMyInd].r = mynormalimage[GIMxInd][GIMyInd].g = mynormalimage[GIMxInd][GIMyInd].b = 1.0;

    if(GIMyInd == 63)
    {
        myimage[GIMxInd][GIMyInd+1].r = myimage[GIMxInd][GIMyInd+1].g = myimage[GIMxInd][GIMyInd+1].b = 1.0;
        mynormalimage[GIMxInd][GIMyInd+1].r = mynormalimage[GIMxInd][GIMyInd+1].g = mynormalimage[GIMxInd][GIMyInd+1].b = 1.0;
    }

    if(GIMxInd < 63)
    {
        GIMxInd++;
    }
    else
    {
        myimage[GIMxInd+1][GIMyInd].r = myimage[GIMxInd+1][GIMyInd].g = myimage[GIMxInd+1][GIMyInd].b = 1.0;
        mynormalimage[GIMxInd+1][GIMyInd].r = mynormalimage[GIMxInd+1][GIMyInd].g = mynormalimage[GIMxInd+1][GIMyInd].b = 1.0;

        GIMxInd = 0;
        GIMyInd++;
        if(GIMyInd > 63)
        {
            myimage[64][64].r = myimage[64][64].g = myimage[64][64].b = 1.0;
            mynormalimage[64][64].r = mynormalimage[64][64].g = mynormalimage[64][64].b = 1.0;

            recreatedMesh = true;
            printf("Finished reconstruction!\n");

            // Build the quad tree
            quadTree = createQuadTree(65, 0, 0, myimage2);
            printf("Built quad tree\n");
        }
    }
    //calcAllNewVertexNormals();
}
*/

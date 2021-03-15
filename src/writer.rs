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

        // println!("part 1.0, cut edges and triangles: begins");
        // {
        //     let tri_idx = rand::thread_rng().gen_range(0..mesh.num_face_indices.len());
        //     let mut triangle = tri(
        //         mesh.indices[tri_idx * 3 + 0],
        //         mesh.indices[tri_idx * 3 + 1],
        //         mesh.indices[tri_idx * 3 + 2],
        //     );
        //
        //     if !tris.contains(&triangle) {
        //         println!("random seed triangle not in the mesh");
        //     }
        //
        //     let mut inst_tris = Vec::<(Triangle, f32, Edge)>::with_capacity(3);
        //
        //     println!("{} total edges before removal", edges.len());
        //
        //     let mut removable_edges = Vec::<Edge>::new();
        //
        //     'part_one: loop {
        //         if !tris.remove(&triangle) {
        //             println!("this is not possible, loop triangle not in the mesh");
        //         }
        //
        //         // find the adjacent tris
        //         {
        //             let edg = edge(triangle[0], triangle[1]);
        //             if let Some(value) = edges.get(&edg) {
        //                 let idx = if value[0] == triangle[2] { 1 } else { 0 };
        //                 let push_trig = tri(triangle[0], triangle[1], value[idx]);
        //                 if tris.contains(&push_trig) {
        //                     // this distance is the distance between the points not on the edge itself
        //                     inst_tris.push((push_trig, dist(triangle[2], value[idx]), edg));
        //                 }
        //             }
        //         };
        //
        //         {
        //             let edg = edge(triangle[0], triangle[2]);
        //             if let Some(value) = edges.get(&edg) {
        //                 let idx = if value[0] == triangle[1] { 1 } else { 0 };
        //                 let push_trig = tri(triangle[0], triangle[2], value[idx]);
        //                 if tris.contains(&push_trig) {
        //                     inst_tris.push((push_trig, dist(triangle[1], value[idx]), edg));
        //                 }
        //             }
        //         };
        //
        //         {
        //             let edg = edge(triangle[1], triangle[2]);
        //             if let Some(value) = edges.get(&edg) {
        //                 let idx = if value[0] == triangle[0] { 1 } else { 0 };
        //                 let push_trig = tri(triangle[1], triangle[2], value[idx]);
        //                 if tris.contains(&push_trig) {
        //                     inst_tris.push((push_trig, dist(triangle[0], value[idx]), edg));
        //                 }
        //             }
        //         };
        //
        //         if inst_tris.len() == 0 {
        //             'stack: while let Some(edge) = removable_edges.pop() {
        //                 // remove edge and the triangle adjacent to it
        //
        //                 // check if there is another triangle attached to this edge
        //                 // if the edge is not in the simplicial complex, it has already been removed and we can move to the next on the stack
        //                 if let Some(value) = edges.get(&edge) {
        //                     let idx = if triangle[0] != value[0] && triangle[1] != value[0] && triangle[2] != value[0] {
        //                         0
        //                     } else if triangle[0] != value[1] && triangle[1] != value[1] && triangle[2] != value[1] {
        //                         1
        //                     } else {
        //                         panic!("this is not possible");
        //                     };
        //                     let push_trig = tri(edge[0], edge[1], value[idx]);
        //                     if tris.contains(&push_trig) {
        //                         // if there is: remove edge and set the triangle for the next loop iteration
        //                         edges.remove(&edge);
        //                         triangle = push_trig;
        //                         continue 'part_one;
        //                     }
        //                 }
        //             }
        //
        //             {
        //                 println!("part 1.0, cut edges and triangles: complete");
        //                 break 'part_one;
        //             }
        //         } else {
        //             // three possible edges to remove from
        //             // sort by distance to the current tri and choose the closest in the simplicial complex
        //             inst_tris.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        //
        //             // remove the corresponding edge
        //             edges.remove(&inst_tris[0].2);
        //
        //             // and set the triangle for the next loop iteration
        //             triangle = inst_tris[0].0;
        //
        //             if inst_tris.len() == 2 {
        //                 removable_edges.push(inst_tris[1].2);
        //             } else if inst_tris.len() == 3 {
        //                 removable_edges.push(inst_tris[1].2);
        //                 removable_edges.push(inst_tris[2].2);
        //             }
        //         }
        //
        //         // empty the temp to prep for next tri
        //         inst_tris.clear();
        //     }
        // }

        // Intermediate step to clean up remaining triangles
        // println!("part 1.5, cut triangles: begins");
        // {
        //     let mut tri_vec: Vec<Triangle> = tris.iter().map(|t| t.clone()).collect();
        //     let mut inst_tris = Vec::<(Triangle, f32, Edge)>::with_capacity(3);
        //     let mut temp = Vec::<Triangle>::with_capacity(tri_vec.len());
        //     'part_one_point_five: loop {
        //         tri_vec.append(&mut temp);
        //         temp.clear();
        //         while let Some(triangle) = tri_vec.pop() {
        //             {
        //                 let edg = edge(triangle[0], triangle[1]);
        //                 if let Some(value) = edges.get(&edg) {
        //                     let idx = if value[0] == triangle[2] { 1 } else { 0 };
        //                     let push_trig = tri(triangle[0], triangle[1], value[idx]);
        //                     if value[idx] != u32::MAX && !tris.contains(&push_trig) {
        //                         // this distance is the distance between the points not on the edge itself
        //                         inst_tris.push((push_trig, dist(triangle[2], value[idx]), edg));
        //                     }
        //                 }
        //             };
        //
        //             {
        //                 let edg = edge(triangle[0], triangle[2]);
        //                 if let Some(value) = edges.get(&edg) {
        //                     let idx = if value[0] == triangle[1] { 1 } else { 0 };
        //                     let push_trig = tri(triangle[0], triangle[2], value[idx]);
        //                     if value[idx] != u32::MAX && !tris.contains(&push_trig) {
        //                         inst_tris.push((push_trig, dist(triangle[1], value[idx]), edg));
        //                     }
        //                 }
        //             };
        //
        //             {
        //                 let edg = edge(triangle[1], triangle[2]);
        //                 if let Some(value) = edges.get(&edg) {
        //                     let idx = if value[0] == triangle[0] { 1 } else { 0 };
        //                     let push_trig = tri(triangle[1], triangle[2], value[idx]);
        //                     if value[idx] != u32::MAX && !tris.contains(&push_trig) {
        //                         inst_tris.push((push_trig, dist(triangle[0], value[idx]), edg));
        //                     }
        //                 }
        //             };
        //
        //             if inst_tris.len() != 0 {
        //                 // three possible edges to remove from
        //                 // sort by distance to the current tri and choose the closest in the simplicial complex
        //                 inst_tris.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        //
        //                 // remove the corresponding edge
        //                 edges.remove(&inst_tris[0].2);
        //
        //             // tris.remove(triangle);
        //             } else {
        //                 temp.push(triangle);
        //             }
        //
        //             // empty the temp to prep for next tri
        //             inst_tris.clear();
        //
        //             tris.remove(&triangle);
        //         }
        //
        //         if temp.is_empty() {
        //             println!("part 1.5, cut triangles: complete");
        //             break 'part_one_point_five;
        //         }
        //     }
        // }

        // let mut cut_vertices: HashMap<u32, Vec<u32>> = HashMap::new();
        // for e in &edges {
        //     cut_vertices.entry(e.0[0]).or_insert(Vec::new()).push(e.0[1]);
        //     cut_vertices.entry(e.0[1]).or_insert(Vec::new()).push(e.0[0]);
        // }

        // second half of the algorithm, simplifying p
        // println!("part 2.0, cut vertices: begins");
        // 'part_two: loop {
        //     let mut removed = false;
        //     for e in &edges {
        //         if cut_vertices.contains_key(&e.0[0]) && cut_vertices.contains_key(&e.0[1]) {
        //             if cut_vertices.get(&e.0[0]).unwrap().len() == 1 {
        //                 let value = cut_vertices.remove(&e.0[0]).unwrap();
        //                 cut_vertices.entry(value[0]).and_modify(|v| {
        //                     v.remove(v.iter().position(|&n| n == e.0[0]).unwrap());
        //                 });
        //
        //                 removed = true;
        //             }
        //
        //             if cut_vertices.get(&e.0[1]).unwrap().len() == 1 {
        //                 let value = cut_vertices.remove(&e.0[1]).unwrap();
        //                 cut_vertices.entry(value[0]).and_modify(|v| {
        //                     v.remove(v.iter().position(|&n| n == e.0[1]).unwrap());
        //                 });
        //
        //                 removed = true;
        //             }
        //         }
        //     }
        //
        //     if !removed {
        //         println!("part 2.0, cut vertices: complete");
        //         break 'part_two;
        //     }
        // }

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
        // A vertex v with valence k in ρ is replicated as k vertices in ρ'. Vertices in ρ that have valence k != 2 in the cut are called cut-nodes. (We still refer to these as cut-nodes when replicated in ρ'.)
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
        // SECTION - Used to define a region (See 'Hierarchical anchors')
        // LINK - Used to link to a file that can be opened within the editor (See 'Link Anchors')

        // Find the cut-nodes in p
        let final_cut_path = Vec::<u32>::new();
        let mut edge_path = Vec::<Edge>::new();
        let mut pairs = Vec::<Edge>::new();
        let mut boundary_groups = Vec::<Vec<u32>>::new();
        {
            let mut remaining_boundary_edges = boundary_edges.clone();
            while let Some(seed) = remaining_boundary_edges.pop() {
                let mut path = Vec::<u32>::new();
                path.push(seed[0]);
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

            let mut consumable_boundary_groups = boundary_groups.clone();

            let first_group = consumable_boundary_groups.pop().unwrap();
            let mut current_group = first_group.clone();
            let mut l_dst;
            loop {
                l_dst = (f32::MAX, u32::MAX, u32::MAX, usize::MAX);
                for (i, group) in consumable_boundary_groups.iter().enumerate() {
                    for j in &current_group {
                        for k in group {
                            let dst = dist(*j, *k);
                            if dst < l_dst.0 {
                                l_dst = (dst, *j, *k, i);
                            }
                        }
                    }
                }

                pairs.push([l_dst.1, l_dst.2]);
                current_group = consumable_boundary_groups.remove(l_dst.3);

                if consumable_boundary_groups.is_empty() {
                    break;
                }
            }

            {
                for j in &current_group {
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

            let mut f = File::create("debug.txt").unwrap();
            for pair in pairs {
                f.write_all(format!("starting on pair: {} => {}\n", pair[0], pair[1]).as_bytes())
                    .unwrap();
                let mut pq = Vec::<(Vec<u32>, f32)>::new();
                let mut traversed = HashSet::<u32>::new();
                pq.push((vec![pair[0]], f32::MAX));
                'djikstra: loop {
                    let current_vertex = pq.pop().unwrap();
                    f.write_all(b"dst: ").unwrap();
                    for v in &vertices[current_vertex.0.last().unwrap()] {
                        if !traversed.contains(v) {
                            let dst = dist(*v, pair[1]);
                            let mut path = current_vertex.0.clone();
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

                    f.write_all(format!("shortest: {} \n", pq.last().unwrap().1).as_bytes())
                        .unwrap();
                    f.flush().unwrap();
                }
                println!("done a node");
                let final_path = pq.pop().unwrap();
                for i in 1..final_path.0.len() {
                    edge_path.push([final_path.0[i - 1], final_path.0[i]]);
                }
            }
            f.flush().unwrap();
        };

        // unit circle thing
        {
            // edge_path
            // pairs
            // boundary_groups
            //    => final_cut_path

            // start at a random point (from pairs?)
            // proceed down edge_path until you hit the other pair
            // proceed down the respective element of "boundary_groups"
        };

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
            for e in edge_path {
                f.write_all(format!("l {} {}\n", e[0] + 1, e[1] + 1).as_bytes())
                    .unwrap();
            }
            f.flush().unwrap();
        };

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
        // "Finally, we find that placing a valence-1 cut-node at a corner of D results in poor geometric behavior, so if this occurs we rotate the boundary parametrization."
        // not sure what that means

        // edges is what is the cut p
        {
            // add edges to a vec with their length
            // map the edges to a ?square? boundary
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

// fn cut() {
// function Cut(mesh M)
//    Remove seed triangle.
//    while there remains an edge e adjacent to only one triangle t
//    Remove e and t.
//    while there remains a vertex v adjacent to only one edge e
//    Remove v and e.
//    Cut p := remaining edges and vertices.
//    if only a single vertex remains in p then
//    Add back two adjacent edges to p.
// }

/* c++

void removeSeedTriangle()
{
    int randTri = ((float)rand())/RAND_MAX * Triangles.size();
    int i = 0, j;
    vector<Tri *>::iterator triItr;

    for(triItr = Triangles.begin(); triItr != Triangles.end(); triItr++)
    {
        if(i == randTri)
        {
            for(j = 0; j < cutPathEdges.size(); j++)
            {
                if((cutPathEdges[j]->v1 == Triangles[i]->v1
                    && cutPathEdges[j]->v2 == Triangles[i]->v2)
                    || (cutPathEdges[j]->v1 == Triangles[i]->v2
                    && cutPathEdges[j]->v2 == Triangles[i]->v1))
                {
                    cutPathEdges[j]->color.x = 0;
                    cutPathEdges[j]->color.y = 1;
                }

                if((cutPathEdges[j]->v1 == Triangles[i]->v2
                    && cutPathEdges[j]->v2 == Triangles[i]->v3)
                    || (cutPathEdges[j]->v1 == Triangles[i]->v3
                    && cutPathEdges[j]->v2 == Triangles[i]->v2))
                {
                    cutPathEdges[j]->color.x = 0;
                    cutPathEdges[j]->color.y = 1;
                }

                if((cutPathEdges[j]->v1 == Triangles[i]->v1
                    && cutPathEdges[j]->v2 == Triangles[i]->v3)
                    || (cutPathEdges[j]->v1 == Triangles[i]->v3
                    && cutPathEdges[j]->v2 == Triangles[i]->v1))
                {
                    cutPathEdges[j]->color.x = 0;
                    cutPathEdges[j]->color.y = 1;
                }
            }
            Triangles.erase(triItr);
            break;
        }
        i++;
    }
}

void createInitialCutPart1()
{
    bool inMoreThanOneTriangle = false;
    bool foundEdge = false;
    int triInd, edgeInd;

    // IF there remains an edge e adjacent to only one triangle t
    //     remove e and t
    for(int i = 0; i < cutPathEdges.size(); i++)
    {
        foundEdge = false;
        for(int j = 0; j < Triangles.size(); j++)
        {
            if((cutPathEdges[i]->v1 == Triangles[j]->v1
                && cutPathEdges[i]->v2 == Triangles[j]->v2)
                || (cutPathEdges[i]->v1 == Triangles[j]->v2
                && cutPathEdges[i]->v2 == Triangles[j]->v1))
            {
                if(!foundEdge)
                {
                    triInd = j;
                    foundEdge = true;
                }
                else
                {
                    inMoreThanOneTriangle = true;
                    break;
                }
            }
            else if((cutPathEdges[i]->v1 == Triangles[j]->v2
                && cutPathEdges[i]->v2 == Triangles[j]->v3)
                || (cutPathEdges[i]->v1 == Triangles[j]->v3
                && cutPathEdges[i]->v2 == Triangles[j]->v2))
            {
                if(!foundEdge)
                {
                    triInd = j;
                    foundEdge = true;
                }
                else
                {
                    inMoreThanOneTriangle = true;
                    break;
                }
            }
            else if((cutPathEdges[i]->v1 == Triangles[j]->v1
                && cutPathEdges[i]->v2 == Triangles[j]->v3)
                || (cutPathEdges[i]->v1 == Triangles[j]->v3
                && cutPathEdges[i]->v2 == Triangles[j]->v1))
            {
                if(!foundEdge)
                {
                    triInd = j;
                    foundEdge = true;
                }
                else
                {
                    inMoreThanOneTriangle = true;
                    break;
                }
            }
            else
            {
                // Do nothing
            }
        }

        // IF one edge found
        if(foundEdge && !inMoreThanOneTriangle)
        {
            edgeInd = i;
            break;
        }
        else
        {
            inMoreThanOneTriangle = false;
        }
    }

    // IF found one edge, remove that edge and triangle
    if(foundEdge && !inMoreThanOneTriangle)
    {
        vector<Tri *>::iterator triItr;
        vector<Edge *>::iterator edgeItr;
        int edgeNum = 0, triNum = 0;

        for(triItr = Triangles.begin(); triItr != Triangles.end(); triItr++)
        {
            // IF this is the triangle to remove
            if(triInd == triNum)
            {
                // Search for each pair of vertices in cutPathEdges
                for(int i = 0; i < cutPathEdges.size(); i++)
                {
                    if((cutPathEdges[i]->v1 == Triangles[triInd]->v1
                        && cutPathEdges[i]->v2 == Triangles[triInd]->v2)
                        || (cutPathEdges[i]->v1 == Triangles[triInd]->v2
                        && cutPathEdges[i]->v2 == Triangles[triInd]->v1))
                    {
                        cutPathEdges[i]->color.x = 0;
                        cutPathEdges[i]->color.y = 1;
                    }
                    else if((cutPathEdges[i]->v1 == Triangles[triInd]->v2
                        && cutPathEdges[i]->v2 == Triangles[triInd]->v3)
                        || (cutPathEdges[i]->v1 == Triangles[triInd]->v3
                        && cutPathEdges[i]->v2 == Triangles[triInd]->v2))
                    {
                        cutPathEdges[i]->color.x = 0;
                        cutPathEdges[i]->color.y = 1;
                    }
                    else if((cutPathEdges[i]->v1 == Triangles[triInd]->v1
                        && cutPathEdges[i]->v2 == Triangles[triInd]->v3)
                        || (cutPathEdges[i]->v1 == Triangles[triInd]->v3
                        && cutPathEdges[i]->v2 == Triangles[triInd]->v1))
                    {
                        cutPathEdges[i]->color.x = 0;
                        cutPathEdges[i]->color.y = 1;
                    }
                    else
                    {
                        // Do nothing
                    }
                }
                Triangles.erase(triItr);
                break;
            }
            triNum++;
        }

        for(edgeItr = cutPathEdges.begin(); edgeItr != cutPathEdges.end(); edgeItr++)
        {
            // IF this is the edge to remove
            if(edgeInd == edgeNum)
            {
                cutPathEdges.erase(edgeItr);
                break;
            }
            edgeNum++;
        }
    }
    else
    {
        removedEandT = true;
    }
}

void createInitialCutPart2()
{
    // IF there remains a vertex v adjacent to only one edge e
    //     remove v and e
    bool inMoreThanOneEdge = false;
    bool foundVertex = false;
    int vertInd, edgeInd;
    int count = 0;

    // IF there remains a vertex v adjacent to only one edge e
    //     remove v and e
    for(int i = 0; i < Vertices.size(); i++)
    {
        foundVertex = false;
        for(int j = 0; j < cutPathEdges.size(); j++)
        {
            if((i == cutPathEdges[j]->v1) || (i == cutPathEdges[j]->v2))
            {
                if(!foundVertex)
                {
                    edgeInd = j;
                    foundVertex = true;
                }
                else
                {
                    inMoreThanOneEdge = true;
                    break;
                }
            }
            else
            {
                // Do nothing
            }
        }

        // IF one vertex found
        if(foundVertex && !inMoreThanOneEdge)
        {
            vertInd = i;
            break;
        }
        else
        {
            inMoreThanOneEdge = false;
        }
    }

    // IF found one edge, remove that edge and triangle
    if(foundVertex && !inMoreThanOneEdge)
    {
        vector<Vector3 *>::iterator vertItr;
        vector<Edge *>::iterator edgeItr;
        int edgeNum = 0, vertNum = 0;

        if(cutPathEdges.size() <= 2)
        {
            removedVandE = true;

            for(int i = 0; i < cutPathEdges.size(); i++)
            {
                cutPathEdges[i]->color.x = 0;
                cutPathEdges[i]->color.y = 1;
                cutPathEdges[i]->color.z = 0;
            }

            return;
        }

        for(edgeItr = cutPathEdges.begin(); edgeItr != cutPathEdges.end(); edgeItr++)
        {
            // IF this is the edge to remove
            if(edgeInd == edgeNum)
            {
                cutPathEdges.erase(edgeItr);
                break;
            }
            edgeNum++;
        }
    }
    else
    {
        removedVandE = true;
    }


    for(int i = 0; i < cutPathEdges.size(); i++)
    {
        bool v1flag = false;
        bool v2flag = false;

        for(int j = 0; j < cutPathEdges.size(); j++)
        {
            if(i != j)
            {
                if(cutPathEdges[i]->v1 == cutPathEdges[j]->v1)
                {
                    v1flag = true;
                }
                else if(cutPathEdges[i]->v1 == cutPathEdges[j]->v2)
                {
                    v1flag = true;
                }
                else if(cutPathEdges[i]->v2 == cutPathEdges[j]->v1)
                {
                    v2flag = true;
                }
                else if(cutPathEdges[i]->v2 == cutPathEdges[j]->v2)
                {
                    v2flag = true;
                }
                else
                {
                    // Do nothing
                }
            }
        }

        if(v1flag && v2flag)
        {

        }
        else
        {
            cutPathEdges[i]->color.x = 1;
            cutPathEdges[i]->color.y = .6;
            cutPathEdges[i]->color.z = .6;
        }
    }
}

void recreateMesh()
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

/* tri-force
   /\
 \/__\/
 /\  /\
/__\/__\
    |
*/

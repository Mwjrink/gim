// create clusters of 128 triangles
// choose a seed triangle, choose the next vertex corresponding to an edge that is closest to the seed triangle

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::io::{BufWriter, Write};
use std::{fs::File, io::prelude::*};

use rand::Rng;

use crate::cluster::*;
use crate::ctree::*;
use crate::mesh::*;

const TRIS_IN_CLUSTER: usize = 128;

pub fn write(mesh: &Mesh) -> CTree {
    // create data structure that can be referenced to find which indices are included in which triangles

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

    let (shared_edges, mut ctree) = split_mesh(&mesh, TRIS_IN_CLUSTER);

    let write_file = File::create("./logs/log").unwrap();
    let mut writer = BufWriter::new(&write_file);
    write!(writer, "shared_edges: \n{:?} \n", shared_edges).unwrap();

    println!("{} clusters generated", ctree.len());

    // panic!();

    // TODO rewrite this, use a different algorithm to find the best combination of clusters to combine, take into account the layer/lod level
    // !SECTION - Combine the clusters into a graph
    {
        let mut clusters: Vec<u32> = (0..ctree.len() as u32).collect();
        let mut shared_edges: Vec<SharedEdges> = shared_edges
            .into_iter()
            .map(|v| {
                // let mut connections = v.1;
                // connections.sort_by(|a, b| b.1.cmp(&a.1));
                SharedEdges {
                    id: v.0,
                    lod: 0,
                    connections: v.1,
                }
            })
            .collect();

        // find 4 clusters that share the most edges with eachother (closest anchors?)
        while clusters.len() > 3 {
            // TODO add weighting to clusters on the same LOD level as this? largest triangle area?
            // ^ this is where the multivariate optimization problem comes into play
            // grab the highest value of sum_of_two_highest_shared_edges and grab the highest index that is in both of the other connecting vecs
            shared_edges.sort_by(|element1, element2| {
                let cmp = element1.lod.cmp(&element2.lod);

                // ensure connections is sorted
                if cmp == Ordering::Equal {
                    return element1
                        .connections
                        .get(0)
                        .cmp(&element2.connections.get(0));
                }

                cmp

                //                 let val1 = {
                //                     let first = element1.connections.get(0);
                //                     let second = element1.connections.get(1);
                //                     if first == None || second == None {
                //                         i32::MIN
                //                     } else {
                //                         first.unwrap().1 + second.unwrap().1
                //                     }
                //                 };
                //
                //                 let val2 = {
                //                     let first = element2.connections.get(0);
                //                     let second = element2.connections.get(1);
                //                     if first == None || second == None {
                //                         i32::MIN
                //                     } else {
                //                         first.unwrap().1 + second.unwrap().1
                //                     }
                //                 };
                //
                //                 val1.cmp(&val2)
            });

            if shared_edges.is_empty() {
                println!(
                    "shared edges is empty but clusters has {} elements left.",
                    clusters.len()
                );
            }

            // ? TODO determine if any clusters are abandoned by combining these, include that one, ie ring abandons tip of ear on bunny
            // ? Is this, ^, actually an issue? it will just get picked up by a later cluster, it would be unoptimal but it wouldn't happen very often at all

            // get the cluster with the largest adjacency, add all the adjacents to a HashMap<u32, u32>
            // choose the highest value, choose it as the next cluster, add all the adjacents of that
            // cluster to the HashMap, if the value exists just add the values together. continue until
            // you have 4 values

            let mut adjacency_map = HashMap::new();

            // TODO the clusters are not adjacent, they are not connected
            let cluster0_id = shared_edges.last().unwrap().id;
            let cluster0_idx = shared_edges.len() - 1;
            let cluster0 = shared_edges.last().unwrap();
            for adj in &cluster0.connections {
                adjacency_map.insert(adj.0, adj.1);
            }

            let cluster1_id = cluster0.connections.get(0).unwrap().0;
            let cluster1_idx = shared_edges
                .iter()
                .position(|e| e.id == cluster1_id)
                .unwrap();
            let cluster1 = shared_edges.get(cluster1_idx).unwrap();
            for adj in &cluster1.connections {
                let old = adjacency_map.insert(adj.0, adj.1);
                if let Some(value) = old {
                    *adjacency_map.get_mut(&adj.0).unwrap() += value;
                }
            }

            let mut adjacency_vec: Vec<(u32, i32)> =
                adjacency_map.iter().map(|v| (*v.0, *v.1)).collect();
            adjacency_vec.sort_by(|a, b| a.1.cmp(&b.1));

            let mut cluster2_id = u32::MAX;
            for adj in adjacency_vec.iter().rev() {
                if adj.0 != cluster0_id && adj.0 != cluster1_id {
                    cluster2_id = adj.0;
                    break;
                }
            }

            let cluster2_idx = shared_edges
                .iter()
                .position(|e| e.id == cluster2_id)
                .unwrap();
            let cluster2 = shared_edges.get(cluster2_idx).unwrap();
            for adj in &cluster2.connections {
                let old = adjacency_map.insert(adj.0, adj.1);
                if let Some(value) = old {
                    *adjacency_map.get_mut(&adj.0).unwrap() += value;
                }
            }

            let mut adjacency_vec: Vec<(u32, i32)> =
                adjacency_map.iter().map(|v| (*v.0, *v.1)).collect();
            adjacency_vec.sort_by(|a, b| a.1.cmp(&b.1));

            let mut cluster3_id = u32::MAX;
            for adj in adjacency_vec.iter().rev() {
                if adj.0 != cluster0_id && adj.0 != cluster1_id && adj.0 != cluster2_id {
                    cluster3_id = adj.0;
                    break;
                }
            }
            let cluster3_idx = shared_edges
                .iter()
                .position(|e| e.id == cluster3_id)
                .unwrap();
            let cluster3 = shared_edges.get(cluster3_idx).unwrap();
            // };

            println!("cluster0 id: {}", cluster0_id);
            println!("cluster0 adjacency: {:?}", cluster0.connections);
            println!();

            println!("cluster1 id: {}", cluster1_id);
            println!("cluster1 adjacency: {:?}", cluster1.connections);
            println!();

            println!("cluster2 id: {}", cluster2_id);
            println!("cluster2 adjacency: {:?}", cluster2.connections);
            println!();

            println!("cluster3 id: {}", cluster3_id);
            println!("cluster3 adjacency: {:?}", cluster3.connections);
            println!();

            // TODO new algorithm for choosing what to merge
            // sort the list by lod level with the a + b as a tie breaker?

            // combine the clusters into one mesh
            // let mut triangles = Vec::<Triangle>::with_capacity(
            //     clusters[cluster0_id].tris.len()
            //         + clusters[cluster1_id].tris.len()
            //         + clusters[cluster2_id].tris.len()
            //         + clusters[cluster3_id].tris.len(),
            // );
            // triangles.extend_from_slice(&clusters[cluster0_id].tris);
            // triangles.extend_from_slice(&clusters[cluster1_id].tris);
            // triangles.extend_from_slice(&clusters[cluster2_id].tris);
            // triangles.extend_from_slice(&clusters[cluster3_id].tris);

            // idx += clusters[cluster3_id].positions.len();

            // let border = clusters[starting_point.id].cut
            //     .symmetric_difference(&clusters[cluster1_id].cut).map(|value| { value.clone() }).collect::<HashSet<[u32; 2]>>()
            //     .symmetric_difference(&clusters[cluster2_id].cut).map(|value| { value.clone() }).collect::<HashSet<[u32; 2]>>()
            //     .symmetric_difference(&clusters[cluster3_id].cut).map(|value| { value.clone() }).collect::<HashSet<[u32; 2]>>();

            let (mesh_indices, mesh_positions) = {
                let mut indices = Vec::with_capacity(
                    ctree.get(&cluster0_id).mesh.indices.len()
                        + ctree.get(&cluster1_id).mesh.indices.len()
                        + ctree.get(&cluster2_id).mesh.indices.len()
                        + ctree.get(&cluster3_id).mesh.indices.len(),
                );

                let mut positions = Vec::with_capacity(
                    ctree.get(&cluster0_id).mesh.positions.len()
                        + ctree.get(&cluster1_id).mesh.positions.len()
                        + ctree.get(&cluster2_id).mesh.positions.len()
                        + ctree.get(&cluster3_id).mesh.positions.len(),
                );

                // might be faster with a .extend_from_slice(clusters[cluster0_id].indices.map(|a| { a + idx })) etc
                let mut idx = 0;
                for i in &ctree.get(&cluster0_id).mesh.indices {
                    indices.push(i + idx);
                }
                positions.extend_from_slice(&ctree.get(&cluster0_id).mesh.positions);
                idx += ctree.get(&cluster0_id).mesh.positions.len() as u32 / 3;
                for i in &ctree.get(&cluster1_id).mesh.indices {
                    indices.push(i + idx);
                }
                positions.extend_from_slice(&ctree.get(&cluster1_id).mesh.positions);
                idx += ctree.get(&cluster1_id).mesh.positions.len() as u32 / 3;
                for i in &ctree.get(&cluster2_id).mesh.indices {
                    indices.push(i + idx);
                }
                positions.extend_from_slice(&ctree.get(&cluster2_id).mesh.positions);
                idx += ctree.get(&cluster2_id).mesh.positions.len() as u32 / 3;
                for i in &ctree.get(&cluster3_id).mesh.indices {
                    indices.push(i + idx);
                }
                positions.extend_from_slice(&ctree.get(&cluster3_id).mesh.positions);

                // TODO I have duplicated vertices with different indices in this list because they have their own indices in their index lists
                remove_duplicated_vertices(&mut indices, &mut positions);

                simplify_indices_positions(&indices, &positions)
            };

            // println!("max index: {}, positions len: {}", mesh_indices.iter().max().unwrap(), mesh_positions.len());

            let combined = Mesh {
                positions: mesh_positions,
                indices: mesh_indices,
            };

            // mesh simplification
            let cluster = simplify(&combined, (TRIS_IN_CLUSTER * 2) as u32);

            println!("writing cluster0 id: {}", &cluster0_id);
            write_mesh_to_file("cluster0", &ctree.get(&cluster0.id).mesh);
            println!("writing cluster1 id: {}", &cluster1_id);
            write_mesh_to_file("cluster1", &ctree.get(&cluster1.id).mesh);
            println!("writing cluster2 id: {}", &cluster2_id);
            write_mesh_to_file("cluster2", &ctree.get(&cluster2.id).mesh);
            println!("writing cluster3 id: {}", &cluster3_id);
            write_mesh_to_file("cluster3", &ctree.get(&cluster3.id).mesh);
            // println!(
            //     "writing new_indices with {} tris",
            //     cluster.indices.len() / 3
            // );
            // write_mesh_to_file("new_indices", &cluster.positions, &cluster.indices);

            if cluster.indices.len() > TRIS_IN_CLUSTER * 3 * 2 {
                panic!("unable to simplify the cluster enough");
            }

            // calculate new cut
            // let (_, tris, _) = generate_data_structures(&cluster.indices);

            let combined_anchor = cluster.calculate_anchor();
            // split cluster into 2, choose any triangle on the border as a seed, same algorithm as above

            // TODO, some free floating triangles
            let mut furthest = ([0.0f32, 0.0, 0.0], 0.0 as f32);
            // let mut furthest_tri: Triangle = [0, 0, 0];
            for idx in (0..cluster.indices.len()).step_by(3) {
                let mid_point = cluster.mid_point_3(
                    cluster.indices[idx + 0],
                    cluster.indices[idx + 1],
                    cluster.indices[idx + 2],
                );
                let dist = cluster.sqr_dist_w_2custom(
                    cluster.mid_point_3(
                        cluster.indices[idx + 0],
                        cluster.indices[idx + 1],
                        cluster.indices[idx + 2],
                    ),
                    combined_anchor,
                );
                if dist > furthest.1 {
                    furthest = (mid_point, dist);
                    // furthest_tri = *tri;
                }
            }

            let mut cluster1_idxs = Vec::with_capacity(TRIS_IN_CLUSTER);
            let mut cluster1_edges = HashMap::<Edge, [u32; 2]>::with_capacity(TRIS_IN_CLUSTER * 3);

            let mut cluster2_idxs = Vec::with_capacity(TRIS_IN_CLUSTER);
            let mut cluster2_edges = HashMap::<Edge, [u32; 2]>::with_capacity(TRIS_IN_CLUSTER * 3);

            let mut sorted_tris = Vec::new();
            for idx in (0..cluster.indices.len()).step_by(3) {
                let mid_point = cluster.mid_point_3(
                    cluster.indices[idx + 0],
                    cluster.indices[idx + 1],
                    cluster.indices[idx + 2],
                );
                let dist1 = cluster.sqr_dist_w_2custom(furthest.0, mid_point);

                sorted_tris.push((
                    dist1,
                    [
                        cluster.indices[idx + 0],
                        cluster.indices[idx + 1],
                        cluster.indices[idx + 2],
                    ],
                ));
            }
            sorted_tris.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

            // TODO also have to make sure all these are connected, could have the two sides of an ear which would be closer but also disconnected
            // ? ^ is disconnected an issue at all? does it matter for rendering? it might make a difference for culling?
            for idx in 0..sorted_tris.len() {
                if idx < TRIS_IN_CLUSTER {
                    cluster1_idxs.push(sorted_tris[idx].1[0]);
                    cluster1_idxs.push(sorted_tris[idx].1[1]);
                    cluster1_idxs.push(sorted_tris[idx].1[2]);

                    {
                        let key: Edge = new_edge(sorted_tris[idx].1[0], sorted_tris[idx].1[1]);
                        let value = cluster1_edges.insert(key, [sorted_tris[idx].1[2], u32::MAX]);
                        if let Some(v) = value {
                            if v[1] != u32::MAX {
                                println!(
                                    "this edge has appeared more than twice: {:?} => {:?} + {}",
                                    key, v, sorted_tris[idx].1[2]
                                );
                            }
                            cluster1_edges.insert(key, [sorted_tris[idx].1[2], v[0]]);
                        }
                    }
                    {
                        let key: Edge = new_edge(sorted_tris[idx].1[0], sorted_tris[idx].1[2]);
                        let value = cluster1_edges.insert(key, [sorted_tris[idx].1[1], u32::MAX]);
                        if let Some(v) = value {
                            if v[1] != u32::MAX {
                                println!(
                                    "this edge has appeared more than twice: {:?} => {:?} + {}",
                                    key, v, sorted_tris[idx].1[1]
                                );
                            }
                            cluster1_edges.insert(key, [sorted_tris[idx].1[1], v[0]]);
                        }
                    }
                    {
                        let key: Edge = new_edge(sorted_tris[idx].1[1], sorted_tris[idx].1[2]);
                        let value = cluster1_edges.insert(key, [sorted_tris[idx].1[0], u32::MAX]);
                        if let Some(v) = value {
                            if v[1] != u32::MAX {
                                println!(
                                    "this edge has appeared more than twice: {:?} => {:?} + {}",
                                    key, v, sorted_tris[idx].1[0]
                                );
                            }
                            cluster1_edges.insert(key, [sorted_tris[idx].1[0], v[0]]);
                        }
                    }
                } else {
                    cluster2_idxs.push(sorted_tris[idx].1[0]);
                    cluster2_idxs.push(sorted_tris[idx].1[1]);
                    cluster2_idxs.push(sorted_tris[idx].1[2]);

                    {
                        let key: Edge = new_edge(sorted_tris[idx].1[0], sorted_tris[idx].1[1]);
                        let value = cluster2_edges.insert(key, [sorted_tris[idx].1[2], u32::MAX]);
                        if let Some(v) = value {
                            if v[1] != u32::MAX {
                                println!(
                                    "this edge has appeared more than twice: {:?} => {:?} + {}",
                                    key, v, sorted_tris[idx].1[2]
                                );
                            }
                            cluster2_edges.insert(key, [sorted_tris[idx].1[2], v[0]]);
                        }
                    }
                    {
                        let key: Edge = new_edge(sorted_tris[idx].1[0], sorted_tris[idx].1[2]);
                        let value = cluster2_edges.insert(key, [sorted_tris[idx].1[1], u32::MAX]);
                        if let Some(v) = value {
                            if v[1] != u32::MAX {
                                println!(
                                    "this edge has appeared more than twice: {:?} => {:?} + {}",
                                    key, v, sorted_tris[idx].1[1]
                                );
                            }
                            cluster2_edges.insert(key, [sorted_tris[idx].1[1], v[0]]);
                        }
                    }
                    {
                        let key: Edge = new_edge(sorted_tris[idx].1[1], sorted_tris[idx].1[2]);
                        let value = cluster2_edges.insert(key, [sorted_tris[idx].1[0], u32::MAX]);
                        if let Some(v) = value {
                            if v[1] != u32::MAX {
                                println!(
                                    "this edge has appeared more than twice: {:?} => {:?} + {}",
                                    key, v, sorted_tris[idx].1[0]
                                );
                            }
                            cluster2_edges.insert(key, [sorted_tris[idx].1[0], v[0]]);
                        }
                    }
                }
            }

            let combined1_anchor = cluster.calculate_anchor();
            let combined2_anchor = cluster.calculate_anchor();

            // calculate cut
            let combined1_cut = {
                let mut cut = HashSet::<Edge>::new();
                for kv in cluster1_edges {
                    if kv.1[0] == u32::MAX || kv.1[1] == u32::MAX {
                        cut.insert(kv.0);
                    }
                }

                cut
            };
            let combined2_cut = {
                let mut cut = HashSet::<Edge>::new();
                for kv in cluster2_edges {
                    if kv.1[0] == u32::MAX || kv.1[1] == u32::MAX {
                        cut.insert(kv.0);
                    }
                }

                cut
            };

            let (indices, positions) =
                simplify_indices_positions(&cluster1_idxs, &cluster.positions.clone());
            let combined1 = Cluster {
                mesh: Mesh { positions, indices },
                cut: combined1_cut,
                anchor: combined1_anchor,
            };
            let (indices, positions) =
                simplify_indices_positions(&cluster2_idxs, &cluster.positions.clone());
            let combined2 = Cluster {
                mesh: Mesh { positions, indices },
                cut: combined2_cut,
                anchor: combined2_anchor,
            };

            let idx1 = ctree.insert(
                combined1,
                &[u32::MAX, u32::MAX],
                &[cluster0_id, cluster1_id, cluster2_id, cluster3_id],
            );
            clusters.push(idx1);

            let idx2 = ctree.insert(
                combined2,
                &[u32::MAX, u32::MAX],
                &[cluster0_id, cluster1_id, cluster2_id, cluster3_id],
            );
            clusters.push(idx2);

            let removals = vec![cluster0_id, cluster1_id, cluster2_id, cluster3_id];
            let mut remove_idxs: Vec<usize> = (0..clusters.len())
                .filter(|idx| removals.contains(&clusters[*idx]))
                .collect();
            remove_idxs.sort();
            for r in remove_idxs.iter().rev() {
                clusters.remove(*r as usize);
            }

            // println!("cluster0_index: {}", cluster0_idx);
            // println!("cluster1_index: {}", cluster1_idx);
            // println!("cluster2_index: {}", cluster2_idx);
            // println!("cluster3_index: {}", cluster3_idx);
            let mut remove_idxs = [cluster0_idx, cluster1_idx, cluster2_idx, cluster3_idx];
            remove_idxs.sort();
            let mut old_shared = Vec::with_capacity(4);
            for r in remove_idxs.iter().rev() {
                old_shared.push(shared_edges.remove(*r as usize));
            }

            // fix the overlapping edges on the other clusters in the lists
            let (mut connections1, mut connections2) = recalculate_shared_edges(
                &ctree,
                [cluster0_id, cluster1_id, cluster2_id, cluster3_id],
                [idx1, idx2],
                &clusters,
                &old_shared
                    .iter()
                    .flat_map(|value| value.connections.clone())
                    .collect(),
                &mut shared_edges,
            );

            let mut lod_levels: Vec<u32> = old_shared.iter().map(|e| e.lod).collect();
            lod_levels.sort();
            let lod = *lod_levels.last().unwrap() + 1;
            connections1.sort_by(|a, b| b.1.cmp(&a.1));
            shared_edges.push(SharedEdges {
                id: idx1,
                lod,
                connections: connections1,
            });
            connections2.sort_by(|a, b| b.1.cmp(&a.1));
            shared_edges.push(SharedEdges {
                id: idx2,
                lod,
                connections: connections2,
            });
        }
        // simplify the large clusters triangles
        // split the cluster into two clusters of TRIS_IN_CLUSTER size
        // add the clusters back to the list of clusters to be combined
        // end when there is only one cluster leftover
    }
    // !SECTION - Combine the clusters into a graph

    // !SECTION - Generate debug.obj mesh
    //     println!("Generate debug.obj mesh");
    //     {
    //         let mut f = File::create("debug.obj").unwrap();
    //         f.write_all(b"mtllib bunny.mtl\no bun_zipper\n").unwrap();
    //         for v_idx in (0..mesh.positions.len()).step_by(3) {
    //             f.write_all(
    //                 format!(
    //                     "v {} {} {}\n",
    //                     mesh.positions[v_idx + 0],
    //                     mesh.positions[v_idx + 1],
    //                     mesh.positions[v_idx + 2]
    //                 )
    //                 .as_bytes(),
    //             )
    //             .unwrap();
    //         }
    //         f.write_all(b"usemtl None\ns off\n").unwrap();
    //
    //         // for t in tris {
    //         //     f.write_all(format!("f {} {} {}\n", t[0] + 1, t[1] + 1, t[2] + 1).as_bytes())
    //         //         .unwrap();
    //         // }
    //         // for verts in cut_vertices {
    //         //     for v in verts.1 {
    //         //         f.write_all(format!("l {} {}\n", verts.0 + 1, v + 1).as_bytes())
    //         //             .unwrap();
    //         //     }
    //         // }
    //         // for e in edges {
    //         //     f.write_all(format!("l {} {}\n", e.0[0] + 1, e.0[1] + 1).as_bytes())
    //         //         .unwrap();
    //         // }
    //         // for p in pairs {
    //         //     f.write_all(
    //         //         format!("l {} {}\n", boundary_groups[p[0]][0] + 1, boundary_groups[p[1]][0] + 1).as_bytes(),
    //         //     )
    //         //     .unwrap();
    //         // }
    //         // {
    //         //     // for idx in 1..final_cut_path.len() {}
    //         //     let path = final_cut_path
    //         //         .iter()
    //         //         .map(|pt| format!("{}", pt + 1))
    //         //         .collect::<Vec<String>>()
    //         //         .join(" ");
    //         //     f.write_all(format!("l {}\n", path).as_bytes()).unwrap();
    //         // }
    //
    //         // !Most recent vvv
    //         // for e in edge_path {
    //         //     f.write_all(format!("l {} {}\n", e[0] + 1, e[1] + 1).as_bytes())
    //         //         .unwrap();
    //         // }
    //         // for group in &boundary_groups {
    //         //     for idx in 1..group.len() {
    //         //         f.write_all(
    //         //             format!("l {} {}\n", group[idx - 1] + 1, group[idx] + 1).as_bytes(),
    //         //         )
    //         //         .unwrap();
    //         //     }
    //         // }
    //
    //         //             f.write_all(format!("o prev_cluster_cut_1\n").as_bytes())
    //         //                 .unwrap();
    //         //             f.write_all(prev_cluster_cut.as_bytes()).unwrap();
    //         //
    //         //             f.write_all(format!("o last_cluster_cut_2\n").as_bytes())
    //         //                 .unwrap();
    //         //             f.write_all(last_cluster_cut.as_bytes()).unwrap();
    //         //
    //         //             f.write_all(format!("o prev_overall_cut_1\n").as_bytes())
    //         //                 .unwrap();
    //         //             f.write_all(prev_overall_cut.as_bytes()).unwrap();
    //         //
    //         //             f.write_all(format!("o last_overall_cut_2\n").as_bytes())
    //         //                 .unwrap();
    //         //             f.write_all(last_overall_cut.as_bytes()).unwrap();
    //
    //         for (idx, cluster) in clusters.iter().enumerate() {
    //             f.write_all(format!("o Cluster{}\n", idx + 1).as_bytes())
    //                 .unwrap();
    //             for idx in (0..cluster.mesh.indices.len()).step_by(3) {
    //                 f.write_all(
    //                     format!(
    //                         "f {} {} {}\n",
    //                         cluster.mesh.indices[idx + 0] + 1,
    //                         cluster.mesh.indices[idx + 1] + 1,
    //                         cluster.mesh.indices[idx + 2] + 1
    //                     )
    //                     .as_bytes(),
    //                 )
    //                 .unwrap();
    //             }
    //             f.write_all("\n".as_bytes()).unwrap();
    //         }
    //
    //         f.flush().unwrap();
    //     };
    // !SECTION - Generate debug.obj mesh

    write_clusters(&ctree.clusters());

    panic!();

    ctree
}

// fn add_tri(
//     tri: &Triangle,
//     tris: &mut HashSet<Triangle>,
//     edges: &mut HashMap<Edge, [u32; 2]>,
//     overall_cluster_cut: &mut HashSet<Edge>,
//     cur_cut: &mut HashSet<Edge>,
// ) {
//     tris.insert(*tri);
//
//     let edge0 = new_edge(tri[0], tri[1]);
//     let edge1 = new_edge(tri[0], tri[2]);
//     let edge2 = new_edge(tri[1], tri[2]);
//
//     let mut add0 = false;
//     let mut add1 = false;
//     let mut add2 = false;
//
//     if let Some(value) = edges.get_mut(&edge0) {
//         if value[0] == u32::MAX {
//             value[0] = tri[2];
//         } else if value[1] == u32::MAX {
//             value[1] = tri[2];
//         }
//     } else {
//         add0 = true;
//     }
//
//     if let Some(value) = edges.get_mut(&edge1) {
//         if value[0] == u32::MAX {
//             value[0] = tri[1];
//         } else if value[1] == u32::MAX {
//             value[1] = tri[1];
//         }
//     } else {
//         add1 = true;
//     }
//
//     if let Some(value) = edges.get_mut(&edge2) {
//         if value[0] == u32::MAX {
//             value[0] = tri[0];
//         } else if value[1] == u32::MAX {
//             value[1] = tri[0];
//         }
//     } else {
//         add2 = true;
//     }
//
//     if add0 {
//         if overall_cluster_cut.contains(&edge0) {
//             // cur_cut.insert(edge0);
//             overall_cluster_cut.remove(&edge0);
//         } else {
//             // cur_cut.remove(&edge0);
//             overall_cluster_cut.insert(edge0);
//         }
//         edges.insert(edge0, new_edge(tri[2], u32::MAX));
//     }
//     // else {
//     //     cur_cut.insert(edge0);
//     // }
//
//     if add1 {
//         if !overall_cluster_cut.contains(&edge1) {
//             // cur_cut.insert(edge1);
//             overall_cluster_cut.remove(&edge1);
//         } else {
//             // cur_cut.insert(edge1);
//             overall_cluster_cut.insert(edge1);
//         }
//         edges.insert(edge1, new_edge(tri[1], u32::MAX));
//     }
//     // else {
//     //     cur_cut.insert(edge1);
//     // }
//
//     if add2 {
//         if !overall_cluster_cut.contains(&edge2) {
//             // cur_cut.insert(edge2);
//             overall_cluster_cut.remove(&edge2);
//         } else {
//             // cur_cut.insert(edge2);
//             overall_cluster_cut.insert(edge2);
//         }
//         edges.insert(edge2, new_edge(tri[0], u32::MAX));
//     }
//     // else {
//     //     cur_cut.insert(edge2);
//     // }
//
//     if !cur_cut.contains(&edge0) {
//         cur_cut.insert(edge0);
//     } else {
//         cur_cut.remove(&edge0);
//     }
//
//     if !cur_cut.contains(&edge1) {
//         cur_cut.insert(edge1);
//     } else {
//         cur_cut.remove(&edge1);
//     }
//
//     if !cur_cut.contains(&edge2) {
//         cur_cut.insert(edge2);
//     } else {
//         cur_cut.remove(&edge2);
//     }
// }

fn recalculate_shared_edges(
    ctree: &CTree,
    combined: [u32; 4],
    new: [u32; 2],
    clusters: &Vec<u32>,
    connections: &Vec<(u32, i32)>,
    shared_edges: &mut Vec<SharedEdges>,
) -> (Vec<(u32, i32)>, Vec<(u32, i32)>) {
    let connections: HashSet<u32> = connections.iter().map(|x| x.0).collect();
    let mut connections1 = Vec::new();
    let mut connections2 = Vec::new();
    for connection in connections {
        // TODO is this right? it can only be on this level?
        if !clusters.contains(&connection) {
            continue;
        }

        // test first
        {
            let connected = ctree.get(&connection);
            let count = ctree.get(&new[0]).cut.intersection(&connected.cut).count();
            if count > 0 {
                connections1.push((connection, count as i32));
            }
        };
        // test second
        {
            let connected = ctree.get(&connection);
            let count = ctree.get(&new[1]).cut.intersection(&connected.cut).count();
            if count > 0 {
                connections2.push((connection, count as i32));
            }
        };

        // shared_edges[connection].remove(,position(combined))

        //         for id in 0..clusters.len() as u32 {
        //             let mut connections = Vec::with_capacity(clusters.len());
        //             for j in 0..clusters.len() as u32 {
        //                 if id == j {
        //                     continue;
        //                 }
        //
        //                 let count = ctree.get(&id).cut.intersection(&ctree.get(&j).cut).count();
        //                 if count > 0 {
        //                     connections.push((j, count as i32));
        //                 }
        //             }
        //             shared_edges.push(SharedEdges { id, connections });
        //         }
    }

    for sedge_idx in 0..shared_edges.len() {
        // if connections
        //     .iter()
        //     .any(|e| e.0 == shared_edges[sedge_idx].id)
        {
            let mut removals = Vec::new();
            for idx in 0..shared_edges[sedge_idx].connections.len() {
                if combined.contains(&shared_edges[sedge_idx].connections[idx].0) {
                    removals.push(idx);
                }
            }
            for r in removals.iter().rev() {
                shared_edges[sedge_idx].connections.remove(*r);
            }
        }

        for c in &connections1 {
            if c.0 == shared_edges[sedge_idx].id {
                shared_edges[sedge_idx].connections.push((new[0], c.1));
                break;
            }
        }

        for c in &connections2 {
            if c.0 == shared_edges[sedge_idx].id {
                shared_edges[sedge_idx].connections.push((new[1], c.1));
                break;
            }
        }
    }

    (connections1, connections2)
}

fn write_clusters(clusters: &Vec<Cluster>) {
    let mut rng = rand::thread_rng();

    let mut f = File::create("clusters.obj").unwrap();
    f.write_all(b"mtllib bunny.mtl\no bun_zipper\n").unwrap();

    let mut offsets = Vec::new();
    let mut offset = 0;
    for cluster in clusters {
        // f.write_all(format!("o Cluster{}\n", idx + 1).as_bytes())
        //     .unwrap();
        let color: [f32; 3] = [rng.gen(), rng.gen(), rng.gen()];
        offsets.push(offset);
        for pos in (0..cluster.mesh.positions.len()).step_by(3) {
            f.write_all(
                format!(
                    "v {} {} {} {} {} {}\n",
                    cluster.mesh.positions[pos + 0],
                    cluster.mesh.positions[pos + 1],
                    cluster.mesh.positions[pos + 2],
                    color[0],
                    color[1],
                    color[2]
                )
                .as_bytes(),
            )
            .unwrap();
        }
        f.write_all("\n".as_bytes()).unwrap();
        offset += cluster.mesh.positions.len() / 3;
    }

    f.write_all(b"usemtl None\ns off\n").unwrap();

    for (idx, cluster) in clusters.iter().enumerate() {
        f.write_all(format!("o Cluster{}\n", idx + 1).as_bytes())
            .unwrap();
        let offset = offsets[idx];
        for idx in (0..cluster.mesh.indices.len()).step_by(3) {
            f.write_all(
                format!(
                    "f {} {} {}\n",
                    cluster.mesh.indices[idx + 0] + 1 + offset as u32,
                    cluster.mesh.indices[idx + 1] + 1 + offset as u32,
                    cluster.mesh.indices[idx + 2] + 1 + offset as u32
                )
                .as_bytes(),
            )
            .unwrap();
        }
        f.write_all("\n".as_bytes()).unwrap();
    }

    f.flush().unwrap();
}

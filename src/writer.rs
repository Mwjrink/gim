use core::panic;
use rand::Rng;
use std::collections::{HashMap, HashSet};
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

        println!("num_face_indices: {}", mesh.num_face_indices.len());

        let mut idx = 0;
        for verts in &mesh.num_face_indices {
            {
                let key: Edge = edge(mesh.indices[idx + 0], mesh.indices[idx + 1]);
                let value = edges.insert(key, [mesh.indices[idx + 2], 0]);
                if let Some(v) = value {
                    if v[1] != 0 {
                        println!("this edge has appeared more than twice: {:?}", key);
                    }
                    edges.insert(key, [mesh.indices[idx + 2], v[0]]);
                }
            }
            {
                let key: Edge = edge(mesh.indices[idx + 0], mesh.indices[idx + 2]);
                let value = edges.insert(key, [mesh.indices[idx + 1], 0]);
                if let Some(v) = value {
                    if v[1] != 0 {
                        println!("this edge has appeared more than twice: {:?}", key);
                    }
                    edges.insert(key, [mesh.indices[idx + 1], v[0]]);
                }
            }
            {
                let key: Edge = edge(mesh.indices[idx + 1], mesh.indices[idx + 2]);
                let value = edges.insert(key, [mesh.indices[idx + 0], 0]);
                if let Some(v) = value {
                    if v[1] != 0 {
                        println!("this edge has appeared more than twice: {:?}", key);
                    }
                    edges.insert(key, [mesh.indices[idx + 0], v[0]]);
                }
            }

            tris.insert(tri(mesh.indices[idx + 0], mesh.indices[idx + 1], mesh.indices[idx + 2]));

            idx += *verts as usize;
        }

        let tri_idx = rand::thread_rng().gen_range(0..mesh.num_face_indices.len());
        let mut triangle = tri(
            mesh.indices[tri_idx * 3 + 0],
            mesh.indices[tri_idx * 3 + 1],
            mesh.indices[tri_idx * 3 + 2],
        );

        if !tris.contains(&triangle) {
            println!("random seed triangle not in the mesh");
        }

        let mut inst_tris = Vec::<(Triangle, f32, Edge)>::with_capacity(3);

        let dist = |a: u32, b: u32| -> f32 {
            let p1 = &mesh.positions[a as usize..(a + 2) as usize];
            let p2 = &mesh.positions[b as usize..(b + 2) as usize];

            f32::sqrt(f32::powf(p1[0] - p2[0], 2.0) + f32::powf(p1[0] - p2[0], 2.0) + f32::powf(p1[0] - p2[0], 2.0))
        };

        println!("{} total edges before removal", edges.len());

        let mut removable_edges = Vec::<Edge>::new();

        'part_one: loop {
            if !tris.remove(&triangle) {
                println!("this is not possible, loop triangle not in the mesh");
            }

            // find the adjacent tris
            {
                let edg = edge(triangle[0], triangle[1]);
                if let Some(value) = edges.get(&edg) {
                    let idx = if value[0] == triangle[2] { 1 } else { 0 };
                    let push_trig = tri(triangle[0], triangle[1], value[idx]);
                    if tris.contains(&push_trig) {
                        // this distance is the distance between the points not on the edge itself
                        inst_tris.push((push_trig, dist(triangle[2], value[idx]), edg));
                    }
                }
            };

            {
                let edg = edge(triangle[0], triangle[2]);
                if let Some(value) = edges.get(&edg) {
                    let idx = if value[0] == triangle[1] { 1 } else { 0 };
                    let push_trig = tri(triangle[0], triangle[2], value[idx]);
                    if tris.contains(&push_trig) {
                        inst_tris.push((push_trig, dist(triangle[1], value[idx]), edg));
                    }
                }
            };

            {
                let edg = edge(triangle[1], triangle[2]);
                if let Some(value) = edges.get(&edg) {
                    let idx = if value[0] == triangle[0] { 1 } else { 0 };
                    let push_trig = tri(triangle[1], triangle[2], value[idx]);
                    if tris.contains(&push_trig) {
                        inst_tris.push((push_trig, dist(triangle[0], value[idx]), edg));
                    }
                }
            };

            if inst_tris.len() == 0 {
                'stack: while let Some(edge) = removable_edges.pop() {
                    // remove edge and the triangle adjacent to it

                    // check if there is another triangle attached to this edge
                    // if the edge is not in the simplicial complex, it has already been removed and we can move to the next on the stack
                    if let Some(value) = edges.get(&edge) {
                        let idx = if triangle[0] != value[0] && triangle[1] != value[0] && triangle[2] != value[0] {
                            0
                        } else if triangle[0] != value[1] && triangle[1] != value[1] && triangle[2] != value[1] {
                            1
                        } else {
                            panic!("this is not possible");
                        };
                        let push_trig = tri(edge[0], edge[1], value[idx]);
                        if tris.contains(&push_trig) {
                            // if there is: remove edge and set the triangle for the next loop iteration
                            edges.remove(&edge);
                            triangle = push_trig;
                            continue 'part_one;
                        }
                    }
                }

                {
                    println!("There are no adjacent triangles to this triangle AND no removables in the queue... fuck");
                    println!(
                        "triggles: {}, edges: {}, removables: {}",
                        tris.len(),
                        edges.len(),
                        removable_edges.len()
                    );
                    break 'part_one;
                }
            } else {
                // three possible edges to remove from
                // sort by distance to the current tri and choose the closest in the simplicial complex
                inst_tris.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

                // remove the corresponding edge
                edges.remove(&inst_tris[0].2);

                // and set the triangle for the next loop iteration
                triangle = inst_tris[0].0;

                if inst_tris.len() == 2 {
                    removable_edges.push(inst_tris[1].2);
                } else if inst_tris.len() == 3 {
                    removable_edges.push(inst_tris[1].2);
                    removable_edges.push(inst_tris[2].2);
                }
            }

            // empty the temp to prep for next tri
            inst_tris.clear();
        }

        // DEBUG
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
            for t in tris {
                f.write_all(format!("f {} {} {}\n", t[0] + 1, t[1] + 1, t[2] + 1).as_bytes())
                    .unwrap();
            }
            f.flush().unwrap();
        };

        // second half of the algorithm, simplifying p
        // loop {
        //     //
        //     break;
        // }
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

/*

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

/*
#[derive(Clone)]
        struct Node {
            tris: Vec<Triangle>,
        }
        let mut nodes = vec![Node { tris: Vec::new() }; mesh.indices.len() / 3];
        let mut tris = Vec::<Triangle>::with_capacity(mesh.num_face_indices.len());
        let mut current_index = 0;
        for f in 0..mesh.num_face_indices.len() {
            let tri = [
                mesh.indices[current_index],
                mesh.indices[current_index + 1],
                mesh.indices[current_index + 2],
                (current_index / 3) as _,
            ];

            nodes[mesh.indices[current_index] as usize].tris.push(tri);
            nodes[mesh.indices[current_index + 1] as usize].tris.push(tri);
            nodes[mesh.indices[current_index + 2] as usize].tris.push(tri);

            tris.push(tri);

            current_index += mesh.num_face_indices[f] as usize;
        }

        // indices x 2 is likely too big but is decent as a worst case scenario and should prevent reallocation in the middle of the loop (not a huge deal but might as well)
        let mut idxs = Vec::<u32>::with_capacity(mesh.indices.len() * 2);
        let mut consumed = HashSet::new();

        fn findTri(tri: &[u32; 4], node: &Node, consumed: &HashSet<usize>) -> (u32, u32) {
            for t in &node.tris {
                let index = t[3] as usize;
                if consumed.contains(&index) {
                    continue;
                }

                if tri[1] == t[0] {
                    if t[1] == tri[0] {
                        return (t[2], t[3]);
                    } else {
                        return (t[1], t[3]);
                    }
                }
                if tri[1] == t[1] {
                    if t[0] == tri[0] {
                        return (t[2], t[3]);
                    } else {
                        return (t[0], t[3]);
                    }
                }
                if tri[1] == t[2] {
                    if t[1] == tri[0] {
                        return (t[0], t[3]);
                    } else {
                        return (t[1], t[3]);
                    }
                }
                if tri[2] == t[0] {
                    if t[1] == tri[0] {
                        return (t[2], t[3]);
                    } else {
                        return (t[1], t[3]);
                    }
                }
                if tri[2] == t[1] {
                    if t[0] == tri[0] {
                        return (t[2], t[3]);
                    } else {
                        return (t[0], t[3]);
                    }
                }
                if tri[2] == t[2] {
                    if t[1] == tri[0] {
                        return (t[0], t[3]);
                    } else {
                        return (t[1], t[3]);
                    }
                }
            }

            (u32::MAX, u32::MAX)
        }

        let mut current = 0;
        {
            let tri = &tris[current];
            idxs.extend_from_slice(tri);
            consumed.insert(current);

            let node = &nodes[tri[0] as usize];

            // try to find a triangle sharing an edge with this triangle
            let (v, c) = findTri(tri, node, &consumed);

            current = c as _;
            idxs.push(v);
            consumed.insert(current);
        }

        loop {
            let tri = &tris[current];

            // not all
            let node = &nodes[tri[0] as usize];

            // try to find a triangle sharing an edge with this triangle
            let (v, c) = findTri(tri, node, &consumed);

            if v == u32::MAX {
                println! {"AHH {} / {}", idxs.len(), mesh.indices.len()};
                break;
            }

            current = c as _;
            idxs.push(v);
            consumed.insert(current);
        }
    }

    println!("complete.");

    */

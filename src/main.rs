use gim::mesh::*;
// 3.15 GB (3,390,337,024 bytes)

// So, say an object that is 10 feet tall is 100 feet away. If I hold up a ruler 3 feet away, then the object in the
// distance would correspond to about how many inches? => x/3 = 10/100

fn main() {
    let obj_file = "./test_assets/bunny.obj".to_string();

    let (models, materials) = tobj::load_obj(&obj_file, true).expect("Failed to load file");

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

        // TODO probably bad cause these could be HUGE, do something else, ...
        // TODO... but that may have to wait for a better asset importer or a switch to gltf2 or usda
        let mut mesh = Mesh {
            positions: mesh.positions.clone(),
            indices: mesh.indices.clone(),
        };

        let _output = gim::writer::write(&mut mesh);
    }

    return ();
}

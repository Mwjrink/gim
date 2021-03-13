// use exr::prelude::simple_image::*;
use std::{fs::File, io::prelude::*};

use gim::reader::Model;

// 3.15 GB (3,390,337,024 bytes)

// So, say an object that is 10 feet tall is 100 feet away. If I hold up a ruler 3 feet away, then the object in the
// distance would correspond to about how many inches? => x/3 = 10/100

// TODO, make it so that the files have the number of bytes as the first 4
// bytes and then read that to read the rest of the file, should be faster.
// can also have albedo, x, y, z etc in the same file

fn main() {
    let obj_file = "./test_assets/bunny.obj".to_string();

    gim::writer::write(&obj_file);

    return ();

    // }

    // load files
    // let model = Model::load("./test_assets/chalet");

    // println!("vertices(x) {}", model.x.len());
    // println!("vertices(y) {}", model.y.len());
    // println!("vertices(z) {}", model.z.len());
    // println!("vertices(albedo) {}", model.albedo.raw_pixels.len());

    // // TODO: test everything not just x
    // // now test to see if the files are the same

    // let mut i = 0;
    // let mut next_face = 0;
    // for f in 0..models[0].mesh.num_face_indices.len() {
    //     let end = next_face + models[0].mesh.num_face_indices[f] as usize;

    //     for t in next_face..end {
    //         let findex: usize = (models[0].mesh.indices[t]) as usize;
    //         if models[0].mesh.positions[findex * 3] != model.x[i] {
    //             println!(
    //                 "unequal x at {}, with: {} vs {}",
    //                 i,
    //                 models[0].mesh.positions[findex * 3],
    //                 model.x[i]
    //             );
    //             return ();
    //         }

    //         if models[0].mesh.positions[findex * 3 + 1] != model.y[i] {
    //             println!(
    //                 "unequal y at {}, with: {} vs {}",
    //                 i,
    //                 models[0].mesh.positions[findex * 3 + 1],
    //                 model.y[i]
    //             );
    //             return ();
    //         }

    //         if models[0].mesh.positions[findex * 3 + 2] != model.z[i] {
    //             println!(
    //                 "unequal z at {}, with: {} vs {}",
    //                 i,
    //                 models[0].mesh.positions[findex * 3 + 2],
    //                 model.z[i]
    //             );
    //             return ();
    //         }

    //         if models[0].mesh.texcoords[findex * 2] != model.u[i] {
    //             println!(
    //                 "unequal u at {}, with: {} vs {}",
    //                 i,
    //                 models[0].mesh.positions[findex * 2],
    //                 model.u[i]
    //             );
    //             return ();
    //         }

    //         if models[0].mesh.texcoords[findex * 2 + 1] != model.v[i] {
    //             println!(
    //                 "unequal v at {}, with: {} vs {}",
    //                 i,
    //                 models[0].mesh.texcoords[findex * 2 + 1],
    //                 model.v[i]
    //             );
    //             return ();
    //         }

    //         i += 1;
    //     }

    //     next_face = end;
    // }

    // println!("all vertex values are equal");

    // {
    //     let albedo = image::open("./test_assets/chalet.jpg").unwrap();

    //     let rgba = albedo.to_rgba();

    //     // write width and height
    //     if model.albedo.width != rgba.width() {
    //         println!("width not equal on albedo");
    //         return ();
    //     }

    //     if model.albedo.height != rgba.height() {
    //         println!("height not equal on albedo");
    //         return ();
    //     }

    //     let raw = rgba.into_raw();
    //     for i in 0..raw.len() {
    //         if raw[i] != model.albedo.raw_pixels[i] {
    //             println!("albedo not equal");
    //             return ();
    //         }
    //     }
    // }

    // println!("all texture values are equal");
}

// function Cut(mesh M)
//    Remove seed triangle.
//    while there remains an edge e adjacent to only one triangle t
//    Remove e and t.
//    while there remains a vertex v adjacent to only one edge e
//    Remove v and e.
//    Cut p := remaining edges and vertices.
//    if only a single vertex remains in p then
//    Add back two adjacent edges to p.

// 3.15 GB (3,390,337,024 bytes)

// So, say an object that is 10 feet tall is 100 feet away. If I hold up a ruler 3 feet away, then the object in the
// distance would correspond to about how many inches? => x/3 = 10/100

fn main() {
    let obj_file = "./test_assets/bunny.obj".to_string();

    gim::writer::write(&obj_file);

    return ();
}

use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::mem;
use std::os::windows::prelude::FileExt;

// enum ClusterCollection {
//     map(HashMap<u32, Cluster>),
//     vec(Vec<Cluster>),
// }

pub struct CTree {
    pub nodes: Vec<Node>,
    pub clusters: Vec<Cluster>,
    pub parents: [u32; 2],
    file: File,
}

pub struct Cluster {
    pub pos: Vec<f32>,
    pub idx: Vec<u32>,
}

#[repr(C)]
pub struct Node {
    pub offset: u32,
    pub children: [u32; 4],
}

pub fn read(file_name: String) -> CTree {
    let mut read_file = File::open(file_name).unwrap();

    let mut u32_buff = [0 as u8; mem::size_of::<u32>()];
    read_file.read_exact(&mut u32_buff).unwrap();
    let nodes_length = u32::from_le_bytes(u32_buff) as usize;

    assert!(mem::size_of::<Node>() == mem::size_of::<u32>() * 5);

    let buffer = Vec::<u8>::with_capacity(mem::size_of::<Node>() * nodes_length as usize);
    read_file.read_exact(&mut buffer).unwrap();
    let ptr = buffer.as_mut_ptr();
    let cap = buffer.capacity();
    let len = buffer.len();
    mem::forget(buffer); // Avoid calling the destructor!

    let nodes: Vec<Node> =
        unsafe { Vec::from_raw_parts(ptr as *mut Node, nodes_length, nodes_length) };

    let mut u32_buff = [0 as u8; mem::size_of::<u32>()];
    read_file.read_exact(&mut u32_buff).unwrap();
    let clusters_length = u32::from_le_bytes(u32_buff) as usize;

    // SeekFrom::Start(42) for seek()
    read_file.seek_read(buf, offset);

    let parents = [u32; 2];

    CTree {
        nodes,
        clusters,
        parents,
    }
}

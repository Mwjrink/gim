use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::mem;
use std::os::windows::prelude::FileExt;

pub struct CTree {
    pub nodes: Vec<Node>,
    pub cluster_cache: Vec<Cluster>,
    pub parents: [u32; 2],
    file: File,
}

#[derive(Clone)]
pub struct Cluster {
    pub pos: Vec<f32>,
    pub idx: Vec<u32>,
}

#[repr(C)]
pub struct Node {
    pub offset: u32,
    pub children: [u32; 4],
}

pub fn read(file_name: &str) -> CTree {
    let mut read_file = File::open(file_name).unwrap();

    let mut u32_buff = [0 as u8; mem::size_of::<u32>()];
    read_file.read_exact(&mut u32_buff).unwrap();
    let nodes_length = u32::from_le_bytes(u32_buff) as usize;

    read_file
        .seek(std::io::SeekFrom::Start(mem::size_of::<u32>() as u64))
        .unwrap();

    assert!(mem::size_of::<Node>() == mem::size_of::<u32>() * 5);

    let mut buffer = Vec::<u8>::with_capacity(mem::size_of::<Node>() * nodes_length as usize);
    unsafe { buffer.set_len(mem::size_of::<Node>() * nodes_length as usize) };
    read_file.read_exact(&mut buffer).unwrap();
    let ptr = buffer.as_mut_ptr();
    // let cap = buffer.capacity();
    // let len = buffer.len();
    mem::forget(buffer); // Avoid calling the destructor!

    let nodes: Vec<Node> =
        unsafe { Vec::from_raw_parts(ptr as *mut Node, nodes_length, nodes_length) };

    // let mut u32_buff = [0 as u8; mem::size_of::<u32>()];
    // read_file.read_exact(&mut u32_buff).unwrap();
    // let clusters_length = u32::from_le_bytes(u32_buff) as usize;

    println!("clusters_length: {}", nodes_length);

    let mut cluster_cache = Vec::with_capacity(nodes_length);
    cluster_cache.resize(
        nodes_length,
        Cluster {
            pos: Vec::with_capacity(0),
            idx: vec![],
        },
    );

    // TODO find the upper most parents?
    let parents = [0 as u32; 2];

    CTree {
        nodes,
        cluster_cache,
        parents,
        file: read_file,
    }
}

impl CTree {
    pub fn get_cluster(&mut self, id: &u32) -> Option<&Cluster> {
        println!("cluster_cache len: {}", self.cluster_cache.len());
        if *id < self.cluster_cache.len() as u32 {
            if !self.cluster_cache[*id as usize].pos.is_empty() {
                Some(&self.cluster_cache[*id as usize])
            } else {
                let mut offset = self.nodes[*id as usize].offset as u64;
                let positions = {
                    let positions_length = read_u32_offset(&mut self.file, &mut offset) as usize;
                    println!("positions_length: {}", positions_length);

                    let buff_size = mem::size_of::<f32>() * positions_length;
                    let mut buffer = Vec::<u8>::with_capacity(buff_size);
                    self.file.seek_read(&mut buffer, offset).unwrap();
                    offset += buff_size as u64;
                    let ptr = buffer.as_mut_ptr();
                    // let cap = buffer.capacity();
                    // let len = buffer.len();
                    mem::forget(buffer); // Avoid calling the destructor!

                    unsafe {
                        Vec::from_raw_parts(ptr as *mut f32, positions_length, positions_length)
                    }
                };

                let indices = {
                    let indices_length = read_u32_offset(&mut self.file, &mut offset) as usize;
                    println!("indices_length: {}", indices_length);

                    let buff_size = mem::size_of::<u32>() * indices_length;
                    let mut buffer = Vec::<u8>::with_capacity(buff_size);
                    self.file.seek_read(&mut buffer, offset).unwrap();
                    // offset += buff_size as u64;
                    let ptr = buffer.as_mut_ptr();
                    // let cap = buffer.capacity();
                    // let len = buffer.len();
                    mem::forget(buffer); // Avoid calling the destructor!

                    unsafe { Vec::from_raw_parts(ptr as *mut u32, indices_length, indices_length) }
                };

                let read_cluster = Cluster {
                    pos: positions,
                    idx: indices,
                };

                self.cluster_cache[*id as usize] = read_cluster;
                Some(&self.cluster_cache[*id as usize])
            }
        } else {
            None
        }
    }

    pub fn get_node(&self, id: &usize) -> &Node {
        &self.nodes[*id]
    }

    pub fn get_cluster_w_offset(&mut self, id: &u32) -> Option<(&Cluster, u32)> {
        // println!("cluster_cache len: {}", self.cluster_cache.len());
        if *id < self.cluster_cache.len() as u32 {
            if !self.cluster_cache[*id as usize].pos.is_empty() {
                Some((
                    &self.cluster_cache[*id as usize],
                    self.nodes[*id as usize].offset,
                ))
            } else {
                let mut offset = self.nodes[*id as usize].offset as u64;
                // println!("offset before: {}", offset);
                let positions = {
                    let positions_length = read_u32_offset(&mut self.file, &mut offset) as usize;
                    // println!("positions_length: {}", positions_length);

                    let buff_size = mem::size_of::<f32>() * positions_length;
                    let mut buffer = Vec::<u8>::with_capacity(buff_size);
                    unsafe { buffer.set_len(buff_size) };
                    self.file.seek(SeekFrom::Start(offset)).unwrap();
                    self.file.read_exact(&mut buffer).unwrap();
                    offset += buff_size as u64;
                    let ptr = buffer.as_mut_ptr();
                    mem::forget(buffer); // Avoid calling the destructor!

                    unsafe {
                        Vec::from_raw_parts(ptr as *mut f32, positions_length, positions_length)
                    }
                };

                let indices = {
                    let indices_length = read_u32_offset(&mut self.file, &mut offset) as usize;
                    // println!("indices_length: {}", indices_length);

                    let buff_size = mem::size_of::<u32>() * indices_length;
                    let mut buffer = Vec::<u8>::with_capacity(buff_size);
                    unsafe { buffer.set_len(buff_size) };
                    self.file.seek(SeekFrom::Start(offset)).unwrap();
                    self.file.read_exact(&mut buffer).unwrap();
                    // offset += buff_size as u64;
                    let ptr = buffer.as_mut_ptr();
                    mem::forget(buffer); // Avoid calling the destructor!

                    unsafe { Vec::from_raw_parts(ptr as *mut u32, indices_length, indices_length) }
                };

                // println!("offset after: {}", offset);

                let read_cluster = Cluster {
                    pos: positions,
                    idx: indices,
                };

                self.cluster_cache[*id as usize] = read_cluster;
                Some((
                    &self.cluster_cache[*id as usize],
                    self.nodes[*id as usize].offset,
                ))
            }
        } else {
            None
        }
    }

    pub fn len(&self) -> usize {
        self.nodes.len()
    }
}

fn read_u32(file: &mut File) -> u32 {
    let mut u32_buff = [0 as u8; mem::size_of::<u32>()];
    file.read_exact(&mut u32_buff).unwrap();
    u32::from_le_bytes(u32_buff)
}

fn read_u32_offset(file: &mut File, offset: &mut u64) -> u32 {
    let mut u32_buff = [0 as u8; mem::size_of::<u32>()];
    file.seek_read(&mut u32_buff, *offset).unwrap();
    *offset += mem::size_of::<u32>() as u64;
    u32::from_le_bytes(u32_buff)
}

// fn read_vec<T>(file: &File, length: usize) {
//     let buffer = Vec::<u8>::with_capacity(mem::size_of::<T>() * length as usize);
//     file.read_exact(&mut buffer).unwrap();
//     let ptr = buffer.as_mut_ptr();
//     let cap = buffer.capacity();
//     let len = buffer.len();
//     mem::forget(buffer); // Avoid calling the destructor!
//
//     let nodes: Vec<T> = unsafe { Vec::from_raw_parts(ptr as *mut T, length, length) };
// }

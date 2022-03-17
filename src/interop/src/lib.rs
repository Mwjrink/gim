#[derive(Clone)]
#[repr(C)]
pub struct Cluster {
    pub pos: Vec<f32>,
    pub idx: Vec<u32>,
}

#[repr(C)]
pub struct Node {
    pub offset: u32,
    pub children: [u32; 4],
}

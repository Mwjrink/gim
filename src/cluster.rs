use crate::mesh::*;
use std::collections::HashSet;

#[derive(Clone)]
pub struct Cluster {
    // TODO make this have an internal Mesh object
    pub positions: Vec<f32>,
    pub indices: Vec<u32>,
    // -------------------------------------------
    pub cut: HashSet<Edge>,
    pub anchor: Point,
}

// TODO pub(crate) instead of all pub
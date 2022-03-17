use std::io::Write;
use std::{fs::File, mem};

use crate::{
    cluster::*,
    mesh::{simplify_indices_positions_cut, Mesh},
};

pub struct CTree {
    pub nodes: Vec<CTreeNode>,
    pub offsets: Vec<usize>,
    pub clusters: Vec<Cluster>,
}

pub struct CTreeNode {
    pub cluster: u32,
    pub parents: [u32; 2],
    pub children: [u32; 4],
    // is this needed in the final tree or just for TempCTree, should there be a TempCTreeNode?
    pub lod: u32,

    // TEMP
    pub offset: u32,
}

impl CTree {
    pub fn new() -> CTree {
        CTree {
            nodes: Vec::new(),
            clusters: Vec::new(),
            offsets: vec![],
        }
    }

    pub fn from_raw(nodes: Vec<CTreeNode>, clusters: Vec<Cluster>) -> CTree {
        CTree {
            nodes,
            clusters,
            offsets: vec![],
        }
    }

    /// You can assume the id of the cluster is its index in the provided vector
    pub fn from_clusters(clusters: Vec<Cluster>) -> CTree {
        let mut nodes = Vec::with_capacity(clusters.len());
        for idx in 0..clusters.len() {
            nodes.push(CTreeNode {
                cluster: idx as u32,
                parents: [u32::MAX, u32::MAX],
                children: [u32::MAX, u32::MAX, u32::MAX, u32::MAX],
                lod: 0,
                offset: 0,
            });
        }

        CTree {
            nodes,
            clusters,
            offsets: vec![],
        }
    }

    pub fn insert(
        &mut self,
        cluster: Cluster,
        parent_ids: &[u32; 2],
        children_ids: &[u32; 4],
    ) -> u32 {
        let id = self.clusters.len() as u32;
        self.clusters.push(cluster);
        let mut lod_map: Vec<u32> = children_ids
            .iter()
            .map(|idx| self.nodes[*idx as usize].lod)
            .collect();
        lod_map.sort();
        self.nodes.push(CTreeNode {
            cluster: id,
            parents: *parent_ids,
            children: *children_ids,
            lod: *lod_map.last().unwrap_or(&0),
            offset: 0,
        });

        for child_id in children_ids {
            if *child_id != u32::MAX {
                let parents = &mut self.nodes[*child_id as usize].parents;
                if parents[0] == u32::MAX {
                    parents[0] = id;
                } else if parents[1] == u32::MAX {
                    parents[1] = id;
                } else {
                    panic!("more than one parent for child: {}", child_id);
                }
            }
        }

        id
    }

    // TODO return an option and clusters.get(id) instead
    pub fn get(&self, id: &u32) -> &Cluster {
        &self.clusters[*id as usize]
    }

    pub fn get_w_offset(&self, id: &u32) -> (&Cluster, usize) {
        assert!(self.nodes[*id as usize].cluster == *id);
        (&self.clusters[*id as usize], self.offsets[*id as usize])
    }

    pub fn get_lod(&self, id: &u32) -> u32 {
        self.nodes[*id as usize].lod
    }

    // TODO return an option and clusters.get(id) instead
    pub fn get_mut(&mut self, id: &u32) -> &mut Cluster {
        &mut self.clusters[*id as usize]
    }

    pub fn len(&self) -> usize {
        self.clusters.len()
    }

    pub fn clusters(&self) -> &Vec<Cluster> {
        &self.clusters
    }

    pub fn get_node(&self, id: &usize) -> &CTreeNode {
        &self.nodes[*id]
    }

    pub fn write_to_file(&mut self, file_name: &str) {
        let mut cluster_buffer = Vec::with_capacity(
            self.clusters.len()
                * (mem::size_of::<u32>() * 128 * 3 + mem::size_of::<f32>() * 128 * 2),
        );
        let mut node_buffer =
            Vec::with_capacity(self.nodes.len() * mem::size_of::<interop::Node>());
        let mut offsets = Vec::with_capacity(self.clusters.len());

        let mut offset = 0;
        // write the clusters themselves
        for cluster in &self.clusters {
            /*
            cluster: u32,
            parents: [u32; 2], // dont write this?
            children: [u32; 4],
            lod: u32, // dont write this?
            */
            // positions list can be a raw list of f32, dump the raw bits from the vec
            cluster_buffer.extend_from_slice(&(cluster.mesh.positions.len() as u32).to_le_bytes());
            for pos in &cluster.mesh.positions {
                cluster_buffer.extend_from_slice(&pos.to_le_bytes());
            }

            // indices list can be 14bits for the first value then 5/5 bits? for the other two?
            cluster_buffer.extend_from_slice(&(cluster.mesh.indices.len() as u32).to_le_bytes());
            for idx in &cluster.mesh.indices {
                cluster_buffer.extend_from_slice(&idx.to_le_bytes());
            }

            offsets.push(offset);
            offset = cluster_buffer.len();
        }

        // TODO should these be sorted so the highest parents are the first ones then you can split from there?
        let cluster_offset =
            self.nodes.len() * mem::size_of::<interop::Node>() + mem::size_of::<u32>();
        for node in &mut self.nodes {
            /*
            cluster: u32,
            parents: [u32; 2], // dont write this?
            children: [u32; 4],
            lod: u32, // dont write this?
            */
            let offset = (offsets[node.cluster as usize] + cluster_offset) as u32;
            node.offset = offset;
            node_buffer.extend_from_slice(&offset.to_le_bytes());
            // write the nodes themselves, write the indices as file pointers/offsets?
            // TODO write the clusters to a temp buffer first and keep track of the offsets of each cluster, then use the offsets as the id here/pointer/fseek-offset for reading
            node_buffer.extend_from_slice(&node.children[0].to_le_bytes());
            node_buffer.extend_from_slice(&node.children[1].to_le_bytes());
            node_buffer.extend_from_slice(&node.children[2].to_le_bytes());
            node_buffer.extend_from_slice(&node.children[3].to_le_bytes());

            // raw parts?
        }

        self.offsets = offsets;

        let mut write_file = File::create(format!("./output/{}", file_name)).unwrap();
        
        // write number of nodes
        write_file
            .write_all(&(self.nodes.len() as u32).to_le_bytes())
            .unwrap();

        write_file.write_all(&node_buffer).unwrap();

        assert!(cluster_offset == node_buffer.len() + mem::size_of::<u32>());

        // write the number of clusters
        // write_file
        //     .write_all(&(self.clusters.len() as u32).to_le_bytes())
        //     .unwrap();

        write_file.write_all(&cluster_buffer).unwrap();
    }
}

pub struct TempCTree {
    pub nodes: Vec<CTreeNode>,
    pub clusters: Vec<TempCluster>,
}

impl TempCTree {
    pub fn new() -> TempCTree {
        TempCTree {
            nodes: Vec::new(),
            clusters: Vec::new(),
        }
    }

    /// You can assume the id of the cluster is its index in the provided vector
    pub fn from_clusters(clusters: Vec<TempCluster>) -> TempCTree {
        let mut nodes = Vec::with_capacity(clusters.len());
        for idx in 0..clusters.len() {
            nodes.push(CTreeNode {
                cluster: idx as u32,
                parents: [u32::MAX, u32::MAX],
                children: [u32::MAX, u32::MAX, u32::MAX, u32::MAX],
                lod: 0,
                offset: 0,
            });
        }

        TempCTree { nodes, clusters }
    }

    pub fn insert(
        &mut self,
        cluster: TempCluster,
        parent_ids: &[u32; 2],
        children_ids: &[u32; 4],
    ) -> u32 {
        let id = self.clusters.len() as u32;
        self.clusters.push(cluster);
        let mut lod_map: Vec<u32> = children_ids
            .iter()
            .map(|idx| self.nodes[*idx as usize].lod)
            .collect();
        lod_map.sort();
        let lod = if let Some(lod_level) = lod_map.last() {
            lod_level + 1
        } else {
            0
        };
        self.nodes.push(CTreeNode {
            cluster: id,
            parents: *parent_ids,
            children: *children_ids,
            lod,
            offset: 0,
        });

        for child_id in children_ids {
            if *child_id != u32::MAX {
                let parents = &mut self.nodes[*child_id as usize].parents;
                if parents[0] == u32::MAX {
                    parents[0] = id;
                } else if parents[1] == u32::MAX {
                    parents[1] = id;
                } else {
                    panic!("more than one parent for child: {}", child_id);
                }
            }
        }

        id
    }

    // TODO return an option and clusters.get(id) instead
    pub fn get(&self, id: &u32) -> &TempCluster {
        &self.clusters[*id as usize]
    }

    // TODO return an option and clusters.get(id) instead
    pub fn get_mut(&mut self, id: &u32) -> &mut TempCluster {
        &mut self.clusters[*id as usize]
    }

    pub fn len(&self) -> usize {
        self.clusters.len()
    }

    pub fn clusters(self) -> Vec<TempCluster> {
        self.clusters
    }

    pub fn to_official(self, positions: &Vec<f32>) -> CTree {
        let mut clusters = Vec::with_capacity(self.clusters.len());
        for cluster in self.clusters {
            let (indices, positions, cut) =
                simplify_indices_positions_cut(&cluster.mesh.indices, positions, &cluster.cut);
            clusters.push(Cluster {
                mesh: Mesh { positions, indices },
                cut,
                anchor: cluster.anchor,
            });
        }

        CTree::from_raw(self.nodes, clusters)
    }

    pub fn get_lod(&self, id: &u32) -> u32 {
        self.nodes[*id as usize].lod
    }

    pub fn has_parent(&self, id: &u32) -> bool {
        let parents = self.nodes[*id as usize].parents;
        parents[0] != u32::MAX || parents[1] != u32::MAX
    }
}

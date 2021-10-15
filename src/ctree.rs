use crate::{
    cluster::*,
    mesh::{simplify_indices_positions_cut, Mesh},
};

pub struct CTree {
    pub nodes: Vec<CTreeNode>,
    pub clusters: Vec<Cluster>,
}

pub struct CTreeNode {
    pub cluster: u32,
    pub parents: [u32; 2],
    pub children: [u32; 4],
}

impl CTree {
    pub fn new() -> CTree {
        CTree {
            nodes: Vec::new(),
            clusters: Vec::new(),
        }
    }

    pub fn from_raw(nodes: Vec<CTreeNode>, clusters: Vec<Cluster>) -> CTree {
        CTree { nodes, clusters }
    }

    /// You can assume the id of the cluster is its index in the provided vector
    pub fn from_clusters(clusters: Vec<Cluster>) -> CTree {
        let mut nodes = Vec::with_capacity(clusters.len());
        for idx in 0..clusters.len() {
            nodes.push(CTreeNode {
                cluster: idx as u32,
                parents: [u32::MAX, u32::MAX],
                children: [u32::MAX, u32::MAX, u32::MAX, u32::MAX],
            });
        }

        CTree { nodes, clusters }
    }

    pub fn insert(
        &mut self,
        cluster: Cluster,
        parent_ids: &[u32; 2],
        children_ids: &[u32; 4],
    ) -> u32 {
        let id = self.clusters.len() as u32;
        self.clusters.push(cluster);
        self.nodes.push(CTreeNode {
            cluster: id,
            parents: *parent_ids,
            children: *children_ids,
        });

        id
    }

    // TODO return an option and clusters.get(id) instead
    pub fn get(&self, id: &u32) -> &Cluster {
        &self.clusters[*id as usize]
    }

    // TODO return an option and clusters.get(id) instead
    pub fn get_mut(&mut self, id: &u32) -> &mut Cluster {
        &mut self.clusters[*id as usize]
    }

    pub fn len(&self) -> usize {
        self.clusters.len()
    }

    pub fn clusters(self) -> Vec<Cluster> {
        self.clusters
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
        self.nodes.push(CTreeNode {
            cluster: id,
            parents: *parent_ids,
            children: *children_ids,
        });

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
}

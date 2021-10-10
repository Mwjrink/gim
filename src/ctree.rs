use crate::cluster::*;

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

    pub fn get(&self, id: &u32) -> &Cluster {
        &self.clusters[*id as usize]
    }
}

# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os

class NetworkDWPC:

    def __init__(self, G: nx.MultiGraph, damping_factor: float=-0.4):
        """
        Initialize DWPC calculator for a NetworkX MultiGraph.

        Args:
            G: NetworkX MultiGraph with 'node_type' node attribute and 'predicate' edge attribute
            damping_factor: Power to raise degrees to (default -0.4)
        """
        self.G = G
        self.damping_factor = damping_factor
        self.metapath_finder = MetapathFinder(G)
        self.node_to_idx = {node: idx for idx, node in enumerate(G.nodes())}
        self.idx_to_node = {idx: node for node, idx in self.node_to_idx.items()}
        self.edge_types = set((data['predicate'] for _, _, data in G.edges(data=True)))
        print(f'Found {len(self.edge_types)} edge types: {self.edge_types}')
        self.edge_matrices = self._create_edge_matrices()

    def find_metapaths(self, source_type: str, target_type: str, max_length: int=2) -> List[List[str]]:
        """Use optimized metapath finder."""
        return self.metapath_finder.find_metapaths(source_type, target_type, max_length)

    def _create_edge_matrices(self) -> Dict[str, sparse.csr_matrix]:
        """Create sparse adjacency matrix for each edge type."""
        n = len(self.G)
        matrices = {}
        for edge_type in self.edge_types:
            matrices[edge_type] = sparse.lil_matrix((n, n))
        for u, v, data in self.G.edges(data=True):
            edge_type = data['predicate']
            u_idx = self.node_to_idx[u]
            v_idx = self.node_to_idx[v]
            matrices[edge_type][u_idx, v_idx] = 1
            if not self.G.is_directed():
                matrices[edge_type][v_idx, u_idx] = 1
        return {k: v.tocsr() for k, v in matrices.items()}

    def get_node_type(self, node: str) -> str:
        """Get type of a node."""
        return self.G.nodes[node]['node_type']

    def get_edge_types_between_nodetypes(self, source_type: str, target_type: str) -> Set[str]:
        """Find all edge types that connect nodes of given types."""
        edge_types = set()
        for u, v, data in self.G.edges(data=True):
            if self.get_node_type(u) == source_type and self.get_node_type(v) == target_type:
                edge_types.add(data['predicate'])
        return edge_types

    def get_node_degree(self, node: str, edge_type: str, direction: str='out') -> int:
        """
        Get degree of a node for specific edge type and direction.

        Args:
            node: Node identifier
            edge_type: Type of edge to count
            direction: 'out' for outgoing edges, 'in' for incoming edges
        """
        count = 0
        if direction == 'out':
            for _, _, data in self.G.edges(node, data=True):
                if data['predicate'] == edge_type:
                    count += 1
        else:
            for u, v, data in self.G.edges(data=True):
                if v == node and data['predicate'] == edge_type:
                    count += 1
        return count

    def _create_degree_matrices(self, metapath: List[str]) -> List[Tuple[sparse.dia_matrix, sparse.dia_matrix]]:
        """Create degree matrices for each step in the metapath."""
        n = len(self.G)
        degree_matrices = []
        for edge_type in metapath:
            adj_matrix = self.edge_matrices[edge_type]
            out_degrees = np.array(adj_matrix.sum(axis=1)).flatten()
            in_degrees = np.array(adj_matrix.sum(axis=0)).flatten()
            out_degrees[out_degrees == 0] = 1
            in_degrees[in_degrees == 0] = 1
            damped_out = np.power(out_degrees, self.damping_factor)
            damped_in = np.power(in_degrees, self.damping_factor)
            degree_matrices.append((sparse.diags(damped_out), sparse.diags(damped_in)))
        return degree_matrices

    def compute_dwpc(self, source: str, target: str, metapath: List[str]) -> float:
        """
        Compute DWPC for a given source, target, and metapath.

        Args:
            source: Source node identifier
            target: Target node identifier
            metapath: List of edge types defining the path

        Returns:
            DWPC value
        """
        if source not in self.node_to_idx or target not in self.node_to_idx:
            return 0.0
        source_idx = self.node_to_idx[source]
        target_idx = self.node_to_idx[target]
        degree_matrices = self._create_degree_matrices(metapath)
        result = self.edge_matrices[metapath[0]]
        result = degree_matrices[0][0].dot(result)
        for i in range(len(metapath) - 1):
            result = result.dot(degree_matrices[i][1])
            result = result.dot(degree_matrices[i + 1][0])
            result = result.dot(self.edge_matrices[metapath[i + 1]])
        result = result.dot(degree_matrices[-1][1])
        return float(result[source_idx, target_idx])

    def compute_simple_path_count(self, source: str, target: str, metapath: List[str]) -> int:
        """
        Compute simple path count (without degree weighting) for a given metapath.
        """
        if source not in self.node_to_idx or target not in self.node_to_idx:
            return 0
        source_idx = self.node_to_idx[source]
        target_idx = self.node_to_idx[target]
        result = self.edge_matrices[metapath[0]]
        for edge_type in metapath[1:]:
            result = result.dot(self.edge_matrices[edge_type])
        return int(result[source_idx, target_idx])

    def get_node_degree(self, node: str, edge_type: str, direction: str='out') -> int:
        """
        Get degree of a node for specific edge type and direction.

        Args:
            node: Node identifier
            edge_type: Type of edge to count
            direction: 'out' for outgoing edges, 'in' for incoming edges
        """
        if node not in self.node_to_idx:
            return 0
        idx = self.node_to_idx[node]
        matrix = self.edge_matrices[edge_type]
        if direction == 'out':
            return int(matrix[idx].sum())
        else:
            return int(matrix[:, idx].sum())

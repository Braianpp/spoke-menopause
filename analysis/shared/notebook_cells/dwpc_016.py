# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os

class MetapathFinder:

    def __init__(self, G: nx.MultiGraph):
        """Initialize with a NetworkX MultiGraph."""
        self.G = G
        self._cache_network_structure()

    def _cache_network_structure(self):
        """Cache network structure for faster metapath finding."""
        self.node_types = set(nx.get_node_attributes(self.G, 'node_type').values())
        self.valid_transitions = defaultdict(set)
        self.edge_types = set()
        for u, v, data in self.G.edges(data=True):
            edge_type = data['predicate']
            u_type = self.G.nodes[u]['node_type']
            v_type = self.G.nodes[v]['node_type']
            self.valid_transitions[u_type, v_type].add(edge_type)
            self.edge_types.add(edge_type)
            if not self.G.is_directed():
                self.valid_transitions[v_type, u_type].add(edge_type)

    def _extend_metapath(self, current_path: List[str], current_type: str, target_type: str, length_left: int, metapaths: List[List[str]], visited_types: Set[str]=None):
        """
        Recursively extend metapath to find all valid paths.

        Args:
            current_path: Current metapath being built
            current_type: Current node type
            target_type: Target node type we want to reach
            length_left: How many more edges we can add
            metapaths: List to collect valid metapaths
            visited_types: Set of node types already in the current path
        """
        if visited_types is None:
            visited_types = {current_type}
        if current_type == target_type and length_left == 0:
            metapaths.append(list(current_path))
            return
        if length_left == 0:
            return
        for next_type in self.node_types:
            edge_types = self.valid_transitions.get((current_type, next_type), set())
            for edge_type in edge_types:
                current_path.append(edge_type)
                self._extend_metapath(current_path, next_type, target_type, length_left - 1, metapaths, visited_types | {next_type})
                current_path.pop()

    def find_metapaths(self, source_type: str, target_type: str, max_length: int=2) -> List[List[str]]:
        """
        Find all possible metapaths between given node types up to max_length.

        Args:
            source_type: Starting node type
            target_type: Ending node type
            max_length: Maximum length of metapaths to find

        Returns:
            List of metapaths (each metapath is a list of edge types)
        """
        if max_length < 1:
            return []
        metapaths = []
        for length in range(1, max_length + 1):
            self._extend_metapath(current_path=[], current_type=source_type, target_type=target_type, length_left=length, metapaths=metapaths)
        return metapaths

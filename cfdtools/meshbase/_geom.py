import cfdtools.api as api
import numpy as np
import scipy.spatial as spspa


class Nodes():
    """class of nodes cloud ans associated methods
    may be used to implement distributed data in the future
    """
    def __init__(self, nodes: np.ndarray) -> None:
        assert len(nodes.shape)==2
        self._nodes = nodes
        self._tree = None
    
    def center(self):
        return self._nodes.mean(axis=0)
    
    def __iadd__(self, other: np.ndarray):
        self._nodes += other
        return self
    
    def kdtree_init(self):
        self._tree = spspa.KDTree(self._nodes)

    def kdtree_query(self, othernodes, **kwargs):
        if self._tree is None:
            self.kdtree_init()
        return self._tree.query(othernodes._nodes, p=2, **kwargs)
    
    def __str__(self):
        return f"Nodes: {self._nodes}"
    

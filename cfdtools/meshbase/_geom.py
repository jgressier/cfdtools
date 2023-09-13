import cfdtools.api as api
import numpy as np
import scipy.spatial as spspa
from cfdtools.utils._dev import lazyprop

class Nodes():
    """class of nodes cloud ans associated methods
    may be used to implement distributed data in the future
    """
    def __init__(self, nodes: np.ndarray) -> None:
        assert len(nodes.shape)==2
        self._nodes = nodes
        self._tree = None

    @lazyprop    
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

    def rotate(self, axis: np.ndarray, angle: float):
        """perform a rotation along given axis and angle (in degrees)"""
        rot = spspa.transform.Rotation.from_rotvec(axis*angle, degrees=True)
        self._nodes = rot.apply(self._nodes)

    def estimate_rotation(self, othernodes, method='unordered'):
        """estimate of rotation self to othernodes using Kabsch algorithm
        (scipy implementation needs the same order of nodes set)

        Args:
            othernodes (_type_): _description_

        Returns:
            Rotation: (scipy) rotation matrix
            float: error related to non matching points
        """
        x = self._nodes-self.center
        y = othernodes._nodes-othernodes.center
        if method == 'ordered':
            rotm, error = spspa.transform.Rotation.align_vectors(y, x)
        elif method == 'unordered':
            _, Sx, Vhx = np.linalg.svd(x)
            _, Sy, Vhy = np.linalg.svd(y)
            # U @ diag(S) = U * S
            rotm, error = spspa.transform.Rotation.align_vectors(Vhy*Sy, Vhx*Sx) 
        else:
            api.error_stop(f"unknown method '{method}' in Nodes.estimate_rotation")
        return rotm, error

    def __str__(self):
        return f"Nodes: {self._nodes}"
    

if __name__ == "__main__":
    n1 = Nodes(np.array([[5., -8., 10.], [1., 2., 10.], [3., 2., 10.], ]))
    n2 = Nodes(np.array([[-2., 1., 10.], [-2., 3., 10.], [8., 5., 10.]]))
    rotm, error = n1.estimate_rotation(n2)
    print(error, rotm.as_rotvec(degrees=True))

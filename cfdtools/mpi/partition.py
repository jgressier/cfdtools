"""Local mesh plus the global identifiers required for distributed work."""

from dataclasses import dataclass

import numpy as np


def _as_global_ids(values, expected_size, name):
    array = np.asarray(values)
    if array.ndim != 1 or array.size != expected_size:
        raise ValueError(
            "{} must be a one-dimensional array of length {}, got shape {}".format(
                name, expected_size, array.shape
            )
        )
    if not np.issubdtype(array.dtype, np.integer):
        raise TypeError("{} must contain integer global identifiers".format(name))
    if np.any(array < 0):
        raise ValueError("{} must contain non-negative identifiers".format(name))
    if np.unique(array).size != array.size:
        raise ValueError("{} must be unique on one rank".format(name))
    return array.astype(np.int64, copy=False)


@dataclass
class PartitionedMesh:
    """A local mesh view with immutable global node and cell identifiers.

    A cell belongs to exactly one rank. Nodes may be present as halos on
    multiple ranks, but every local occurrence uses its global node ID.
    """

    local_mesh: object
    global_node_ids: np.ndarray
    global_cell_ids: np.ndarray
    comm: object = None
    owned_node_mask: np.ndarray = None

    def __post_init__(self):
        from cfdtools.mpi.comm import get_comm

        self.comm = get_comm(self.comm)
        self.global_node_ids = _as_global_ids(
            self.global_node_ids, self.local_mesh.nnode, "global_node_ids"
        )
        self.global_cell_ids = _as_global_ids(
            self.global_cell_ids, self.local_mesh.ncell, "global_cell_ids"
        )
        if self.owned_node_mask is None:
            self.owned_node_mask = np.ones(self.local_mesh.nnode, dtype=bool)
        else:
            mask = np.asarray(self.owned_node_mask, dtype=bool)
            if mask.ndim != 1 or mask.size != self.local_mesh.nnode:
                raise ValueError(
                    "owned_node_mask must be one-dimensional with length {}".format(
                        self.local_mesh.nnode
                    )
                )
            self.owned_node_mask = mask

    @property
    def rank(self):
        return self.comm.rank

    @property
    def size(self):
        return self.comm.size

    @property
    def n_owned_nodes(self):
        return int(np.count_nonzero(self.owned_node_mask))

    @classmethod
    def from_serial_mesh(cls, mesh, comm=None):
        """Wrap an unpartitioned mesh; valid only for a one-rank communicator."""

        from cfdtools.mpi.comm import get_comm

        communicator = get_comm(comm)
        if communicator.size != 1:
            raise ValueError(
                "from_serial_mesh() cannot invent global IDs for multiple ranks; "
                "provide an explicit PartitionedMesh instead"
            )
        return cls(
            local_mesh=mesh,
            global_node_ids=np.arange(mesh.nnode, dtype=np.int64),
            global_cell_ids=np.arange(mesh.ncell, dtype=np.int64),
            comm=communicator,
        )

import numpy as np
import pytest

from cfdtools.meshbase.simple import Cube
from cfdtools.mpi import (
    FaceCandidate,
    PartitionedMesh,
    SerialComm,
    canonical_face_key,
    face_owner,
    get_comm,
)


def test_serial_comm_collectives():
    comm = SerialComm()
    assert (comm.rank, comm.size) == (0, 1)
    assert comm.bcast("mesh") == "mesh"
    assert comm.gather("mesh") == ["mesh"]
    assert comm.scatter(["mesh"]) == "mesh"
    assert comm.allgather("mesh") == ["mesh"]
    assert comm.allreduce(3) == 3
    assert comm.alltoall(["mesh"]) == ["mesh"]


def test_get_comm_keeps_explicit_communicator():
    comm = SerialComm()
    assert get_comm(comm) is comm


def test_partitioned_mesh_wraps_serial_mesh():
    mesh = Cube(2, 1, 1).export_mesh()
    distributed = PartitionedMesh.from_serial_mesh(mesh, comm=SerialComm())
    assert distributed.n_owned_nodes == mesh.nnode
    assert np.array_equal(distributed.global_node_ids, np.arange(mesh.nnode))
    assert np.array_equal(distributed.global_cell_ids, np.arange(mesh.ncell))


def test_partitioned_mesh_rejects_duplicate_local_global_ids():
    mesh = Cube(1, 1, 1).export_mesh()
    with pytest.raises(ValueError, match="global_node_ids must be unique"):
        PartitionedMesh(
            mesh,
            global_node_ids=np.zeros(mesh.nnode, dtype=np.int64),
            global_cell_ids=np.arange(mesh.ncell, dtype=np.int64),
            comm=SerialComm(),
        )


def test_face_key_is_orientation_independent_and_owner_is_stable():
    key = canonical_face_key("quad4", [9, 2, 7, 4])
    assert key == canonical_face_key("quad4", [4, 7, 2, 9])
    assert face_owner(key, 4) == face_owner(key, 4)
    candidate = FaceCandidate("quad4", (9, 2, 7, 4), 12, 3, (9, 2, 7, 4))
    assert candidate.key == key

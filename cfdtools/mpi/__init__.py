"""MPI-ready building blocks with a dependency-free serial fallback.

The package deliberately does not distribute a meshbase.Mesh yet. It provides
the communication and global-ID contracts required before doing so.
"""

from cfdtools.mpi.comm import SerialComm, get_comm, mpi_available, require_mpi
from cfdtools.mpi.connectivity import FaceCandidate, canonical_face_key, face_owner
from cfdtools.mpi.partition import PartitionedMesh

__all__ = [
    "FaceCandidate",
    "PartitionedMesh",
    "SerialComm",
    "canonical_face_key",
    "face_owner",
    "get_comm",
    "mpi_available",
    "require_mpi",
]

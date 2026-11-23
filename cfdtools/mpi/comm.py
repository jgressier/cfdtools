"""Small communication abstraction for optional MPI support.

Algorithms receive a communicator explicitly. When mpi4py is unavailable, the
serial implementation keeps the same code path usable in ordinary tests.
"""


class SerialComm:
    """Subset of the mpi4py communicator API for one process."""

    rank = 0
    size = 1

    def barrier(self):
        """Synchronize the single participant."""

    def bcast(self, value, root=0):
        self._check_root(root)
        return value

    def gather(self, value, root=0):
        self._check_root(root)
        return [value]

    def scatter(self, values, root=0):
        self._check_root(root)
        if len(values) != 1:
            raise ValueError("SerialComm.scatter() expects exactly one value")
        return values[0]

    def allgather(self, value):
        return [value]

    def allreduce(self, value, op=None):
        return value

    def alltoall(self, values):
        if len(values) != 1:
            raise ValueError("SerialComm.alltoall() expects exactly one value")
        return [values[0]]

    @staticmethod
    def _check_root(root):
        if root != 0:
            raise ValueError("SerialComm has only root 0")


def _mpi_module():
    try:
        from mpi4py import MPI
    except ImportError:
        return None
    return MPI


def mpi_available():
    """Return whether mpi4py can be imported in this interpreter."""

    return _mpi_module() is not None


def require_mpi():
    """Return mpi4py.MPI or explain how to install the optional dependency."""

    mpi = _mpi_module()
    if mpi is None:
        raise RuntimeError(
            "MPI support requires mpi4py. Install the optional dependency with "
            "'python -m pip install cfdtools[mpi]'."
        )
    return mpi


def get_comm(comm=None):
    """Return an explicit communicator, MPI.COMM_WORLD, or a serial fallback."""

    if comm is not None:
        return comm
    mpi = _mpi_module()
    return SerialComm() if mpi is None else mpi.COMM_WORLD

try:
    from mpi4py import MPI

    _comm = MPI.COMM_WORLD
except:
    _comm = None


def get_mpirank():
    return 0 if _comm is None else _comm.rank

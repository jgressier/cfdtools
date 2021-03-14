import cfdtools.api as api

class mesh():
    """versatile mesh object
    """
    def __init__(self, ncell, nnode):
        self.ncell = ncell
        self.nnode = nnode

    
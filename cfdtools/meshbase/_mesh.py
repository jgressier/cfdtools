import cfdtools.api as api

class mesh():
    """versatile mesh object
    """
    def __init__(self, ncell, nnode):
        self.ncell = ncell
        self.nnode = nnode
        self.nface = 0
        self._params = {}
        self._cell2node = {}
        self._face2node = {}
        self._face2cell = {}
        self._cell2node = {}
    
    def set_cell2node(self, cell2node):

    def set_face2cell(self, face2cell):

    def set_face2node(self, face2node):

    def printinfo(self):

    def check(self):
        
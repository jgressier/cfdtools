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
        self._nodes = {}
        self._bocos = {}
        self._celldata = {}
        self._nodedata = {}
        self._facedata = {}
        self._cellprop = {}
        
    def set_nodescoord_nd(self, xyz):
        for i,c in enumerate(['x', 'y', 'z']):
            self._nodes[c] = xyz[:,i]
    
    def set_cell2node(self, cell2node):
        self._cell2node = cell2node

    def set_face2cell(self, face2cell):
        self._face2cell = face2cell

    def set_face2node(self, face2node):
        self._face2node = face2node

    def set_bocos(self, bocos):
        self._bocos = bocos

    def set_params(self, params):
        self._params = params

    def update_params(self, params):
        self._params.update(params)

    def set_celldata(self, celldata):
        self._celldata = celldata

    def set_facedata(self, facedata):
        self._facedata = facedata

    def set_nodedata(self, nodedata):
        self._nodedata = nodedata

    def set_partition(self, partition):
        self._cellprop['partition'] = partition

    def printinfo(self):
        print("ncell:",self.ncell)
        print("nnode:",self.nnode)
        print("nface:",self.nface)
        print("params:",self._params)

    def check(self):
        # check cell2node and cell numbers
        return True


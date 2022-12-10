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
    
    def set_nodescoord_xyz(self, x, y, z):
        self._nodes['x'] = x
        self._nodes['y'] = y
        self._nodes['z'] = z
    
    def set_cell2node(self, cell2node: dict):
        """set cell to node connectivity as a dict

        Args:
            cell2node (dict): dict of ndarray
        """
        self._cell2node = cell2node
        self._check_cell2node()

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

    def pop_nodedata(self, name):
        return self._nodedata.pop(name) if name in self._nodedata.keys() else None

    def pop_celldata(self, name):        
        return self._celldata.pop(name) if name in self._celldata.keys() else None

    def printinfo(self):
        print("ncell:",self.ncell)
        if self._cell2node:
            for celltype, elemco in self._cell2node.items():
                print(f"  {celltype}: {elemco.shape}")
        else:
            print("  no cell/node connectivity")
        print("nnode:",self.nnode)
        print("nface:",self.nface)
        if self._face2node:
            print(f"  {self._face2node}")
        else:
            print("  no face/node connectivity")
        if self._face2cell:
            print(f"  {self._face2cell}")
        else:
            print("  no face/cell connectivity")
        print(f"bocos: {self._bocos.keys()}")
        for name, boco in self._bocos.items():
            print(f"  {name}: {boco}")
        print("params:",self._params)

    def _check_cell2node(self):
        assert isinstance(self._cell2node, dict)

    def check(self):
        # check cell2node and cell numbers
        api.io.print('std','ckeck: at least cell and node numbers')
        assert(self.ncell >0)
        assert(self.nnode >0)
        api.io.print('std','ckeck: at least cell/node or face/node face/cell connectivity')
        assert(not self._cell2node or (not self._face2node and not self._face2cell))
        self._check_cell2node()
        api.io.print('std','ckeck: done')
        return True

    def morph(self, fmorph):
        '''
        change x, y, z position with a function of x, y, z, returning new x, y, z
        '''
        newx, newy, newz = fmorph(self._nodes['x'], self._nodes['y'], self._nodes['z'])
        self.set_nodescoord_xyz(newx, newy, newz)
import cfdtools.api as api
import cfdtools.meshbase._connectivity as conn
import cfdtools.meshbase._elements as ele

class submeshmark():
    def __init__(self, name):
        self._name = name
        self.properties = {}

    @property
    def name(self):
        return self._name

    @property
    def geodim(self):
        return self._geodim

    @geodim.setter
    def set_geodim(self, geodim):
        assert geodim in ele.geomdim.keys()
        self._geodim = geodim

    @property
    def index(self):
        return self._index

    @index.setter
    def index(self, index: conn.indexlist):
        self._index = index

    @property
    def properties(self):
        return self._properties

    @properties.setter
    def properties(self, properties: dict):
        self._properties = properties

class mesh():
    """versatile mesh object
    """
    __available_facetypes = ( 'mixed', 'internal', 'boundaries')

    def __init__(self, ncell, nnode):
        self.ncell = ncell
        self.nnode = nnode
        self.nface = 0
        self._params = {}
        self._cell2node = {}
        self._faces = {}
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
    
    def set_cell2node(self, cell2node):
        """set cell to node connectivity as a dict

        Args:
            cell2node (dict): dict of ndarray
        """
        #print("cell2node : ", cell2node)
        self._cell2node = cell2node
        self._check_cell2node()

    def set_faces(self, facetype: str, face2node: conn.elem_connectivity, face2cell: conn.indexindirection = None):
        """set faces connectivity with face type et optional face/cell connectivity

        Args:
            facetype (str): _description_
            face2node (conn.elem_connectivity): _description_
            face2cell (conn.indexindirection, optional): _description_. Defaults to None.
        """
        if facetype in self.__available_facetypes:
            self._faces[facetype] = {'face2node' : face2node, 'face2cell': face2cell}
        else:
            api.io.error_stop(f"bad face type: {facetype} since {self.__available_facetypes} expected")

    def add_boco(self, boco: submeshmark):
        self._bocos[boco.name] = boco

    def remove_boco(self, name):
        self._bocos.pop(name)

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
        api.io.print("std", f"ncell: {self.ncell}")
        if self._cell2node:
            self._cell2node.print()
        else:
            api.io.print("std", "  no cell/node connectivity")
        api.io.print("nnode:",self.nnode)
        api.io.print("nface:",self.nface)
        if self._faces:
            for t, facedict in self._faces.items():
                api.io.print("std", f"  type {t}")
        else:
            api.io.print("std", "  no face/node connectivity")
        api.io.print(f"bocos: {self._bocos.keys()}")
        for name, boco in self._bocos.items():
            if isinstance(boco, dict):
                api.io.print("std", f"  BC {name}: {boco.keys()}")
            else:
                api.io.print("std", f"  BC {name}: {boco}")
        api.io.print("std", "params:",self._params)

    def _check_cell2node(self):
        assert isinstance(self._cell2node, conn.elem_connectivity)
        #for etype, conn in self._cell2node.items():
        #    assert etype in ele.elem2faces.keys()
        return True

    def _make_face_connectivity(self):
        faces_elem = self._cell2node.create_faces_from_elems()
        intfaces, intf2c, boundfaces, boundf2c = conn.find_duplicates(faces_elem)
        #print(boundfaces, boundf2c.conn)

    def check(self):
        # check cell2node and cell numbers
        #api.io.print('std','ckeck: at least cell and node numbers')
        assert(self.ncell >0)
        assert(self.nnode >0)
        #api.io.print('std','ckeck: at least cell/node or face/node face/cell connectivity')
        #assert(not self._cell2node or (not self._face2node and not self._face2cell))
        #assert self._check_cell2node() # not compulsory
        #api.io.print('std','ckeck: done')
        return True

    def morph(self, fmorph):
        '''
        change x, y, z position with a function of x, y, z, returning new x, y, z
        '''
        newx, newy, newz = fmorph(self._nodes['x'], self._nodes['y'], self._nodes['z'])
        self.set_nodescoord_xyz(newx, newy, newz)

    
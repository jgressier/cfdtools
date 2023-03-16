import cfdtools.api as api
import cfdtools.meshbase._connectivity as conn
import cfdtools.meshbase._elements as ele
from cfdtools.utils.math import minavgmax
import itertools
import numpy as np

class submeshmark():

    # authorized geomdim type and actual dimension
    _available_geomdim = (
        'node', 'intnode', 'bdnode',
        'face', 'intface', 'bdface',
        'cell' )

    _available_types = (
        'internal',
        'boundary',
        'perio_cart',
        'perio_cylx',
        'perio_cyly',
        'perio_cylz',
    )

    def __init__(self, name):
        self._name = name
        self._geodim = None
        self._properties = {}

    @property
    def name(self):
        return self._name

    @property
    def geodim(self):
        return self._geodim

    @geodim.setter
    def geodim(self, geodim):
        assert geodim in self._available_geomdim
        self._geodim = geodim

    @property
    def index(self):
        return self._index

    @index.setter
    def index(self, index: conn.indexlist):
        assert self.geodim in self._available_geomdim
        self._index = index

    @property
    def type(self):
        return self._properties['type']

    @type.setter
    def type(self, type):
        assert type in self._available_types
        self._properties['type'] = type

    def nodebased(self):
        return self._geodim in {'node', 'bdnode', 'intnode'}
    
    def facebased(self):
        return self._geodim in {'face', 'bdface', 'intface'}
    
    @property
    def properties(self):
        return self._properties

    @properties.setter
    def properties(self, properties: dict):
        self._properties = properties

    def __str__(self):
        return f"{self.name} ({self.geodim}): {self.index}"

class mesh():
    """versatile mesh object

    nodes are compulsory for a valid mesh

    cells (to nodes) connectivity is not compulsory but may be necessary
    for many operations

    faces (to nodes) connectivity is not compulsory ; it can exists either
        - `mixed` internal and boundaries, with mixed index
        - separated `internal` and `boundary` with separated index

    
    """
    __available_facetypes = ( 'mixed', 'internal', 'boundary')

    def __init__(self, ncell=0, nnode=0):
        self._ncell = ncell
        self._nnode = nnode
        self.nface = 0
        self._params = {}
        self._cell2node = None
        self._faces = {}
        self._nodes = {}
        self._bocos = {}
        self._celldata = {}
        self._nodedata = {}
        self._facedata = {}
        self._cellprop = {}

    @property
    def ncell(self):
        return self._ncell

    @property
    def nnode(self):
        return self._nnode

    def set_nodescoord_nd(self, xyz: np.ndarray):
        assert xyz.shape[0] == self.nnode
        for i,c in enumerate(['x', 'y', 'z']):
            self._nodes[c] = xyz[:,i]
    
    def set_nodescoord_xyz(self, x, y, z):
        for c in [x, y, z]:
            assert len(c) == self.nnode
        self._nodes['x'] = x
        self._nodes['y'] = y
        self._nodes['z'] = z

    def nodescoord(self, ndarray=False):
        coords = (self._nodes[c] for c in ['x', 'y', 'z'])
        return np.column_stack(coords) if ndarray else coords
    
    def set_cell2node(self, cell2node: conn.elem_connectivity):
        """set cell to node connectivity as a dict

        Args:
            cell2node (dict): dict of ndarray
        """
        self._cell2node = cell2node
        self._ncell = self._cell2node.nelem
        self._check_cell2node()

    def add_faces(self, facetype: str, face2node: conn.elem_connectivity, face2cell: conn.indexindirection = None):
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
        self.nface = np.sum([fcon['face2node'].nelem for _, fcon in self._faces.items()])

    def pop_faces(self, facetype: str):
        if facetype in self.__available_facetypes:
            if facetype in self._faces.keys():
                for key, item in self._faces[facetype].items():
                    del item
                self._faces.pop(facetype)

    def export_mixedfaces(self):
        mixedfaces_con = conn.elem_connectivity()
        mixedfaces_con.importfrom_merge((self._faces['boundary']['face2node'], 
                                         self._faces['internal']['face2node']))
        face2cell = conn.indexindirection()
        face2cell.conn = np.concatenate((self._faces['boundary']['face2cell'].conn, 
                                         self._faces['internal']['face2cell'].conn), axis=0)
        #print('merge',face2cell.conn)
        return mixedfaces_con, face2cell
                                                  

    def add_boco(self, boco: submeshmark):
        self._bocos[boco.name] = boco

    def remove_boco(self, name):
        self._bocos.pop(name)

    def bocomarks_set_node_to_face(self):
        """transform node list into face list
        face/node connectivity must exist and must be splitted into
        'internal'/'boundary' instead of 'mixed'
        """
        def face_in_nodelist(face, nodelist):
            return all(map(lambda n: n in nodelist, face))
        assert 'boundary' in self._faces.keys()
        index_face_tuples = self._faces['boundary']['face2node'].index_elem_tuples()
        for _, boco in self._bocos.items():
            if boco.nodebased():
                nodeset = set(boco.index.list())
                listface_index = [i for i,_ in 
                                filter(lambda t: face_in_nodelist(t[1], nodeset), 
                                        index_face_tuples)]
                boco.geodim = 'bdface'
                boco.index = conn.indexlist(list=listface_index)

    def seekmark(self, name: str)->submeshmark:
        """look for diffent marks set to find mark name"""
        # only _bocos for now
        return self._bocos[name]
    
    def exportmark_asmesh(self, name):
        meshmark = self.seekmark(name)
        newmesh = mesh()
 
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
    
    def reindex_boundaryfaces(self):
        assert 'boundary' in self._faces.keys(), "can only reindex faces according to boco if separated in 'boundary' list"
        for _,boco in self._bocos.items():
            assert boco.geodim in ('face', 'bdface'), "boco marks must be faces index"
        oldindex = list(itertools.chain(*[ boco.index.list() for _,boco in self._bocos.items() ]))
        assert np.all(np.unique(oldindex)==sorted(oldindex)), "some faces are marked by several boundary marks"
        assert min(oldindex) == 0, "first face index (0) is not marked as a boundary"
        assert max(oldindex) == len(oldindex)-1, "boundary faces must be indexed first before reindexing"
        newindex = np.full_like(oldindex, -1)
        newindex[oldindex] = np.arange(len(oldindex))
        assert min(newindex) == 0, "inconsistency: there must not be -1 index"
        # reindex boco
        for _,boco in self._bocos.items():
            boco.index = conn.indexlist(list=newindex[boco.index.list()])
            boco.index.compress() # try to (and must) make it a range
        # reindex boundary faces
        for _, fdict in self._faces['boundary']['face2node'].items():
            fdict['index'] = conn.indexlist(list=newindex[fdict['index'].list()])
            #fdict['index'].compress() # not expected
        if 'face2cell' in self._faces['boundary']:
            self._faces['boundary']['face2cell'].conn = self._faces['boundary']['face2cell'].conn[oldindex,:]

    def printinfo(self, detailed=False):
        api.io.print("std", f"nnode: {self.nnode}")
        for c in ('x', 'y', 'z'):
            api.io.print("std", "  {} min:avg:max = {:.3f}:{:.3f}:{:.3f}".format(c, *minavgmax(self._nodes[c])))
            
        api.io.print("std", f"ncell: {self.ncell}")
        if self._cell2node:
            self._cell2node.print()
        else:
            api.io.print("std", "  no cell/node connectivity")
        api.io.print('std', "nnode:",self.nnode)
        api.io.print('std', "nface:",self.nface)
        if self._faces:
            for t, facedict in self._faces.items():
                api.io.print("std", f"  type {t}: {' '.join(facedict['face2node'].elems())}")
                facedict['face2node'].print(prefix='  . ', detailed=detailed)
        else:
            api.io.print("std", "  no face/node connectivity")
        api.io.print('std', f"bocos: {' '.join(self._bocos.keys())}")
        for name, boco in self._bocos.items():
            api.io.print("std", f"  BC {name}: {boco}")
        api.io.print("std", "params:",self._params)

    def _check_cell2node(self):
        if self._cell2node is not None:
            assert isinstance(self._cell2node, conn.elem_connectivity), "cell2node connecitivity is not the expected class"
            assert self.ncell == self._cell2node.nelem, f"inconsistent size of cells {self.ncell} and {self._cell2node.nelem}"
        #for etype, conn in self._cell2node.items():
        #    assert etype in ele.elem2faces.keys()
        return True

    def make_face_connectivity(self):
        intfaces, intf2c, boundfaces, boundf2c = self._cell2node.create_faces_from_elems()
        self.pop_faces('mixed') # remove if it exists
        self.add_faces('internal', intfaces, intf2c)
        self.add_faces('boundary', boundfaces, boundf2c)

    def check(self):
        # check cell2node and cell numbers
        assert(self.ncell >0)
        assert(self.nnode >0)
        self._check_cell2node()
        #api.io.print('std','ckeck: at least cell/node or face/node face/cell connectivity')
        #assert(not self._cell2node or (not self._face2node and not self._face2cell))
        #assert self._check_cell2node() # not compulsory
        return True

    def morph(self, fmorph):
        '''
        change x, y, z position with a function of x, y, z, returning new x, y, z
        '''
        newx, newy, newz = fmorph(self._nodes['x'], self._nodes['y'], self._nodes['z'])
        self.set_nodescoord_xyz(newx, newy, newz)
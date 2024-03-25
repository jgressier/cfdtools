import itertools
import logging

import numpy as np

import cfdtools.api as api
import cfdtools.data as _data
import cfdtools.meshbase._connectivity as _conn
from cfdtools.utils.maths import minavgmax

log = logging.getLogger(__name__)


class submeshmark:
    # authorized geomdim type and actual dimension
    _available_geodim = (
        'node',
        'intnode',
        'bdnode',
        'face',
        'intface',
        'bdface',
        'cell',
    )

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
        assert geodim in self._available_geodim
        self._geodim = geodim

    @property
    def index(self):
        return self._index

    @index.setter
    def index(self, index: _conn.indexlist):
        assert self.geodim in self._available_geodim
        self._index = index

    @property
    def type(self):
        return self._properties['type']

    @type.setter
    def type(self, mtype):
        assert mtype in self._available_types
        self._properties['type'] = mtype

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
        return f"{self.name:12} ({self.geodim}): {self.index}"


class Mesh:
    """versatile mesh object

    nodes are compulsory for a valid mesh

    cells (to nodes) connectivity is not compulsory but may be necessary
    for many operations

    faces (to nodes) connectivity is not compulsory ; it can exists either
        - `mixed` internal and boundaries, with mixed index
        - separated `internal` and `boundary` with separated index

    """

    __available_facetypes = ('mixed', 'internal', 'boundary')

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
        self._facedata = None
        self._cellprop = {}

    @property
    def ncell(self):
        return self._ncell

    @property
    def nnode(self):
        return self._nnode

    def set_nodescoord_nd(self, xyz: np.ndarray):
        assert xyz.shape[0] == self.nnode
        for i, c in enumerate(['x', 'y', 'z']):
            self._nodes[c] = xyz[:, i]

    def set_nodescoord_xyz(self, x, y, z):
        for c in [x, y, z]:
            assert len(c) == self.nnode
        self._nodes['x'] = x
        self._nodes['y'] = y
        self._nodes['z'] = z

    def nodescoord(self, ndarray=False):
        coords = tuple(self._nodes[c] for c in ['x', 'y', 'z'])
        return np.column_stack(coords) if ndarray else coords

    def set_cell2node(self, cell2node: _conn.elem_connectivity):
        """set cell to node connectivity

        Args:
            cell2node (elem_connectivity): connectivity of cells
        """
        self._cell2node = cell2node
        self._ncell = self._cell2node.nelem
        self._check_cell2node()

    def add_faces(
        self,
        facetype: str,
        face2node: _conn.elem_connectivity,
        face2cell: _conn.indexindirection = None,
    ):
        """set faces connectivity with face type et optional face/cell connectivity

        Args:
            facetype (str): _description_
            face2node (_conn.elem_connectivity): _description_
            face2cell (_conn.indexindirection, optional): _description_. Defaults to None.
        """
        if facetype in self.__available_facetypes:
            self._faces[facetype] = {'face2node': face2node, 'face2cell': face2cell}
        else:
            api.error_stop(f"bad face type: {facetype} since {self.__available_facetypes} expected")
        self.nface = np.sum([fcon['face2node'].nelem for _, fcon in self._faces.items()])

    def pop_faces(self, facetype: str):
        if facetype in self.__available_facetypes:
            if facetype in self._faces.keys():
                for key, item in self._faces[facetype].items():
                    del item
                self._faces.pop(facetype)

    def export_mixedfaces(self):
        mixedfaces_con = _conn.elem_connectivity()
        mixedfaces_con.importfrom_merge(
            (self._faces['boundary']['face2node'], self._faces['internal']['face2node'])
        )
        face2cell = _conn.indexindirection()
        face2cell.conn = np.concatenate(
            (
                self._faces['boundary']['face2cell'].conn,
                self._faces['internal']['face2cell'].conn,
            ),
            axis=0,
        )
        # print('merge', face2cell.conn)
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
                # get all face index whose nodes are all in nodeset
                listface_index = [
                    i for i, _ in filter(lambda t: face_in_nodelist(t[1], nodeset), index_face_tuples)
                ]
                boco.geodim = 'bdface'
                boco.index = _conn.indexlist(ilist=listface_index)
                # print(boco.name, len(nodeset), len(boco.index.list()))

    def list_boco_index(self):
        return list(itertools.chain(*[boco.index.list() for boco in self._bocos.values()]))

    def make_unmarked_BC(self, name="unmarked_faces"):
        """check all boundaring faces are marked and create a specific boco if not"""
        if 'boundary' in self._faces.keys():
            # number of (boundary) faces in face connectivity
            # nbdface = self._faces['boundary']['face2node'].nelem
            for _, boco in self._bocos.items():
                assert boco.geodim in ('face', 'bdface'), "boco marks must be faces index"
            list_marked = self.list_boco_index()
            list_missing = list(set(self._faces['boundary']['face2node'].all_index()) - set(list_marked))
            if list_missing:
                boco = submeshmark(name)
                boco.geodim = 'bdface'
                boco.type = 'boundary'
                boco.properties['periodic_transform'] = None
                boco.index = _conn.indexlist(ilist=list_missing)
                self.add_boco(boco)
        else:
            list_missing = []
            log.warning("can only reindex faces according to boco if separated in 'boundary' list")
        return list_missing

    def seekmark(self, name: str) -> submeshmark:  # What is this method ??
        """look for diffent marks set to find mark name"""
        # only _bocos for now
        return self._bocos[name]

    def exportmark_asmesh(self, name):  # What is this method ??
        # meshmark = self.seekmark(name)  # Never used. Uncomment if useful...
        newmesh = Mesh()
        return newmesh

    def export_extruded(self, direction=None, extrude=None, domain="fluid"):
        if direction is None:
            direction = np.array([0.0, 0.0, 1.0])
        if extrude is None:
            extrude = [0.0, 1.0]
        extrude_range = np.array(extrude)
        nrange = extrude_range.size
        assert nrange > 1, "extrusion only possible for at least 2 planes"
        newmesh = Mesh(ncell=self.ncell * nrange, nnode=self.nnode * nrange)
        # SHOULD CHECK DIRECTION AND MESH ORIENTATION
        # extrude nodes
        ntotnode = self.nnode
        newcoords = np.tile(self.nodescoord(ndarray=True), (nrange, 1))
        for i, s in enumerate(extrude_range):
            newcoords[i * ntotnode : (i + 1) * ntotnode, :] += s * np.array(direction)
        newmesh.set_nodescoord_nd(newcoords)
        # extrude cells
        newmesh.set_cell2node(self._cell2node.extrude(nrange, ntotnode))
        # extrude faces if any
        # extrude/extend existing marks
        for _, boco in self._bocos.items():
            assert boco.nodebased(), "extrusion only possible with node marks"
            newboco = submeshmark(boco.name)
            newboco.geodim = boco.geodim
            newboco.type = boco.type
            index = np.tile(boco.index.list(), (nrange))
            nbcnode = boco.index.size
            for i in range(nrange):
                index[i * nbcnode : (i + 1) * nbcnode] += i * ntotnode
            newboco.index = _conn.indexlist(ilist=index.tolist())
            newmesh.add_boco(newboco)
        # create initial 2D domain as boco
        newboco = submeshmark(name=domain + '0')
        newboco.geodim = 'bdnode'
        newboco.type = 'boundary'
        index = self._cell2node.nodelist()
        newboco.index = _conn.indexlist(ilist=index)
        newmesh.add_boco(newboco)
        # create extruded 2D domain as boco
        newboco = submeshmark(name=domain + '1')
        newboco.geodim = 'bdnode'
        newboco.type = 'boundary'
        index = (np.array(index) + (nrange - 1) * ntotnode).tolist()
        newboco.index = _conn.indexlist(ilist=index)
        newmesh.add_boco(newboco)
        return newmesh

    def set_params(self, params):
        self._params = params

    def update_params(self, params):
        self._params.update(params)

    def set_celldata(self, celldata: _data.DataSet):
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
        return self._celldata.data[name] if self._celldata else None

    def reindex_boundaryfaces(self):
        assert (
            'boundary' in self._faces.keys()
        ), "can only reindex faces according to boco if separated in 'boundary' list"
        # number of (boundary) faces in face connectivity
        nbdface = self._faces['boundary']['face2node'].nelem
        for _, boco in self._bocos.items():
            assert boco.geodim in ('face', 'bdface'), "boco marks must be faces index"
        oldindex = self.list_boco_index()
        # checks
        c_unique = np.all(np.unique(oldindex) == sorted(oldindex))
        if not c_unique:  # pragma: no cover
            log.error("  some faces are marked by several boundary marks")
        c_min0 = min(oldindex) == 0
        if not c_min0:  # pragma: no cover
            log.error(
                "  first face index (0) is not marked as a boundary\n" "  some boundary faces may be missing",
            )
        c_max = max(oldindex) <= len(oldindex) - 1
        if not c_max:  # pragma: no cover
            log.error(
                "  max face reference is greater than the number of found faces\n"
                "  boundary faces must be indexed first before reindexing",
            )
        c_lengths = len(oldindex) == nbdface
        if not c_lengths:
            log.error(f"  some boundary faces are not marked: {nbdface-len(oldindex)}")
        if not (c_unique and c_min0 and c_max):
            api.error_stop("inconsistent face marks when reordering")
        newindex = np.full_like(oldindex, -1)
        newindex[oldindex] = np.arange(len(oldindex))
        assert min(newindex) >= 0, "inconsistency: there must not be -1 index"
        # reindex boco
        for _, boco in self._bocos.items():
            boco.index = _conn.indexlist(ilist=newindex[boco.index.list()].tolist())
            boco.index.compress()  # try to (and must) make it a range
        # reindex boundary faces
        for _, fdict in self._faces['boundary']['face2node'].items():
            fdict['index'] = _conn.indexlist(ilist=newindex[fdict['index'].list()].tolist())
            # fdict['index'].compress() # not expected
        if 'face2cell' in self._faces['boundary']:
            self._faces['boundary']['face2cell'].conn = self._faces['boundary']['face2cell'].conn[oldindex, :]

    def printinfo(self, detailed=False):
        log.info(f"nnode: {self.nnode}")
        for c in ('x', 'y', 'z'):
            log.info(
                f"  {c} min:avg:max =" + " {:.3f}:{:.3f}:{:.3f}".format(*minavgmax(self._nodes[c])),
            )

        log.info(f"ncell: {self.ncell}")
        if self._cell2node:
            self._cell2node.print()
        else:
            log.info("  no cell/node connectivity")
        log.info(f"nface: {self.nface}")
        if self._faces:
            for t, facedict in self._faces.items():
                log.info(f"  type {t}: {' '.join(facedict['face2node'].elems())}")
                facedict['face2node'].print(prefix='  . ', detailed=detailed)
        else:
            log.info("  no face/node connectivity")
        log.info(f"bocos: {' '.join(self._bocos.keys())}")
        for name, boco in self._bocos.items():
            log.info(f"  BC {boco}")
        log.info(f"params: {self._params}")

    def _check_cell2node(self):
        if self._cell2node is not None:
            assert isinstance(
                self._cell2node, _conn.elem_connectivity
            ), "cell2node connecitivity is not the expected class"
            assert (
                self.ncell == self._cell2node.nelem
            ), f"inconsistent size of cells {self.ncell} and {self._cell2node.nelem}"
        # for etype, econn in self._cell2node.items():
        #    assert etype in _elem.elem2faces.keys()
        return True

    def make_face_connectivity(self):
        (
            intfaces,
            intf2c,
            boundfaces,
            boundf2c,
        ) = self._cell2node.create_faces_from_elems()
        self.pop_faces('mixed')  # remove if it exists
        self.add_faces('internal', intfaces, intf2c)
        self.add_faces('boundary', boundfaces, boundf2c)

    def check(self):
        """check a consistent mesh with adequate boundary conditions"""
        # check cell2node and cell numbers
        assert self.ncell > 0
        assert self.nnode > 0
        self._check_cell2node()
        # log.info('ckeck: at least cell/node or face/node face/cell connectivity')
        # assert(not self._cell2node or (not self._face2node and not self._face2cell))
        # assert self._check_cell2node() # not compulsory
        assert not self.make_unmarked_BC()
        return True

    def morph(self, fmorph):
        '''
        change x, y, z position with a function of x, y, z, returning new x, y, z
        '''
        newx, newy, newz = fmorph(self._nodes['x'], self._nodes['y'], self._nodes['z'])
        self.set_nodescoord_xyz(newx, newy, newz)

    def scale(self, scale_xyz):
        '''
        scale x, y, z positions
        '''
        self.set_nodescoord_xyz(
            self._nodes['x'] * scale_xyz[0],
            self._nodes['y'] * scale_xyz[1],
            self._nodes['z'] * scale_xyz[2],
        )

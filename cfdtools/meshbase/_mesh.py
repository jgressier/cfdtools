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
        "node",
        "intnode",
        "bdnode",
        "face",
        "intface",
        "bdface",
        "cell",
    )

    _available_types = (
        "internal",
        "boundary",
        "perio_cart",
        "perio_cylx",
        "perio_cyly",
        "perio_cylz",
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
        return self._properties["type"]

    @type.setter
    def type(self, mtype):
        assert mtype in self._available_types
        self._properties["type"] = mtype

    def nodebased(self):
        return self._geodim in {"node", "bdnode", "intnode"}

    def facebased(self):
        return self._geodim in {"face", "bdface", "intface"}

    @property
    def properties(self):
        return self._properties

    @properties.setter
    def properties(self, properties: dict):
        self._properties = properties

    def __str__(self):
        return f"{self.name:12} ({self.geodim}): {self.index}"


class Mesh:
    """Versatile mesh object.

    A mesh holds node coordinates and, optionally, cell and face
    connectivity to those nodes:

    - Node coordinates are compulsory for a valid mesh.
    - Cell-to-node connectivity is not compulsory, but many operations
      (e.g. face creation, extrusion) require it.
    - Face-to-node connectivity is not compulsory. It can exist either as
      a `mixed` set (internal and boundary faces with a shared index) or
      as separated `internal` and `boundary` sets with their own indices.

    """

    __available_facetypes = ("mixed", "internal", "boundary")

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
        for i, c in enumerate(["x", "y", "z"]):
            self._nodes[c] = xyz[:, i]

    def set_nodescoord_xyz(self, x, y, z):
        for c in [x, y, z]:
            assert len(c) == self.nnode
        self._nodes["x"] = x
        self._nodes["y"] = y
        self._nodes["z"] = z

    def nodescoord(self, ndarray=False):
        coords = tuple(self._nodes[c] for c in ["x", "y", "z"])
        return np.column_stack(coords) if ndarray else coords

    def set_cell2node(self, cell2node: _conn.elem_connectivity):
        """Set element-node connectivity.

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
        """Set faces connectivity with face type et optional face/cell connectivity.

        Args:
            facetype (str): _description_
            face2node (_conn.elem_connectivity): _description_
            face2cell (_conn.indexindirection, optional): _description_. Defaults to None.

        """
        if facetype not in self.__available_facetypes:
            api.error_stop(
                f"bad face type: {facetype} since {self.__available_facetypes} expected"
            )

        self._faces[facetype] = {"face2node": face2node, "face2cell": face2cell}
        self.nface = sum(fcon["face2node"].nelem for fcon in self._faces.values())

    def pop_faces(self, facetype: str):
        if facetype in self.__available_facetypes:
            if facetype in self._faces.keys():
                for key, item in self._faces[facetype].items():
                    del item
                self._faces.pop(facetype)

    def export_mixedfaces(self):
        mixedfaces_con = _conn.elem_connectivity()
        mixedfaces_con.importfrom_merge(
            (self._faces["boundary"]["face2node"], self._faces["internal"]["face2node"])
        )
        face2cell = _conn.indexindirection()
        face2cell.conn = np.concatenate(
            (
                self._faces["boundary"]["face2cell"].conn,
                self._faces["internal"]["face2cell"].conn,
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
        """Transform node list into face list.

        face/node connectivity must exist and must be splitted into
        'internal'/'boundary' instead of 'mixed'
        """
        assert "boundary" in self._faces.keys()
        # Per boco, test membership on the already-stored face2node arrays
        # (no face_nodes concatenation: that would allocate nbdface x
        # nodes_per_face and spike the RSS peak). For each boundary face, all
        # of its nodes must belong to the boco.
        face2node = self._faces["boundary"]["face2node"]
        for boco in self._bocos.values():
            if boco.nodebased():
                boco_nodes = boco.index.to_list()
                face_indices = []
                for e2n in face2node.values():
                    nodes = e2n["elem2node"]
                    # mask of the faces of this type whose every node is in
                    # the boco; select their indices
                    keep = np.isin(nodes, boco_nodes).all(axis=1)
                    face_indices.append(e2n["index"].to_list()[keep])
                if any(a.size for a in face_indices):
                    boco.geodim = "bdface"
                    boco.index = _conn.indexlist(ilist=np.concatenate(face_indices))

    def list_boco_index(self):
        if not self._bocos:
            return np.array([], dtype=np.int64)
        return np.concatenate([b.index.to_list() for b in self._bocos.values()])

    def make_unmarked_BC(self, name="unmarked_faces"):
        """Check all boundaring faces are marked and create a specific boco if not."""
        if "boundary" in self._faces:
            # number of (boundary) faces in face connectivity
            # nbdface = self._faces['boundary']['face2node'].nelem
            for boco in self._bocos.values():
                assert boco.geodim in ("face", "bdface"), (
                    "boco marks must be faces index"
                )
            list_marked = self.list_boco_index()
            all_idx = self._faces["boundary"]["face2node"].all_index()
            # Which boundary face indices are not marked by any boco?
            list_missing = all_idx[~np.isin(all_idx, list_marked)]
            if list_missing.size:
                boco = submeshmark(name)
                boco.geodim = "bdface"
                boco.type = "boundary"
                boco.properties["periodic_transform"] = None
                boco.index = _conn.indexlist(ilist=list_missing)
                self.add_boco(boco)
            return bool(list_missing.size)
        else:
            log.warning(
                "can only reindex faces according to boco if separated in 'boundary' list"
            )
            return False

    def seekmark(self, name: str) -> submeshmark:  # What is this method ??
        """Look for diffent marks set to find mark name."""
        # only _bocos for now
        return self._bocos[name]

    def exportmark_asmesh(self, name):  # What is this method ??
        # meshmark = self.seekmark(name)  # Never used. Uncomment if useful...
        newmesh = Mesh()
        return newmesh

    def export_extruded(self, direction=None, extrude=None, domain="fluid"):
        import cfdtools.meshbase._elements as _elem

        if direction is None:
            direction = np.array([0.0, 0.0, 1.0])
        else:
            direction = np.array(direction)
        if extrude is None:
            extrude = [0.0, 1.0]
        extrude_range = np.array(extrude)
        nplanes = extrude_range.size
        if nplanes <= 1:
            raise ValueError(
                f"Extrusion only possible for at least 2 planes. Actual: {nplanes}"
            )

        # determine if quadratic
        etype = next(iter(self._cell2node.elems()))
        # number of nodes of the source element
        fnnode = _elem.nnode_elem[etype]
        target_type = _elem.extruded_face[etype]
        target_nnode = _elem.nnode_elem[target_type]
        is_quadratic = target_nnode == 3 * fnnode

        ntotnode = self.nnode
        # total nodes: for quadratic, there's a middle layer between each pair of planes
        if is_quadratic:
            nnewnode = ntotnode * (2 * nplanes - 1)
        else:
            nnewnode = ntotnode * nplanes
        newmesh = Mesh(ncell=self.ncell * (nplanes - 1), nnode=nnewnode)

        # extrude nodes
        coords = self.nodescoord(ndarray=True)
        newcoords = np.zeros((nnewnode, 3))
        if is_quadratic:
            # layout: plane0, mid01, plane1, mid12, plane2, ...
            for i, s in enumerate(extrude_range):
                newcoords[2 * i * ntotnode : (2 * i + 1) * ntotnode, :] = (
                    coords + s * direction
                )
                if i < nplanes - 1:
                    mid_s = 0.5 * (s + extrude_range[i + 1])
                    newcoords[(2 * i + 1) * ntotnode : (2 * i + 2) * ntotnode, :] = (
                        coords + mid_s * direction
                    )
        else:
            for i, s in enumerate(extrude_range):
                newcoords[i * ntotnode : (i + 1) * ntotnode, :] = coords + s * direction
        newmesh.set_nodescoord_nd(newcoords)
        # extrude cells (inodeshift is always ntotnode)
        newmesh.set_cell2node(self._cell2node.extrude(nplanes, ntotnode))
        # extrude faces if any
        # extrude/extend existing marks
        plane_step = 2 * ntotnode if is_quadratic else ntotnode
        for boco in self._bocos.values():
            assert boco.nodebased(), "extrusion only possible with node marks"
            newboco = submeshmark(boco.name)
            newboco.geodim = boco.geodim
            newboco.type = boco.type
            index = np.tile(boco.index.to_list(), (nplanes))
            nbcnode = boco.index.size
            for i in range(nplanes):
                index[i * nbcnode : (i + 1) * nbcnode] += i * plane_step
            newboco.index = _conn.indexlist(ilist=index)
            newmesh.add_boco(newboco)
        # create initial 2D domain as boco
        newboco = submeshmark(name=domain + "0")
        newboco.geodim = "bdnode"
        newboco.type = "boundary"
        index = self._cell2node.nodelist()
        newboco.index = _conn.indexlist(ilist=index)
        newmesh.add_boco(newboco)
        # create extruded 2D domain as boco
        newboco = submeshmark(name=domain + "1")
        newboco.geodim = "bdnode"
        newboco.type = "boundary"
        index = np.array(index) + (nplanes - 1) * plane_step
        newboco.index = _conn.indexlist(ilist=index)
        newmesh.add_boco(newboco)
        return newmesh

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
        self._cellprop["partition"] = partition

    def pop_nodedata(self, name):
        return self._nodedata.pop(name) if name in self._nodedata.keys() else None

    def pop_celldata(self, name):
        return self._celldata.data[name] if self._celldata else None

    def reindex_boundaryfaces(self):
        assert "boundary" in self._faces.keys(), (
            "can only reindex faces according to boco if separated in 'boundary' list"
        )
        # number of (boundary) faces in face connectivity
        nbdface = self._faces["boundary"]["face2node"].nelem
        for boco in self._bocos.values():
            assert boco.geodim in ("face", "bdface"), "boco marks must be faces index"
        # compact numpy array of all boco face indices (a Python list of boxed
        # ints would spike the RSS peak)
        oldindex = np.concatenate([b.index.to_list() for b in self._bocos.values()])
        # checks
        c_unique = np.all(np.unique(oldindex) == np.sort(oldindex))
        if not c_unique:  # pragma: no cover
            log.error("  some faces are marked by several boundary marks")
        c_min0 = oldindex.min() == 0
        if not c_min0:  # pragma: no cover
            log.error(
                "  first face index (0) is not marked as a boundary\n"
                "  some boundary faces may be missing",
            )
        c_max = oldindex.max() <= len(oldindex) - 1
        if not c_max:  # pragma: no cover
            log.error(
                "  max face reference is greater than the number of found faces\n"
                "  boundary faces must be indexed first before reindexing",
            )
        c_lengths = len(oldindex) == nbdface
        if not c_lengths:
            log.error(
                f"  some boundary faces are not marked: {nbdface - len(oldindex)}"
            )
        if not (c_unique and c_min0 and c_max):
            api.error_stop("inconsistent face marks when reordering")
        newindex = np.full_like(oldindex, -1)
        newindex[oldindex] = np.arange(len(oldindex))
        assert min(newindex) >= 0, "inconsistency: there must not be -1 index"
        # reindex boco
        for boco in self._bocos.values():
            boco.index = _conn.indexlist(ilist=newindex[boco.index.to_list()])
            boco.index.compress()  # try to (and must) make it a range
        # reindex boundary faces
        for fdict in self._faces["boundary"]["face2node"].values():
            fdict["index"] = _conn.indexlist(
                ilist=newindex[fdict["index"].to_list(dtype=oldindex.dtype)]
            )
            # fdict['index'].compress() # not expected
        if "face2cell" in self._faces["boundary"]:
            self._faces["boundary"]["face2cell"].conn = self._faces["boundary"][
                "face2cell"
            ].conn[oldindex, :]

    def printinfo(self, detailed=False):
        log.info(f"nnode: {self.nnode}")
        for c in ("x", "y", "z"):
            log.info(
                f"  {c} min:avg:max ="
                + " {:.3f}:{:.3f}:{:.3f}".format(*minavgmax(self._nodes[c])),
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
                facedict["face2node"].print(prefix="  . ", detailed=detailed)
        else:
            log.info("  no face/node connectivity")
        log.info(f"bocos: {' '.join(self._bocos.keys())}")
        for name, boco in self._bocos.items():
            log.info(f"  BC {boco}")
        log.info(f"params: {self._params}")

    def _check_cell2node(self):
        if self._cell2node is not None:
            assert isinstance(self._cell2node, _conn.elem_connectivity), (
                "cell2node connecitivity is not the expected class"
            )
            assert self.ncell == self._cell2node.nelem, (
                f"inconsistent size of cells {self.ncell} and {self._cell2node.nelem}"
            )
        # for etype, econn in self._cell2node.items():
        #    assert etype in _elem.elem2faces.keys()
        return True

    def make_face_connectivity(self, method="current", nblocks=None):
        """Build faces from cells with the selected connectivity implementation.

        Args:
            method: One of ``"current"``, ``"legacy"`` or ``"multiblock"``.
            nblocks: Number of consecutive cell-index blocks. Required when
                ``method="multiblock"``; ignored otherwise.

        """
        if method == "multiblock":
            if nblocks is None:
                raise ValueError("nblocks is required when method='multiblock'")
            (
                intfaces,
                intf2c,
                boundfaces,
                boundf2c,
            ) = self._cell2node.create_faces_from_elems_multiblock(nblocks)
        else:
            methods = {
                "current": self._cell2node.create_faces_from_elems,
                "legacy": self._cell2node._create_faces_from_elems_legacy,
            }
            if method not in methods:
                raise ValueError(f"unknown face connectivity method: {method}")
            (
                intfaces,
                intf2c,
                boundfaces,
                boundf2c,
            ) = methods[method]()
        self.pop_faces("mixed")  # remove if it exists
        self.add_faces("internal", intfaces, intf2c)
        self.add_faces("boundary", boundfaces, boundf2c)

    def check(self):
        """Check a consistent mesh with adequate boundary conditions."""
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
        """Change x, y, z position with a function of x, y, z, returning new x, y, z."""
        newx, newy, newz = fmorph(self._nodes["x"], self._nodes["y"], self._nodes["z"])
        self.set_nodescoord_xyz(newx, newy, newz)

    def scale(self, scale_xyz):
        """Scale x, y, z positions."""
        self.set_nodescoord_xyz(
            self._nodes["x"] * scale_xyz[0],
            self._nodes["y"] * scale_xyz[1],
            self._nodes["z"] * scale_xyz[2],
        )

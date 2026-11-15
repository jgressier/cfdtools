# cgns.py
import logging
from pathlib import Path

try:
    from functools import cache  # python >= 3.9
except ImportError:
    from functools import lru_cache  #  3.6 <= python <= 3.8

    cache = lru_cache(maxsize=None)
    del lru_cache

import numpy as np

from cfdtools.api import error_stop, fileformat_reader  # , memoize
from cfdtools.hdf5 import h5File, h5_str
from cfdtools.meshbase._mesh import Mesh, submeshmark
import cfdtools.data as _data
import cfdtools.meshbase._connectivity as _conn
import cfdtools.meshbase._elements as _elem

log = logging.getLogger(__name__)

cgtype = {}
ele_cgns2local = {2: 'node1', 3: 'bar2', 5: 'tri3', 7: 'quad4', 9: 'quad9', 17: 'hexa8', 19: 'hexa27'}


def cgnstype(obj):
    cgnsdatatype = obj.attrs.get('label')
    return cgnsdatatype


def dict_cgnstype(obj, cgtype):
    return {name: obj for name, obj in obj.items() if cgnstype(obj) == cgtype}


def cg_gridlocation(bc):
    if cgnstype(bc) not in [b"BC_t", b"GridConnectivity_t"]:
        error_stop("Invalid BC type encountered in cg_gridlocation")
    if "GridLocation" in bc.keys():
        bcloc = h5_str(bc["GridLocation/ data"])
    else:
        bcloc = "Vertex"
    return bcloc


class cgnszone:

    def __init__(self, zone, geodim=None) -> None:
        self._zone = zone
        self._zonetype = h5_str(zone["ZoneType/ data"])
        self._geodim = geodim
        if self._zonetype != "Unstructured":
            error_stop("Only Unstructured zone expected")
        # look for Elements
        self._elems = dict_cgnstype(zone, b'Elements_t')
        # look for ZoneBC and BC
        self._zonebc = dict_cgnstype(zone, b'ZoneBC_t')
        self._zonebc.update(dict_cgnstype(zone, b"ZoneGridConnectivity_t"))
        self._BCs = {}
        for zbc in self._zonebc.values():
            self._BCs.update(dict_cgnstype(zbc, b'BC_t'))
            self._BCs.update(dict_cgnstype(zbc, b'GridConnectivity_t'))

    @property
    def nnode(self):
        return self._zone[' data'][0, 0]

    @property
    def ncell(self):
        return self._zone[' data'][1, 0]

    def coords(self):
        x = self._zone['GridCoordinates']['CoordinateX/ data'][:]
        y = self._zone['GridCoordinates']['CoordinateY/ data'][:]
        z = self._zone['GridCoordinates']['CoordinateZ/ data'][:]
        return x, y, z

    def elemcon(self, geodim):
        cellconn = _conn.elem_connectivity()
        for _, elements in self._elems.items():
            cgnstype = elements[" data"][0]
            etype = ele_cgns2local[cgnstype]
            nnode = _elem.nnode_elem[etype]
            # Keep only elements that belong to the requested dimension.
            if _elem.dim_elem[etype] != geodim:
                continue

            # Change the node numbering convention from CGNS to 0-based.
            index = _conn.indexlist(irange=elements["ElementRange/ data"][:] - 1)
            econ = elements["ElementConnectivity/ data"][:].reshape((-1, nnode))
            # Change the node numbering convention from CGNS to 0-based.
            econ -= 1
            cellconn.add_elems(etype, econ, index)
        return cellconn

    def export_cellcon(self):
        return self.elemcon(self._geodim)

    # @memoize
    @cache
    def export_facecon(self):
        return self.elemcon(self._geodim - 1)

    def export_BC(self, BC):
        if "FamilyName" in BC.keys():
            name = h5_str(BC["FamilyName/ data"])
        else:
            name = Path(BC.name).name  # extract final name of
        boco = submeshmark(name)
        boco.geodim = 'node'  # don't know if node, intnode or bdnode
        boco.type = 'boundary'
        boco.properties['BCtype'] = h5_str(BC[" data"])
        boco.properties['periodic_transform'] = None
        if "PointList" in BC.keys():
            indexlist = (BC["PointList/ data"][:] - 1).ravel().tolist()
            gridloc = cg_gridlocation(BC)
        elif "ElementList" in BC.keys():  # not in CGNS norm
            indexlist = (BC["ElementList/ data"][:] - 1).ravel().tolist()
            gridloc = "FaceCenter"
            if len(indexlist) == self.ncell:
                gridloc = "CellCenter"
                boco.type = 'internal'
        else:
            error_stop(f"Unknown indexing of BC mark: {name}")
        # convert to node marks
        if gridloc == "FaceCenter":
            nodelist = self.export_facecon().nodes_of_indexlist(indexlist)
        elif gridloc == "Vertex":
            nodelist = indexlist
            if len(nodelist) == self.nnode:
                boco.type = 'internal'
        elif gridloc == "CellCenter":
            nodelist = indexlist  # cells indeed
            boco.geodim = 'cell'
        else:
            error_stop(f'unknown gridlocation {gridloc}')
        boco.index = _conn.indexlist(ilist=nodelist)  # must start at 0
        return boco


class cgnsfile(h5File):
    def __init__(self, filename: str):
        super().__init__(filename)
        self.open()

    def open(self):
        super().open()
        self._cgnsver = self._h5file['CGNSLibraryVersion'][' data'][0]

    def printinfo(self):
        super().printinfo()

    def list_bases(self):
        return [bname for bname, base in self._h5file.items() if cgnstype(base) == b'CGNSBase_t']


@fileformat_reader('CGNS', '.cgns')
class cgnsMesh:
    def __init__(self, filename) -> None:
        self._filename = filename
        self._ncell = None
        self._exclude_center_points = False
        self._celldata = _data.DataSet("cellaverage")

    @property
    def ncell(self):
        return self._ncell

    def read_data(self, zone=None, exclude_center_points=False):
        log.info(f"> CGNS reader: starts reading {self._filename}")
        self._exclude_center_points = exclude_center_points
        # Use the CGNS node ordering convention for face extraction. This is the
        # module default, but a previous gmsh read may have switched the global.
        _elem.elem2faces = _elem.cgns_elem2faces
        # Check file exists
        if not Path(self._filename).exists():
            error_stop(f"File not found: {self._filename!r}")
        self._file = cgnsfile(self._filename)
        # get BASE list
        self._bases = self._file.list_bases()
        for base in self._bases:
            self._zones = dict_cgnstype(self._file._h5file[base], b'Zone_t')
        # geo dimension from base
        self._geodim = self._file._h5file[self._bases[0]][" data"][0]
        if zone is None:
            if len(self._zones) != 1:
                error_stop("Multiple zones found, must specify which zone to export")
            name = list(self._zones.keys())[0]
        else:
            name = zone
        self._zonename = name
        self._zone = cgnszone(self._zones[name], self._geodim)
        self._ncell = self._zone.ncell

    def printinfo(self):
        # super().printinfo()
        self._file.printinfo()
        log.info(f"CGNS version: {self._file._cgnsver}")
        log.info(f"bases: {self._bases}")
        log.info(f"zones: {list(self._zones.keys())}")
        for zn in self._zones.keys():
            log.info(f"  Zone {zn}")

    def export_mesh(self):
        log.info("> export mesh from CGNS")
        cgzone = self._zone
        log.info(
            f"Parse zone {self._zonename} ({self._geodim}D) ncell: {cgzone.ncell}, nnode: {cgzone.nnode}"
        )
        # get coordinates and cell connectivity
        x, y, z = cgzone.coords()
        cellcon = cgzone.export_cellcon()
        # boundary conditions (collect, filtering out full-domain internal marks)
        bocos = []
        # boundary conditions
        for _, bc in cgzone._BCs.items():
            boco = cgzone.export_BC(bc)
            # filter full domain
            if boco.type in ['internal']:
                log.info(f"  filter internal mark {boco.name}")
            else:
                bocos.append(boco)
        # optionally remove the central (27th) node of each HEXA27 element
        if self._exclude_center_points:
            x, y, z = self.__remove_27th_point_hexa27(cellcon, bocos, x, y, z)

        # assemble mesh
        meshdata = Mesh(ncell=cgzone.ncell, nnode=len(x))
        meshdata.set_nodescoord_xyz(x, y, z)
        meshdata.set_cell2node(cellcon)
        for boco in bocos:
            log.info(f"  add boco {boco}")
            meshdata.add_boco(boco)
        if self._exclude_center_points:
            meshdata.set_celldata(self._celldata)

        # meshdata.check()
        # meshdata.printinfo()
        return meshdata

    def __remove_27th_point_hexa27(self, cellcon, bocos, x, y, z):
        """Remove the 27th local node of each HEXA27 element.

        The HEXA27 connectivity is reduced to its first 26 nodes.
        The coordinates of the removed nodes are deleted and all remaining node
        indices (cell and boundary) are renumbered.
        Removed coordinates are kept as "point27" cell data.

        Parameters
        ----------
        bocos : dict
            Boundary‑face connectivity dictionary. (`{bnd_tag: {elt_type: np.ndarray}}`).
        cellcon : (nel,27) array_like
            HEXA27 connectivity
        x, y, z : (nnode,ndim) array_like
            Coordinates
        Returns
        -------
        x : ndarray
            Coordinates without point27
        """
        if 'hexa27' not in cellcon.elems():
            log.warning("exclude_center_points: no hexa27 element found, nothing to remove")
            return x, y, z
        conn = cellcon['hexa27']

        # Nodes in 27th column
        removed_nodes = conn[:, 26]

        # Store the actual coordinates of those centre points (for celldata)
        removed_coords = np.column_stack((x[removed_nodes], y[removed_nodes], z[removed_nodes]))

        # Keep mask
        initial_nb_nodes = len(x)
        keep = np.ones(initial_nb_nodes, dtype=bool)
        keep[removed_nodes] = False

        # Renumber map old -> new
        new_id = -np.ones(initial_nb_nodes, dtype=int)
        new_id[keep] = np.arange(np.sum(keep))

        # Renumber connectivity
        # reduce HEXA27 to its first 26 nodes and renumber
        cellcon._elem2node['hexa27']['elem2node'] = new_id[conn[:, :26]]
        for boco in bocos:
            if boco.nodebased():
                boco.index = _conn.indexlist(ilist=new_id[np.asarray(boco.index.list())].tolist())
        self._celldata.add_data("point27", removed_coords)
        log.info(f"  excluded {len(removed_nodes)} hexa27 center nodes")
        return x[keep], y[keep], z[keep]

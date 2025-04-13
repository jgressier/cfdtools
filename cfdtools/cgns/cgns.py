# cgns.py
import logging
from pathlib import Path

try:
    from functools import cache  # python >= 3.9
except ImportError:
    from functools import lru_cache  #  3.6 <= python <= 3.8

    cache = lru_cache(maxsize=None)
    del lru_cache
from cfdtools.api import error_stop, fileformat_reader  # , memoize
from cfdtools.hdf5 import h5File, h5_str
from cfdtools.meshbase._mesh import Mesh, submeshmark
import cfdtools.meshbase._connectivity as _conn
import cfdtools.meshbase._elements as _elem

log = logging.getLogger(__name__)

cgtype = {}
ele_cgns2local = {2: 'node1', 3: 'bar2', 5: 'tri3', 7: 'quad4', 17: 'hexa8'}


def cgnstype(obj):
    cgnsdatatype = obj.attrs.get('label')
    return cgnsdatatype


def dict_cgnstype(obj, cgtype):
    return {name: obj for name, obj in obj.items() if cgnstype(obj) == cgtype}


def cg_gridlocation(bc):
    assert cgnstype(bc) in [b"BC_t", b"GridConnectivity_t"]
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
        assert self._zonetype == "Unstructured", "Only Unstructured zone expected"
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
            # extract cell connectivity only
            if _elem.dim_elem[etype] == geodim:
                index = _conn.indexlist(irange=elements["ElementRange/ data"][:] - 1)
                econ = elements["ElementConnectivity/ data"][:].reshape((-1, nnode))
                econ -= 1  # shift node index (starts 0)
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

    @property
    def ncell(self):
        return self._ncell

    def read_data(self, zone=None):
        log.info(f"> CGNS reader: starts reading {self._filename}")
        # Check file exists
        if not Path(self._filename).exists():
            error_stop(f"File not found: {self._filename!r}")
        self._file = cgnsfile(self._filename)
        # get BASE list
        self._bases = self._file.list_bases()
        for base in self._bases:
            # print('base', self._file._h5file[base].name, self._file._h5file[base][" data"][:])
            self._zones = dict_cgnstype(self._file._h5file[base], b'Zone_t')
        # geo dimension from base
        self._geodim = self._file._h5file[self._bases[0]][" data"][0]
        if zone is None:
            assert len(self._zones) == 1, "Multiple zones found, must specify which zone to export"
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
            # for bcn, bc in

    def export_mesh(self):
        # log.info(f"> export mesh ") # printed by parent
        cgzone = self._zone
        log.info(
            f"Parse zone {self._zonename} ({self._geodim}D) ncell: {cgzone.ncell}, nnode: {cgzone.nnode}"
        )
        meshdata = Mesh(ncell=cgzone.ncell, nnode=cgzone.nnode)
        # get coordinates
        meshdata.set_nodescoord_xyz(*cgzone.coords())
        # # cell connectivity
        meshdata.set_cell2node(cgzone.export_cellcon())
        # boundary conditions
        for _, bc in cgzone._BCs.items():
            boco = cgzone.export_BC(bc)
            # filter full domain
            if boco.type in ['internal']:
                log.info(f"  filter internal mark {boco.name}")
            else:
                log.info(f"  add boco {boco}")
                meshdata.add_boco(boco)
        # meshdata.check()
        # meshdata.printinfo()
        return meshdata


# if __name__ == "__main__":
#     f = cgnsMesh(filename="./examples/MESH.nogit/cavity-degen.hdf")
#     f.read_data()
#     f.printinfo()
#     f.export_mesh()

# ElementType_t := Enumeration(
#      ElementTypeNull, ElementTypeUserDefined, NODE, BAR_2, BAR_3,
#      TRI_3, TRI_6, QUAD_4, QUAD_8, QUAD_9,
#      TETRA_4, TETRA_10, PYRA_5, PYRA_14,
#      PENTA_6, PENTA_15, PENTA_18, HEXA_8, HEXA_20, HEXA_27,
#      MIXED, PYRA_13, NGON_n, NFACE_n,
#      BAR_4, TRI_9, TRI_10, QUAD_12, QUAD_16,
#      TETRA_16, TETRA_20, PYRA_21, PYRA_29, PYRA_30,
#      PENTA_24, PENTA_38, PENTA_40, HEXA_32, HEXA_56, HEXA_64 );

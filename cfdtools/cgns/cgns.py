# cgns.py
from cfdtools.api import io, fileformat_reader
from cfdtools.hdf5 import h5file, h5_str
from cfdtools.meshbase._mesh import Mesh, submeshmark
import cfdtools.meshbase._connectivity as conn
import cfdtools.meshbase._elements as ele
import numpy as np

ele_cgns2local = {
    2: 'node1',
    3: 'bar2',
    5: 'tri3',
    7: 'quad4',
    17: 'hexa8'
}

def cgnstype(obj):
    cgnsdatatype = obj.attrs.get('label')
    return cgnsdatatype

def dict_cgnstype(obj, cgtype):
    return { name: obj for name, obj in obj.items() if cgnstype(obj) == cgtype }

class cgnszone():
    def __init__(self, zone, geodim=None) -> None:
        self._zone = zone
        self._zonetype = h5_str(zone["ZoneType/ data"])
        self._geodim = geodim
        assert self._zonetype == "Unstructured", "Only Unstructured zone expected"
        # look for Elements
        self._elems = dict_cgnstype(zone, b'Elements_t')
        # look for ZoneBC and BC
        self._zonebc = dict_cgnstype(zone, b'ZoneBC_t')
        self._BCs = {}
        for zbc in self._zonebc.values():
            self._BCs.update(dict_cgnstype(zbc, b'BC_t'))
        #print(self._BCs)
 
    @property
    def nnode(self):
        return self._zone[' data'][0,0]

    @property
    def ncell(self):
        return self._zone[' data'][1,0]
    
    def coords(self):
        x = self._zone['GridCoordinates']['CoordinateX/ data'][:]
        y = self._zone['GridCoordinates']['CoordinateY/ data'][:]
        z = self._zone['GridCoordinates']['CoordinateZ/ data'][:]
        return x, y, z

    def export_cellcon(self):
        cellconn = conn.elem_connectivity()
        for name, elements in self._elems.items():
            cgnstype = elements[" data"][0]
            etype = ele_cgns2local[cgnstype]
            #print(cgnstype,etype, ele.elem_dim[etype], self._geodim)
            # extract cell connectivity only
            if ele.elem_dim[etype] == self._geodim:
                #print(elements["ElementRange/ data"][:])
                index = conn.indexlist(range=elements["ElementRange/ data"][:]-1)
                econ = elements["ElementConnectivity/ data"][:].reshape((-1,ele.nnode_elem[etype]))-1
                #print(econ)
                cellconn.add_elems(etype, econ, index)
        return cellconn

    def export_facecon(self):
        pass

    def export_BC(self, BC):
        name = h5_str(BC["FamilyName/ data"])
        boco = submeshmark(name)
        boco.geodim = 'node' # don't know if node, intnode or bdnode
        boco.type = 'boundary'
        boco.properties['BCtype'] = h5_str(BC[" data"])
        boco.properties['periodic_transform'] = None
        assert "PointList" in BC.keys(), "only PointList implemented"
        nodelist = (BC["PointList/ data"][:]-1).ravel().tolist()
        boco.index = conn.indexlist(list=nodelist) # must start at 0
        if boco.index.size == self.nnode:
            boco.type = 'internal'
        return boco

class cgnsfile(h5file):

    def __init__(self, filename: str):
        super().__init__(filename)
        self.open()

    def open(self):
        super().open()
        self._cgnsver = self._h5file['CGNSLibraryVersion'][' data'][0]

    def printinfo(self):
        super().printinfo()
    
    def list_bases(self):
        return [ bname for bname, base in self._h5file.items()
                       if cgnstype(base) == b'CGNSBase_t' ]

@fileformat_reader('CGNS', '.cgns')
class cgnsMesh():
    def __init__(self, filename) -> None:
        self._filename = filename

    def read_data(self):
        self._file = cgnsfile(self._filename)
        # get BASE list
        self._bases = self._file.list_bases()
        for base in self._bases:
            #print('base', self._file._h5file[base].name, self._file._h5file[base][" data"][:])
            self._zones = dict_cgnstype(self._file._h5file[base], b'Zone_t')

    def printinfo(self):
        #super().printinfo()
        self._file.printinfo()
        io.print('std', 'CGNS version:', self._file._cgnsver)
        io.print('std', 'bases:',self._bases)
        io.print('std', 'zones:',list(self._zones.keys()))
        for zn, z in self._zones.items():
            io.print('std', f"  Zone {zn}")
            #for bcn, bc in 

    def export_mesh(self, zone=None):
        io.print('std', f"> export mesh ")
        # geo dimension from base
        self._geodim = self._file._h5file[self._bases[0]][" data"][0]
        if zone is None:
            assert len(self._zones)==1, "Multiple zones found, must specify which zone to export"
            name = list(self._zones.keys())[0]
        else:
            name = zone
        cgzone = cgnszone(self._zones[name], self._geodim)
        io.print('std', f"Parse zone {name} ({self._geodim}D) ncell: {cgzone.ncell}, nnode: {cgzone.nnode}")
        meshdata = Mesh(ncell=cgzone.ncell, nnode=cgzone.nnode)
        # get coordinates
        meshdata.set_nodescoord_xyz(*cgzone.coords())
        # # cell connectivity
        meshdata.set_cell2node(cgzone.export_cellcon())
        # boundary conditions
        for name, bc in cgzone._BCs.items():
            boco = cgzone.export_BC(bc)
            # filter full domain
            if boco.type in ['internal']:
                io.print('std', f"  filter internal mark {boco.name}")
            else:
                io.print('std', f"  add boco {boco}")
                meshdata.add_boco(boco)
        #meshdata.check()
        meshdata.printinfo()
        return meshdata

if __name__ == "__main__":
    f = cgnsMesh(filename="./examples/MESH.nogit/cavity-degen.hdf")
    f.read_data()
    f.printinfo()
    f.export_mesh()

# ElementType_t := Enumeration(
#      ElementTypeNull, ElementTypeUserDefined, NODE, BAR_2, BAR_3,
#      TRI_3, TRI_6, QUAD_4, QUAD_8, QUAD_9,
#      TETRA_4, TETRA_10, PYRA_5, PYRA_14,
#      PENTA_6, PENTA_15, PENTA_18, HEXA_8, HEXA_20, HEXA_27,
#      MIXED, PYRA_13, NGON_n, NFACE_n,
#      BAR_4, TRI_9, TRI_10, QUAD_12, QUAD_16,
#      TETRA_16, TETRA_20, PYRA_21, PYRA_29, PYRA_30,
#      PENTA_24, PENTA_38, PENTA_40, HEXA_32, HEXA_56, HEXA_64 );
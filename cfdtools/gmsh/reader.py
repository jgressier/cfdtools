# coding: utf8

# Import modules
import collections
import logging

import numpy as np

import cfdtools.meshbase._mesh as _mesh
import cfdtools.meshbase._connectivity as _conn
import cfdtools.meshbase._elements as _elem
from cfdtools.gmsh._gmsh import gmshtype_elem  # , nodes_per_cell

import cfdtools.api as api

log = logging.getLogger(__name__)


@api.fileformat_reader('GMSH', '.msh')
class reader(api._files):
    """Implementation of the reader to read Gmsh meshes."""

    def __init__(self, filename, cIntegrity=False):
        """Initialization of a GMSH reader.

        param: filename: file name [type string]
        """
        super().__init__(filename)
        self.check_integrity = cIntegrity

    @property
    def ncell(self):
        return self._ncell

    def read_data(self):
        log.info(f"> GMSH reader: starts reading {self.filename}")
        # Check file exists
        if not self.exists():
            api.error_stop(f"File not found: {self.filename!r}")

        # Set the node ordering convention
        _elem.elem2faces = _elem.gmsh_elem2faces

        # Read file.
        filename = self.filename
        fam, bctype, x, y, z, elts = self.__read_sections(filename)

        # fam: dict with index: family name
        # bctype: dict with name: type
        # elts: list of list[index, element type, x, x, x, nodes index ]
        log.info("Analyze...")

        # Check for 3D
        self._maxdim = np.amax([_elem.dim_elem[gmshtype_elem[e[1]]] for e in elts])
        # Define list of elements
        if self._maxdim < 2 or self._maxdim > 3:
            api.error_stop(f"unexpected max element dimension: {self._maxdim}")

        log.info("%dD mesh", self._maxdim)
        mesh_elt = [etype for etype, dim in _elem.dim_elem.items() if dim == self._maxdim]
        bc_elt = [etype for etype, dim in _elem.dim_elem.items() if dim == self._maxdim - 1]

        # Initialize returned variables
        boundaries = {}
        connectivity = collections.OrderedDict({})
        self._nboco = 0
        bnd, bnd2fam = collections.defaultdict(list), {}

        # Volume Connectivity
        self._ncell = 0
        log.info("  extract volume connectivity")
        bnd_connect = 1
        for elt in elts:
            if elt[1] in gmshtype_elem.keys():
                elt_type = gmshtype_elem[elt[1]]
                if elt_type in mesh_elt:
                    self._ncell += 1
                    if elt_type not in connectivity.keys():
                        connectivity[elt_type] = []
                    ntag = elt[2]  # expected to be 2
                    connectivity[elt_type].append(elt[3 + ntag :])
                    bnd_connect += 1
                else:
                    if elt_type not in bc_elt:
                        log.warning("%s not in available types %r.", elt_type, mesh_elt)
            else:
                log.warning("%d not in available types %r.", elt[1], gmshtype_elem.keys())

        # Boundary Connectivity
        log.info("  extract boundaries connectivity")
        connect_bc = {}
        for elt in elts:
            if elt[1] in gmshtype_elem.keys():
                elt_type = gmshtype_elem[elt[1]]
                if elt_type in bc_elt:
                    ntag = elt[2]  # expected to be 2
                    # bnd_tag = elt[2 + tag]
                    # bnd_fam = elt[2 + tag - 1]
                    bnd_tag = elt[4]
                    bnd_fam = elt[3]
                    if bnd_tag not in connect_bc.keys():
                        connect_bc[bnd_tag] = collections.defaultdict(list)
                    # if elt_type in connect_bc[bnd_tag].keys():
                    #     connect_bc[bnd_tag][elt_type].append(elt[3+ntag:])
                    # else:
                    #    connect_bc[bnd_tag][elt_type] = []
                    connect_bc[bnd_tag][elt_type].append(elt[3 + ntag :])
                    # if bnd_tag not in bnd.keys():
                    #     bnd[bnd_tag] = []
                    bnd[bnd_tag].append(bnd_connect)
                    # if bnd_fam not in bnd2fam:
                    bnd2fam[bnd_tag] = bnd_fam
                    bnd_connect += 1
        log.info(f"    tags are {bnd2fam}")
        # Reindex connectivities
        # Element-to-vertex
        for elt_type in connectivity.keys():
            connectivity[elt_type] = np.array(connectivity[elt_type]) - 1
        # Boundary patches element-to-vertex
        for bnd_tag in connect_bc.keys():
            for elt_type in connect_bc[bnd_tag].keys():
                connect_bc[bnd_tag][elt_type] = np.array(connect_bc[bnd_tag][elt_type]) - 1

        # Fill the dictionary of boundary conditions
        self._famm = []
        if fam is not None:  # B.C. are defined
            for bnd_tag in bnd.keys():
                family = fam[str(bnd2fam[bnd_tag])]
                if family != "fluid":
                    self.__create_bnd(
                        boundaries,
                        np.array(bnd[bnd_tag]),
                        family,
                        bctype,
                        connect_bc[bnd_tag],
                    )
        else:  # define B.C.
            for bnd_tag in bnd.keys():
                family = "fluid"
                self.__create_bnd(
                    boundaries,
                    np.array(bnd[bnd_tag]),
                    family,
                    bctype,
                    connect_bc[bnd_tag],
                )

        # Add the fluid boundary condition with the internal type now
        boundaries["int_fluid"] = {
            "slicing": None,
            "type": None,
            "periodic_transform": None,
        }
        boundaries["int_fluid"]["type"] = "internal"
        boundaries["int_fluid"]["periodic_transform"] = np.zeros((16,), dtype=np.float64)
        boundaries["int_fluid"]["slicing"] = []
        for elt_type in connectivity.keys():
            raveled = np.unique(connectivity[elt_type].ravel())
            boundaries["int_fluid"]["slicing"] += raveled.tolist()
        boundaries["int_fluid"]["slicing"] = np.array(boundaries["int_fluid"]["slicing"])

        self._elems = elts
        self._coords = (x, y, z)
        self._cellconnectivity = connectivity
        self._boundaries = boundaries

    def export_mesh(self):
        log.info("> export gmsh mesh to cfdtools mesh data")
        meshdata = _mesh.Mesh(nnode=len(self._coords[0]))
        meshdata.set_nodescoord_xyz(*self._coords)
        # meshdata.set_face2node(self.mesh['connectivity']['noofa'])
        cellconn = _conn.elem_connectivity()
        # extract cell connectivity only
        for etype, econn in self._cellconnectivity.items():
            if _elem.dim_elem[etype] == self._maxdim:
                cellconn.add_elems(etype, econn)
        meshdata.set_cell2node(cellconn)
        for name, bc_dict in self._boundaries.items():
            if bc_dict['type'] == 'boundary':
                boco = _mesh.submeshmark(name)
                boco.geodim = 'bdnode'
                boco.properties['type'] = bc_dict['type']
                boco.properties['periodic_transform'] = bc_dict['periodic_transform']
                boco.index = _conn.indexlist(ilist=bc_dict['slicing'])
                meshdata.add_boco(boco)
        # meshdata.set_celldata(self.variables['cells'])
        # meshdata.set_nodedata(self.variables['nodes'])
        # meshdata.set_facedata(self.variables['faces'])
        # meshdata.set_params(self.mesh['params'])
        return meshdata

    def __create_bnd(self, boundaries, window, family, bctype, connectivity):
        fff = family
        self._nboco += 1

        # Loop to concatenate multiple faces/entities into the same boundary (For complex geometries)

        if fff in self._famm:
            # Temporary list to be merged below
            bc2 = {"slicing": None}
            bc2["slicing"] = []
            # existing array is tranformed into a list and back to array
            bc = boundaries[family]
            bc["type"] = "None"
            bc["periodic_transform"] = "None"
            bc["slicing"] = bc["slicing"].tolist()
            for elt_type in connectivity.keys():
                raveled = np.unique(connectivity[elt_type].ravel())
                bc2["slicing"] += raveled.tolist()
            # log.info(bc["slicing"])
            bc["slicing"] += bc2["slicing"]
            bc["slicing"] = np.unique(bc["slicing"])
            bc["type"] = "boundary"
            bc["periodic_transform"] = np.zeros((16,), dtype=np.float64)

        else:
            boundaries[family] = {
                "slicing": None,
                "type": None,
                "periodic_transform": None,
            }

            bc = boundaries[family]
            bc["slicing"] = []
            for elt_type in connectivity.keys():
                raveled = np.unique(connectivity[elt_type].ravel())
                bc["slicing"] += raveled.tolist()
            bc["slicing"] = np.array(bc["slicing"])
            bc["type"] = "boundary"
            bc["periodic_transform"] = np.zeros((16,), dtype=np.float64)
            self._famm.append(fff)

    def __get_section_from_name(self, msh, sectionName):
        assert sectionName.startswith('$')
        ibeg = msh.index(["$" + sectionName[1:]])
        iend = msh.index(["$End" + sectionName[1:]])
        return msh[ibeg + 1 : iend]

    def __read_sections_V2(self, msh):
        # ------------------------------------------
        # Reading msh file for version 2.0 and below
        # ------------------------------------------

        log.info("Running version 2.0 reader")

        # Find the families.
        log.info("  parse Physical Names")
        if ["$PhysicalNames"] in msh:
            families = self.__get_section_from_name(msh, "$PhysicalNames")
            fam = {}
            bctype = {}
            for i in range(1, int(families[0][0]) + 1):
                fam[families[i][1]] = families[i][2][1:-1]
                bctype[families[i][2][1:-1]] = families[i][0]
            log.info(f"    found Physical names {fam}")
            log.info(f"    found Entities {bctype}")
        else:
            fam = None
            bctype = None

        # Find the coordinates.
        log.info("  parse Nodes")
        coords = self.__get_section_from_name(msh, "$Nodes")
        nnodes = int(coords[0][0])
        x = [float(coords[i][1]) for i in range(1, nnodes + 1)]
        y = [float(coords[i][2]) for i in range(1, nnodes + 1)]
        z = [float(coords[i][3]) for i in range(1, nnodes + 1)]

        # Find the elements.
        log.info("  parse Elements")
        elements = self.__get_section_from_name(msh, "$Elements")
        elts = []
        nelems = int(elements[0][0])
        for i in range(1, nelems + 1):
            elts.append([int(j) for j in elements[i]])

        return fam, bctype, x, y, z, elts

    def __read_sections_V4(self, msh):
        # ------------------------------------------
        # Reading msh file for version 4.0 and above
        # ------------------------------------------

        log.info("Running version 4.x reader")

        # Find the families.
        log.info("  parse Physical Names")
        if ["$PhysicalNames"] in msh:
            families = self.__get_section_from_name(msh, "$PhysicalNames")
            fam = {}
            bctype = {}
            for i in range(1, int(families[0][0]) + 1):
                fam[families[i][1]] = families[i][2][1:-1]
                bctype[families[i][2][1:-1]] = families[i][0]
            log.info(f"    found Physical names {fam}")
            log.info(f"    found Entities {bctype}")
        else:
            fam = None
            bctype = None

        # To find the entity number used for concatenation of the mesh.
        # The returned values corresponds to the physical group
        def find_ent(j):
            entities = self.__get_section_from_name(msh, "$Entities")
            addd = 0  # temporary variable to seek the line and return
            bot = int(entities[0][0]) + int(entities[0][1])
            top = int(entities[0][2]) + bot
            for i in range(bot, top + 1):
                if j == int(entities[i][0]) and len(entities[i]) >= 8:
                    addd = int(entities[i][8])
                    break
            return addd

        # Find the coordinates.
        log.info("  parse Nodes")
        coords = self.__get_section_from_name(msh, "$Nodes")
        nnodes = int(coords[0][3])
        nodes = 1
        count = 1
        counter = 1
        x = [None] * (nnodes)
        y = [None] * (nnodes)
        z = [None] * (nnodes)
        while nodes < nnodes:
            cnt = int(coords[counter][3])
            for i in range(counter + 1, counter + cnt + 1):
                pos = int(coords[i][0]) - 1
                x[pos] = float(coords[i + cnt][0])
                y[pos] = float(coords[i + cnt][1])
                z[pos] = float(coords[i + cnt][2])
            nodes += cnt
            counter = counter + 2 * cnt + 1

        # Find the elements.
        log.info("  parse Elements")
        elements = self.__get_section_from_name(msh, "$Elements")
        # header is: numEntityBlocks, numElements, minIndex, maxIndex
        elts = []
        counter = 1
        count = 1
        neltypes = int(elements[0][0])
        # entityblock: dimEntity EntityIndex ElemType numElements
        # elements: ElementIndex NodesIndex
        for i in range(0, neltypes):
            a = int(elements[count][2])  # dimension
            b = int(2)  # No. of tags
            d = int(elements[count][1])  # Entity group
            c = find_ent(d)  # Physical group
            nxtrange = int(elements[count][3])
            for j in range(count + 1, count + nxtrange + 1):
                elt1 = [int(k) for k in elements[j]]
                elt1.insert(1, a)
                elt1.insert(2, b)
                elt1.insert(3, c)
                elt1.insert(4, d)
                elts.append(elt1)
            count += nxtrange + 1

        return fam, bctype, x, y, z, elts

    def __read_sections(self, filename):
        # Read the entire mesh.
        msh = []
        log.info(f"Reading file {filename!r}...")
        with open(filename) as fid:
            for l in fid:
                msh.append(l.split())
        log.info(" done")

        # Find version of the GMSH used
        version = self.__get_section_from_name(msh, "$MeshFormat")
        self.version = int(float(version[0][0]))
        assert int(version[0][1]) == 0, "only ASCII version is supported"
        assert int(version[0][2]) == 8, "size of float must be 8 (64bits)"

        if self.version <= 2:
            return self.__read_sections_V2(msh)
        elif self.version >= 4:
            return self.__read_sections_V4(msh)
        else:
            api.error_stop(f"unexpected mesh version: {self.version} (expected <= 2 or >= 4)")

# coding: utf8

# Import modules
import collections
import cfdtools.meshbase._mesh as _mesh
import cfdtools.meshbase._connectivity as _conn
from cfdtools.gmsh._gmsh import gmshelt2canelt, nodes_per_cell
#import os

import cfdtools.api as api
import numpy as np

# from operator import itemgetter


@api.fileformat_reader('GMSH', '.msh')
class reader(api._files):
    """Implementation of the reader to read Gmsh meshes."""

    def read_data(self):
        api.io.print('std',f'GMSH reader: starts reading {self.filename}')
        # Check file exists
        if not self.exists():
            print("Fatal error. File %s cannot be found."%(self.filename))
            exit()

        # Read file.
        filename = self.filename
        fam, bctype, x, y, z, elts = self.__read_sections(filename)
        #print(fam, bctype, elts)
        # fam: dict with index: family name
        # bctype: dict with name: type
        # elts: list of list[index, element type, x, x, x, nodes index ]

        # Check for 3D
        dim = "2"
        for elt in elts:
            if elt[1] in gmshelt2canelt.keys():
                elt_type = gmshelt2canelt[elt[1]]
                if (
                    elt_type != "tri"
                    and elt_type != "qua"
                    and elt_type != "tri2"
                    and elt_type != "qua2"
                    and elt_type != "lin"
                    and elt_type != "lin2"
                    and elt_type != "node"
                ):
                    dim = "3"
                    break

        # Define list of element
        if dim == "3":
            mesh_elt = ["tet", "hex", "pri", "pyr", "tet2", "hex2", "pri2", "pyr2"]
            bc_elt = ["tri", "qua", "lin2", "tri2", "qua2"]
        else:
            mesh_elt = ["tri", "qua", "tri2", "qua2"]
            bc_elt = ["lin", "lin2"]

        # Initialize returned variables
        boundaries = {}
        connectivity = collections.OrderedDict({})
        families = {}

        # And some local stuff
        global bnd_idx
        bnd, bnd2fam = {}, {}

        # Volume Connectivity
        bnd_connect = 1
        for elt in elts:
            if elt[1] in gmshelt2canelt.keys():
                elt_type = gmshelt2canelt[elt[1]]
                if elt_type in mesh_elt:
                    if elt_type not in connectivity.keys():
                        connectivity[elt_type] = []
                    tag = elt[2]
                    connectivity[elt_type].append(elt[3 + tag :])
                    bnd_connect += 1
        # Boundary Connectivity
        connect_bc = {}
        for elt in elts:
            if elt[1] in gmshelt2canelt.keys():
                elt_type = gmshelt2canelt[elt[1]]
                if elt_type in bc_elt:
                    bnd_tag = elt[2 + tag]
                    bnd_fam = elt[2 + tag - 1]

                    if bnd_tag not in connect_bc.keys():
                        connect_bc[bnd_tag] = {}

                    if elt_type in connect_bc[bnd_tag].keys():
                        connect_bc[bnd_tag][elt_type].append(elt[3 + tag :])
                    else:
                        connect_bc[bnd_tag][elt_type] = []
                        connect_bc[bnd_tag][elt_type].append(elt[3 + tag :])

                    if bnd_tag not in bnd.keys():
                        bnd[bnd_tag] = []

                    bnd[bnd_tag].append(bnd_connect)
                    # if bnd_fam not in bnd2fam:
                    bnd2fam[bnd_tag] = bnd_fam
                    bnd_connect += 1
        # Reindex connectivities
        # Element-to-vertex
        for elt_type in connectivity.keys():
            connectivity[elt_type] = np.array(connectivity[elt_type]) - 1
        # Boundary patches element-to-vertex
        for bnd_tag in connect_bc.keys():
            for elt_type in connect_bc[bnd_tag].keys():
                connect_bc[bnd_tag][elt_type] = (
                    np.array(connect_bc[bnd_tag][elt_type]) - 1
                )

        # Fill the dictionary of boundary conditions
        global famm
        famm = []
        bnd_idx = 0
        global bound
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
        boundaries["int_fluid"]["periodic_transform"] = np.zeros(
            (16,), dtype=np.float64
        )
        boundaries["int_fluid"]["slicing"] = []
        for elt_type in connectivity.keys():
            raveled = np.unique(connectivity[elt_type].ravel())
            boundaries["int_fluid"]["slicing"] += raveled.tolist()
        boundaries["int_fluid"]["slicing"] = np.array(
            boundaries["int_fluid"]["slicing"]
        )
        # print(connectivity)
        meshdata = _mesh.mesh(len(elts), len(x))
        # np.array(zip(x, y, z)), connectivity, boundaries, None, None, None
        meshdata.set_nodescoord_xyz(x, y, z)
        # meshdata.set_face2cell(self.mesh['connectivity']['cvofa'])
        # meshdata.set_face2node(self.mesh['connectivity']['noofa'])
        meshdata.set_cell2node(connectivity)
        meshdata.set_bocos(boundaries)
        # meshdata.set_celldata(self.variables['cells'])
        # meshdata.set_nodedata(self.variables['nodes'])
        # meshdata.set_facedata(self.variables['faces'])
        # meshdata.set_params(self.mesh['params'])        
        return meshdata

    def __create_bnd(self, boundaries, window, family, bctype, connectivity):

        fff = family
        global bnd_idx
        bnd_idx += 1

        # Loop to concatenate multiple faces/entities into the same boundary (For complex geometries)

        if fff in famm:
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
            api.io.print('std',bc["slicing"])
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
            famm.append(fff)

    def __read_sections(self, filename):
        # Read the entire mesh.
        msh = []
        fid = open(filename)
        for l in fid:
            msh.append(l.split())

        # Find version of the GMSH used
        ibeg = msh.index(["$MeshFormat"])  # To be safe, use these indicators
        iend = msh.index(["$EndMeshFormat"])
        version = msh[ibeg + 1 : iend]
        global ver
        ver = int(float(version[0][0]))

        # -------------------------------------
        # Reading the msh file for version 2.0
        # -------------------------------------

        if ver <= 2:
            api.io.print('std',"--Running version 2.0 reader--")
            # Find the families.
            if ["$PhysicalNames"] in msh:
                ibeg = msh.index(["$PhysicalNames"])
                iend = msh.index(["$EndPhysicalNames"])
                families = msh[ibeg + 1 : iend]
                fam = {}
                bctype = {}
                for i in range(1, int(families[0][0]) + 1):
                    fam[families[i][1]] = families[i][2][1:-1]
                    bctype[families[i][2][1:-1]] = families[i][0]

            else:
                fam = None
                bctype = None
            # Find the coordinates.
            ibeg = msh.index(["$Nodes"])
            iend = msh.index(["$EndNodes"])
            coordinates = msh[ibeg + 1 : iend]
            x, y, z = [], [], []
            for i in range(1, int(coordinates[0][0]) + 1):
                x.append(float(coordinates[i][1]))
                y.append(float(coordinates[i][2]))
                z.append(float(coordinates[i][3]))
            # api.io.print('std',x)
            # Find the elements.
            ibeg = msh.index(["$Elements"])
            iend = msh.index(["$EndElements"])
            elements = msh[ibeg + 1 : iend]
            elts = []
            for i in range(1, int(elements[0][0]) + 1):
                elts.append([int(j) for j in elements[i]])

            fid.close()
            return fam, bctype, x, y, z, elts

        # ---------------------------------------
        # msh reading for version 4.0 and above
        # ---------------------------------------

        elif ver >= 4:
            api.io.print('std',"--Running 4.0 reader--")
            # Find the families.
            if ["$PhysicalNames"] in msh:
                ibeg = msh.index(["$PhysicalNames"])
                iend = msh.index(["$EndPhysicalNames"])
                families = msh[ibeg + 1 : iend]
                fam = {}
                bctype = {}
                for i in range(1, int(families[0][0]) + 1):
                    fam[families[i][1]] = families[i][2][1:-1]
                    bctype[families[i][2][1:-1]] = families[i][0]

            else:
                fam = None
                bctype = None
            # api.io.print('std',(fam, bctype))
            # To find the entitie number used for concantination of the
            # mesh. The returned values corresponds to the physical group
            def find_ent(j):
                ibeg = msh.index(["$Entities"])
                iend = msh.index(["$EndEntities"])
                entities = msh[ibeg + 1 : iend]
                addd = 0  # temporary variable to seek the line and return
                low = int(entities[0][0]) + int(entities[0][1])
                up = int(entities[0][2]) + low
                for i in range(low, up + 1):
                    if j == int(entities[i][0]) and len(entities[i]) >= 8:
                        addd = int(entities[i][8])
                        break
                return addd

                # Find the coordinates.

            ibeg = msh.index(["$Nodes"])
            iend = msh.index(["$EndNodes"])
            coordinates = msh[ibeg + 1 : iend]
            counter = 1
            maxnodes = int(coordinates[0][3])
            nodes = 1
            rng = []
            count = 1
            x = [None] * (maxnodes)
            y = [None] * (maxnodes)
            z = [None] * (maxnodes)
            while nodes < maxnodes:
                cnt = int(coordinates[counter][3])
                for i in range(counter + 1, counter + cnt + 1):
                    pos = int(coordinates[i][0]) - 1
                    x[pos] = float(coordinates[i + cnt][0])
                    y[pos] = float(coordinates[i + cnt][1])
                    z[pos] = float(coordinates[i + cnt][2])
                nodes += cnt
                counter = counter + 2 * cnt + 1
                # Find the elements.
            ibeg = msh.index(["$Elements"])
            iend = msh.index(["$EndElements"])
            elements = msh[ibeg + 1 : iend]
            elts = []
            # a, b, c, d = [], [], [], []
            counter = 1
            count = 1
            global eltts
            eltts = []
            totcomp = int(elements[0][0])
            # api.io.print('std',totcomp)
            totelm = int(elements[0][1])
            for i in range(0, totcomp):
                a = int(elements[count][2])  # dimension
                b = int(2)  # No. of tags
                d = int(elements[count][1])  # Entity group
                c = find_ent(d)  # Physical group
                nxtrange = int(elements[count][3])
                for j in range(count + 1, count + nxtrange + 1):
                    elt1 = [int(k) for k in elements[j]]
                    # elt2=[a, b, c, d]
                    # eltadd = elt2+elt1
                    # elt1.insert(0,y)
                    elt1.insert(1, a)
                    elt1.insert(2, b)
                    elt1.insert(3, c)
                    elt1.insert(4, d)
                    elts.append(elt1)

                count = count + nxtrange + 1
            elts = np.array(elts)
            # api.io.print('std',elts)
            fid.close()
            return fam, bctype, x, y, z, elts


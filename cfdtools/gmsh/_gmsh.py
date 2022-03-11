# coding: utf8

# Import modules
import collections
import cfdtools.meshbase._mesh as _mesh
import cfdtools.meshbase._connectivity as _conn
#import os

import cfdtools.api as api
import numpy as np

# from operator import itemgetter

# Gmsh element types to canonical types see description below after ReaderGmsh class object
gmshelt2canelt = {
    1: "lin",
    2: "tri",
    3: "qua",
    4: "tet",
    5: "hex",
    6: "pri",
    7: "pyr",
    8: "lin2",
    9: "tri2",
    10: "qua2",
    11: "tet2",
    12: "hex2",
    13: "pri2",
    14: "pyr2",
    15: "node",
}
# Actual number of vertices for a given cell type
nodes_per_cell = {
    "bi": 2,
    "tri": 3,
    "qua": 4,
    "tet": 4,
    "hex": 8,
    "pri": 6,
    "pyr": 5,
}

@api.fileformat_reader('GMSH', '.msh')
class reader(api._files):
    """Implementation of the reader to read Gmsh meshes."""

    def read_data(self):
        api.io.print('std',f'GMSH reader: starts reading {self.filename}')
        # Check file exists
        if not self._exists:
            print("Fatal error. File %s cannot be found."%(self.filename))
            exit()

        # Read file.
        filename = self.filename
        fam, bctype, x, y, z, elts = self.__read_sections(filename)

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

        meshdata = _mesh.mesh(len(elts), len(x))
        # np.array(zip(x, y, z)), connectivity, boundaries, None, None, None
        # meshdata.set_nodescoord_nd(self.mesh['coordinates'])
        # meshdata.set_face2cell(self.mesh['connectivity']['cvofa'])
        # meshdata.set_face2node(self.mesh['connectivity']['noofa'])
        # meshdata.set_bocos(self.mesh['bocos'])
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


# #================================================================================================
# #================================================================================================

# From GMSH doc -
# 1  : 2-node line.
# 2  : 3-node triangle.
# 3  : 4-node quadrangle.
# 4  : 4-node tetrahedron.
# 5  : 8-node hexahedron.
# 6  : 6-node prism.
# 7  : 5-node pyramid.
# 8  : 3-node second order line (2 nodes associated with the vertices and 1 with the edge).
# 9  : 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
# 10 : 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
# 11 : 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges).
# 12 : 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).
# 13 : 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).
# 14 : 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).
# 15 : 1-node point.
# 16 : 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).
# 17 : 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges).
# 18 : 15-node second order prism (6 nodes associated with the vertices and 9 with the edges).
# 19 : 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges).
# 20 : 9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)
# 21 : 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
# 22 : 12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)
# 23 : 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)
# 24 : 15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)
# 25 : 21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)
# 26 : 4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)
# 27 : 5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)
# 28 : 6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)
# 29 : 20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)
# 30 : 35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)
# 31 : 56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)
# 92 : 64-node third order hexahedron (8 nodes associated with the vertices, 24 with the edges, 24 with the faces, 8 in the volume)
# 93 : 125-node fourth order hexahedron (8 nodes associated with the vertices, 36 with the edges, 54 with the faces, 27 in the volume)

# Line:                   Line3:           Line4:

# 0----------1 --> u      0-----2----1     0----2----3----1

# Triangle:               Triangle6:          Triangle9/10:          Triangle12/15:

# v
# ^                                                                   2
# |                                                                   | \
# 2                       2                    2                      9   8
# |`\                     |`\                  | \                    |     \
# |  `\                   |  `\                7   6                 10 (14)  7
# |    `\                 5    `4              |     \                |         \
# |      `\               |      `\            8  (9)  5             11 (12) (13) 6
# |        `\             |        `\          |         \            |             \
# 0----------1 --> u      0-----3----1         0---3---4---1          0---3---4---5---1

# Quadrangle:            Quadrangle8:            Quadrangle9:

#       v
#       ^
#       |
# 3-----------2          3-----6-----2           3-----6-----2
# |     |     |          |           |           |           |
# |     |     |          |           |           |           |
# |     +---- | --> u    7           5           7     8     5
# |           |          |           |           |           |
# |           |          |           |           |           |
# 0-----------1          0-----4-----1           0-----4-----1

# Tetrahedron:                          Tetrahedron10:

#                    v
#                  .
#                ,/
#               /
#            2                                     2
#          ,/|`\                                 ,/|`\
#        ,/  |  `\                             ,/  |  `\
#      ,/    '.   `\                         ,6    '.   `5
#    ,/       |     `\                     ,/       8     `\
#  ,/         |       `\                 ,/         |       `\
# 0-----------'.--------1 --> u         0--------4--'.--------1
#  `\.         |      ,/                 `\.         |      ,/
#     `\.      |    ,/                      `\.      |    ,9
#        `\.   '. ,/                           `7.   '. ,/
#           `\. |/                                `\. |/
#              `3                                    `3
#                 `\.
#                    ` w
# Hexahedron:             Hexahedron20:          Hexahedron27:

#        v
# 3----------2            3----13----2           3----13----2
# |\     ^   |\           |\         |\          |\         |\
# | \    |   | \          | 15       | 14        |15    24  | 14
# |  \   |   |  \         9  \       11 \        9  \ 20    11 \
# |   7------+---6        |   7----19+---6       |   7----19+---6
# |   |  +-- |-- | -> u   |   |      |   |       |22 |  26  | 23|
# 0---+---\--1   |        0---+-8----1   |       0---+-8----1   |
#  \  |    \  \  |         \  17      \  18       \ 17    25 \  18
#   \ |     \  \ |         10 |        12|        10 |  21    12|
#    \|      w  \|           \|         \|          \|         \|
#     4----------5            4----16----5           4----16----5

# Prism:                      Prism15:               Prism18:

#            w
#            ^
#            |
#            3                       3                      3
#          ,/|`\                   ,/|`\                  ,/|`\
#        ,/  |  `\               12  |  13              12  |  13
#      ,/    |    `\           ,/    |    `\          ,/    |    `\
#     4------+------5         4------14-----5        4------14-----5
#     |      |      |         |      8      |        |      8      |
#     |    ,/|`\    |         |      |      |        |    ,/|`\    |
#     |  ,/  |  `\  |         |      |      |        |  15  |  16  |
#     |,/    |    `\|         |      |      |        |,/    |    `\|
#    ,|      |      |\        10     |      11       10-----17-----11
#  ,/ |      0      | `\      |      0      |        |      0      |
# u   |    ,/ `\    |    v    |    ,/ `\    |        |    ,/ `\    |
#     |  ,/     `\  |         |  ,6     `7  |        |  ,6     `7  |
#     |,/         `\|         |,/         `\|        |,/         `\|
#     1-------------2         1------9------2        1------9------2

# Pyramid:                     Pyramid13:                   Pyramid14:

#                4                            4                            4
#              ,/|\                         ,/|\                         ,/|\
#            ,/ .'|\                      ,/ .'|\                      ,/ .'|\
#          ,/   | | \                   ,/   | | \                   ,/   | | \
#        ,/    .' | `.                ,/    .' | `.                ,/    .' | `.
#      ,/      |  '.  \             ,7      |  12  \             ,7      |  12  \
#    ,/       .' w |   \          ,/       .'   |   \          ,/       .'   |   \
#  ,/         |  ^ |    \       ,/         9    |    11      ,/         9    |    11
# 0----------.'--|-3    `.     0--------6-.'----3    `.     0--------6-.'----3    `.
#  `\        |   |  `\    \      `\        |      `\    \     `\        |      `\    \
#    `\     .'   +----`\ - \ -> v  `5     .'        10   \      `5     .' 13     10   \
#      `\   |    `\     `\  \        `\   |           `\  \       `\   |           `\  \
#        `\.'      `\     `\`          `\.'             `\`         `\.'             `\`
#           1----------------2            1--------8-------2           1--------8-------2
#                     `\
#                        u

###################################################################################################

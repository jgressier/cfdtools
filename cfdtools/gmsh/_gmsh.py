#!/usr/bin/env python

# Import modules
import collections
import numpy as np
import os
#from operator import itemgetter

# Gmsh element types to canonical types see description below after ReaderGmsh class object
gmshelt2canelt = {1: 'lin',
                  2: 'tri',
                  3: 'qua',
                  4: 'tet',
                  5: 'hex',
                  6: 'pri',
                  7: 'pyr',
                  8: 'lin2',
                  9: 'tri2',
                  10: 'qua2',
                  11: 'tet2',
                  12: 'hex2',
                  13: 'pri2',
                  14: 'pyr2',
                  15: 'node',}
# Actual number of vertices for a given cell type
nodes_per_cell = {'bi': 2,
                  'tri': 3,
                  'qua': 4,
                  'tet': 4,
                  'hex': 8,
                  'pri': 6,
                  'pyr': 5,}

class ReaderGmsh():
    '''Implementation of the reader to read Gmsh meshes.'''

    def __init__(self, filename):

      print
      print ":: READER GMSH ::"
      print

      self.filename = filename

    def read_data(self):

      # Check file exists
      if not os.path.isfile(self.filename):
        print "Fatal error. File %s cannot be found."%(self.filename)
        exit()

      # Read file.
      filename = self.filename
      fam, bctype, x, y, z, elts = self.__read_sections(filename)

	
	 
      #Check for 3D
      dim = '2'
      for elt in elts:
        if elt[1] in gmshelt2canelt.keys():
          elt_type = gmshelt2canelt[elt[1]]
          if elt_type != 'tri' and elt_type != 'qua' and \
            elt_type != 'tri2' and elt_type != 'qua2' and \
            elt_type != 'lin'  and elt_type != 'lin2' and \
            elt_type != 'node':
             dim = '3'
             break
	 
      # Define list of element
      if dim == '3':
        mesh_elt = ['tet', 'hex', 'pri', 'pyr',
                    'tet2', 'hex2', 'pri2', 'pyr2']
        bc_elt = ['tri', 'qua', 'lin2', 'tri2', 'qua2']
      else:
        mesh_elt = ['tri', 'qua', 'tri2', 'qua2']
        bc_elt = ['lin', 'lin2']

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
            connectivity[elt_type].append(elt[3+tag:])
            bnd_connect += 1
        # Boundary Connectivity
      connect_bc = {}
      for elt in elts:
        if elt[1] in gmshelt2canelt.keys():
            elt_type = gmshelt2canelt[elt[1]]
            if elt_type in bc_elt:
                bnd_tag = elt[2+tag]
                bnd_fam = elt[2+tag-1]

                if bnd_tag not in connect_bc.keys():
                    connect_bc[bnd_tag] = {}

                if elt_type in connect_bc[bnd_tag].keys():
                    connect_bc[bnd_tag][elt_type].append(elt[3+tag:])
                else:
                    connect_bc[bnd_tag][elt_type] = []
                    connect_bc[bnd_tag][elt_type].append(elt[3+tag:])

                if bnd_tag not in bnd.keys():
                    bnd[bnd_tag] = []

                bnd[bnd_tag].append(bnd_connect)
	        #if bnd_fam not in bnd2fam:
                bnd2fam[bnd_tag] = bnd_fam
                bnd_connect += 1
      # Reindex connectivities
      # Element-to-vertex
      for elt_type in connectivity.keys():
          connectivity[elt_type] = np.array(connectivity[elt_type])-1
      # Boundary patches element-to-vertex
      for bnd_tag in connect_bc.keys():
          for elt_type in connect_bc[bnd_tag].keys():
              connect_bc[bnd_tag][elt_type] = np.array(connect_bc[bnd_tag][elt_type])-1

      # Fill the dictionary of boundary conditions
      global famm
      famm = []
      bnd_idx = 0
      global bound
      if fam is not None: # B.C. are defined
          for bnd_tag in bnd.keys():
              family = fam[str(bnd2fam[bnd_tag])]
              if family != 'fluid':
                self.__create_bnd(boundaries, np.array(bnd[bnd_tag]), family,
                                  bctype, connect_bc[bnd_tag])
      else:  # define B.C.
          for bnd_tag in bnd.keys():
              family = 'fluid'
              self.__create_bnd(boundaries, np.array(bnd[bnd_tag]), family,
                                bctype, connect_bc[bnd_tag])

      # Add the fluid boundary condition with the internal type now
      boundaries["int_fluid"] = {'slicing':None, 'type':None, 'periodic_transform':None}
      boundaries["int_fluid"]["type"] = "internal"
      boundaries["int_fluid"]["periodic_transform"] = np.zeros((16,), dtype=np.float64)
      boundaries["int_fluid"]["slicing"] = []
      for elt_type in connectivity.keys():
          raveled = np.unique(connectivity[elt_type].ravel())
          boundaries["int_fluid"]["slicing"] += raveled.tolist()
      boundaries["int_fluid"]["slicing"] = np.array(boundaries["int_fluid"]["slicing"])

      return np.array(zip(x,y,z)), connectivity, boundaries, None, None, None

    def __create_bnd(self, boundaries, window, family, bctype, connectivity):
      
      fff = family
      global bnd_idx
      bnd_idx += 1

      #Loop to concatenate multiple faces/entities into the same boundary (For complex geometries)  
      
      if fff in famm:
	  #Temporary list to be merged below
	  bc2 = {'slicing':None}
	  bc2["slicing"] = []
	  #existing array is tranformed into a list and back to array
	  bc = boundaries[family]
	  bc['type'] = 'None'
          bc['periodic_transform'] = 'None'
	  bc["slicing"]=bc["slicing"].tolist()
	  for elt_type in connectivity.keys():
              raveled = np.unique(connectivity[elt_type].ravel())
              bc2["slicing"] += raveled.tolist()
	  print bc["slicing"]
          bc["slicing"] +=(bc2["slicing"])
	  bc["slicing"] = np.unique(bc["slicing"])
	  bc['type'] = 'boundary'
          bc['periodic_transform'] = np.zeros((16,), dtype=np.float64)

      else:
          boundaries[family] = {'slicing':None, 'type':None, 'periodic_transform':None}
	
          bc = boundaries[family]
          bc['slicing'] = []
          for elt_type in connectivity.keys():
              raveled = np.unique(connectivity[elt_type].ravel())
              bc["slicing"] += raveled.tolist()
          bc["slicing"] = np.array(bc["slicing"])
          bc['type'] = 'boundary'
          bc['periodic_transform'] = np.zeros((16,), dtype=np.float64)
          famm.append(fff)

	
    def __read_sections(self, filename):
        # Read the entire mesh.
        msh = []
        fid = open(filename)
        for l in fid:
            msh.append(l.split())
	
	    # Find version of the GMSH used
	ibeg = msh.index(['$MeshFormat']) # To be safe, use these indicators
        iend = msh.index(['$EndMeshFormat'])
	version = msh[ibeg+1:iend]
	global ver 
	ver = int(float(version[0][0]))
	    
	
	#-------------------------------------
 	# Reading the msh file for version 2.0
	#-------------------------------------

	if ver <= 2:
	  print "--Running version 2.0 reader--"
          # Find the families.
          if ['$PhysicalNames'] in msh:
             ibeg = msh.index(['$PhysicalNames'])
             iend = msh.index(['$EndPhysicalNames'])
             families = msh[ibeg+1:iend]
             fam = {}
             bctype = {}
             for i in range(1, int(families[0][0])+1):
                fam[families[i][1]] = families[i][2][1:-1]
                bctype[families[i][2][1:-1]] = families[i][0]

          else:
              fam = None
              bctype = None
         # Find the coordinates.
          ibeg = msh.index(['$Nodes'])
          iend = msh.index(['$EndNodes'])
          coordinates = msh[ibeg+1:iend]
          x, y, z = [], [], []	
          for i in range(1, int(coordinates[0][0])+1):	
              x.append(float(coordinates[i][1]))
              y.append(float(coordinates[i][2]))
              z.append(float(coordinates[i][3]))
	  #print x
         # Find the elements.
          ibeg = msh.index(['$Elements'])
          iend = msh.index(['$EndElements'])
          elements = msh[ibeg+1:iend]
          elts = []
          for i in range(1, int(elements[0][0])+1):
              elts.append([int(j) for j in elements[i]])

	  fid.close()
          return fam, bctype, x, y, z, elts


	#--------------------------------------- 
	# msh reading for version 4.0 and above  
	#---------------------------------------

	elif ver >= 4:
	  print "--Running 4.0 reader--"
          # Find the families.
          if ['$PhysicalNames'] in msh:
             ibeg = msh.index(['$PhysicalNames'])
             iend = msh.index(['$EndPhysicalNames'])
             families = msh[ibeg+1:iend]
             fam = {}
             bctype = {}
             for i in range(1, int(families[0][0])+1):
                fam[families[i][1]] = families[i][2][1:-1]
                bctype[families[i][2][1:-1]] = families[i][0]
             
	  else:	
             fam = None
             bctype = None
	  #print (fam, bctype)
	# To find the entitie number used for concantination of the 
	# mesh. The returned values corresponds to the physical group
	  def find_ent(j):
	      ibeg = msh.index(['$Entities'])
	      iend = msh.index(['$EndEntities'])
	      entities = msh[ibeg+1:iend]
	      addd = 0 #temporary variable to seek the line and return 	
	      low = int(entities[0][0])+int(entities[0][1])
	      up = int(entities[0][2])+low
	      for i in range(low, up+1):
	          if j == int(entities[i][0]) and len(entities[i])>=8:
	             addd = int(entities[i][8])
		     break
	      return addd;
	
        # Find the coordinates.
          ibeg = msh.index(['$Nodes'])
          iend = msh.index(['$EndNodes'])
          coordinates = msh[ibeg+1:iend]
          counter = 1
	  maxnodes = int(coordinates[0][3]) 
          nodes = 1
	  rng = []
	  count = 1
	  x = [None]*(maxnodes)
          y = [None]*(maxnodes)
          z = [None]*(maxnodes)
	  while nodes < maxnodes:
	        cnt = int(coordinates[counter][3])
	        for i in range (counter+1,counter+cnt+1):
		    pos = int(coordinates[i][0])-1
		    x[pos] = float(coordinates[i+cnt][0])
           	    y[pos] = float(coordinates[i+cnt][1])
           	    z[pos] = float(coordinates[i+cnt][2]) 
	        nodes += cnt
                counter = counter+2*cnt+1
          # Find the elements.
          ibeg = msh.index(['$Elements'])
          iend = msh.index(['$EndElements'])
          elements = msh[ibeg+1:iend]
          elts = []
	  #a, b, c, d = [], [], [], []
	  counter = 1
	  count = 1
	  global eltts
	  eltts = []
	  totcomp = int(elements[0][0])
	  #print totcomp 
	  totelm = int(elements[0][1])
	  for i in range(0,totcomp):
	      a = int(elements[count][2]) # dimension  
	      b = int(2)			  # No. of tags
	      d = int(elements[count][1]) # Entity group
	      c = find_ent(d)             # Physical group
	      nxtrange=int(elements[count][3])
	      for j in range(count+1,count+nxtrange+1):	
		  elt1=([int(k) for k in elements[j]])
		 # elt2=[a, b, c, d]
		 # eltadd = elt2+elt1
		  #elt1.insert(0,y)
		  elt1.insert(1,a)
                  elt1.insert(2,b)
                  elt1.insert(3,c)
                  elt1.insert(4,d)
		  elts.append(elt1)
		
	      count=count+nxtrange+1
          elts = np.array(elts)
	  #print elts
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

if __name__ == "__main__":
    '''
    The script is supposed to be used with command line arguments
    but if it is not, it runs a test on a pre-defined file name.
    '''

    # Module import for script use
    import copy
    import os
    import sys
    import vtk as _vtk
    from vtk.util import numpy_support

    # Parse command line arguments seeking a restart file name
    defaultFName = "gmsh_test.msh"
    if len(sys.argv) == 1:
        fName = defaultFName
        print "Use requires the name of a gmsh file to be input in argument."
        print "Using a default name, %s."%(defaultFName)
    else:
        fName = defaultFName
        for arg in sys.argv[1:]:
            if '--gmshname' in arg:
                fName = arg.split('=')[-1]
    if fName == defaultFName:
        print "Use requires the name of a gmsh file to be input in argument."
        print "Using a default name, %s."%(defaultFName)

    # Check given file name
    if not os.path.isfile(fName):
        print "Fatal error: %s cannot be found. Exiting."%(fName)
        exit()

    # Else, proceed
    readr = ReaderGmsh(fName)
    xyz, co, bocos, fam = readr.read_data()

    ###################################################################################################

    # Now, we can export the data under a more human-friendly format
    #
    #- IN VTK
    #
    # Build point coordinates for VTK
    nb_point = xyz.shape[0]
    if hasattr(_vtk, 'vtkSOADataArrayTemplate'):
        # VTK 8.1.0+: use SOA coordinate array data structure (zero-copy)
        points_vtk = _vtk.vtkSOADataArrayTemplate[np.float64]()  # shall avoid copy
        points_vtk.SetNumberOfComponents(3)
        points_vtk.charlesx_arrays = []
        for index in range(3):
            np_array = copy.deepcopy(xyz[:,index])
            points_vtk.SetArray(index, np_array, nb_point, True, True)
            points_vtk.charlesx_arrays.append(np_array)  # record array ref so it won't be gc'ed
    else:
        # VTK 8.0.1-: build a new AOS coordinate array
        points_vtk = numpy_support.numpy_to_vtk(xyz, deep=True)
    vtkpoints = _vtk.vtkPoints()
    vtkpoints.SetData(points_vtk)
    # Build connectivity for VTK
    vtkprimitive = {
        'bi': _vtk.vtkLine(),
        'tri': _vtk.vtkTriangle(),
        'qua': _vtk.vtkQuad(),
        'tet': _vtk.vtkTetra(),
        'hex': _vtk.vtkHexahedron(),
        'pri': _vtk.vtkWedge(),
        'pyr': _vtk.vtkPyramid(),
    }
    cells = np.empty((0, ), dtype=np.int64)
    cell_types = np.empty((0, ), dtype=np.int64)
    offsets = np.empty((0, ), dtype=np.int64)
    offset_start = 0
    total_nb_cells = 0

    for uns_type, connect in co.items():

        nb_cells = connect.shape[0]

        # record cell-type for each individual cell
        if uns_type in vtkprimitive:
            cell_type = vtkprimitive[uns_type].GetCellType()
        else:
            raise ValueError(str_error('Unknown cell type'))
        cell_types = np.append(cell_types, np.tile(cell_type, (nb_cells, 1)))

        # put number of vertices before each cell
        cells = np.append(cells, np.concatenate((np.tile(nodes_per_cell[uns_type], (nb_cells, 1)), connect),
                                                axis=1).flat)

        # start offset of each cell in connectivity array
        offset_stop = offset_start + nb_cells * (nodes_per_cell[uns_type] + 1)
        offsets = np.append(offsets, np.arange(offset_start, offset_stop, nodes_per_cell[uns_type] + 1))
        offset_start = offset_stop
        total_nb_cells += nb_cells
    idtype_vtk = _vtk.vtkIdTypeArray().GetDataType()
    uchartype_vtk = _vtk.vtkUnsignedCharArray().GetDataType()

    cells_vtk = numpy_support.numpy_to_vtk(cells, deep=True, array_type=idtype_vtk)
    cell_array = _vtk.vtkCellArray()
    cell_array.SetCells(total_nb_cells, cells_vtk)
    # Build VTK unstructured Grid
    if len(cell_types) == 1:
        vtk_obj = _vtk.vtkUnstructuredGrid()
        vtk_obj.SetPoints(vtkpoints)
        vtk_obj.SetCells(cell_types[0], cell_array)
    else:
        cell_types_vtk = numpy_support.numpy_to_vtk(cell_types, deep=True, array_type=uchartype_vtk)
        offsets_vtk = numpy_support.numpy_to_vtk(offsets, deep=True, array_type=idtype_vtk)
        vtk_obj = _vtk.vtkUnstructuredGrid()
        vtk_obj.SetPoints(vtkpoints)
        vtk_obj.SetCells(cell_types_vtk, offsets_vtk, cell_array)

    elemdata = {
        "nodes": vtk_obj.GetPointData(),
        "cells": vtk_obj.GetCellData(),
    }
    # And add the variables
    for idx, var in enumerate(['x','y','z']):
        np_array = eval(var)
        # cast without copy if possible
        if np_array.flags.contiguous:
            np_array = np_array.astype(np.float64, copy=False)
        else:
            np_array = np_array.astype(np.float64)
        # zero-copy if contiguous array
        vtkarray = numpy_support.numpy_to_vtk(np_array, deep=False)
        vtkarray.charlesx_array = np_array  # avoid garbage collection
        vtkarray.SetName(var)
        #
        elemdata["nodes"].AddArray(vtkarray)
    # Finally, write the file
    writer = _vtk.vtkXMLUnstructuredGridWriter()
    if _vtk.vtkVersion.GetVTKMajorVersion() >= 6:
        writer.SetInputData(vtk_obj)
    else:
        writer.SetInput(vtk_obj)
    filename = "from_gmsh.00000"
    writer.SetFileName(filename + '.vtu')
    writer.SetDataModeToBinary()
    writer.Write()

import numpy as np

geodim = {'point': 0, 'line': 1, 'surface': 2, 'volume': 3}

# Declare all element types here.
# itype, name of the element, number of nodes, geometry dimension, element given by the extrusion of the element
elem_properties = [
    (0, '   node1', 1, 0, '    bar2'),
    (0, '    bar2', 2, 1, '   quad4'),
    (0, '    tri3', 3, 2, '  penta6'),
    (0, '   quad4', 4, 2, '   hexa8'),
    (0, '  tetra4', 4, 3, '    none'),
    (0, '   pyra5', 5, 3, '    none'),
    (0, '  prism6', 6, 3, '    none'),
    (0, '   hexa8', 8, 3, '    none'),
]
# remove leading spaces of strings
elem_properties = [tuple(i.lstrip() if isinstance(i, str) else i for i in u) for u in elem_properties]

dim_elem = {e: d for _, e, _, d, _ in elem_properties}
nnode_elem = {e: n for _, e, n, _, _ in elem_properties}
face_from_nnode = {n: e for _, e, n, d, _ in elem_properties if d <= 2}
extruded_face = {e: x for _, e, _, d, x in elem_properties if d <= 2}

del elem_properties

###############################################################################
# CGNS convention uses the convention of outward normal for element faces.
# see https://cgns.github.io/CGNS_docs_current/sids/conv.html#unstructgrid.
# Copy the definition of the list of faces from an element type.
# Change the starting index of the ordering indices.
# Permute nodes to get the chosen normal orientation of faces (optional).


def set_starting_index_to_zero(face_nodes):
    # face_nodes gives nodes associated to a face as in the CGNS convention
    # take from cgns and remove 1 to indices to start at 0
    return [[int(x.replace('N', '')) - 1 for x in node] for node in face_nodes.values()]


def set_inward_normal(connectivity):
    # set inward normal: face connectivity respects a CGNS convention
    # permute second and fourth elements for quad elements.
    # permute second and third elements for tri elements.
    return [
        [elt[0], elt[2], elt[1]] if len(elt) == 3 else [elt[0], elt[3], elt[2], elt[1]]
        for elt in connectivity
    ]


def group_faces_by_type(connectivity):
    # only tri3 and quad4
    # no error check
    return {
        'tri3': [elt for elt in connectivity if len(elt) == 3],
        'quad4': [elt for elt in connectivity if len(elt) == 4],
    }


###############################################################################
# From CGNS
# see https://cgns.github.io/CGNS_docs_current/sids/conv.html#unst_quad
# QUAD 4 2D
cgns_elem2faces = {
    'quad4': {'bar2': [[0, 1], [1, 2], [2, 3], [3, 0]]},
}

###############################################################################
# From CGNS
# see https://cgns.github.io/CGNS_docs_current/sids/conv.html#unst_hexa
# HEXA 8 3D
face_corners = {
    # Face  Corner Nodes
    'F1': ['N1', 'N4', 'N3', 'N2'],
    'F2': ['N1', 'N2', 'N6', 'N5'],
    'F3': ['N2', 'N3', 'N7', 'N6'],
    'F4': ['N3', 'N4', 'N8', 'N7'],
    'F5': ['N1', 'N5', 'N8', 'N4'],
    'F6': ['N5', 'N6', 'N7', 'N8'],
}

connectivity_hexa8 = set_starting_index_to_zero(face_corners)
connectivity_hexa8 = set_inward_normal(connectivity_hexa8)

# List of quad4 faces for a hexa8 element.
cgns_elem2faces['hexa8'] = group_faces_by_type(connectivity_hexa8)

###############################################################################
# HEXA 27 3D
face_midedge = {
    # Face  Mid-Edge Nodes
    'F1': ['N12', 'N11', 'N10', 'N9'],
    'F2': ['N9', 'N14', 'N17', 'N13'],
    'F3': ['N10', 'N15', 'N18', 'N14'],
    'F4': ['N11', 'N16', 'N19', 'N15'],
    'F5': ['N13', 'N20', 'N16', 'N12'],
    'F6': ['N17', 'N18', 'N19', 'N20'],
}

connectivity_midedge_hexa27 = set_starting_index_to_zero(face_midedge)
# connectivity_midedge_hexa27 = set_inward_normal(connectivity_midedge_hexa27)

face_midface = {
    # Mid-Face Node
    'F1': ['N21'],
    'F2': ['N22'],
    'F3': ['N23'],
    'F4': ['N24'],
    'F5': ['N25'],
    'F6': ['N26'],
}

connectivity_midface_hexa27 = set_starting_index_to_zero(face_midface)

connectivity_hexa27 = [
    [*a, *b, *c]
    for a, b, c in zip(connectivity_hexa8, connectivity_midedge_hexa27, connectivity_midface_hexa27)
]

cgns_elem2faces['hexa27'] = {'quad9': connectivity_hexa27}

###############################################################################
# From CGNS
# see https://cgns.github.io/CGNS_docs_current/sids/conv.html#unst_penta
# PENTA 6 3D
face_corners = {
    # Face  Corner Nodes
    'F1': ['N1', 'N2', 'N5', 'N4'],
    'F2': ['N2', 'N3', 'N6', 'N5'],
    'F3': ['N3', 'N1', 'N4', 'N6'],
    'F4': ['N1', 'N3', 'N2'],
    'F5': ['N4', 'N5', 'N6'],
}

connectivity_prism6 = set_starting_index_to_zero(face_corners)
connectivity_prism6 = set_inward_normal(connectivity_prism6)

# List of quad4 and tri3 faces for a prism6 element.
cgns_elem2faces['prism6'] = group_faces_by_type(connectivity_prism6)

###############################################################################
# From CGNS
# see https://cgns.github.io/CGNS_docs_current/sids/conv.html#unst_tetra
# TETRA 4 3D
face_corners = {
    # Face  Corner Nodes
    'F1': ['N1', 'N3', 'N2'],
    'F2': ['N1', 'N2', 'N4'],
    'F3': ['N2', 'N3', 'N4'],
    'F4': ['N3', 'N1', 'N4'],
}

connectivity_tetra4 = set_starting_index_to_zero(face_corners)
connectivity_tetra4 = set_inward_normal(connectivity_tetra4)

# List of tri3 faces for a tetra4 element.
cgns_elem2faces['tetra4'] = group_faces_by_type(connectivity_tetra4)

###############################################################################
# From CGNS
# https://cgns.github.io/CGNS_docs_current/sids/conv.html#unst_pyramid
# PYRA 5 3D
face_corners = {
    # Face  Corner Nodes
    'F1': ['N1', 'N4', 'N3', 'N2'],
    'F2': ['N1', 'N2', 'N5'],
    'F3': ['N2', 'N3', 'N5'],
    'F4': ['N3', 'N4', 'N5'],
    'F5': ['N4', 'N1', 'N5'],
}

connectivity_pyra5 = set_starting_index_to_zero(face_corners)
connectivity_pyra5 = set_inward_normal(connectivity_pyra5)

# List of quad4 and tri3 faces for a pyra5 element.
cgns_elem2faces['pyra5'] = group_faces_by_type(connectivity_pyra5)

###############################################################################
# Nodes are ordered with the CGNS convention
# Transformation to go from the CGNS convention to the GMSH convention
# https://gmsh.info/doc/texinfo/gmsh.html#Low-order-elements
hexa8_cgns2gmsh = [1, 2, 3, 4, 5, 6, 7, 8]
# fmt: off
hexa27_cgns2gmsh = [1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 14, 10, 11, 13, 15, 16, 17, 19, 20, 18, 21, 22, 24, 25, 23, 26, 27]
prism6_cgns2gmsh = [1, 2, 3, 4, 5, 6]
tetra4_cgns2gmsh = [1, 2, 3, 4]
pyra5_cgns2gmsh = [1, 2, 3, 4, 5]
###############################################################################
# Transformation to go from the GMSH convention to the CGNS convention
gmsh2cgns = {
    'tri3': [1, 2, 3],
    'quad4': [1, 2, 3, 4],
    'quad9': [1, 2, 3, 4, 5, 6, 7, 8, 9],
    'tetra4': [1, 2, 3, 4],
    'prism6': [1, 2, 3, 4, 5, 6],
    'pyra5': [1, 2, 3, 4, 5],
    'hexa8': [1, 2, 3, 4, 5, 6, 7, 8],
    'hexa27': [1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 10, 14, 11, 15, 16, 17, 20, 18, 19, 21, 22, 25, 23, 24, 26, 27]
}


def cgns2gmsh(element_cgns2gmsh, face_connectivity):
    # starting index of face_connectivity must be zero
    # no error check
    return [np.array([i - 1 for i in element_cgns2gmsh])[face].tolist() for face in face_connectivity]

###############################################################################
# fmt: on


gmsh_connectivity_hexa8 = cgns2gmsh(hexa8_cgns2gmsh, connectivity_hexa8)
gmsh_connectivity_hexa27 = cgns2gmsh(hexa27_cgns2gmsh, connectivity_hexa27)
gmsh_connectivity_prism6 = cgns2gmsh(prism6_cgns2gmsh, connectivity_prism6)
gmsh_connectivity_tetra4 = cgns2gmsh(tetra4_cgns2gmsh, connectivity_tetra4)
gmsh_connectivity_pyra5 = cgns2gmsh(pyra5_cgns2gmsh, connectivity_pyra5)

gmsh_elem2faces = {
    'tri3': {'bar2': [[0, 1], [1, 2], [2, 0]]},
    'quad4': {'bar2': [[0, 1], [1, 2], [2, 3], [3, 0]]},
    'hexa8': group_faces_by_type(gmsh_connectivity_hexa8),
    'hexa27': {'quad9': gmsh_connectivity_hexa27},
    'prism6': group_faces_by_type(gmsh_connectivity_prism6),
    'pyra5': group_faces_by_type(gmsh_connectivity_pyra5),
    'tetra4': group_faces_by_type(gmsh_connectivity_tetra4),
}

# either CGNS or GMSH to be set in the reader. Default is CGNS.
elem2faces = cgns_elem2faces

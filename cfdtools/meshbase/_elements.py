# import cfdtools.api as api

geodim = {'point': 0, 'line': 1, 'surface': 2, 'volume': 3}

elem_properties = [  # itype, name, nnode, geodim, extruded
    (0, '   node1', 1, 0, '    bar2'),
    (0, '    bar2', 2, 1, '   quad4'),
    (0, '    tri3', 3, 2, '  penta6'),
    (0, '   quad4', 4, 2, '   hexa8'),
    (0, '  tetra4', 4, 3, '    none'),
    (0, '   hexa8', 8, 3, '    none'),
]

elem_properties = [tuple(i.lstrip() if isinstance(i, str) else i for i in u) for u in elem_properties]

elem_dim = {e: d for _, e, _, d, _ in elem_properties}
nnode_elem = {e: n for _, e, n, _, _ in elem_properties}
face_from_nnode = {n: e for _, e, n, d, _ in elem_properties if d <= 2}
extruded_face = {e: x for _, e, _, d, x in elem_properties if d <= 2}

del elem_properties

# define list of faces from an element type
#   faces are defined with inward normal
# see https://cgns.github.io/CGNS_docs_current/sids/conv.html#unst_hexa
elem2faces = {
    'quad4': {'bar2': [[0, 1], [1, 2], [2, 3], [3, 0]]},
    'hexa8': {'quad4': [[0, 1, 2, 3], [4, 7, 6, 5], [0, 4, 5, 1], [2, 6, 7, 3], [0, 3, 7, 4], [1, 5, 6, 2]]},
}

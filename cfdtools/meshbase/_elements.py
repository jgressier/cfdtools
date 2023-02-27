import cfdtools.api as api

geodim = { 'point': 0, 'line':1, 'surface':2, 'volume':3}

# TODO: set properties (dim, nnodes) for all elems and build reverse map

elem_dim = { e: d for d, elems in [
    (0, ['node1']), (1, ['bar2']),
    (2, ['tri3', 'quad4']),
    (3, ['tetra4', 'hexa8']) ] for e in elems }

# define list of faces from an element type
#   faces are defined with inward normal 
# see https://cgns.github.io/CGNS_docs_current/sids/conv.html#unst_hexa
elem2faces = {
    'quad4' : {
        'bar2' : [ [0, 1], [1, 2], [2, 3], [3, 0]]
    },
    'hexa8' : {
        'quad4' : [ [0, 1, 2, 3], [4, 7, 6, 5], 
                    [0, 4, 5, 1], [2, 6, 7, 3],
                    [0, 3, 7, 4], [1, 5, 6, 2] ]
    }
}
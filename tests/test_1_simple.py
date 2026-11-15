# import pytest
import cfdtools.meshbase.simple as simplemesh
from itertools import product
import numpy as np
from pathlib import Path
# import filecmp

# _datadir = Path("./tests/data")
# _builddir = Path("./tests/build")


def test_cube_mini():
    cube = simplemesh.Cube(2, 2, 2)
    rmesh = cube.export_mesh()
    rmesh.printinfo()
    assert rmesh.check()


def test_cube_large():
    cube = simplemesh.Cube(100, 100, 100)
    rmesh = cube.export_mesh()
    rmesh.printinfo()
    assert rmesh.check()


def test_cube_vectorized_connectivity_keeps_structured_order():
    mesh = simplemesh.Cube(2, 1, 1).export_mesh()
    assert np.array_equal(
        mesh._cell2node['hexa8'],
        [[0, 4, 6, 2, 1, 5, 7, 3], [4, 8, 10, 6, 5, 9, 11, 7]],
    )


def test_cube_vectorized_generation_matches_index_formulas():
    cube = simplemesh.Cube(2, 3, 4)
    mesh = cube.export_mesh()
    element_ijk = np.array(list(product(range(2), range(3), range(4))))
    ei, ej, ek = element_ijk.T
    expected_cells = np.column_stack(
        [
            cube.nodeglobindex_ijk(ei, ej, ek),
            cube.nodeglobindex_ijk(ei + 1, ej, ek),
            cube.nodeglobindex_ijk(ei + 1, ej + 1, ek),
            cube.nodeglobindex_ijk(ei, ej + 1, ek),
            cube.nodeglobindex_ijk(ei, ej, ek + 1),
            cube.nodeglobindex_ijk(ei + 1, ej, ek + 1),
            cube.nodeglobindex_ijk(ei + 1, ej + 1, ek + 1),
            cube.nodeglobindex_ijk(ei, ej + 1, ek + 1),
        ]
    )
    assert np.array_equal(mesh._cell2node['hexa8'], expected_cells)

    node_ijk = np.array(list(product(range(3), range(4), range(5))))
    expected_coordinates = np.column_stack(cube.nodes_xyz(*node_ijk.T))
    assert np.allclose(mesh.nodescoord(ndarray=True), expected_coordinates)


def test_cube_vectorized_boundary_marks_keep_historical_order():
    cube = simplemesh.Cube(2, 3, 4)
    mesh = cube.export_mesh()
    boundary_ijk = {
        'imin': product([0], range(4), range(5)),
        'imax': product([2], range(4), range(5)),
        'jmin': product(range(3), [0], range(5)),
        'jmax': product(range(3), [3], range(5)),
        'kmin': product(range(3), range(4), [0]),
        'kmax': product(range(3), range(4), [4]),
    }
    for name, ijk in boundary_ijk.items():
        ni, nj, nk = np.array(list(ijk)).T
        expected = cube.nodeglobindex_ijk(ni, nj, nk).tolist()
        mark = mesh.seekmark(name)
        assert mark.nodebased()
        assert mark.index.list() == expected

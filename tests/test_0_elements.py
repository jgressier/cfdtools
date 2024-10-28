import cfdtools.meshbase._elements as _elem
import numpy as np
import pytest


def test_starting_index_to_zero():
    face_corners = {
        'F1': ['N1', 'N4', 'N3', 'N2'],
        'F2': ['N1', 'N2', 'N6'],
        'F3': ['N2', 'N7'],
        'F4': ['N8'],
    }
    res = _elem.set_starting_index_to_zero(face_corners)
    assert res == [[0, 3, 2, 1], [0, 1, 5], [1, 6], [7]]


def test_inward_normal():
    connectivity = [[1, 2, 3], [4, 3, 2, 1]]
    res = _elem.set_inward_normal(connectivity)
    assert res == [[1, 3, 2], [4, 1, 2, 3]]


def test_group_faces_by_type():
    connectivity = [[1, 2, 3], [4, 3, 2, 1]]
    res = _elem.group_faces_by_type(connectivity)
    assert res == {'tri3': [[1, 2, 3]], 'quad4': [[4, 3, 2, 1]]}


def test_cgns2gmsh():
    res = _elem.cgns2gmsh([2, 4, 3, 1], [[0, 1, 2, 3], [3, 2, 1, 0]])
    assert res == [[1, 3, 2, 0], [0, 2, 3, 1]]

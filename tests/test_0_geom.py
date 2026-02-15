import cfdtools.meshbase._geom as geom
import numpy as np
import pytest


def test_Nodes_2d():
    nset = geom.Nodes(np.array([[1., 2.], [3., 2.], [5., 8.]]))
    assert np.allclose(nset.center, [3., 4.])
    # # translate    
    # nset += np.array([10., 10., 10]) # error: array shapes
    # assert np.allclose(nset.center, [13., 14.])


def test_Nodes_3d():
    nset = geom.Nodes(np.array([[1., 2., 10.], [3., 2., 10.], [5., 8., 10.]]))
    assert np.allclose(nset.center, [3., 4., 10.])
    # # translate    
    # nset += 10.  # error: does not add
    # assert np.allclose(nset.center, [13., 14., 20.])


@pytest.mark.parametrize('method', ['ordered']) #, 'unordered'])
def test_Nodes_rotation_ordered(method):
    n1 = geom.Nodes(np.array([[1., 2., 10.], [3., 2., 10.], [5., -8., 10.]]))
    n2 = geom.Nodes(np.array([[-2., 1., 10.], [-2., 3., 10.], [8., 5., 10.]]))
    rotm, error = n1.estimate_rotation(n2, method=method)
    assert error < 1.E-6
    assert np.allclose(rotm.as_rotvec(degrees=True), [0., 0., 90.])


# TO BE DEVELOPPED
# @pytest.mark.parametrize('method', ['unordered'])
# def test_Nodes_rotation_unordered(method):
#     n1 = geom.Nodes(np.array([[5., -8., 10.], [1., 2., 10.], [3., 2., 10.], ]))
#     n2 = geom.Nodes(np.array([[-2., 1., 10.], [-2., 3., 10.], [8., 5., 10.]]))
#     rotm, error = n1.estimate_rotation(n2, method=method)
#     print(error, rotm.as_rotvec(degrees=True))

import cfdtools.utils.polybase as poly
import numpy as np
import pytest

def test_initlegendre():
    order = 3
    b = poly.poly1d(order, series='legendre')
    assert b._name == 'legendre'
    assert len(b._basis) == order

def test_initpolynomial():
    order = 3
    b = poly.poly1d(3, series='polynomial')
    assert b._name == 'polynomial'
    assert len(b._basis) == order

@pytest.mark.parametrize('series', ['legendre', 'polynomial'])
def test_consistency(series):
    for order in range(2, 10):
        b = poly.poly1d(3, series=series)
        for pos in [-1, 0., .5, 1.]:
            assert sum(b.interp_weights(pos)) == pytest.approx(1.)
            assert sum(b.diff_weights(pos)) == pytest.approx(0.)

def test_compare_interp():
    order = 5
    b1 = poly.poly1d(3, series='legendre')
    b2 = poly.poly1d(3, series='polynomial')
    assert np.allclose(b1.interp_weights(-1), b2.interp_weights(-1))

def test_compare_diff():
    order = 5
    b1 = poly.poly1d(3, series='legendre')
    b2 = poly.poly1d(3, series='polynomial')
    # interpolation coefficients must be the same, non dependant of poly basis
    assert np.allclose(b1.diff_weights(-1), b2.diff_weights(-1))
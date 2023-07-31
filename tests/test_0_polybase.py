import cfdtools.utils.polybase as poly
import numpy as np
import pytest

@pytest.mark.parametrize('series', ['legendre', 'polynomial'])
def test_init(series):
    order = 3
    b = poly.poly1d(order, series=series)
    assert b._name == series
    assert len(b._basis) == order

@pytest.mark.parametrize('series', ['legendre', 'polynomial'])
def test_consistency(series):
    for order in range(2, 10):
        b = poly.poly1d(order, series=series)
        for pos in [-1.0, 0.0, 0.5, 1.0]:
            assert sum(b.value_weights(pos)) == pytest.approx(1.0, rel=1e-6)
            assert sum(b.gradv_weights(pos)) == pytest.approx(0.0, abs=1e-12)

def test_compare_value():
    order = 5
    b1 = poly.poly1d(order, series='legendre')
    b2 = poly.poly1d(order, series='polynomial')
    # interpolation coefficients must be the same, non dependant of poly basis
    assert np.allclose(b1.value_weights(-1), b2.value_weights(-1))
    assert b1.value_weights(-1) == pytest.approx(b2.value_weights(-1))

def test_compare_gradv():
    order = 5
    b1 = poly.poly1d(order, series='legendre')
    b2 = poly.poly1d(order, series='polynomial')
    # differentiation coefficients must be the same, non dependant of poly basis
    assert np.allclose(b1.gradv_weights(-1), b2.gradv_weights(-1))
    assert b1.gradv_weights(-1) == pytest.approx(b2.gradv_weights(-1))

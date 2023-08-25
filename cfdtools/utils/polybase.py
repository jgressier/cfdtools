try:
    from ..utils._dev import lazyprop
except ImportError:
    from _dev import lazyprop
import numpy as np
import numpy.linalg as lin
import numpy.polynomial.legendre as polyL
import numpy.polynomial.polynomial as poly
from scipy.special import roots_legendre  # , eval_legendre

polybasis_dict = {'legendre': polyL.Legendre, 'polynomial': poly.Polynomial}


class poly1d:
    def __init__(self, order: int, series='legendre') -> None:
        self._name = series
        self._order = order
        self._roots, self._wghts = roots_legendre(order)
        self._basis = [(0,) * i + (1,) for i in range(order)]
        self._polyeval = polybasis_dict[series]
        self._gradbasis = [self._polyeval(ibas).deriv().coef for ibas in self._basis]

    @lazyprop
    def vdm(self):
        """vandermonde matrix for base to Lagrange"""
        return np.array([self._polyeval(ibas)(self._roots) for ibas in self._basis])

    @lazyprop
    def invvdm(self):
        return lin.inv(self.vdm)

    def interp_weights(self, x):
        return self.invvdm @ np.array([self._polyeval(ibas)(x) for ibas in self._basis])

    value_weights = interp_weights

    def diff_weights(self, x):
        return self.invvdm @ np.array([self._polyeval(ibas)(x) for ibas in self._gradbasis])

    gradv_weights = diff_weights


if __name__ == "__main__":
    with np.printoptions(precision=12, floatmode='fixed', sign=' '):
        # b = poly1d(3, series='legendre')
        # b = poly1d(3, series='polynomial')
        for p in polybasis_dict:
            b = poly1d(4, series=p)
            print(f"    basis:  {b._basis}")
            print(f"    roots:  {b._roots}")
            print(f"    wghts:  {b._wghts}")
            print(f"(-1)value:  {b.value_weights(-1)}")
            print(f"(+1)value:  {b.value_weights(+1)}")
            print(f"(-1)gradv:  {b.gradv_weights(-1)}")
            print(f"(+1)gradv:  {b.gradv_weights(+1)}")

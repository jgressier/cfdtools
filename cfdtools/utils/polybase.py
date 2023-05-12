from ..utils._dev import lazyprop
import numpy as np
import numpy.linalg as lin
import numpy.polynomial.legendre as polyL
import numpy.polynomial.polynomial as poly
from scipy.special import roots_legendre, eval_legendre

polybasis_dict = {
    'legendre' : polyL.Legendre,
    'polynomial' : poly.Polynomial
}

class poly1d():
    def __init__(self, order: int, series='legendre') -> None:
        self._name = series
        self._order = order
        self._roots, self._weights = roots_legendre(order)
        self._basis = [ (0,)*i+(1,) for i in range(order) ]
        self._polyeval = polybasis_dict[series]
        self._gradbasis = [
            self._polyeval(ibas).deriv().coef
                for ibas in self._basis
        ]

    @lazyprop
    def vdm(self):
        """vandermonde matrix for base to Lagrange
        
        """
        return np.array([
            self._polyeval(ibas)(self._roots)
                for ibas in self._basis
        ])

    @lazyprop
    def invvdm(self):
        return lin.inv(self.vdm)

    def interp_weights(self, x):
        return self.invvdm @ np.array([
            self._polyeval(ibas)(x)
                for ibas in self._basis
        ])

    def diff_weights(self, x):
        return self.invvdm @ np.array([
            self._polyeval(ibas)(x)
                for ibas in self._gradbasis
        ])


if __name__ == "__main__":
    b = poly1d(3, series='legendre')
    print(b._roots)
    #print(b._basis)
    print("interp", b.interp_weights(-1))
    print("diff", b.diff_weights(-1))

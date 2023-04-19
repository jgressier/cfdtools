import cfdtools.api as api



# class indexeddata:
#     def __init__(self, type='nodal', index='direct'):
#         self._type = type
#         self._index = index

_available_Xrep = ('nodal', 'cellaverage', 'spectralcell')
_available_Trep = ('instant', 'timeevol', 'fourier', 'pod')


class DataSet():

    def __init__(self, Xrep='cellaverage', ndof=1, Trep='instant'):
        self._ndof = ndof
        self._data = dict()

    @property
    def ndof(self):
        return self._ndof

    @ndof.setter
    def ndof(self, ndof):
        assert ndof == 1 or self._Xrep == 'spectralcell'
        self._ndof = ndof

    @property
    def Xrep(self):
        return self._Xrep

    @Xrep.setter
    def Xrep(self, Xrep):
        self._Xrep = Xrep

    @property
    def Trep(self):
        return self._Trep

    @Trep.setter
    def Trep(self, Trep):
        self._Trep = Trep

    def add_data(self, name, data, time=None):
        self._data[name] = data

    def data(self):
        return self._data
    
    def items(self):
        return self._data.items()

    def __getitem__(self, name):
        return self._data[name]


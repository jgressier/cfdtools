import cfdtools.api as api
import numpy as np

_available_Xrep = ('nodal', 'cellaverage', 'spectralcell')
_available_Trep = ('instant', 'timeevol', 'fourier', 'pod')

class DataSet():

    def __init__(self, Xrep='cellaverage', ndof=1, Trep='instant'):
        self._ndof = ndof
        self._data = dict()
        self.Trep = Trep
        self.Xrep = Xrep
        self._properties = dict()

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
        assert Xrep in _available_Xrep
        self._Xrep = Xrep

    @property
    def Trep(self):
        return self._Trep

    @Trep.setter
    def Trep(self, Trep):
        assert Trep in ('instant')
        self._Trep = Trep

    def add_data(self, name, data, time=None):
        if time:
            self._properties['time'] = time
        self._data[name] = data

    def data(self):
        return self._data
    
    def keys(self):
        return self._data.keys()

    def items(self):
        return self._data.items()

    def __getitem__(self, name):
        return self._data[name]

# class manufactured_DataSet(DataSet):

#     def __init__(self, datafunctions, Xrep='cellaverage', ndof=1, Trep='instant'):
#         super().__init__(Xrep, ndof, Trep)

class DataSetList(DataSet):

    def __init__(self, ndataset, Xrep='cellaverage', ndof=1, Trep='timeevol'):
        super().__init__(Xrep, ndof, Trep)
        self._ndataset = ndataset

    @property
    def Trep(self):
        return self._Trep

    @Trep.setter
    def Trep(self, Trep):
        assert Trep in _available_Trep
        self._Trep = Trep


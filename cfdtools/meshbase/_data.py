import cfdtools.api as api
import h5py
import numpy as np


class DataSetBase():

    _available_Xrep = ('nodal', 'cellaverage')
    _available_Trep = ('instant')

    def __init__(self, Xrep='cellaverage', ndof=1, Trep='instant'):
        self._ndof = ndof
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
        assert Xrep in self._available_Xrep
        self._Xrep = Xrep

    @property
    def Trep(self):
        return self._Trep

    @Trep.setter
    def Trep(self, Trep):
        # self._available_Xrep is essential to self-adapt to class
        assert Trep in self._available_Trep 
        self._Trep = Trep

class DataSet(DataSetBase):

    _available_Xrep = ('nodal', 'cellaverage', 'spectralcell')
    _available_Trep = ('instant', 'timeevol', 'fourier', 'pod')

    def __init__(self, Xrep='cellaverage', ndof=1, Trep='instant'):
        super().__init__(Xrep, ndof, Trep)
        self._data = dict()

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

class DataSetList(DataSetBase):

    _available_Xrep = ('nodal', 'cellaverage', 'spectralcell')
    _available_Trep = ('instant', 'timeevol', 'pod')

    def __init__(self, ndataset, Xrep='cellaverage', ndof=1, Trep='timeevol'):
        super().__init__(Xrep, ndof, Trep)
        self._ndataset = ndataset
        self._datalist = []

    def add_datalist(self, datalist: dict, time=None):
        if time:
            datalist['time'] = time
        self._datalist.append(datalist)

    def dumphdf(self, hgroup: h5py.Group,  options={}):
        n = len(self._datalist)
        for i, datadict in enumerate(self._datalist):
            datagroup = hgroup.create_group(f"i{i:06}")
            for vname, var in datadict.items():
                datagroup.create_dataset(vname, data=var, **options)
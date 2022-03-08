import cfdtools.api as api

_data_representation = []

class indexeddata():
    def __init__(self, type='nodal', index='direct'):
        self._type = type
        self._index = index

class celldata(indexeddata):
    def __init__(self, type='cellaverage', index='direct', ndof=1):
        super().__init__(type, index)
        self._ndof = ndof

    def set_data(self, data):
        self._data = data

    def data(self):
        return self._data

    def ndof(self):
        return self._ndof

    def __str__(self):
        return "type: {}\nndof: {}".format(self._type, self._ndof)
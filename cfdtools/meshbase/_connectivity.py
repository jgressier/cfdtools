import cfdtools.api as api

class genericindex():
    __available_types = ['list', 'range']
    
    def __init__(self):
        self._type = None

    def _delete(self):
        self._list = []
        self._range = []

    def set_list(self, ilist):
        self._delete()
        self._type = 'list'
        self._list = ilist

    def set_range(self, irange):
        self._delete()
        self._type = 'range'
        self._range = irange

    def range(self):
        if self._type == 'range':
            return self._range
        else:
            api.error_stop("unable to get range from list connectivity")

    def list(self):
        if self._type == 'range':
            return list(arange(self._range))
        elif self._type == 'list':
            return self._list
        else:
            api.error_stop("unknown type")


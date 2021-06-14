import os

_fileformat_map = {}

def fileformat_reader(name, extension):
    """decorator to register fileformat properties for given name in api._fileformat_map
    """
    def decorator(thisclass):
        properties = { 'reader': thisclass, 'ext': extension }
        if name in _fileformat_map.keys():
            _fileformat_map[name].update(properties)
        else: 
            _fileformat_map[name] = properties
        return thisclass
    return decorator

def fileformat_writer(name, extension):
    """decorator to register fileformat properties for given name  in api._fileformat_map
    """
    def decorator(thisclass):
        properties = { 'writer': thisclass, 'ext': extension }
        if name in _fileformat_map.keys():
            _fileformat_map[name].update(properties)
        else: 
            _fileformat_map[name] = properties
        return thisclass
    return decorator

class api_output():
    """class to handle library outputs
    """
    _available = ['internal', 'error', 'warning', 'std', 'debug']
    _default = ['internal', 'error', 'warning', 'std' ]

    def __init__(self, list=None):
        self._api_output = []
        if list==None:
            self.set_default()
        else:
            self.set_modes(list)

    def set_modes(self, modes):
        modelist = modes if type(modes) is list else [modes]
        check = all([mode in self._available for mode in modelist])
        if check:
            self._api_output = modelist
        else:
            self.print('internal','some output modes are unknown')
        return check

    def get_modes(self):
        return self._api_output

    def set_default(self):
        self.set_modes(self._default)

    def print(self, mode, *args):
        if mode in self._api_output:
            print(mode+': ',*args)

io = api_output()
#print(io.get_modes())

class _files():

    def __init__(self, filename):
      self.filename = filename
      self._exists  = os.path.isfile(self.filename)

    def __str__(self):
        s = '  filename: '+self.filename
        return s

    def printinfo(self):
        print(self)
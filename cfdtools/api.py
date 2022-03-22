from pathlib import Path 
import numpy as np

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

def _printreadable(string, value):
    if isinstance(value, (int, float, str, np.int32, np.int64)):
        print(string+':',value)
    elif isinstance(value, np.ndarray):
        if value.size <= 10:
            print(string+': ndarray',value.shape, value)
        else:
            print(string+': ndarray',value.shape)
    else:
        print(string+': '+str(type(value)))

class api_output():
    """class to handle library outputs
    """
    _available = ['internal', 'error', 'warning', 'std', 'debug']
    _default = ['internal', 'error', 'warning', 'std' ]

    def __init__(self, list=None):
        self._api_output = []
        if list is None:
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

    def __init__(self, filename: str):
        self._path = Path(filename)

    @property
    def filename(self):
        return str(self._path)

    @filename.setter
    def set_filename(self, filename):
        self._path = Path(filename)

    def exists(self):
        return self._path.exists()

    def __str__(self):
        s = '  filename: '+self.filename
        return s

    def remove_dir(self):
        self._path = Path(self._path.name)

    def change_suffix(self, ext: str):
        self._path = self._path.with_suffix(ext)

    def find_safe_newfile(self):
        """Returns safe destination, adding (n) if needed"""
        safepath = self._path
        # components strings
        dir = safepath.parent
        stem = safepath.stem
        suff = safepath.suffix
        i = 0
        while safepath.exists():
            i += 1
            #safepath = safepath.with_stem(stem+f'({i})') # only python >= 3.9
            #print(dir,stem,f'({i})',suff)
            safepath = Path(dir/(stem+f'({i})'+suff))
        self._path = safepath
        return i>0
            
    def printinfo(self):
        print(self)
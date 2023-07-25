from functools import wraps
from pathlib import Path
import numpy as np
import time

_fileformat_map = {}


def fileformat_reader(name, extension):
    """decorator to register fileformat properties for given name in api._fileformat_map

    a reader is a class which is initialized with a filename
    and has the following functions
    - read_data()
    - export_mesh() which returns a meshbase._mesh.mesh class
    """

    def decorator(thisclass):
        properties = {'reader': thisclass, 'ext': extension}
        if name in _fileformat_map:
            _fileformat_map[name].update(properties)
        else:
            _fileformat_map[name] = properties
        return thisclass

    return decorator


def fileformat_writer(name, extension):
    """decorator to register fileformat properties for given name  in api._fileformat_map"""

    def decorator(thisclass):
        properties = {'writer': thisclass, 'ext': extension}
        if name in _fileformat_map:
            _fileformat_map[name].update(properties)
        else:
            _fileformat_map[name] = properties
        return thisclass

    return decorator


def _printreadable(string, value):
    if isinstance(value, (int, float, str, np.int32, np.int64)):
        print(string + ':', value)
    elif isinstance(value, np.ndarray):
        if value.size <= 10:
            print(string + ': ndarray', value.shape, value)
        else:
            print(string + ': ndarray', value.shape)
    else:
        print(string + ': ' + str(type(value)))


class api_output:
    """class to handle library outputs"""

    _prefix = {
        'internal': 'int:',
        'error': 'ERROR:',
        'warning': 'WARNING:',
        'std': '',
        'debug': 'debug:',
    }
    _available = list(_prefix.keys())
    _default = ['internal', 'error', 'warning', 'std']

    def __init__(self, iolist=None):
        self._api_output = []
        if iolist is None:
            self.set_default()
        else:
            self.set_modes(iolist)

    def set_modes(self, modes):
        modelist = modes if type(modes) is list else [modes]
        check = all([mode in self._available for mode in modelist])
        if check:
            self._api_output = modelist
        else:
            self.print('internal', 'some output modes are unknown')
        return check

    def get_modes(self):
        return self._api_output

    def set_default(self):
        self.set_modes(self._default)

    def print(self, mode, *args, **kwargs):
        if mode in self._api_output:
            print(self._prefix[mode], *args, **kwargs)

    def printstd(self, *args, **kwargs):
        self.print('std', *args, **kwargs)

    def printdebug(self, *args, **kwargs):
        self.print('debug', *args, **kwargs)

    def warning(self, *args, **kwargs):
        self.print('warning', *args, **kwargs)


io = api_output()
# print(io.get_modes())


def error_stop(msg):
    # io.print('error', msg)
    raise RuntimeError(msg)


class _files:
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
        s = '  filename: ' + self.filename
        return s

    def change_dir(self, newdir):
        self._path = Path(newdir) / Path(self._path.name)

    def remove_dir(self):
        self._path = Path(self._path.name)

    def change_suffix(self, ext: str):
        self._path = self._path.with_suffix(ext)

    def find_safe_newfile(self):
        """Returns safe destination, adding (n) if needed"""
        safepath = self._path
        # components strings
        folder = safepath.parent
        stem = safepath.stem
        suff = safepath.suffix
        i = 0
        while safepath.exists():
            i += 1
            # safepath = safepath.with_stem(stem+f'({i})') # only python >= 3.9
            # print(dir,stem,f'({i})',suff)
            safepath = Path(folder / (stem + f'({i})' + suff))
        self._path = safepath
        return i > 0

    def printinfo(self):
        print(self)


class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""


class Timer:  # from https://realpython.com/python-timer/
    default_tab = 60
    default_msg = ""

    def __init__(self, task="", msg=default_msg, nelem=None, tab=default_tab):
        self.reset()
        self._nelem = nelem
        self._tab = tab
        self._task = task
        self._msg = msg

    def reset(self):
        self._start_time = None
        self._task = ""
        self._col = 0
        self._elapsed = 0.

    @property
    def elapsed(self):
        return self._elapsed
        
    def start(self):
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")
        self._start_time = time.perf_counter()

    def pause(self):
        """Stop the timer, add elapsed and do NOT report"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")
        self._elapsed += time.perf_counter() - self._start_time
        self._start_time = None

    def stop(self, nelem=None):
        """Stop the timer, and report the elapsed time"""
        if nelem is not None:
            self._nelem = nelem
        self.pause()
        normalized_time_ms = (
            0.0 if self._nelem is None else 1e6 * self._elapsed / self._nelem
        )
        if self._nelem is None:
            io.printstd(f"{' '*(self._tab-self._col)}wtime: {self._elapsed:0.4f}s")
        else:
            io.printstd(f"{' '*(self._tab-self._col)}wtime: {self._elapsed:0.4f}s | {normalized_time_ms:0.4f}µs/elem",)
        # reset
        self.reset()

    def __enter__(self):
        io.print('std', self._task, end='')
        self._col = len(self._task) + 1
        self.start()

    def __exit__(self, *exitoptions):
        self.stop()


def memoize(f):
    cache = {}

    @wraps(f)
    def wrapper(*args):
        if not args in cache:
            cache[args] = f(*args)
        # Warning: You may wish to do a deepcopy here if returning objects
        return cache[args]

    return wrapper

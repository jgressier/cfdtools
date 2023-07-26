from cfdtools import __version__
from cfdtools.api import io, _files, error_stop

try:
    import h5py
    import_h5py = True
except ImportError:
    import_h5py = False


from h5py import Group # to be available in _hdf5

_available_types = ('external', 'dataset', 'datalist', 'probes', 'cfdmesh')


def h5_str(obj):
    str = ""
    for c in map(chr, obj[:]):
        str += c
    return str


class h5File(_files):
    def __init__(self, filename: str):
        super().__init__(filename)

    def open(self, mode='r', datatype=None):
        try:
            self._h5file = h5py.File(self._path, mode=mode)
        except:
            io.print('error', "h5py could not be imported")
            raise
        if mode == 'w':
            assert datatype in _available_types
            self._h5file.attrs.update({'cfdtools_version': __version__, 'cfd_datatype': datatype})
        elif mode == 'r':
            if self.exists():
                # only for cgns file (I guess)
                self._h5ver = self._h5file.get(' hdf5version', None)
                self._datatype = self._h5file.attrs.get('cfd_datatype', None)
            else:
                error_stop(f"unable to find {self.filename} (mode r)")
        else:
            error_stop("unknown mode for opening h5 file")

    def close(self):
        return self._h5file.close()
    
    @property
    def datatype(self):
        return self._datatype

    def __getitem__(self, item):
        return self._h5file[item]

    def printinfo(self):
        super().printinfo()
        io.printstd(h5_str(self._h5ver))


# if __name__ == "__main__":
#     box = Cube(2, 1, 1)
#     mesh = box.export_mesh()
#     mesh.printinfo()

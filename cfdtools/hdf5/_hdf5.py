from cfdtools.api import io, _files

try:
    import h5py
    import_h5py = True
except ImportError:
    import_h5py = False


def h5_str(obj):
    str = ""
    for c in map(chr, obj[:]):
        str += c
    return str


class h5file(_files):
    def __init__(self, filename: str):
        super().__init__(filename)

    def open(self, mode='r'):
        self._h5file = h5py.File(self._path, mode=mode)
        if self.exists():
            self._h5ver = self._h5file[' hdf5version']

    def __getitem__(self, item):
        return self._h5file[item]

    def printinfo(self):
        super().printinfo()
        io.printstd(h5_str(self._h5ver))


# if __name__ == "__main__":
#     box = Cube(2, 1, 1)
#     mesh = box.export_mesh()
#     mesh.printinfo()

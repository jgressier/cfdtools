from cfdtools.api import io, _files
try:
    import h5py
    import_h5py = True
except:
    import_h5py = False

def h5_str(obj):
    str = ""
    for c in map(chr, obj[:]):
        str += c
    return str

class h5file(_files):
    def __init__(self, filename: str):
        super().__init__(filename)
        self._h5file = h5py.File(self._path)

    def open(self):
        if self.exists():
            self._h5ver = self._h5file[' hdf5version']

    def __getitem__(self, item):
        return 
    
    def printinfo(self):
        super().printinfo()
        io.print('std', h5_str(self._h5ver))

# if __name__ == "__main__":
#     box = Cube(2, 1, 1)
#     mesh = box.export_mesh()
#     mesh.printinfo()
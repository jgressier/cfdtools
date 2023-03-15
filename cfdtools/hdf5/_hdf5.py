from cfdtools.api import io, _files
try:
    import h5py
    import_h5py = True
except:
    import_h5py = False

class h5file(_files):
    def read(self):
        self._h5file = h5py.File(self._path)
        self._h5ver = self._h5file['hdf5version']

    def printinfo(self):
        super().printinfo()
        io.print('std', f"HDF5 version {self._h5ver}")

# if __name__ == "__main__":
#     box = Cube(2, 1, 1)
#     mesh = box.export_mesh()
#     mesh.printinfo()
import logging

import numpy as np

from cfdtools import __version__
from cfdtools.api import _files, error_stop

try:
    import h5py
    from h5py import Group  # to be available in _hdf5

    import_h5py = True
except ImportError:  # pragma: no cover
    import_h5py = False

_available_types = ('external', 'dataset', 'datalist', 'probes', 'cfdmesh')

log = logging.getLogger(__name__)


def h5_str(obj):
    return "".join(map(chr, obj))


class h5File(_files):
    def __init__(self, filename: str):
        super().__init__(filename)
        self._openedmode = None
        self._datatype = None
        self._cfdtools_version = None

    def open(self, mode='r', datatype=None, version=None):
        """open hdf5 file and parse some version info

        Args:
            mode (str, default: 'r'): could be 'r', 'r+', 'w', 'w-', 'x', 'a' (see h5py.File)
            datatype (_type_, optional): _description_. Defaults to None.
            version (_type_, optional): _description_. Defaults to None.
        """
        try:
            self._h5file = h5py.File(self._path, mode=mode)
        except NameError:
            raise NameError(f"{self._path} can not be open in {mode} mode.")
        if mode in ('w', 'w-', 'x'):
            assert datatype in _available_types
            self._h5file.attrs.update({'cfdtools_version': __version__, 'cfd_datatype': datatype})
            if datatype is not None and (version is None):
                error_stop('an existing cfdtools datatype must provide a version number for writing')
            if datatype is not None:
                self._h5file.attrs.update({'data_version': version})
        elif mode in ('r', 'r+'):
            # only for cgns file (I guess)
            self._h5ver = self._h5file.get(' hdf5version', None)
            self._datatype = self._h5file.attrs.get('cfd_datatype', None)
            self._cfdtools_version = self._h5file.attrs.get('cfdtools_version', None)
            # if datatype, read or set to 0 (backward compatibility) else ignore
            self._version = self._h5file.attrs.get('data_version', 0 if self._datatype else None)
        else:
            error_stop("unknown mode for opening h5 file")
        self._openedmode = mode

    def close(self):
        return self._h5file.close()

    @property
    def datatype(self):
        return self._datatype

    def __getitem__(self, item):
        return self._h5file[item]

    def write_unsmesh(self, cellnode: dict, nodes: np.ndarray, path="/mesh", meshtype="unsmesh", **options):
        hgroup = self._h5file.create_group(path)
        hgroup.attrs['meshtype'] = meshtype
        hgroup.create_dataset("nodes", data=nodes, **options)
        hcells = hgroup.create_group("cells")
        for ctype, cellco in cellnode.items():
            hcells.create_dataset(ctype, data=cellco, **options)

    def printinfo(self):
        super().printinfo()
        if self._openedmode:
            if self._h5ver:
                log.info(f"         hdf5 version: {h5_str(self._h5ver)}")
            if self._cfdtools_version:
                log.info(f"     cfdtools version: {self._cfdtools_version}")
            if self._datatype:
                log.info(f"        cfdtools type: {self._datatype}")
            if self._version:
                log.info(f"cfdtools data version: {self._version}")


# if __name__ == "__main__":
#     box = Cube(2, 1, 1)
#     mesh = box.export_mesh()
#     mesh.printinfo()

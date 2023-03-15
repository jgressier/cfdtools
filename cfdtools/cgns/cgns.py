# cgns.py
from cfdtools.api import io
from cfdtools.hdf5 import h5file
import numpy as np

class cgnsfile(h5file):
    def __init__(self, filename: str):
        super().__init__(filename)

    def read(self):
        super().read()
        self._cgnsver = self._h5file['CGNSLibraryVersion']

class cgnsmesh():
    def __init__(self) -> None:
        pass

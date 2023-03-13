# Import modules
import numpy as np
from cfdtools.vtk._vtk import map_ele
import cfdtools.meshbase._mesh as _mesh
import cfdtools.api as api
try:
    import pyvista as pv
    from pyvista import CellType
    importpyvista=True
except:
    importpyvista=False

class writer():
    ''' Implementation of the writer to write binary Vtk files using pyvista'''

    def __init__(self, mesh):
        self._mesh = mesh

    def set_mesh(self):
        self._coords = self._mesh.nodescoord(ndarray=True)
        self._celldict = { 
            map_ele[etype]: elem2node['elem2node']
            for etype, elem2node in self._mesh._cell2node.items() }
        self._grid = pv.UnstructuredGrid(self._celldict, self._coords)

    def write_mesh(self, filename):
        self._grid.save(filename)
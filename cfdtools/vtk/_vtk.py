try:
    import pyvista as pv
    from pyvista import CellType
    importpyvista = True
except ImportError:
    importpyvista = False
import cfdtools.meshbase._elements as ele
import cfdtools.api as api
import cfdtools.utils.maths as maths
from pathlib import Path
import numpy as np
import scipy.spatial as spspa


map_ele = {'bar2': CellType.LINE, 'quad4': CellType.QUAD, 'hexa8': CellType.HEXAHEDRON}


@api.fileformat_writer("VTK", '.vtu')
class vtkMesh:
    '''Implementation of the writer to write binary Vtk files using pyvista'''

    def __init__(self, mesh=None):
        self._mesh = mesh
        if mesh:
            self.set_mesh(mesh)

    def set_mesh(self, mesh):
        self._coords = self._mesh.nodescoord(ndarray=True)
        self._celldict = {
            map_ele[etype]: elem2node['elem2node']
            for etype, elem2node in self._mesh._cell2node.items()
        }
        self._grid = pv.UnstructuredGrid(self._celldict, self._coords)

    def write_data(self, filename):
        self._grid.save(filename)

    def read(self, filename):
        self.__init__()
        self._grid = pv.read(filename)

    def pyvista_grid(self):
        return self._grid

    def plot(self, background='white', show_edges=True, *args, **kwargs):
        self.pyvista_grid().plot(background='white', show_edges=True, *args, **kwargs)

class vtkList():

    def __init__(self, filelist, verbose=False) -> None:
        self._list = filelist
        self._verbose = verbose

    @property
    def nfile(self):
        return len(self._list)
    
    def exist(self):
        return all(Path(file).exists() for file in self._list)

    def check_order(self, pos='cellcenter', tol=1e-10):
        mappos = { 'cellcenter' : lambda m: m.cell_centers().points,
                   'node' : lambda m: m.points }        
        count = 0
        ref = mappos[pos](pv.read(self._list[0]))
        for name in self._list:
            mesh = pv.read(name)
            d = np.max(maths.distance(ref, mappos[pos](mesh)))
            if d > tol:
                count += 1
                if self._verbose:
                    api.io.printstd(f"  . {name}: {d}")
        if self._verbose:
            api.io.printstd(f"  {count}/{self.nfile} grids are not {pos}-coincident")
        return count > 0

    def read(self, filterdata=None, reorder=False):
        assert self.nfile > 0
        self._mesh = pv.read(self._list[0])

try:
    import pyvista as pv
    from pyvista import CellType
    importpyvista=True
except:
    importpyvista=False
import cfdtools.meshbase._elements as ele
import cfdtools.api as api


map_ele = {
    'hexa8' : CellType.HEXAHEDRON
}

@api.fileformat_writer("VTK", '.vtu')
class vtkMesh():
    ''' Implementation of the writer to write binary Vtk files using pyvista'''

    def __init__(self, mesh):
        self._mesh = mesh
        self.set_mesh()

    def set_mesh(self):
        self._coords = self._mesh.nodescoord(ndarray=True)
        self._celldict = { 
            map_ele[etype]: elem2node['elem2node']
            for etype, elem2node in self._mesh._cell2node.items() }
        self._grid = pv.UnstructuredGrid(self._celldict, self._coords)

    def write_mesh(self, filename):
        self._grid.save(filename)

    def pyvista_grid(self):
        return self._grid
    
    def plot(self, background='white', show_edges=True, *args, **kwargs):
        self.pyvista_grid().plot(background='white', show_edges=True, *args, **kwargs)
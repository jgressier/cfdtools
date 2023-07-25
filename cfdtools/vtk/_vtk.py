try:
    import pyvista as pv
    from pyvista import CellType
    importpyvista = True
except ImportError:
    importpyvista = False
import cfdtools.api as api
import cfdtools.utils.maths as maths
import cfdtools.meshbase._elements as ele
import cfdtools.meshbase._mesh as cfdmesh
from cfdtools.data import DataSetList
import cfdtools.hdf5 as hdf5
from pathlib import Path
import numpy as np
import scipy.spatial as spspa


vtktype_ele = {'bar2': CellType.LINE, 'quad4': CellType.QUAD, 'hexa8': CellType.HEXAHEDRON}
ele_vtktype = { i: etype for etype, i in vtktype_ele.items() }

@api.fileformat_writer("VTK", '.vtu')
class vtkMesh:
    '''Implementation of the writer to write binary Vtk files using pyvista'''

    def __init__(self, mesh=None, pvmesh=None):
        self._mesh = mesh
        if mesh and pvmesh:
            api.error_stop("vtkMesh cannot be initialized by both structures")
        if mesh:
            self.set_mesh(mesh)
        if pvmesh:
            self.set_pvmesh(pvmesh)

    def set_mesh(self, mesh: cfdmesh.Mesh):
        self._mesh = mesh
        coords = self._mesh.nodescoord(ndarray=True)
        try:
            self._celldict = {
                vtktype_ele[etype]: elem2node['elem2node']
                for etype, elem2node in self._mesh._cell2node.items()
            }
            self._grid = pv.UnstructuredGrid(self._celldict, coords)
        except:
            api.error_stop("pyvista (with CellType) could not be imported")

    def set_pvmesh(self, pvmesh: cfdmesh.Mesh):
        self._grid = pvmesh

    def write_data(self, filename):
        self._grid.save(filename)

    def read(self, filename):
        self.__init__()
        self._grid = pv.read(filename)

    @property
    def pyvista_grid(self):
        return self._grid

    def brief(self):
        if self._grid:
            m = self._grid
            api.io.printstd(f"pyvista object: {m}")
            api.io.printstd("> mesh")
            api.io.printstd(f"  bounds : {m.bounds}")
            api.io.printstd(f"  ncells : {m.n_cells}")
            api.io.printstd(f"  npoints: {m.n_points}")
            api.io.printstd("> data")
            api.io.printstd(f"  cell  data names: {m.cell_data.keys()}")
            api.io.printstd(f"  point data names: {m.point_data.keys()}")
            api.io.printstd(f"  field data names: {m.field_data.keys()}")
            #api.io.printstd("> properties")
        else:
            api.io.warning("  no pyvista mesh available")

    def plot(self, background='white', show_edges=True, *args, **kwargs):
        self.pyvista_grid().plot(background=background, show_edges=show_edges, *args, **kwargs)

    def dumhdfgroup(self, hgroup: hdf5.Group, **options):
        hgroup.attrs['meshtype'] = 'unsvtk'
        hgroup.create_dataset("nodes", data=self._grid.points, **options)
        hcells = hgroup.create_group("cells")
        for itype, cellco in self._grid.cells_dict.items():
            hcells.create_dataset(ele_vtktype[itype], data=cellco, **options)



class vtkList():

    def __init__(self, filelist, verbose=False) -> None:
        self._list = filelist
        self._verbose = verbose

    @property
    def nfile(self):
        return len(self._list)
    
    def allexist(self):
        return all(Path(file).exists() for file in self._list)

    def check_order(self, pos='cellcenter', tol=1e-10):
        mappos = { 'cellcenter' : lambda m: m.cell_centers().points,
                   'node' : lambda m: m.points }        
        count = 0
        assert self.nfile > 0
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
        return count == 0

    def read(self, filterdata=None, reorder=False, tol=1e-10):
        count = 0
        assert self.nfile > 0
        Tread = api.Timer()
        Tcomp = api.Timer()
        Tsort = api.Timer()
        Tread.start()
        self._mesh = pv.read(self._list[0])
        self._ncell = self._mesh.n_cells
        ctrRef = self._mesh.cell_centers().points
        Tread.pause()
        Tsort.start()
        tree = spspa.KDTree(ctrRef)
        Tsort.pause()
        #
        self._data = DataSetList(self.nfile, Xrep='cellaverage', Trep='instant')
        # may add alive-progress or other
        for name in self._list:
            Tread.start()
            vtk = pv.read(name)
            Tread.pause()
            namelist = filterdata if filterdata else vtk.cell_data.keys()
            Tcomp.start()
            ctr = vtk.cell_centers().points
            d = np.max(maths.distance(ctrRef, ctr))
            Tcomp.pause()
            if d > tol:
                count += 1
                #if self._verbose:
                #    api.io.printstd(f"  . {name}: {d}")
                Tsort.start()
                dfinal, index = tree.query(ctr, p=2)
                # reverse indexing to sort new arrays
                rindex = index.copy()
                rindex[index] = np.arange(index.size) 
                assert np.max(dfinal) <= tol
                # automatically deals with differnt shapes
                datalist = {
                    name: vtk.cell_data[name][rindex] for name in namelist
                }
                Tsort.pause()
                self._data.add_datalist(datalist)
        if self._verbose:
            api.io.printstd(f"  {count}/{self.nfile} grids were not coincident")
            api.io.printstd(f"       file reading: {Tread.elapsed:.2f}s")
            api.io.printstd(f"    grid comparison: {Tcomp.elapsed:.2f}s")
            api.io.printstd(f"    data reordering: {Tsort.elapsed:.2f}s")

    def dumphdf(self, filename, **options):
        file = hdf5.h5File(filename)
        file.find_safe_newfile()
        file.open(mode="w", datatype='datalist')
        hmesh = file._h5file.create_group("mesh")
        vtkmesh = vtkMesh(pvmesh=self._mesh)
        vtkmesh.dumhdfgroup(hmesh, **options)
        #
        hdata = file._h5file.create_group("datalist")
        self._data._dumphdfgroup(hdata, **options)
        file.close()
        

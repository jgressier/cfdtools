import logging
import numpy as np
import scipy.spatial as spatial
from typing import Callable

try:
    import pyvista as pv
    from pyvista import CellType

    importpyvista = True
except ImportError:
    importpyvista = False

import cfdtools.api as api
import cfdtools.meshbase._connectivity as _conn
import cfdtools.meshbase._mesh as cfdmesh
import cfdtools.utils.maths as maths
import cfdtools.hdf5 as hdf5

log = logging.getLogger(__name__)

if importpyvista:
    vtktype_ele = {
        'bar2': CellType.LINE,
        'tri3': CellType.TRIANGLE,
        'quad4': CellType.QUAD,
        'tetra4': CellType.TETRA,
        'pyra5': CellType.PYRAMID,
        'prism6': CellType.WEDGE,
        'hexa8': CellType.HEXAHEDRON,
    }
    ele_vtktype = {i: etype for etype, i in vtktype_ele.items()}
    PYVISTA2XDMF = {
        CellType.TRIANGLE: "Triangle",
        CellType.QUAD: "Quadrilateral",
        CellType.HEXAHEDRON: "Hexahedron",
        CellType.TETRA: "Tetrahedron",
        CellType.PYRAMID: "Pyramid",
        CellType.WEDGE: "Prism",
    }
else:
    log.warning("pyvista missing or failing at import: some features will be missing.")


@api.fileformat_writer("VTK", '.vtu')
class vtkMesh:
    """Implementation of the writer to write binary Vtk files using pyvista"""

    _version = 1

    def __init__(self, mesh=None, pvmesh=None):
        """Initialisation of vtkMesh, encapsulation of pyvista vtk mesh

        Args:
            mesh (meshbase, optional): cfdtools meshbase object. Defaults to None.
            pvmesh (pyvista, optional): direct pyvista mesh. Defaults to None.
        """
        if mesh and pvmesh:
            api.error_stop("vtkMesh cannot be initialized by both structures")
        self._reset()
        if mesh:
            self.set_mesh(mesh)
        if pvmesh:
            self.set_pvmesh(pvmesh)

    def _reset(self):
        self._mesh = None
        self._grid = None
        self._volume = None

    def set_mesh(self, mesh: cfdmesh.Mesh):
        self._reset()
        self._mesh = mesh
        coords = self._mesh.nodescoord(ndarray=True)
        self._celldict = {
            vtktype_ele[etype]: elem2node['elem2node'] for etype, elem2node in self._mesh._cell2node.items()
        }
        self._grid = pv.UnstructuredGrid(self._celldict, coords)

    def set_pvmesh(self, pvmesh: cfdmesh.Mesh):
        self._reset()
        self._grid = pvmesh

    def write_data(self, filename):
        self._grid.save(filename)

    def read(self, filename):
        self.__init__()
        self.set_pvmesh(pv.read(filename))
        return self

    def export_mesh(self):
        """generates a cfdtools meshbase from vtk connectivity, boundary conditions are missing

        Returns:
            meshbase.Mesh: cell to node connectivity and node positions
        """
        cellnode = _conn.elem_connectivity()
        for itype, cellco in self._grid.cells_dict.items():
            cellnode.add_elems(ele_vtktype[itype], cellco)
        mesh = cfdmesh.Mesh(ncell=self.ncell, nnode=self.nnode)
        mesh.set_cell2node(cellnode)
        mesh.set_nodescoord_nd(self._grid.points)
        return mesh

    @property
    def pyvista_grid(self):
        return self._grid

    @property
    def ncell(self):
        return self._grid.n_cells

    @property
    def nnode(self):
        return self._grid.n_points

    def brief(self):
        if self._grid:
            m = self._grid
            log.info(f"pyvista object: {m.DATA_TYPE_NAME}")
            log.info("> mesh")
            log.info(f"    ncells : {m.n_cells}")
            log.info(f"    npoints: {m.n_points}")
            for i, c in enumerate(['X', 'Y', 'Z']):
                cmin, cmax = m.bounds[i * 2 : i * 2 + 2]
                log.info(f"  {c} bounds:  {cmin} - {cmax}")
            log.info("> data")
            log.info(f"  cell  data names: {m.cell_data.keys()}")
            log.info(f"  point data names: {m.point_data.keys()}")
            log.info(f"  field data names: {m.field_data.keys()}")
            # log.info("> properties")
        else:
            log.warning("  no pyvista mesh available")

    def plot(self, background='white', show_edges=True, *args, **kwargs):
        self.pyvista_grid.plot(background=background, show_edges=show_edges, *args, **kwargs)

    def importhdfgroup(self, hgroup: hdf5.Group, verbose=False):
        assert hgroup.attrs['meshtype'] in ('unsmesh',)
        points = np.array(hgroup["nodes"])
        celldict = {vtktype_ele[etype]: np.array(elem2node) for etype, elem2node in hgroup["cells"].items()}
        self.set_pvmesh(pv.UnstructuredGrid(celldict, points))

    def xdmf_content(self, filename):
        """Create the XDMF content associated to the mesh.

        :param str filename: Name of the output HDF5 file.
        :return: The XDMF content.
        :rtype: list(str)
        """
        dim2type = {1: "X", 2: "XY", 3: "XYZ"}
        geomtype = dim2type[self._grid.points.shape[1]]
        lines = [f'<Geometry GeometryType="{geomtype}">']
        lines += [f'<DataItem Dimensions="{np.prod(self._grid.points.shape)}" Format="HDF">']
        lines += [f"{filename}:/mesh/nodes"]
        lines += ["</DataItem>"]
        lines += ["</Geometry>"]

        for itype, cellco in self._grid.cells_dict.items():
            lines += [
                f'<Topology NumberOfElements="{cellco.shape[0]:d}" TopologyType="{PYVISTA2XDMF[itype]}">'
            ]
            lines += [f'<DataItem Dimensions="{np.prod(cellco.shape)}" Format="HDF">']
            lines += [f"{filename}:/mesh/cells/{ele_vtktype[itype]}"]
            lines += ["</DataItem>"]
            lines += ["</Topology>"]

        return lines

    def volumes(self):
        if not self._volume:
            self._volume = self._grid.compute_cell_sizes().cell_data["Volume"]
        return self._volume

    def dumphdf(self, filename, overwrite=False, **options):
        """dump pyvista mesh and (future) data to cfdtools hdf5 file
        (in the future) must be converted to dataset and dump done by dataset

        Args:
            filename (string or pathlib): file name
            overwrite (bool): overwrite file or find new safe name
            **options: passed to dumpdfgroup and h5py.createdataset

        Returns:
            string: file name of the actual (safe) saved file
        """
        file = hdf5.h5File(filename)
        if not overwrite:
            file.find_safe_newfile()
        file.open(mode="w", datatype='dataset', version=self._version)
        # map element types from vtk type to cfdtools type
        celldict = {ele_vtktype[itype]: cellco for itype, cellco in self._grid.cells_dict.items()}
        file.write_unsmesh(celldict, self._grid.points, **options)
        file.close()
        return file.filename

    def vtk_extract_cells(self, selected_cell_ids: list):
        """re-implementation of extract_cell because vtk/pyvista does not keep the order
        celltypes is the list a cell types codes (vtk based)
        offsets is the list of offsets for each cell
        cells is a coompound defintion of cells with concatenated (size, point ids)
        points is the list of points, shape (:,3)
        """
        with api.Timer("vtu extract slice"):
            mesh = self._grid
            orig_cells = mesh.cells
            celltypes = mesh.celltypes
            offsets = mesh.offset
            points = mesh.points
            #
            new_celltypes = [ celltypes[cid] for cid in selected_cell_ids ]
            extract_vtex = np.concatenate([orig_cells[offsets[cid]+cid+1:offsets[cid+1]+cid+1] for cid in selected_cell_ids])
            nvtex = np.array([orig_cells[offsets[cid]+cid] for cid in selected_cell_ids], dtype=np.uint8)
            used_point_ids = list(set(extract_vtex))
            # Create mapping for point ID reindexing
            point_id_map = {old: new for new, old in enumerate(used_point_ids)}
            new_vtex = np.array(list(map(lambda x: point_id_map[x], extract_vtex)), dtype=np.int64)
            # Rebuild VTK-style cell array
            vtk_cells = np.empty(len(new_vtex)+len(new_celltypes), dtype=np.int64)
            i = int(0) ; j = 0
            for size in nvtex:
                vtk_cells[i] = size
                vtk_cells[i+1:i+1+size] = new_vtex[j:j+size]
                i += np.int64(size) + 1 # size + 1, size is uint8
                j += np.int64(size)
            new_celltypes = np.array(new_celltypes, dtype=np.uint8) # cell type is small
            new_points = points[used_point_ids, :]
            # create extracted unstructured grid
            new_grid = pv.UnstructuredGrid(vtk_cells, new_celltypes, new_points)
            new_mesh = vtkMesh()
            new_mesh.set_pvmesh(new_grid)
        return new_mesh

    def vtk_zconvolution(self, splitcoords: Callable = None, tol = 1.e-6, nmode=0, snapshot=False, rms=False, phase=False, select: dict = None):
        """average the mesh data along a splitcoords and return a new vtkMesh object

        Options can be replaced by a dictionary of options which keys are celldata names
        e.g. select = {'Q': 'avg', 'U': ['no_avg', 'rms', {'nmode': 2}, 'phase' ]} will compute the average of Q and the average, RMS, 2 fourier modes of U
        (phase included), all other data will be ignored.

        Args:
            splitcoords (Callable, optional): function to split the coordinates. Defaults to None.
            tol (float, optional): tolerance for distance between planes. Defaults to 1.e-6.
            nmode (int, optional): number of modes to keep. Defaults to 0.
            snapshot (bool, optional): keep the snapshot data. Defaults to False.
            rms (bool, optional): compute the RMS of the data. Defaults to False.
            select (dict, optional): dictionary of variables/options selection to apply. Defaults to None.
        Returns:
            vtkMesh: averaged mesh
        """
        # parse dict if necessary and apply general options
        current_options = {'avg': True, 'rms': rms, 'phase': phase, 'nmode' : nmode, 'snapshot': snapshot }
        if select is None:
            directives = { keyname: current_options.copy() for keyname in self._grid.cell_data.keys() }
        else:
            # dict value can be only one (string) key or a list of a dict
            directives = {}
            for key, value in select.items():
                if key not in self._grid.cell_data.keys():
                    log.error(f"  unknown name {key} in VTK cell data")
                    raise ValueError(f"  unknown key {key} used in select option of vtk_zconvolution()")
                directives[key] = {**current_options}
                if isinstance(value, str): value = [value]
                if isinstance(value, list):
                    for v in value:
                        if isinstance(v, str):
                            if v == 'no_avg':
                                directives[key]['avg'] = False
                            else:
                                directives[key][v] = True
                        elif isinstance(v, dict):
                            directives[key].update(v)
                        else: # pragma: no cover
                            raise ValueError(f"unknown type {type(v)} for {key}")
                else: # pragma: no cover
                    raise ValueError(f"unknown type {type(value)} for {key}")
        #
        log.info("> computes coords and compute chunks for planes")
        with api.Timer():
            # splitcoords
            if splitcoords is None:
                splitcoords = lambda x: (x[:,0:2], x[:, 2]) # default to z splitcoords
            if not callable(splitcoords): # pragma: no cover
                raise ValueError("splitcoords must be a callable function")
            # compute cell centers and positions in slice (default xy) and transverse (z)
            cell_centers = self._grid.cell_centers().points
            xypos, zpos = splitcoords(cell_centers)
            # sort
            isort = np.argsort(zpos)
            ichunk = np.where(np.abs(np.diff(zpos[isort])) > tol)[0] + 1
            # check chunk size
            nchunk = len(ichunk)+1
            if nchunk <= 1: # pragma: no cover
                log.error("  no chunks found, please check the direction function")
                raise ValueError("  no chunks found")
            chunks = np.split(isort, ichunk)
            lens = list(map(len, chunks))
            log.info(f"  {len(chunks)} chunks found with min:max sizes {min(lens)}:{max(lens)}")
            if min(lens) != max(lens): # pragma: no cover
                log.error("  chunks have different sizes, please check the direction function")
                raise ValueError("  chunks have different sizes")
            chunksize = lens[0]
            # check distance between planes
            zavg = [np.average(zpos[chunk]) for chunk in chunks]
            zdiff = np.diff(zavg)
            zdiffavg, zdiffstd = (func(zdiff) for func in (np.average, np.std))
            log.info(f"  average distance (stddev) between planes: {zdiffavg:.2e} ({zdiffstd:.2e})")
            if zdiffstd > tol:
                log.error(f"  max distance between planes is above tolerance {zdiffstd:.2e}")
                raise ValueError("  max distance between planes is above tolerance")
        # compute all index with kdtree localization
        log.info("> sort and localize cells in planes")
        with api.Timer():
            ikplanes = np.zeros((nchunk, chunksize), dtype=np.int64)
            tree = spatial.KDTree(xypos[chunks[0]])
            for i, chunk in enumerate(chunks):
                if i == 0:
                    ikplanes[0,:] = chunk
                else:
                    dfinal, index = tree.query(xypos[chunk], p=2) # p=2 is the norm
                    dmin, davg, dmax = maths.minavgmax(dfinal)
                    if dmax > tol:
                        log.error(f"  max distance is above tolerance {tol}: {np.count_nonzero(dfinal > tol)} cells incriminated")
                        log.error(f"  min:avg:max = {dmin:.2e}:{davg:.2e}:{dmax:.2e}")
                        raise ValueError("  max distance is above tolerance")
                    ikplanes[i, index] = chunk
        # check
        # shape of xypos[ikplanes] is (nchunk, chunksize, 2)
        stddev = np.std(xypos[ikplanes], axis=0)
        if np.any(stddev > tol):
            log.error("  stddev of points in planes is above tolerance")
            raise ValueError("  stddev of points in planes is above tolerance")
        # create one slice of mesh
        log.info("> create new sliced mesh with FFT data")
        new_mesh = self.vtk_extract_cells(chunks[0])
        # and copy/create data
        for name, options in directives.items():
            data = self._grid.cell_data[name][ikplanes]
            var_nmode = options.get('nmode')
            log.info(f"  {name:<12}: {var_nmode} mode(s) of {data.shape} + {', '.join([k for k, v in options.items() if v])} ")
            if options['snapshot']:
                new_mesh._grid.cell_data[name] = data[0,...]
            if options['avg']:
                new_mesh._grid.cell_data[f"{name}_avg"] = np.average(data, axis=0)
            if var_nmode > 0:
                datafft = np.fft.fft(data, norm='ortho', axis=0)[:var_nmode+1]
                for i in range(1, var_nmode+1):
                    new_mesh._grid.cell_data[f"{name}_k{i}"] = np.abs(datafft[i,...])
                    if options['phase']:
                        new_mesh._grid.cell_data[f"{name}_k{i}_phase"] = np.angle(datafft[i,...])
            if options['rms']:
                new_mesh._grid.cell_data[f"{name}_rms"] = np.std(data, axis=0)
        return new_mesh
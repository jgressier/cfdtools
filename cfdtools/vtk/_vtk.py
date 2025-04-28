import logging

import numpy as np
import scipy.spatial as spatial
from typing import Callable
import sys

out_is_tty = sys.stdout.isatty()
# CHK # # switch for check+debug counters
# CHK # do_iter_check = False
if out_is_tty:
    from tqdm import tqdm
    # CHK # # imports for check+debug counters
    # CHK # if do_iter_check:
    # CHK #     import matplotlib.pyplot as plt
    # CHK #     from time import time
else:
    tqdm = lambda x: x

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
        else:
            log.warning("  no pyvista mesh available")
            return False

        log.info(f"pyvista object:")
        mstr_list = f"{m}".split(chr(10))
        for l in mstr_list[:1]:
            log.info(f"  {l}")
        log.info(f"  Mesh")
        for l in mstr_list[1:]:
            log.info(f"    | {l[2:]}")
        # # To get all possible attributes:
        # for w in dir(m):
        #     print(f"{w}:", (str(m.__getattr__(w)).split(chr(10)))[0])
        # return
        # log.info(f"pyvista object: {m.DATA_TYPE_NAME}") # Attribute is a method
        # log.info(f"pyvista object: {m.DATA_TYPE_NAME()}") # Wrong output
        log.info(f"pyvista object:")
        log.info(f"  {type(m).__name__} ({hex(id(m))}) [{m.__vtkname__}]")
        log.info(f"  Mesh")
        log.info(f"    | N Cells:  {m.n_cells}")
        log.info(f"    | N Points: {m.n_points}")
        for i, c in enumerate('XYZ'):
            cmin, cmax = m.bounds[i * 2 : i * 2 + 2]
            log.info(f"    |   {c} bounds:  {cmin:10.3e},{cmax:10.3e}")
        if hasattr(m, 'n_arrays'):
            log.info(f"    | N Arrays: {m.n_arrays}")
        log.info(f"  Data")
        log.info(f"    | cell  data names: {m.cell_data.keys()}")
        log.info(f"    | point data names: {m.point_data.keys()}")
        log.info(f"    | field data names: {m.field_data.keys()}")
        # log.info("  Properties")

        return True

    def diag(self):
        if self._grid:
            m = self._grid
        else:
            log.warning("  no pyvista mesh available")
            return False

        if not self.brief():
            return False

        ### Compute cell volumes
        ###=====================

        if out_is_tty:
            log.info(f"  ...Computing volumes of {m.n_cells} cells...")
        volume = self.volumes()
        volmin, volmax = min(volume), max(volume)
        voltot, volava = sum(volume), np.mean(volume)
        vol_nb, volavg = len(volume), np.exp(np.mean(np.log(volume)))

        log.info(f"  Volumes")
        log.info(f"    | vol nb / total vol :    {vol_nb:12d} / {voltot:12.6e}")
        log.info(f"    | (arith / geom) avg :    {volava:12.6e} / {volavg:12.6e}")
        log.info(f"    | min vol / max vol  :    {volmin:12.6e} / {volmax:12.6e}")
        if out_is_tty:
            sys.stdout.flush()

        ### Compute volume ratios from vtk cells and points
        ###================================================

        # m.cells: 1-D (flat) array containing the concatenation of, for each cell:
        #           [n==nb_pts_of_cell, 1st_pt_of_cell, (...), nth_pt_of_cell]

        ### Fill lists: cells_to_points, points_to_cells
        #
        # cells_to_points: list indexed with cell number containing, for each cell:
        #               [1st_pt_of_cell, (...), nth_pt_of_cell]
        # points_to_cells: list indexed with point number containing, for each point:
        #               [1st_cell_of_pt, (...), nth_cell_of_pt]

        ### Build cells_to_points: list for each cell of points of cell
        #
        cells_to_points = [[] for _ in range(m.n_cells)]
        # Initialize counter for m.cells
        idx_mc = 0
        # Initialize counter for cells_to_points
        idx_pc = 0
        # Loop on all elements in m.cells
        while idx_mc < len(m.cells):
            # Number of points of next cell
            c_npt = m.cells[idx_mc]
            # Number of indices for next cell (== 1 + nb_points)
            c_nptp1 = c_npt + 1
            # Number of next cells with the same nb_points
            nc_c_npt = np.argmax(m.cells[idx_mc::c_nptp1] != c_npt)
            # If zero:
            if nc_c_npt == 0:
                # All remaining cells
                nc_c_npt = (len(m.cells) - idx_mc) // c_nptp1
            # Compute index increment
            n_idx_mc = nc_c_npt * c_nptp1
            # Compute cell sublist and extent_base
            cells_sublist = m.cells[idx_mc : idx_mc + n_idx_mc]
            cl_subl_table = cells_sublist.reshape((nc_c_npt, c_nptp1))
            # Complement cells_to_points list
            cells_to_points[idx_pc : idx_pc + nc_c_npt] = (cp[1:] for cp in cl_subl_table)
            # Increment counter by nb_points * nb_indices
            idx_mc += n_idx_mc
            idx_pc += nc_c_npt

        ### Build points_to_cells: list for each point of cells of point
        #
        points_to_cells = [[] for _ in range(m.n_points)]
        # Loop on all cells in cells_to_points
        for ic, c in enumerate(cells_to_points):
            # For every point in cell
            for p in c:
                # Add cell for point in points_to_cells
                points_to_cells[p].append(ic)

        ### Build cell_cells: list for each cell of dict for each neighbour cell of common points
        #
        # Initialize cell_cells: list of nbcells empty dicts
        cell_cells = [{} for _ in range(len(cells_to_points))]
        # Loop on each point and list of cells thereof
        for p, pc in enumerate(points_to_cells):
            # Loop on all cells of points
            for ca in pc:
                # Loop on all cells of points again
                for cb in pc:
                    # Ignore if cellb == cella
                    if cb == ca:
                        continue
                    # If cellb not yet neighbour of cella
                    if cb not in cell_cells[ca]:
                        # Initialize key cellb for cella
                        cell_cells[ca][cb] = []
                    # Add common point
                    cell_cells[ca][cb] += [p]

        ### Activate ONLY if VERY FEW cells
        ### for ica, ca in enumerate(cell_cells):
        ###     print(f"{ica}:", '{',
        ###           ', '.join([f"{i}: {str(x):16}" for i, x in ca.items()]),
        ###           '}')

        node_cells = [[] for _ in range(len(m.points))]
        edge_cells = {}
        face_cells = {}
        for ca, ca_cb_pts in enumerate(cell_cells):
            for cb, pts in ca_cb_pts.items():
                if len(pts) == 1:
                    node_cells[pts[0]] += [ca, cb]
                elif len(pts) == 2:
                    tpts = tuple(pts)
                    if tpts not in edge_cells:
                        edge_cells[tpts] = []
                    edge_cells[tpts] += [ca, cb]
                else:  # len(pts) > 2
                    tpts = tuple(pts)
                    if tpts not in face_cells:
                        face_cells[tpts] = []
                    face_cells[tpts] += [ca, cb]
                del cell_cells[cb][ca]

        ### Activate ONLY if VERY FEW cells
        ### print("Edges:")
        ### for e in edge_cells.items(): print(f"{e[0]}: {e[1]}")
        ### print("Faces:")
        ### for f in face_cells.items(): print(f"{f[0]}: {f[1]}")
        ### for c in cells_to_points: print(c)
        ### for p in points_to_cells: print(p)

        ### #print(type(cells_to_points), [type(c) for c in cells_to_points])
        ### #print(type(points_to_cells), [type(c) for c in points_to_cells])
        ### cells_to_points = np.array([np.array(plist) for plist in cells_to_points], dtype=object)
        ### points_to_cells = np.array([np.array(clist) for clist in points_to_cells], dtype=object)
        ### #print(type(cells_to_points), [type(c) for c in cells_to_points])
        ### #print(type(points_to_cells), [type(c) for c in points_to_cells])
        ### print(cells_to_points.shape)
        ### print(points_to_cells.shape)

        ### Compute ratios
        #
        # Backup point_to_cells
        bk_points_to_cells = [[_ for _ in pc] for pc in points_to_cells]
        # Initialize ratios
        rat_node, rat_edge, rat_face = [1.0, 1.0, 1.0]
        # Initialize counters
        count_rat_all = 0
        count_rat_node = 0
        count_rat_edge = 0
        count_rat_face = 0
        # CHK # # initialize check+debug counters
        # CHK # if do_iter_check:
        # CHK #     nb_iter_ca = nb_iter_pi = nb_iter_cb = 0
        # loop on all cells: ca
        # for ca in range(len(cells_to_points)):
        # CHK # # initialize time distribution
        # CHK # if do_iter_check and out_is_tty:
        # CHK #     time0 = time()
        # CHK #     mytime = [0 for i in range(m.n_cells)]
        for ca in tqdm(range(len(cells_to_points))):
            # CHK # # record time distribution
            # CHK # if do_iter_check and out_is_tty:
            # CHK #     mytime[ca] = time() - time0
            # CHK # if do_iter_check:
            # CHK #     nb_iter_ca += 1
            # neighbour cells already processed: initialize
            ca_done = []
            # loop on those points of ca: pi
            for pi in cells_to_points[ca]:
                # remove cell ca for point pi (to not be processed again)
                points_to_cells[pi].remove(ca)
                # CHK # # increment check+debug counters
                # CHK # if do_iter_check:
                # CHK #     nb_iter_pi += 1
                # cells of pi not already processed (as a previous ca, or as cb for this ca)
                pi_cells = [c for c in points_to_cells[pi] if c not in ca_done]
                # loop on those cells of pi: cb
                for cb in pi_cells:
                    # CHK # # increment check+debug counters
                    # CHK # if do_iter_check:
                    # CHK #     nb_iter_cb += 1
                    # neighbour cells already processed: update (for the next pi)
                    ca_done += [cb]
                    # Compute the absolute volume ratio log of ca/cb
                    rat = volume[ca] / volume[cb]
                    if rat < 1.0:
                        rat = 1 / rat
                    count_rat_all += 1
                    # Common points of ca and cb
                    pcom = tuple(sorted(p for p in cells_to_points[ca] if p in cells_to_points[cb]))
                    # if len(pcom) > 2:  # more than two common points for a face
                    if len(pcom) == 1:  # Only one common point
                        rat_node = max(rat_node, rat)
                        count_rat_node += 1
                    elif len(pcom) == 2:  # Two common points for an edge
                        rat_edge = max(rat_edge, rat)
                        count_rat_edge += 1
                    else:  # len(pcom) > 2:  # more than two common points for a face
                        rat_face = max(rat_face, rat)
                        count_rat_face += 1
        # CHK # # plot time distribution
        # CHK # if do_iter_check and out_is_tty:
        # CHK #     plt.plot(mytime[::10])
        # CHK #     plt.show()

        # Recover point_to_cells
        points_to_cells = bk_points_to_cells

        # MAN # # 2D manual example with 2x2 quad cells
        # MAN # class mm:
        # MAN #     points = [i for i in range(9)]
        # MAN # cellpoints = [[0,1,3,4], [1,2,4,5], [3,4,6,7], [4,5,7,8]]
        # MAN # pointcells = [[0], [0,1], [1], [0,2], [0,1,2,3], [1,3], [2], [2,3], [3]]
        # MAN # cells_to_points = { c: lp for c, lp in enumerate(cellpoints)}
        # MAN # points_to_cells = { p: lc for p, lc in enumerate(pointcells)}
        # MAN # cp_copy = {c: [p for p in lp] for c, lp in cells_to_points.items()}
        # MAN # for pi in mm.points:
        # MAN #     print("point", pi)
        # MAN #     pci = points_to_cells[pi]
        # MAN #     for ca in pci[:-1]:
        # MAN #         for cb in pci[1:]:
        # MAN #             pcom = tuple(sorted(p for p in cp_copy[ca] if p in cp_copy[cb]))
        # MAN #             if len(pcom) == 1:
        # MAN #                 print("single point:", pi, (ca,cb))
        # MAN #             elif len(pcom) == 2:
        # MAN #                 print(" edge points:", pcom, (ca,cb))
        # MAN #             else:
        # MAN #                 print(" face points:", pcom, (ca,cb))
        # MAN #         pci = pci[1:]
        # MAN #         cp_copy[ca].remove(pi)
        # MAN #         print("rmove", pi, "from", ca)

        log.info(f"  Computed ratios: {count_rat_all:8d}")
        log.info(f"    |   node ratios: {count_rat_node:8d}")
        log.info(f"    |   edge ratios: {count_rat_edge:8d}")
        log.info(f"    |   face ratios: {count_rat_face:8d}")
        log.info(f"    | max vol ratio /nodes  : {rat_node:12.6e}")
        log.info(f"    | max vol ratio /edges  : {rat_edge:12.6e}")
        log.info(f"    | max vol ratio /faces  : {rat_face:12.6e}")
        # CHK # # display check+debug counters
        # CHK # if do_iter_check:
        # CHK #     log.info(f"      cell__iter = {nb_iter_ca:6}")
        # CHK #     log.info(f"      c_pnt_iter = {nb_iter_pi:6}")
        # CHK #     log.info(f"      com_c_iter = {nb_iter_cb:6}")

        return True

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
            # self._grid.compute_cell_sizes(progress_bar=out_is_tty) is bugged
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

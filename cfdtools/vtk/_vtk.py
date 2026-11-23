import logging
from collections import defaultdict
from itertools import combinations

import numpy as np
import scipy.spatial as spatial
from typing import Callable
import sys

out_is_tty = sys.stdout.isatty()
idem = lambda x: x # pass #
try:
    if not out_is_tty:
        raise RuntimeError('not_a_tty')
    from tqdm import tqdm # pass #
except Exception as e:
    del e
    tqdm = idem # pass #

try:
    import pyvista as pv
    from pyvista import CellType

    importpyvista = True
    # PyVista 0.49 introduced accessors for the unpadded cell arrays.
    _HAS_CELL_ARRAY_ACCESSORS = all(
        hasattr(pv.UnstructuredGrid, name)
        for name in ('cell_offsets', 'cell_connectivity')
    )
except ImportError:
    importpyvista = False
    _HAS_CELL_ARRAY_ACCESSORS = False

import cfdtools.api as api
import cfdtools.meshbase._connectivity as _conn
import cfdtools.meshbase._elements as _elem
import cfdtools.meshbase._mesh as cfdmesh
import cfdtools.utils.maths as maths
import cfdtools.hdf5 as hdf5

log = logging.getLogger(__name__)


def _cell_point_ids(mesh, cell_id):
    """Return one cell's point ids across supported PyVista APIs."""
    if _HAS_CELL_ARRAY_ACCESSORS:
        offsets = mesh.cell_offsets
        return mesh.cell_connectivity[offsets[cell_id] : offsets[cell_id + 1]]
    return np.asarray(mesh.get_cell(int(cell_id)).point_ids, dtype=np.int64)


def _relations_from_incidence(incidence, excluded_pairs=()):
    """Build unique cell-pair relations from entity-to-cell incidence."""
    excluded_pairs = set(excluded_pairs)
    pairs = set()
    relations = []
    for entity_nodes, incident_cells in incidence.items():
        for pair in combinations(sorted(set(incident_cells)), 2):
            if pair not in excluded_pairs and pair not in pairs:
                pairs.add(pair)
                relations.append((entity_nodes, *pair))
    return relations, pairs


def _topological_relations(mesh):
    """Return exclusive face, edge and node cell relations.

    Relations are built from the explicit topology of the supported linear
    element types instead of inferring geometry from a shared-node count.
    """
    face_to_cells = defaultdict(list)
    edge_to_cells = defaultdict(list)
    node_to_cells = [[] for _ in range(mesh.n_points)]

    for cell_id, vtk_type in enumerate(mesh.celltypes):
        try:
            elemtype = ele_vtktype[vtk_type]
            element_faces = [
                face
                for faces in _elem.elem2faces[elemtype].values()
                for face in faces
            ]
        except KeyError as error:
            raise ValueError(
                f"unsupported VTK cell type {int(vtk_type)} in vtkMesh.diag()"
            ) from error
        if not element_faces:
            raise ValueError(
                f"no face topology available for element type {elemtype}"
            )

        point_ids = np.asarray(_cell_point_ids(mesh, cell_id), dtype=np.int64)
        for point_id in point_ids:
            node_to_cells[point_id].append(cell_id)

        for local_face in element_faces:
            face_nodes = tuple(sorted(point_ids[local_face].tolist()))
            face_to_cells[face_nodes].append(cell_id)

        if _elem.dim_elem[elemtype] == 3:
            local_edges = {
                tuple(sorted((face[i - 1], face[i])))
                for face in element_faces
                for i in range(len(face))
            }
            for local_edge in local_edges:
                edge_nodes = tuple(sorted(point_ids[list(local_edge)].tolist()))
                edge_to_cells[edge_nodes].append(cell_id)

    face_relations, face_pairs = _relations_from_incidence(face_to_cells)
    edge_relations, edge_pairs = _relations_from_incidence(
        edge_to_cells,
        excluded_pairs=face_pairs,
    )
    node_incidence = {
        (point_id,): incident_cells
        for point_id, incident_cells in enumerate(node_to_cells)
        if incident_cells
    }
    node_relations, _ = _relations_from_incidence(
        node_incidence,
        excluded_pairs=face_pairs | edge_pairs,
    )
    return {
        'Face': face_relations,
        'Edge': edge_relations,
        'Node': node_relations,
    }, node_to_cells


def _validate_cell_volumes(cell_volume):
    """Reject empty, non-finite and non-positive cell-volume arrays."""
    if cell_volume.size == 0:
        raise ValueError("cannot diagnose a mesh without cells")
    invalid = ~np.isfinite(cell_volume) | (cell_volume <= 0.0)
    if np.any(invalid):
        indices = np.flatnonzero(invalid)
        preview = ', '.join(map(str, indices[:10]))
        raise ValueError(
            f"{indices.size} cells have non-positive or non-finite volumes"
            f" (first indices: {preview})"
        )


def _volume_ratio_fields(mesh, cell_volume):
    """Compute cell and point volume-ratio fields by topological relation."""
    _validate_cell_volumes(cell_volume)
    relations, node_to_cells = _topological_relations(mesh)
    point_ratio = {
        name: np.zeros(mesh.n_points, dtype=float)
        for name in ('Face', 'Edge', 'Node')
    }
    cell_ratio = {
        name: np.ones(mesh.n_cells, dtype=float)
        for name in ('Face', 'Edge', 'Node')
    }
    maxima = {}
    counts = {}

    for name, entity_relations in relations.items():
        maximum = 1.0
        for entity_nodes, cell_a, cell_b in entity_relations:
            ratio = max(
                cell_volume[cell_a] / cell_volume[cell_b],
                cell_volume[cell_b] / cell_volume[cell_a],
            )
            cell_ratio[name][[cell_a, cell_b]] = np.maximum(
                cell_ratio[name][[cell_a, cell_b]],
                ratio,
            )
            point_ratio[name][list(entity_nodes)] = np.maximum(
                point_ratio[name][list(entity_nodes)],
                ratio,
            )
            maximum = max(maximum, ratio)
        maxima[name] = maximum
        counts[name] = len(entity_relations)

    # Give points without a direct relation the maximum value of their
    # incident cells. Isolated points remain NaN and are therefore visible.
    for point_id, incident_cells in enumerate(node_to_cells):
        for name in ('Face', 'Edge', 'Node'):
            if point_ratio[name][point_id] == 0.0:
                point_ratio[name][point_id] = (
                    np.max(cell_ratio[name][incident_cells])
                    if incident_cells
                    else np.nan
                )

    return point_ratio, cell_ratio, maxima, counts


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
            api.error_stop("vtkMesh cannot be initialized by both mesh and pvmesh structures")
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
            pvm = grid = self._grid
        else:
            log.warning("  no pyvista mesh available")
            return False

        # Using pyvista internel representation
        # log.info(f"pyvista object:")
        #
        # mstr_list = f"{pvm}".split(chr(10))
        # mstr_list[1:1] = ["Mesh"]
        # for l in mstr_list:
        #     log.info(f"  {l}")
        #
        # log.info(f"pyvista object: {pvm.DATA_TYPE_NAME}") # Attribute is a method
        # log.info(f"pyvista object: {pvm.DATA_TYPE_NAME()}") # Wrong output
        #
        log.info(f"pyvista object:")
        log.info(f"  {type(pvm).__name__} ({hex(id(pvm))}) [{pvm.__vtkname__}]")
        log.info(f"  Mesh")
        log.info(f"    | N Cells:  {pvm.n_cells}")
        log.info(f"    | N Points: {pvm.n_points}")
        for i, c in enumerate('XYZ'):
            cmin, cmax = pvm.bounds[i * 2 : i * 2 + 2]
            log.info(f"    |   {c} bounds:  {cmin:10.3e},{cmax:10.3e}")
        if hasattr(pvm, 'n_arrays'):
            log.info(f"    | N Arrays: {pvm.n_arrays}")
        log.info(f"  Data")
        log.info(f"    | cell  data names: {pvm.cell_data.keys()}")
        log.info(f"    | point data names: {pvm.point_data.keys()}")
        log.info(f"    | field data names: {pvm.field_data.keys()}")
        # log.info("  Properties")

        return True

    def diag(self, data=False):
        """Diagnose cell volumes and topological volume-ratio variations.

        When ``data`` is true, return a new :class:`vtkMesh` containing the
        volume and node-, edge- and face-neighbour ratio fields. Otherwise,
        log a summary and return ``True``. A mesh with invalid cell volumes is
        rejected because its ratios and geometric mean are undefined.
        """
        return_data = data

        if self._grid is None:
            log.warning("  no pyvista mesh available")
            return False
        pvm = self._grid

        # Print brief and return False if error
        if not self.brief():
            return False

        log.info("")
        full_timer = api.Timer("Computing volumes and size ratios")
        full_timer.enter()

        ### Compute cell volumes
        ###=====================

        with api.Timer(f"  Computing volumes of {pvm.n_cells} cells"):
            volume = CellVolume = self.volumes()
        _validate_cell_volumes(CellVolume)
        volmin, volmax = min(volume), max(volume)
        voltot, volava = sum(volume), np.mean(volume)
        vol_nb, volavg = len(volume), np.exp(np.mean(np.log(volume)))

        log.info(f"  Volumes")
        log.info(f"    | vol nb / total vol :    {vol_nb:12d} / {voltot:12.6e}")
        log.info(f"    | (arith / geom) avg :    {volava:12.6e} / {volavg:12.6e}")
        log.info(f"    | min vol / max vol  :    {volmin:12.6e} / {volmax:12.6e}")

        with api.Timer("  Building topology and computing cell size ratios"):
            PointVolumeRatio, CellVolumeRatio, maxima, counts = (
                _volume_ratio_fields(pvm, CellVolume)
            )

        count_rat_all = sum(counts.values())
        full_timer.stop(show=True)

        log.info(f"  Computed ratios: {count_rat_all:8d}")
        for name in ('Node', 'Edge', 'Face'):
            log.info(f"    | {name.lower():>6} ratios: {counts[name]:8d}")
        for name in ('Node', 'Edge', 'Face'):
            log.info(
                f"    | max vol ratio /{name.lower():<5}: {maxima[name]:12.6e}"
            )

        # Build and display histogram of ratios
        #
        min_rat_face = min(CellVolumeRatio['Face'])
        max_rat_face = maxima['Face']
        if np.isclose(max_rat_face, min_rat_face):
            rw = 2
            hseg = np.array([min_rat_face, max_rat_face])
            hist = (np.array([pvm.n_cells]), hseg)
        else:
            lrat = np.log(max_rat_face/min_rat_face)
            l10l = np.log10(lrat)
            rw = 2 - max(-6, (np.floor(l10l - .3) + 0.5).astype(int))
            nseg = max(2, 16 * min(1, lrat * 1E6))
            hseg = np.logspace(
                        np.log10(min_rat_face),
                        np.log10(max_rat_face), nseg + 1)
            hist = np.histogram(CellVolumeRatio['Face'], hseg)
        hmax = max(hist[0])
        hw_0 = (np.ceil(np.log10(hmax)) + 0.5).astype(int)
        wmax = 56
        wada = wmax - hw_0 - 2 * rw
        hwid = (hist[0] / hmax * wada + 0.5).astype(int)
        hada = hist[0] / pvm.n_cells

        log.info("")
        log.info("  Histogram of volume ratios /faces")
        log.info("  " + (wmax + 23) * "-")
        for u, v1, v2, h, ha in zip(*hist, hseg[1:], hwid, hada, ):
            log.info(f"  {v1:{rw+2}.{rw}f} "
                     f"- {v2:{rw+2}.{rw}f} : {u:{hw_0}}"
                     f"{' ' if ha < .99995 else ''}"
                     f" ({ha:6.2%}) |{h*'@'}|")
        log.info("  " + (wmax + 23) * "-")

        if return_data:
            new_pvm = pvm.copy(deep=False)
            # Remove all existing data
            for data_mode in ['cell', 'point', 'field']:
                dataset = getattr(new_pvm, f"{data_mode}_data")
                for k in dataset:
                    dataset.remove(k)
            new_vtkm = vtkMesh(pvmesh=new_pvm)
            new_vtkm._grid.cell_data['C_Volume'] = CellVolume
            new_vtkm._grid.point_data['P_FaceVR'] = PointVolumeRatio['Face']
            new_vtkm._grid.point_data['P_EdgeVR'] = PointVolumeRatio['Edge']
            new_vtkm._grid.point_data['P_NodeVR'] = PointVolumeRatio['Node']
            new_vtkm._grid.cell_data['C_FaceVR'] = CellVolumeRatio['Face']
            new_vtkm._grid.cell_data['C_EdgeVR'] = CellVolumeRatio['Edge']
            new_vtkm._grid.cell_data['C_NodeVR'] = CellVolumeRatio['Node']
            return new_vtkm

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
        if self._volume is None:
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
        cell_offsets locates each cell in the unpadded cell_connectivity array
        points is the list of points, shape (:,3)
        """
        with api.Timer("vtu extract slice"):
            mesh = self._grid
            selected_cell_ids = np.asarray(selected_cell_ids, dtype=np.int64)
            celltypes = mesh.celltypes
            points = mesh.points

            if _HAS_CELL_ARRAY_ACCESSORS:
                offsets = mesh.cell_offsets
                connectivity = mesh.cell_connectivity
                selected_vertices = [
                    connectivity[offsets[cid] : offsets[cid + 1]]
                    for cid in selected_cell_ids
                ]
            else:  # pragma: no cover - compatibility with PyVista < 0.49
                selected_vertices = [
                    np.asarray(mesh.get_cell(int(cid)).point_ids, dtype=np.int64)
                    for cid in selected_cell_ids
                ]

            new_celltypes = celltypes[selected_cell_ids]
            extract_vtex = np.concatenate(selected_vertices)
            nvtex = np.fromiter(
                (vertices.size for vertices in selected_vertices),
                dtype=np.int64,
                count=selected_cell_ids.size,
            )

            # Keep points in first-occurrence order and build their new indices.
            used_point_ids, first, inverse = np.unique(
                extract_vtex,
                return_index=True,
                return_inverse=True,
            )
            occurrence_order = np.argsort(first, kind='stable')
            sorted_to_occurrence = np.empty_like(occurrence_order)
            sorted_to_occurrence[occurrence_order] = np.arange(occurrence_order.size)
            new_vtex = sorted_to_occurrence[inverse]
            used_point_ids = used_point_ids[occurrence_order]

            # Rebuild VTK-style cell array
            vtk_cells = np.empty(len(new_vtex)+len(new_celltypes), dtype=np.int64)
            i = int(0) ; j = 0
            for size in nvtex:
                vtk_cells[i] = size
                vtk_cells[i+1:i+1+size] = new_vtex[j:j+size]
                i += np.int64(size) + 1 # size + 1, size is uint8
                j += np.int64(size)
            new_celltypes = np.asarray(new_celltypes, dtype=np.uint8) # cell type is small
            new_points = points[used_point_ids, :]
            # create extracted unstructured grid
            new_grid = pv.UnstructuredGrid(vtk_cells, new_celltypes, new_points)
            new_mesh = vtkMesh()
            new_mesh.set_pvmesh(new_grid)
        return new_mesh

    def vtk_zconvolution(self,
            splitcoords: Callable = None,
            tol = 1.e-6,
            nmode=0,
            snapshot=False, rms=False, phase=False,
            select: dict = None):
        """average the mesh data along a splitcoords and return a new vtkMesh object

        Options can be replaced by a dictionary of options which keys are celldata names
        e.g. select = {'Q': 'avg', 'U': ['no_avg', 'rms', {'nmode': 2}, 'phase' ]}
            will compute:
                the average of Q, and
                the average (??), RMS, 2 fourier modes of U (phase included),
            all other data will be ignored.

        Args:
            splitcoords (Callable, optional): function to split the coordinates.        Defaults to None.
            tol         (float,    optional): tolerance for distance between planes.    Defaults to 1.e-6.
            nmode       (int,      optional): number of modes to keep.                  Defaults to 0.
            snapshot    (bool,     optional): keep the snapshot data.                   Defaults to False.
            rms         (bool,     optional): compute the RMS of the data.              Defaults to False.
            select      (dict,     optional): dictionary of variables/options selection to apply. Defaults to None.
        Returns:
            vtkMesh: averaged mesh
        """
        # parse dict if necessary and apply general options
        current_options = {'avg': True, 'rms': rms, 'phase': phase, 'nmode' : nmode, 'snapshot': snapshot }
        if select is None:
            directives = { keyname: current_options.copy() for keyname in self._grid.cell_data.keys() }
        else:
            # dict value can be only one (string) key or a list of a dict (??)
            directives = {}
            unknown_keys = [key for key in select if key not in self._grid.cell_data.keys()]
            if unknown_keys:
                uk, ns = unknown_keys, 's'
                if len(uk) == 1:
                    uk, ns = uk[0], ''
                log.error(f"  unknown name{ns} {uk} in VTK cell data")
                raise ValueError(f"  unknown key{ns} {uk} used in select option of vtk_zconvolution()")
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

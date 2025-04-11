import logging

import numpy as np

try:
    import pyvista as pv
    from pyvista import CellType

    importpyvista = True
except ImportError:
    importpyvista = False

import cfdtools.api as api
import cfdtools.meshbase._connectivity as _conn
import cfdtools.meshbase._mesh as cfdmesh
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
        if mesh:
            self.set_mesh(mesh)
        if pvmesh:
            self.set_pvmesh(pvmesh)

    def _reset(self):
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
            return

        mm = f"{m}"
        log.info(f"pyvista object: {m}")
        log.info(f"pyvista object:")
        [log.info(f"|   {l}") for l in f"{m}".split(chr(10))]
        log.info(f"> mesh")
        log.info(f"   |  mesh bounds : {m.bounds}")
        log.info(f"   |  mesh ncells : {m.n_cells}")
        log.info(f"   |  mesh npoints: {m.n_points}")
        log.info(f"> data")
        log.info(f"   |  cell  data names: {m.cell_data.keys()}")
        log.info(f"   |  point data names: {m.point_data.keys()}")
        log.info(f"   |  field data names: {m.field_data.keys()}")
        # log.info("> properties")

    def diag(self):
        if self._grid:
            m = self._grid
        else:
            log.warning("  no pyvista mesh available")
            return

        log.info(f"pyvista object:")
        #self.brief()
        volume = self.volumes()
        volmin, volmax = min(volume), max(volume)
        voltot, volava = sum(volume), np.mean(volume)
        vol_nb, volavg = len(volume), np.exp(np.mean(np.log(volume)))

        log.info(f"> volumes")
        log.info(f"  vol nb / total vol :    {vol_nb:12d} / {voltot:12.6e}")
        log.info(f"  (arith / geom) avg :    {volava:12.6e} / {volavg:12.6e}")
        log.info(f"  min vol / max vol  :    {volmin:12.6e} / {volmax:12.6e}")

        cell_points = []
        point_cells = [[] for _ in range(m.n_points)]
        idxc = 0
        ic = 0
        while idxc < len(m.cells):
            nbp = m.cells[idxc]
            idxc += 1
            plist = m.cells[idxc:idxc + nbp]
            idxc += nbp
            # cell_points[ic] = [list of points of cell ic]
            cell_points += [ plist ]
            for ip in plist:
            # point_cells[ip] = [list of cells of point ip]
                point_cells[ip] += [ic]
            ic += 1
        lrat_node, lrat_edge, lrat_face = [0.0, 0.0, 0.0]
        c_done = []
        nb_log_comp = 0
        nb_node_count = 0
        nb_edge_count = 0
        nb_face_count = 0
        pdbg = False
        psum = False
        if pdbg:
            pnt_com_dict = {}
            edg_com_dict = {}
            fac_com_dict = {}
            print("processing")
        # loop on all cells: ca
        for ca in range(len(cell_points)):
            c_done += [ca]
            ca_done = []
            if pdbg:
                print(f"cell {ca:2}:")
            # points of ca not already processed (as points of a previous ca)
            ca_points = cell_points[ca]
            # loop on those points of ca: pi
            for pi in ca_points:
                if pdbg:
                    print(f"    point {pi:2}:") #, end='')
                # cells of pi not already processed (as a previous ca, or as cb for this ca)
                pi_cells = [c for c in point_cells[pi] if c not in c_done + ca_done]
                # loop on those cells of pi: cb
                for cb in pi_cells:
                    ca_done += [cb]
                    if pdbg:
                        print(f"        {cb:2}", end='')
                    # Compute the absolute volume ratio log of ca/cb
                    lrat = abs(np.log(volume[ca]/volume[cb]))
                    nb_log_comp += 1
                    # Common points of ca and cb
                    pcom = tuple(sorted(p for p in cell_points[ca] if p in cell_points[cb]))
                    if len(pcom) == 1: # Only one common point
                        lrat_node = max(lrat_node, lrat)
                        nb_node_count += 1
                        if pdbg:
                            print("   ", "pcom:", end='')
                            for p in pcom:
                                if p not in pnt_com_dict:
                                    pnt_com_dict[p] = []
                                pnt_com_dict[p] += [[ca, cb]]
                    elif len(pcom) == 2: # Two common points for an edge
                        lrat_edge = max(lrat_edge, lrat)
                        nb_edge_count += 1
                        if pdbg:
                            print("   ", "ecom:", end='')
                            for p in pcom:
                                if p == pi: continue
                                k = (pi, p)
                                if k not in edg_com_dict:
                                    edg_com_dict[k] = []
                                edg_com_dict[k] += [[ca, cb]]
                    else: # len(pcom) > 2: more than two common points for a face
                        lrat_face = max(lrat_face, lrat)
                        nb_face_count += 1
                        if pdbg:
                            print("   ", "fcom:", end='')
                            fac_com_dict |= {pcom: [[ca, cb]]}
                    if pdbg:
                        print("", pcom, [ca,cb])
        if pdbg:
            print("processing end")
        if psum:
            s_pnt_com_keys = sorted(pnt_com_dict.keys())
            s_edg_com_keys = sorted(edg_com_dict.keys())
            s_fac_com_keys = sorted(fac_com_dict.keys())
            print("pnt_com_dict:")
            pprev = None
            for p in [_ for _ in s_pnt_com_keys]:
                if p == pprev: continue
                pprev = p
                print("   ", {k: v for k, v in pnt_com_dict.items() if k == p})
            print("edg_com_dict:")
            eprev = None
            for e in [_ for _ in s_edg_com_keys]:
                if e == eprev: continue
                eprev = e
                print("   ", {k: v for k, v in edg_com_dict.items() if k == e})
            print("fac_com_dict:")
            fprev = None
            for f in [_ for _ in s_fac_com_keys]:
                if f == fprev: continue
                fprev = f
                print("   ", {k: v for k, v in fac_com_dict.items() if k == f})
        log.info(f"computed ratios: {nb_log_comp:8d}")
        log.info(f"    node ratios: {nb_node_count:8d}")
        log.info(f"    edge ratios: {nb_edge_count:8d}")
        log.info(f"    face ratios: {nb_face_count:8d}")
        rat_node = np.exp(lrat_node)
        rat_edge = np.exp(lrat_edge)
        rat_face = np.exp(lrat_face)
        log.info(f"  max vol ratio /nodes  : {rat_node:12.6e}")
        log.info(f"  max vol ratio /edges  : {rat_edge:12.6e}")
        log.info(f"  max vol ratio /faces  : {rat_face:12.6e}")


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

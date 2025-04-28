import logging
import os
from pathlib import Path

import numpy as np
import scipy.spatial as spatial

try:
    import pyvista as pv
    importpyvista = True
except ImportError:  # pragma: no cover
    importpyvista = False

import cfdtools.api as api
import cfdtools.utils.maths as maths
from cfdtools.vtk import vtkMesh
from cfdtools.data import DataSetList

# import cfdtools.hdf5 as hdf5

log = logging.getLogger(__name__)


class vtkList:
    def __init__(self, filelist, verbose=False) -> None:
        self._list = filelist
        self._verbose = verbose

    @property
    def nfile(self):
        return len(self._list)

    def allexist(self):
        return all(Path(file).exists() for file in self._list)

    def check_order(self, pos='cellcenter', tol=1e-10):
        mappos = {
            'cellcenter': lambda m: m.cell_centers().points,
            'node': lambda m: m.points,
        }
        count = 0
        assert self.nfile > 0
        ref = mappos[pos](pv.read(self._list[0]))
        for name in self._list:
            mesh = pv.read(name)
            d = np.max(maths.distance(ref, mappos[pos](mesh)))
            if d > tol:
                count += 1
                if self._verbose:
                    log.info(f"  . {name}: {d}")
        if self._verbose:
            log.info(f"  {count}/{self.nfile} grids are not {pos}-coincident")
        return count == 0

    def read(self, filterdata=None, reorder=False, tol=1e-10, tolverbose=1e-10):
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
        if self._verbose:
            log.info("> build kd-tree")
        Tsort.start()
        tree = spatial.KDTree(ctrRef)
        Tsort.pause()
        #
        self._data = DataSetList(self.nfile, Xrep='cellaverage', Trep='instant')
        # may add alive-progress or other
        if self._verbose:
            log.info("> read all files, get only CELL data")
            if filterdata:
                log.info(f"  select only {', '.join(filterdata)}")
        for name in self._list:
            if self._verbose:
                log.info(f"  - {name}")
            Tread.start()
            vtk = pv.read(name)
            Tread.pause()
            namelist = filterdata if filterdata else vtk.cell_data.keys()
            Tcomp.start()
            ctr = vtk.cell_centers().points
            d = np.max(maths.distance(ctrRef, ctr))
            Tcomp.pause()
            if d > tol:  # should be sized by domain size
                count += 1
                # if self._verbose:
                #     log.info(f"  . {name}: {d}")
                Tsort.start()
                dfinal, index = tree.query(ctr, p=2) # p=2 is the norm
                # reverse indexing to sort new arrays
                rindex = index.copy()
                rindex[index] = np.arange(index.size)
                dmin, davg, dmax = maths.minavgmax(dfinal)
                if dmax > tolverbose:
                    log.warning(f"  max distance if above tolerance {tolverbose}: {np.count_nonzero(dfinal > tolverbose)} cells incriminated")
                    log.warning(f"  min:avg:max = {dmin:.2e}:{davg:.2e}:{dmax:.2e}")
                if dmax > tol:
                    log.error(f"  max distance is above tolerance {tol}: {np.count_nonzero(dfinal > tol)} cells incriminated")
                    log.error(f"  min:avg:max = {dmin:.2e}:{davg:.2e}:{dmax:.2e}")
                    raise ValueError("  max distance is above tolerance")
                Tsort.pause()
                # automatically deals with different shapes
                datalist = {name: vtk.cell_data[name][rindex] for name in namelist}
            else:
                datalist = {name: vtk.cell_data[name] for name in namelist}
            time = vtk.field_data.get("TimeValue", None)
            self._data.add_datalist(datalist, time=time)
        if self._verbose:
            log.info(f"  {count}/{self.nfile} grids were not coincident")
            log.info(f"       file reading: {Tread.elapsed:.2f}s")
            log.info(f"    grid comparison: {Tcomp.elapsed:.2f}s")
            log.info(f"    data reordering: {Tsort.elapsed:.2f}s")

    def dumphdf(self, filename: str, overwrite: bool = False, xdmf: bool = False, **options):
        """Convert a list of VTK files into a single HDF5 file.

        :param str filename: Path of the output HDF5 file.
        :param xdmf: Output a XMF file as well if True.
        :param options: Options for h5py.
        :type options: dict
        :return: The effective output HDF5 filename.
        :rtype: str
        """
        vtkmesh = vtkMesh(pvmesh=self._mesh)
        self._data.set_mesh(vtkmesh.export_mesh())
        actual_filename = self._data.dumphdf(filename, overwrite=overwrite, **options)
        if xdmf:
            self._dumpxdmf(actual_filename, vtkmesh, self._data)
        return actual_filename

    def _dumpxdmf(self, filepath: str, vtkmesh: vtkMesh, data: DataSetList):
        """Create the XMF file associated to the single HDF5 file created by dumphdf()."""

        def xdmf_domain(content):
            lines = ['<Xdmf Version="3.0">']
            lines += ["<Domain>"]
            lines += content
            lines += ["</Domain>"]
            lines += ["</Xdmf>"]
            return lines

        # get the XDMF content for the mesh geometry and topology
        geometry_content = vtkmesh.xdmf_content(filepath)
        # get the XDMF content for the time data series
        geom_dat_content = data.xdmf_content(filepath, geometry_content)
        content = xdmf_domain(geom_dat_content)

        # write the XMF file
        with open(os.path.splitext(filepath)[0] + ".xmf", "w") as fid:
            fid.write('\n'.join(content))

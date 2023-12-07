# -*- coding: utf-8 -*-
"""Module for writing ic3 restart files in `HDF5`_ format.

Todo:
    * rework inheritance

.. _HDF5:
   https://www.hdfgroup.org/

"""
# Standard library imports

# Third party imports
import h5py

# cfdtools imports
import cfdtools.api as api
from cfdtools import __version__
import cfdtools.meshbase._mesh as _mesh

from cfdtools.ic3._ic3 import type2zonekind
from cfdtools.ic3 import writerV2


@api.fileformat_writer('IC3', '.h5')
class writer(writerV2.writer):
    """Writer of ic3 restart files in HDF5 format."""

    __version__ = "4"

    def __init__(self, mesh: _mesh.Mesh, endian='native'):
        """Initialization of a ic3 restart file writer."""
        super().__init__(mesh, endian)

    def write_data(self, mesh_filename: str, solution_filename: str = None):
        """Write the ic3 restart file in the HDF5 format."""
        api.io.printstd(f"> Writing file {mesh_filename!r}")
        # Open the mesh file
        with h5py.File(mesh_filename, "w") as fid:
            api.io.printstd("> Writing mesh attributes")
            self._write_root_attributes(fid)

            api.io.printstd("> Writing coordinates")
            self._write_coordinates(fid)

            api.io.printstd("> Writing connectivities")
            self._write_connectivities(fid)

            api.io.printstd("> Writing boundaries")
            self._write_boundaries(fid)

            api.io.printstd("> Writing informative values")
            self._write_mesh_parameters()

            # Not needed up to now
            # api.io.printstd("> Writing variables")
            # self._write_mesh_vars(fid)

            api.io.printstd("> End of Mesh File")

        if solution_filename is None:
            return

        # Open the solution file
        with h5py.File(solution_filename, "w") as fid:
            api.io.printstd("> Writing solution attributes")
            self._write_root_solution_attributes(fid)

            api.io.printstd("> Writing variables")
            self._write_solution_vars(fid)

            api.io.printstd("> Writing statistics")
            self._write_solution_stats(fid)

            api.io.printstd("> End of Solution File")

    def _write_root_attributes(self, fid: h5py.File):
        """Write the root attributes of a mesh file."""
        fid.attrs.update(
            {
                'cfdtools_version': __version__,
                'hdf5_version': h5py.version.hdf5_version,
                'ic3_format_version': 4,
                'meshtype': "unstructured",
                'no_count': self._mesh.nnode,
                "fa_count": self._mesh.nface,
                "cv_count": self._mesh.ncell,
            }
        )
        api.io.printdebug("  " + ', '.join(f'{key}:{value}' for key, value in fid.attrs.items()))

    def _write_coordinates(self, fid: h5py.File):
        """Write the coordinates of a mesh file."""
        coordinates = self._mesh.nodescoord(ndarray=True)
        api.io.printdebug('  Coordinates (IC3 UGP_IO_X_NO)')
        fid.create_dataset("coordinates", data=coordinates.ravel())
        fid["coordinates"].attrs['system'] = 'cartesian'

    def _write_connectivities(self, fid: h5py.File):
        """Write the connectivities of a mesh file."""
        api.io.printdebug('  Group Connectivities')
        conn_group = fid.create_group("/Connectivities")

        if self._mesh._cell2node is not None:
            api.io.printdebug('  Element-Node Connectivity')
            cell_group = conn_group.create_group("Cell")
            for ctype, cellco in dict(self._mesh._cell2node).items():
                cell_group.create_dataset(ctype, data=cellco)

        api.io.printdebug('  Face-Node Connectivity')
        face_group = conn_group.create_group("Face")
        api.io.printdebug('  Face-Element Connectivity (CGNS ParentElements, IC3 CVOFA)')
        face_group.create_dataset('cvofa', data=self.f2e.ravel().tolist(), dtype='i4')
        face_group["cvofa"].attrs['synonyms'] = 'CGNS NFACE_n ParentElements'

        api.io.printdebug('  Face-Node Connectivity (CGNS NGON_n, IC3 NOOFA_I_AND_V)')
        # Node count per face
        # ElementStartOffset to be computed?
        face_group.create_dataset('noofa_i', data=self.f2v["noofa"], dtype='i4')
        face_group["noofa_i"].attrs['synonyms'] = 'CGNS NGON_n ElementStartOffset'
        face_group.create_dataset('noofa_v', data=self.f2v["noofa_v"], dtype='i4')
        face_group["noofa_v"].attrs['synonyms'] = 'CGNS NGON_n'

        return

    def _write_boundaries(self, fid: h5py.File):
        """Write the boundaries of a mesh file."""
        # Face zones, a.k.a boundary condition patches
        api.io.printdebug("  Markers of faces / Face zones / boundary condition patches")
        bnd_group = fid.create_group("Boundaries")
        last_boco = 0
        labels = []
        kinds = []
        max_range = []
        min_range = []
        for key, boco in self.bocos.items():
            assert key == boco.name
            assert boco.geodim in ('face', 'bdface'), "boco marks must be faces index"
            assert boco.type in type2zonekind.keys(), f"unsupported type of boco for IC3 output: {boco.type}"
            assert boco.index.type == 'range', "indexing must be a range and may need reordering"

            labels.append(key)
            ifmin, ifmax = boco.index.range()
            min_range.append(ifmin)
            max_range.append(ifmax)
            api.io.printdebug(f"  - ({boco.type}) {boco.name}: {ifmin}-{ifmax}")
            kinds.append(type2zonekind[boco.type])

            last_boco = max(last_boco, ifmax)
            if "periodic_transform" in boco.properties:
                if boco.properties["periodic_transform"] is not None:
                    per_group = bnd_group.create_group(boco.name)
                    per_group.create_dataset(
                        "periodic_transform",
                        3,
                        dtype="f8",
                        data=boco.properties["periodic_transform"][:3],
                    )

        labels.append('internal-domain')
        kinds.append(type2zonekind['internal'])
        ifmin, ifmax = last_boco + 1, self.params['fa_count']
        api.io.printdebug(f"   - additional mark (FA_ZONE) for internal faces: {ifmin}-{ifmax}")
        min_range.append(ifmin)
        max_range.append(ifmax)

        min_range, labels, kinds = zip(*sorted(zip(min_range, labels, kinds), key=lambda x: x[0]))
        min_range = list(min_range)
        min_range.append(ifmax)

        dt = h5py.string_dtype(encoding='utf-8', length=128)
        bnd_group.create_dataset('labels', len(labels), dtype=dt, data=labels)
        bnd_group.create_dataset('offsets', len(min_range), dtype='i4', data=min_range)
        bnd_group.create_dataset('kinds', len(kinds), dtype='i4', data=kinds)

        return

    def _write_mesh_parameters(self):
        """Write information about the mesh file."""
        return

    def _write_mesh_vars(self, fid: h5py.File):
        """Write all the variables into a mesh file.

        Scalars, vectors and tensors all together.
        """
        # Start with the vertex variables Scalar Vector Tensor
        nno = self.params["no_count"]

        nd_group = fid.create_group("NodeData")
        for ndname, nddata in self.vars["nodes"].items():
            # Scalar
            api.io.printdebug('  Node Data')
            if len(nddata.shape) == 1:
                api.io.printdebug('  ' + ndname)
                nd_group.create_dataset('node_gb_index', data=nddata, dtype='i8')

        # Start with the cell variables Scalar Vector Tensor
        return

    def _write_root_solution_attributes(self, fid: h5py.File):
        """Write the root attributes of a solution file."""
        fid.attrs.update(
            {
                "step": self._mesh._params["step"],
                "time": self._mesh._params["time"],
            }
        )
        api.io.printdebug("  " + ', '.join(f'{key}:{value}' for key, value in fid.attrs.items()))

    def _write_solution_vars(self, fid: h5py.File):
        """Write all the variables into a solution file.

        Scalars, vectors and tensors all together.
        """
        # Start with the vertex variables Scalar Vector Tensor

        # Start with the cell variables Scalar Vector Tensor
        if all(x in self.vars["cells"].keys() for x in ['RHO', 'RHOU', 'RHOE']):
            api.io.printdebug('  Cell Data')
            for cvname in ['RHO', 'RHOU', 'RHOE']:
                cvdata = self.vars["cells"][cvname]
                api.io.printdebug('  ' + cvname)
                fid.create_dataset(cvname.lower(), data=cvdata.ravel(order='C'), dtype='f8')
        return

    def _write_solution_stats(self, fid: h5py.File):
        """Write all the statistics into a solution file.

        Scalars and vectors all together.
        """
        # Vertex variables

        # Cell variables
        stats_var = [
            var for var in self.vars["cells"].keys() if any(x in var for x in ['_AVG', '_RMS', '_REY'])
        ]
        if not stats_var:
            return

        volume_group = fid.create_group("VolumeStats")

        weights = {name: value for name, value in self._mesh._params.items() if name.endswith("_wgt")}
        volume_group.attrs.update(weights)

        for cvname in stats_var:
            cvdata = self.vars["cells"][cvname]
            api.io.printdebug('  ' + cvname)
            volume_group.create_dataset(cvname.lower(), data=cvdata.ravel(order='C'), dtype='f8')

        return

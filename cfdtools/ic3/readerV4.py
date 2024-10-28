# coding: utf8
"""Module for reading ic3 restart files in `HDF5`_ format.

Todo:
    * rework inheritance

.. _HDF5:
   https://www.hdfgroup.org/

"""
# Standard library imports
import logging
from pathlib import Path

# Third party imports
import h5py
import numpy as np

# cfdtools imports
import cfdtools.api as api
import cfdtools.meshbase._mesh as _mesh
import cfdtools.data as _data
import cfdtools.meshbase._connectivity as _conn
from cfdtools.ic3._ic3 import zonekind2type

log = logging.getLogger(__name__)


@api.fileformat_reader('IC3V4', '.h5')
class reader:
    """Reader of ic3 restart files in HDF5 format."""

    __version__ = "4"

    def __init__(self, mesh_filename: str, solution_filename: str = None):
        """Initialization of a ic3 restart file reader."""
        self._mesh_filename = Path(mesh_filename)
        if not self._mesh_filename.exists():
            raise FileExistsError()
        self._solution_filename = solution_filename
        if solution_filename is not None:
            self._solution_filename = Path(solution_filename)
            if not self._solution_filename.exists():
                raise FileExistsError()
        # Initialize the mesh
        self.mesh = {
            "params": {
                "no_count": 0,
                "fa_count": 0,
                "cv_count": 0,
                "noofa_count": 0,
                "nboco": 0,
            },
            "connectivity": {"noofa": {}, "cvofa": {}, "nkeys": 0},
            "coordinates": None,
            "bocos": {"nfa_b": 0, "nfa_bp": 0},
            "partition": None,
        }
        # Initialize the state dictionary
        self.simulation_state = {"step": 0, "time": 0}
        # Initialize the variable dictionary
        self.variables = {"nodes": {}, "cells": {}, "faces": {}}
        self.celldata = _data.DataSet('cellaverage')
        self.nodedata = _data.DataSet('nodal')
        self.facedata = _data.DataSet('nodal')

    @property
    def ncell(self):
        return self._ncell

    def __str__(self):
        s = '%22s' % 'mesh filename: ' + str(self._mesh_filename) + '\n'
        s += '%22s' % 'solution filename: ' + str(self._solution_filename) + '\n'
        s += '%22s' % 'simulation: ' + str(list(self.simulation_state.keys())) + '\n'
        s += '%22s' % 'mesh keys: ' + str(list(self.mesh.keys())) + '\n'
        s += '%22s' % 'variable keys: ' + str(list(self.variables.keys()))
        return s

    def printinfo(self):
        print(self)
        log.info("- mesh properties")
        for key, item in self.mesh.items():
            if isinstance(item, dict):
                for key2, item2 in item.items():
                    if isinstance(item2, dict):
                        for key3, item3 in item2.items():
                            api._printreadable('  mesh.' + key + '.' + key2 + '.' + key3, item3)
                    else:
                        api._printreadable('  mesh.' + key + '.' + key2, item2)
            else:
                api._printreadable('  mesh.' + key, item)
        log.info("- variable properties")
        for key, item in self.variables.items():
            for key2, value in item.items():
                api._printreadable('variables.' + key + '.' + key2, value)
                # for key3,value3 in value.items():
                #     print('variables.'+key+'.'+key2+'.'+key3+':'+str(type(value3)))

    def read_data(self):
        """Read the ic3 restart file in the HDF5 format."""
        log.info("Reader restart IC3 V4")

        # Open the file for binary reading
        log.debug(f"Opening file {self._mesh_filename!r}")
        with h5py.File(self._mesh_filename, "r") as fid:
            # log.info("Reading version")
            # self._read_version()

            log.info("Reading connectivity...")
            self._read_connectivity(fid)

        if self._solution_filename is None:
            return

        log.debug(f"Opening file {self._solution_filename!r}")
        with h5py.File(self._solution_filename, "r") as fid:
            log.info("Reading variables...")
            self._read_variables(fid)

            log.info("Reading informative values...")
            self._read_informative_values(fid)

            self._ncell = self.mesh['params']['cv_count']  # for generic writer and timer

    def export_mesh(self):
        """Export Mesh."""
        meshdata = _mesh.Mesh(self.mesh['params']['cv_count'], self.mesh['params']['no_count'])
        meshdata.set_nodescoord_nd(self.mesh['coordinates'])
        face2cell = _conn.indexindirection(self.mesh['connectivity']['cvofa']['cvofa'])
        face2node = _conn.elem_connectivity()
        face2node.importfrom_compressedindex(self.mesh['connectivity']['noofa'])
        meshdata.add_faces('mixed', face2node, face2cell)
        for boco in self.mesh['bocos']:
            meshdata.add_boco(boco)
        meshdata.set_celldata(self.celldata)
        meshdata.set_nodedata(self.nodedata)
        meshdata.set_facedata(self.facedata)
        meshdata.set_params(self.mesh['params'])
        meshdata.update_params(self.simulation_state)
        if 'partition' in self.mesh.keys():
            meshdata.set_partition(self.mesh['partition'])
        return meshdata

    def _read_connectivity(self, fid: h5py.File):
        """Read the connectivities of a mesh file."""
        # Store the size informations at the right place
        # and set the local shortcut names
        # number of nodes/vertices
        no_count = self.mesh["params"]["no_count"] = fid.attrs["no_count"]
        # number of faces
        fa_count = self.mesh["params"]["fa_count"] = fid.attrs["fa_count"]
        # number of cells/volume elements
        cv_count = self.mesh["params"]["cv_count"] = fid.attrs["cv_count"]
        log.info(
            f"mesh with {cv_count} cells, {fa_count} faces and {no_count} nodes",
        )
        # Face-based connectivity
        #
        # - First, NOOFA
        #
        log.info("Parsing face to node connectivity...")

        # Get the node count per face
        nno_per_face = fid["/Connectivities/Face/noofa_i"][()]
        # number of nodes from faces
        noofa_count = self.mesh["params"]["noofa_count"] = np.sum(nno_per_face)

        uniq, counts = np.unique(nno_per_face, return_counts=True)
        for facesize, nfacesize in zip(uniq, counts):
            log.info(f"  {nfacesize} faces of {facesize} nodes")
        # Initialize the proper connectivity arrays in self.mesh
        face2node_index = np.concatenate(([0], np.cumsum(nno_per_face)))
        # np.sum(nno_per_face) == noofa_count
        assert face2node_index[0] == 0
        assert face2node_index[-2] == noofa_count - nno_per_face[-1]
        # store in 8 bytes
        face2node_value = fid["/Connectivities/Face/noofa_v"][()].astype(np.int64)
        zface2node = _conn.compressed_listofindex(face2node_index, face2node_value)
        zface2node.check()
        self.mesh["connectivity"]["noofa"] = zface2node
        log.info("end of face/vertex connectivity")

        del nno_per_face, uniq, counts
        #
        # - Second, CVOFA
        #
        log.info("Parsing face to cell connectivity...")

        # store in 8 bytes
        self.mesh["connectivity"]["cvofa"]["cvofa"] = (
            fid["/Connectivities/Face/cvofa"][()].astype(np.int64).reshape((fa_count, 2))
        )

        log.info("end of face/cell connectivity")

        # Checks and a few associations
        assert self.mesh["connectivity"]["cvofa"]["cvofa"].max() < cv_count
        assert self.mesh["connectivity"]["cvofa"]["cvofa"].max() == cv_count - 1
        uniq, counts = np.unique(self.mesh["connectivity"]["cvofa"]["cvofa"][:, 1], return_counts=True)
        # Number of assigned boundary faces
        try:
            iwhere = np.where(uniq == -1)[0][0]
        except IndexError:
            self.mesh["params"]["nfa_b"] = 0
        else:
            self.mesh["params"]["nfa_b"] = counts[iwhere]
        del uniq, counts
        # Number of periodic boundary faces
        self.mesh["params"]["nfa_bp"] = np.count_nonzero(
            self.mesh["connectivity"]["cvofa"]["cvofa"][:, 1] < -1
        )
        # All periodic boundary faces to -1
        # mask = self.mesh["connectivity"]["cvofa"]["cvofa"][:,1] < -1
        # self.mesh["connectivity"]["cvofa"]["cvofa"][:,1][mask] = -1
        # print("RC",self.mesh["connectivity"]["cvofa"]["cvofa"])
        #
        # The boundary conditions now
        log.info("Parsing boundary conditions...")
        self.mesh['bocos'] = []  # init list of bocos
        if 'Boundaries' in fid:
            offsets = fid['Boundaries']['offsets'][()]
            min_range = offsets[:-1]
            max_range = offsets[1:] - 1
            labels = fid['Boundaries']['labels'][()]
            kinds = fid['Boundaries']['kinds'][()]
            for label, kind, mini, maxi in zip(labels, kinds, min_range, max_range):
                log.info("%s %d %d", label, mini, maxi)
                # 3 ints: kind, face range (begin, start)
                self.mesh["params"]["nboco"] += 1
                boco = _mesh.submeshmark(label.decode())
                boco.type = zonekind2type[kind]
                boco.geodim = 'intface' if boco.type == 'internal' else 'bdface'
                boco.index = _conn.indexlist(irange=[mini, maxi])
                if label in fid['Boundaries']:
                    if 'periodic_transform' in fid['Boundaries'][label]:
                        boco.properties["periodic_transform"] = fid['Boundaries'][label][
                            'periodic_transform'
                        ][()]
                        log.info("  periodic transformation: %s", boco.properties["periodic_transform"])
                famin, famax = boco.index.range()
                sta = face2node_index[famin]
                try:
                    sto = face2node_index[famax + 1]
                except IndexError:
                    sto = face2node_value.size

                boco.properties["slicing"] = np.unique(face2node_value[sta:sto])
                self.mesh['bocos'].append(boco)

        log.info("end of boundary conditions")

        # Parse the header of the partition information
        # log.info("Parsing partitioning information...")
        # h = restartSectionHeader()
        # h.readReqVar(self.fid, self.byte_swap, ["UGP_IO_CV_PART"])
        # self.mesh["partition"] = {}
        # self.mesh["partition"]['npart'] = h.idata[1]
        # self.mesh["partition"]['icvpart'] = np.zeros((cv_count,), dtype=np.int32)
        # s = BinaryRead(self.fid, "%di" % cv_count, self.byte_swap, type2nbytes["int32"] * cv_count)
        # self.mesh["partition"]['icvpart'] = np.asarray(s)
        # log.info("end of partition")

        # The coordinates of the vertices finally
        log.info("Parsing vertices coordinates...")
        self.mesh["coordinates"] = np.ascontiguousarray(fid['coordinates'][()].reshape(no_count, 3))
        log.info("end of node coordinates")

    def _read_informative_values(self, fid: h5py.File):
        """Read all the values also stored in a restart file.

        i.e. the step number, the time, the timestep.
        input:  handle on an open restart file, [type file identifier]
                endianness flag [boolean]
        output: simulation state structure containing informations about the current state
                of the simulation
        """
        for name in fid.attrs:
            self.simulation_state[name] = fid.attrs[name]
        if 'VolumeStats' in fid:
            self._read_informative_values(fid['VolumeStats'])

    def _read_variables(self, fid: h5py.File):
        """Read all the variables from a solution file."""
        # set the local shortcut names
        no_count = self.mesh["params"]["no_count"]
        fa_count = self.mesh["params"]["fa_count"]
        cv_count = self.mesh["params"]["cv_count"]

        def read_vars(fid: h5py.File, var_names):
            for var in var_names:
                var_up = var.upper()
                field = fid[var][()]
                nb_fields = field.size // cv_count
                if nb_fields > 1:
                    field = field.reshape(-1, nb_fields)
                self.celldata.add_data(var_up, field)
                if nb_fields > 1:
                    for component in range(field.shape[1]):
                        print_stats_scalar(var_up + "-%d" % component, field[:, component])
                else:
                    print_stats_scalar(var_up, field)

        read_vars(fid, ['rho', 'rhou', 'rhoe'])

        if 'VolumeStats' in fid:
            read_vars(fid['VolumeStats'], fid['VolumeStats'].keys())


def print_stats_scalar(name, dnp):
    """Print statistics for a scalar variable."""
    log.info(
        "  %s%s:  %+.5e / %+.5e / %+.5e (min/mean/max)."
        % (name, ' ' * (20 - len(name)), dnp.min(), np.mean(dnp), dnp.max())
    )

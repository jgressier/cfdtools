# coding: utf8

import logging
import numpy as np
import sys

import cfdtools.api as api
import cfdtools.meshbase._mesh as _mesh
import cfdtools.data as _data
import cfdtools.meshbase._connectivity as _conn

from cfdtools.ic3._ic3 import (
    binreader,
    type2nbytes,
    restartSectionHeader,
    ic3_restart_codes,
    BinaryRead,
    zonekind2type,
    properties_ugpcode,
)

log = logging.getLogger(__name__)
###################################################################################################


def print_stats_scalar(name, dnp):
    '''Print statistics for a scalar variable.

    :param name: Name of the scalar variable.
    :type name: str
    :param dnp: numpy array of the scalar variable.
    :type dnp: ndarray
    '''
    log.info(
        "  %s%s:  %+.5e / %+.5e / %+.5e (min/mean/max)."
        % (name, ' ' * (20 - len(name)), dnp.min(), np.mean(dnp), dnp.max())
    )


def print_stats_vector(name, dnp):
    '''Print statistics for a vector variable.

    :param name: Name of the vector variable.
    :type name: str
    :param dnp: numpy array of the vector variable.
    :type dnp: ndarray
    '''
    for component in range(dnp.shape[1]):
        log.info(
            "  %s%s%s:  %+.5e / %+.5e / %+.5e (min/mean/max)."
            % (name, "-%d" % component, ' ' * (20 - len(name)), dnp.min(), np.mean(dnp), dnp.max())
        )


@api.fileformat_reader('IC3', '.ic3')
class reader(binreader):
    '''Implementation of the reader to read IC3 restart files.'''

    def __init__(self, filename, cIntegrity=False):
        '''
        Initialization of an IC3 restart reader.
        Just save the filename and a boolean for an integrity check.
        input   : IC3 restart file name [type string]
                  whether to check the integrity of the file beforehand [type boolean]
        '''
        super().__init__(filename)
        self.check_integrity = cIntegrity
        self.ic3_version = -1

    @property
    def ncell(self):
        return self._ncell

    def __str__(self):
        s = '    filename: ' + self.filename
        s += '\n   simulation: ' + str(list(self.simulation_state.keys()))
        s += '\n    mesh keys: ' + str(list(self.mesh.keys()))
        s += '\nvariable keys: ' + str(list(self.variables.keys()))
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
        '''
        Main method of the IC3 restart reader.
        Parses in order the file using sub-methods described below.
        output: the mesh itself
                a lot of variables stored in the restart file
                information on the state of the simulation
        '''
        log.info("READER RESTART IC3")

        if not self.exists():  # pragma: no cover
            api.error_stop("Fatal error. File %s cannot be found." % (self.filename))

        # Open the file for binary reading
        log.debug(f"Opening file {self.filename!r}")
        with open(self.filename, "rb") as self.fid:
            #
            log.info("Reading binary file header")
            self._ReadRestartHeader()
            #
            log.info("Reading connectivity...")
            self._ReadRestartConnectivity()
            #
            log.info("Reading informative values...")
            self._ReadInformativeValues()
            #
            log.info("Reading variables...")
            self.celldata = _data.DataSet('cellaverage')
            self.nodedata = _data.DataSet('nodal')
            self.facedata = _data.DataSet('nodal')
            self._ReadRestartVar()
            #
            self._ncell = self.mesh['params']['cv_count']  # for generic writer and timer

        # file identifier as shortcut only
        del self.fid

    def export_mesh(self):
        # return (
        #     self.mesh["coordinates"],
        #     self.mesh["connectivity"]["e2v"],
        #     self.mesh["bocos"],
        #     self.variables["nodes"],
        #     self.variables["cells"],
        #     (self.simulation_state, self.mesh["params"]),
        # )
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

    def _ReadRestartConnectivity(self):
        '''
        Method reading the first blocks passed the header, containing informations
        on the nodes, the faces, the cells and the connectivity in between those
        input:  handle on an open restart file [type file identifier]
                endianness flag [boolean]
        output: mesh structure containing all the necessary information to build the grid
        '''

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

        # Get the header
        h = restartSectionHeader()
        h.readReqVar(self.fid, self.byte_swap, ["UGP_IO_NO_FA_CV_NOOFA_COUNTS"])
        # Store the size informations at the right place
        # and set the local shortcut names
        # number of nodes/vertices
        no_count = self.mesh["params"]["no_count"] = h.idata[0]
        # number of faces
        fa_count = self.mesh["params"]["fa_count"] = h.idata[1]
        # number of cells/volume elements
        cv_count = self.mesh["params"]["cv_count"] = h.idata[2]
        # number of nodes from faces
        noofa_count = self.mesh["params"]["noofa_count"] = h.idata[3]
        log.info(
            f"mesh with {cv_count} cells, {fa_count} faces and {no_count} nodes",
        )
        # Integrity check
        if self.check_integrity:
            # Check the restart is whole by counting the global ids of no, fa and cv
            # For nodes
            log.info("  Checking nodes integrity ..")
            sys.stdout.flush()
            h = restartSectionHeader()
            if h.readVar(self.fid, self.byte_swap, ["UGP_IO_NO_CHECK"]):
                assert h.idata[0] == no_count
                assert h.id[0] == ic3_restart_codes["UGP_IO_NO_CHECK"]
                s = BinaryRead(self.fid, "%di" % no_count, self.byte_swap, type2nbytes["int32"] * no_count)
                nodes_id = np.asarray(s)
                assert np.allclose(nodes_id, np.arange(no_count))
                log.info("end of node integrity")
                sys.stdout.flush()
                del nodes_id
                # For faces
                log.info("  Checking faces integrity ..")
                sys.stdout.flush()
                h = restartSectionHeader()
                h.readReqVar(self.fid, self.byte_swap, ["UGP_IO_FA_CHECK"])
                s = BinaryRead(self.fid, "%di" % fa_count, self.byte_swap, type2nbytes["int32"] * fa_count)
                faces_id = np.asarray(s)
                assert np.allclose(faces_id, np.arange(fa_count))
                log.info("end of face integrity")
                sys.stdout.flush()
                del faces_id
                # For cells
                log.info("  Checking cells integrity ..")
                sys.stdout.flush()
                h = restartSectionHeader()
                h.readReqVar(self.fid, self.byte_swap, ["UGP_IO_CV_CHECK"])
                assert h.idata[0] == cv_count
                assert h.id[0] == ic3_restart_codes["UGP_IO_CV_CHECK"]
                s = BinaryRead(self.fid, "%di" % cv_count, self.byte_swap, type2nbytes["int32"] * cv_count)
                cells_id = np.asarray(s)
                assert np.allclose(cells_id, np.arange(cv_count))
                log.info("end of cell integrity")
                sys.stdout.flush()
                del cells_id
            else:
                log.warning("check integrity: no check, UGP_IO_NO_CHECK missing")

        # The two connectivities now
        #
        # - First, NOOFA
        #
        log.info("  Parsing face to node connectivity...")
        sys.stdout.flush()
        h = restartSectionHeader()
        h.readReqVar(self.fid, self.byte_swap, ["UGP_IO_NOOFA_I_AND_V"])
        assert h.idata[0] == fa_count
        assert h.idata[1] == noofa_count

        # Get the node count per face
        s = BinaryRead(self.fid, "%di" % fa_count, self.byte_swap, type2nbytes["int32"] * fa_count)
        nno_per_face = np.asarray(s)
        uniq, counts = np.unique(nno_per_face, return_counts=True)
        for facesize, nfacesize in zip(uniq, counts):
            log.info(f"  {nfacesize} faces of {facesize} nodes")
        # Initialize the proper connectivity arrays in self.mesh
        face2node_index = np.concatenate(([0], np.cumsum(nno_per_face)))
        # np.sum(nno_per_face) == noofa_count
        assert face2node_index[0] == 0
        assert face2node_index[-2] == noofa_count - nno_per_face[-1]
        # Now loop on the restart file to fill the connectivities
        s = BinaryRead(self.fid, "%di" % noofa_count, self.byte_swap, type2nbytes["int32"] * noofa_count)
        # store in 8 bytes
        face2node_value = np.asarray(s).astype(np.int64)
        zface2node = _conn.compressed_listofindex(face2node_index, face2node_value)
        zface2node.check()
        self.mesh["connectivity"]["noofa"] = zface2node
        log.info("end of face/vertex connectivity")
        sys.stdout.flush()
        del nno_per_face, uniq, counts
        #
        # - Second, CVOFA
        #
        log.info("  Parsing face to cell connectivity...")
        sys.stdout.flush()
        h = restartSectionHeader()
        h.readReqVar(self.fid, self.byte_swap, ["UGP_IO_CVOFA"])

        assert h.idata[0] == fa_count
        assert h.idata[1] == 2
        # Fill the connectivities
        s = BinaryRead(self.fid, "%di" % (fa_count * 2), self.byte_swap, type2nbytes["int32"] * fa_count * 2)
        # store in 8 bytes
        self.mesh["connectivity"]["cvofa"]["cvofa"] = np.asarray(s).astype(np.int64).reshape((fa_count, 2))

        log.info("end of face/cell connectivity")
        sys.stdout.flush()

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
        log.info("  Parsing boundary conditions...")
        sys.stdout.flush()
        self.mesh['bocos'] = []  # init list of bocos
        while True:
            h = restartSectionHeader()
            if not h.readVar(self.fid, self.byte_swap, ["UGP_IO_FA_ZONE"], reset_offset=False):
                break
            # UGP_IO_FA_ZONE header
            # 3 ints: kind, face range (begin, start)
            self.mesh["params"]["nboco"] += 1
            boco = _mesh.submeshmark(h.name)
            boco.type = zonekind2type[h.idata[0]]
            boco.geodim = 'intface' if boco.type == 'internal' else 'bdface'
            boco.index = _conn.indexlist(irange=[h.idata[1], h.idata[2]])
            boco.properties["periodic_transform"] = h.rdata
            #
            famin, famax = boco.index.range()
            sta = face2node_index[famin]
            try:
                sto = face2node_index[famax + 1]
            except IndexError:
                sto = face2node_value.size
            #
            boco.properties["slicing"] = np.unique(face2node_value[sta:sto])
            self.mesh['bocos'].append(boco)
            if h.idata[0] == 6:
                break
        log.info("end of boundary conditions")
        sys.stdout.flush()

        # Parse the header of the partition information
        log.info("  Parsing partitioning information...")
        sys.stdout.flush()
        h = restartSectionHeader()
        h.readReqVar(self.fid, self.byte_swap, ["UGP_IO_CV_PART"])

        self.mesh["partition"] = {}
        self.mesh["partition"]['npart'] = h.idata[1]
        self.mesh["partition"]['icvpart'] = np.zeros((cv_count,), dtype=np.int32)
        s = BinaryRead(self.fid, "%di" % cv_count, self.byte_swap, type2nbytes["int32"] * cv_count)
        self.mesh["partition"]['icvpart'] = np.asarray(s)
        log.info("end of partition")
        sys.stdout.flush()

        # The coordinates of the vertices finally
        log.info("  Parsing vertices coordinates...")
        sys.stdout.flush()
        h = restartSectionHeader()
        h.readReqVar(self.fid, self.byte_swap, ["UGP_IO_X_NO"])

        assert h.idata[0] == no_count
        assert h.idata[1] == 3
        s = BinaryRead(
            self.fid, "%dd" % (no_count * 3), self.byte_swap, type2nbytes["float64"] * no_count * 3
        )
        self.mesh["coordinates"] = np.ascontiguousarray(np.asarray(s).reshape(no_count, 3))
        log.info("end of node coordinates")
        sys.stdout.flush()

    def _ReadInformativeValues(self):
        """
        Method reading all the values also stored in a restart file,
        i.e. the step number, the time, the timestep.
        input:  handle on an open restart file, [type file identifier]
                endianness flag [boolean]
        output: simulation state structure containing informations about the current state
                of the simulation
        """

        # Initialize the state dictionary
        self.simulation_state = {"step": 0, "dt": 0, "time": 0, "wgt": {}}

        # First, a header saying data to introduce to this block we are parsing now

        # removed as not used at the moment
        # h = restartSectionHeader()
        # if(not readVar(self.fid, self.byte_swap,"UGP_IO_DATA")
        # if (not varfound): raise ("UGP_IO_DATA not found")

        reset_offset = True
        while True:
            h = restartSectionHeader()
            if not h.readVar(self.fid, self.byte_swap, ["UGP_IO_I0", "UGP_IO_D0"], reset_offset=reset_offset):
                break
            reset_offset = False
            if h.id[0] == ic3_restart_codes["UGP_IO_I0"]:
                self.simulation_state[h.name.lower()] = h.idata[0]
            elif h.id[0] == ic3_restart_codes["UGP_IO_D0"]:
                self.simulation_state[h.name.lower()] = h.rdata[0]

    # def get_datacell_properties(self):
    #     return self.variables["cells"]["_info"]

    # def get_ndof(self):
    #     return self.get_datacell_properties["ndof"]

    def _set_ndof_properties(self, intndof: int):
        ndof = 1
        if self.ic3_version < 0:
            raise ValueError("unknown IC3 version number")
        elif self.ic3_version < 3:
            pass
            # ignore ints
            # if intndof != 0: # 0 is expected value for version 1 and 2
            #     log.warning("unexpected non zero value for ndof in IC3 version 1 and 2")
        else:  # version >= 3
            if intndof == 0:  # 0 is NOT expected
                log.warning(
                    "unexpected zero value for ndof in IC3 version 3; set to 1 !",
                )
                ndof = 1
            else:
                ndof = intndof
        # if self.get_datacell_properties().get("ndof", ndof) != ndof:
        #     raise ValueError("Inconsistent cell data size (ndof)")
        # self.get_datacell_properties()["ndof"] = ndof
        return ndof

    def _ReadRestartVar(self):
        '''
        Method reading the variables from the restart file
        input:  handle on an open restart file, [type file identifier]
                endianness flag [boolean]
                mesh structure
        output: structure containing all the variables
        '''

        # Some extra modules
        import copy

        # Initialize the variable dictionary
        self.variables = {"nodes": {}, "cells": {}, "faces": {}}
        # set the local shortcut names
        no_count = self.mesh["params"]["no_count"]
        fa_count = self.mesh["params"]["fa_count"]
        cv_count = self.mesh["params"]["cv_count"]

        # First come the scalars
        log.info("  First the scalars...")
        reset_offset = True
        while True:
            h = restartSectionHeader()
            if not h.readVar(
                self.fid,
                self.byte_swap,
                [
                    "UGP_IO_NO_D1",
                    "UGP_IO_NO_II1",
                    "UGP_IO_FA_D1",
                    "UGP_IO_CV_D1",
                    "UGP_IO_CV_II1",
                ],
                reset_offset=reset_offset,
            ):
                break
            reset_offset = False
            #
            typechar = properties_ugpcode[h.id[0]]['structcode']
            typesize = properties_ugpcode[h.id[0]]['size']
            nptype = properties_ugpcode[h.id[0]]['numpytype']
            #
            if h.idata[0] == no_count:
                s = BinaryRead(self.fid, "%d" % no_count + typechar, self.byte_swap, typesize * no_count)
                scalar = self.variables["nodes"][h.name] = np.asarray(s).astype(nptype)
                self.nodedata.add_data(h.name, scalar)
            elif h.idata[0] == fa_count:
                s = BinaryRead(self.fid, "%d" % fa_count + typechar, self.byte_swap, typesize * fa_count)
                scalar = self.variables["faces"][h.name] = np.asarray(s).astype(nptype)
            elif h.idata[0] == cv_count:
                ndof = self._set_ndof_properties(h.idata[1])
                if ndof > 1:
                    self.celldata.Xrep = 'spectralcell'
                    self.celldata.ndof = ndof
                log.info(
                    "cell variable section of size ncells * ndofs {}x{}".format(cv_count, ndof),
                )
                s = BinaryRead(
                    self.fid, "%d" % (ndof * cv_count) + typechar, self.byte_swap, typesize * ndof * cv_count
                )
                scalar = np.asarray(s).astype(nptype)
                # If multiple connectivities, order the tables correctly
                if self.mesh["connectivity"]["nkeys"] > 1:
                    aux = copy.deepcopy(scalar)
                    scalar = np.empty((0,), dtype=aux.dtype)
                    for uns_type, indices in self.mesh["connectivity"]["cell_indices"].items():
                        scalar = np.concatenate((scalar, aux[indices]))
                self.celldata.add_data(h.name, scalar)
            else:
                api.error_stop(f"Fatal error. Incoherence in dataset {h.name}. Exiting.")
            print_stats_scalar(h.name, scalar)
        log.info("  end of scalars")

        # Then the vectors
        log.info("  Then the vectors...")
        reset_offset = True
        while True:
            h = restartSectionHeader()
            if not h.readVar(
                self.fid,
                self.byte_swap,
                [
                    "UGP_IO_NO_D3",
                    "UGP_IO_FA_D3",
                    "UGP_IO_CV_D3",
                ],
                reset_offset=reset_offset,
            ):
                break
            reset_offset = False
            #
            if h.idata[0] == no_count:
                s = BinaryRead(
                    self.fid, "%dd" % (no_count * 3), self.byte_swap, type2nbytes["float64"] * no_count * 3
                )
                vector = self.variables["nodes"][h.name] = np.asarray(s).reshape((no_count, 3))
            elif h.idata[0] == fa_count:
                self.variables["faces"][h.name] = np.zeros((fa_count, 3), dtype=np.float64)
                s = BinaryRead(
                    self.fid, "%dd" % (fa_count * 3), self.byte_swap, type2nbytes["float64"] * fa_count * 3
                )
                vector = self.variables["faces"][h.name] = np.asarray(s).reshape((fa_count, 3))
            elif h.idata[0] == cv_count:
                ndof = self._set_ndof_properties(h.idata[1])
                if ndof > 1:
                    self.celldata.Xrep = 'spectralcell'
                    self.celldata.ndof = ndof
                s = BinaryRead(
                    self.fid,
                    "%dd" % (ndof * cv_count * 3),
                    self.byte_swap,
                    type2nbytes["float64"] * ndof * cv_count * 3,
                )
                vector = np.asarray(s).reshape((ndof * cv_count, 3))
                # If multiple connectivities, gotta order the tables correctly
                if self.mesh["connectivity"]["nkeys"] > 1:
                    aux = copy.deepcopy(vector)
                    vector = np.empty((0, 3), dtype=aux.dtype)
                    for indices in self.mesh["connectivity"]["cell_indices"].values():
                        vector = np.concatenate((vector, aux[indices, :]), axis=0)
                self.celldata.add_data(h.name, vector)
            else:
                api.error_stop(f"Fatal error. Incoherence in dataset {h.name}. Exiting.")
            print_stats_vector(h.name, vector)
        log.info("  end of vectors")

        # Then the tensors
        log.info("  Then the tensors...")
        reset_offset = True
        while True:
            h = restartSectionHeader()
            if not h.readVar(self.fid, self.byte_swap, ["UGP_IO_CV_D33"], reset_offset=reset_offset):
                break
            reset_offset = False

            if h.idata[0] == cv_count:
                ndof = self._set_ndof_properties(h.idata[1])
                if ndof > 1:
                    self.celldata.Xrep = 'spectralcell'
                    self.celldata.ndof = ndof
                s = BinaryRead(
                    self.fid,
                    "%dd" % (ndof * cv_count * 3 * 3),
                    self.byte_swap,
                    type2nbytes["float64"] * ndof * cv_count * 3 * 3,
                )
                tensor = np.asarray(s).reshape((ndof * cv_count, 3, 3))
                # If multiple connectivities, gotta order the tables correctly
                if self.mesh["connectivity"]["nkeys"] > 1:
                    aux = copy.deepcopy(tensor)
                    tensor = np.empty((0, 3, 3), dtype=aux.dtype)
                    for indices in self.mesh["connectivity"]["cell_indices"].values():
                        tensor = np.concatenate((tensor, aux[indices, :, :]), axis=0)
                self.celldata.add_data(h.name, tensor)
                print_stats_scalar(h.name, tensor)
            else:
                api.error_stop(f"Fatal error. Incoherence in dataset {h.name}. Exiting.")
        log.info("  end of tensors")

    # def __reachedEOF(self):
    #     '''
    #     Method to check whether the reader reached EOF.
    #     It should happen right after reading all the variable blocks.
    #     input   : handle on an open restart file, [type file identifier]
    #               endianness flag [boolean]
    #     output  : check on end-of-file [boolean]
    #     '''

    #     # Read the next header
    #     h = restartSectionHeader()
    #     h.read(self.fid, self.byte_swap)
    #     if h.name != "EOF" or\
    #        h.id != ic3_restart_codes["UGP_IO_EOF"]:
    #        return False

    #     return True


###################################################################################################

# if __name__ == "__main__":
#     '''
#     The script is supposed to be used with command line arguments
#     but if it is not, it runs a test on a pre-defined file name.
#     '''

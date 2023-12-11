# Import modules
import logging

import numpy as np

import cfdtools.api as api

from cfdtools.ic3._ic3 import (
    ic3_restart_codes,
    type2nbytes,
    type2zonekind,
    struct_endian,
    BinaryWrite,
    restartSectionHeader,
)

log = logging.getLogger(__name__)

###################################################################################################


@api.fileformat_writer('IC3', '.ic3')
class writer:
    '''Implementation of the writer to write ic3 restart files'''

    __version__ = "2"

    def __init__(self, mesh, endian='native'):
        """
        Initialization of a ic3 restart file writer.
        """
        log.info("> Initialization of IC3 writer V" + self.__version__)
        if not endian in struct_endian.keys():
            raise ValueError("unknown endian key")
        else:
            self.endian = endian
        # self.vars = {"nodes": {}, "cells": {}}
        self._mesh = mesh
        self.params = {}
        # Initialize the simulation state
        self.set_simstate()
        self.set_mesh()
        self.set_bocos()
        self.set_vars()

    def set_mesh(self):
        """
        Store the coordinates of the nodes and the
        element-to-vertex connectivity. Also triggers the creation
        of a face-to-element and face-to-vertex connectivity
        as per CharlesX implementation.
        """
        timer = api.Timer()
        log.info("> Check connectivity and compute mandatory")
        if not self._mesh._faces:  # empty dict of faces
            with api.Timer(task="  generate all faces"):
                self._mesh.make_face_connectivity()
        if any([boco.nodebased() for _, boco in self._mesh._bocos.items()]):
            with api.Timer(task="  change boco marks (node to face)"):
                self._mesh.bocomarks_set_node_to_face()
        if 'boundary' in self._mesh._faces.keys():
            if self._mesh.make_unmarked_BC(name="unmarked"):
                log.info("  create a specific boco mark for unmarked faces: (unmarked)")
        # self._mesh.printinfo()
        if any([boco.index.type == 'list' for _, boco in self._mesh._bocos.items()]):
            with api.Timer(task="  reindex boundary faces with boco marks and compress"):
                self._mesh.reindex_boundaryfaces()

        log.info("Setting coordinates and connectivity arrays")

        # Nodes
        self.coordinates = np.stack([self._mesh._nodes[c] for c in 'xyz'], axis=1)

        # Compute the number of nodes and elements
        assert self.coordinates.shape[0] == self._mesh.nnode
        self.params["no_count"] = self.coordinates.shape[0]
        # Elements, count throughout all connectivities
        self.params["cv_count"] = self._mesh.ncell

        ### TODO : here most params are directly copied from params after IC3 reading
        ### will need an actual extract of data
        #
        # Store the information properly
        self.params["no_count"] = self._mesh.nnode
        self.params["fa_count"] = self._mesh.nface
        self.params["cv_count"] = self._mesh.ncell

        # check face connectivity
        if 'mixed' in self._mesh._faces.keys():
            zface2node = self._mesh._faces['mixed']['face2node'].exportto_compressedindex()
            self.f2e = self._mesh._faces['mixed']['face2cell'].conn
        else:
            with api.Timer(task="  compressing faces connectivity"):
                mixedfaces, f2cell = self._mesh.export_mixedfaces()
                zface2node = mixedfaces.exportto_compressedindex()
            self.f2e = f2cell.conn
        self.f2v = {}
        self.f2v["noofa"] = zface2node._index[1:] - zface2node._index[:-1]
        self.f2v["noofa_v"] = zface2node._value
        self.params["noofa_count"] = zface2node._value.size

        if not 'partition' in self._mesh._cellprop.keys():
            self.params["partition"] = {}
            self.params["partition"]['npart'] = 1
            self.params["partition"]['icvpart'] = np.zeros((self._mesh.ncell,), dtype=np.int32)
        else:
            self.params['partition'] = self._mesh._cellprop['partition']

    def set_bocos(self, nboco=None):
        """
        Set the boundary conditions dictionary for the restart.
        It will be later written to the file using the
        __WriteRestartConnectivity method.
        The bocos slicing is expressed in terms of faces.
        """
        log.info("Setting boundary conditions")
        self.bocos = {key: boco for key, boco in self._mesh._bocos.items() if not boco.type == 'internal'}
        # self.bocos.pop("nfa_b")
        # self.bocos.pop("nfa_bp")

    def set_simstate(self, state={}):
        """
        Set the simulation state for the restart.
        It will be later written to the file using the
        __WriterInformativeValues method.
        It gives informations to CX about the current iteration number,
        the timestep, and the current time.
        """

        # Initialize the current simulation state
        self.simstate = {"wgt": {}, "step": 0, "dt": 0, "time": 0}
        # update with _mesh._params then state if key exists
        for refdict in (self._mesh._params, state):
            for key in self.simstate.keys():
                if key in refdict.keys():
                    self.simstate[key] = refdict[key]

    def set_vars(self):
        """
        Set the variables for the restart.
        They will be later written to the file using the
        __WriteRestartVar method.
        """
        log.info("Setting variables")
        self.vars = {
            "nodes": self._mesh._nodedata,
            "cells": self._mesh._celldata,
        }
        # Start with the variables stored at the vertices
        # for key, item in self._mesh._nodedata.items():
        #     log.info("  node data: " + key)
        #     self.vars["nodes"][key] = item

        # # Then the variables stored at the cells:
        # for key, item in self._mesh._celldata.items():
        #     log.info("  cell data: " + key)
        #     self.vars["cells"][key] = item

    def write_data(self, filename):
        """
        Main method of the ic3 restart file writer
        """
        log.info(f"> WRITING FILE {filename!r}")
        self.filename = filename
        # Open the file for binary reading
        with open(self.filename, "wb") as self.fid:
            #
            log.info("> check consistency before writing")
            if not self.check():
                raise RuntimeError("Inconsistent data to write")
            #
            log.info("> Writing header")
            self.__WriteRestartHeader()
            #
            log.info("> Writing connectivity")
            self.__WriteRestartConnectivity()
            #
            log.info("> Writing informative values")
            self.__WriteInformativeValues()
            #
            log.info("> Writing variables")
            self.__WriteRestartVar()
            #
            log.info("> End of file")
            header = restartSectionHeader()
            header.name = "EOF"
            header.id = ic3_restart_codes["UGP_IO_EOF"]
            header.skip = header.hsize
            header.write(self.fid, self.endian)
            #
            # self.fid is closed
        del self.fid

    def check(self):
        check_error = True
        # bad field size
        keylist = []
        # for key, cellitem in self.vars["cells"].items():
        #     if cellitem.ndof() != 1:
        #         keylist.append(key)
        # if len(keylist) >= 1:
        #     check_error = False
        #     log.error('wrong size (ndof) of cell data: ' + keylist.__str__()
        #     )
        return check_error

    def __WriteRestartHeader(self):
        """
        Method writing the header of a restart file.
        It is composed of two integers, the "magic number" used as a flag for endianness
        and the CharlesX version number.
        input   : handle on an open restart file, [type file identifier]
        """
        # Write the two integers
        BinaryWrite(self.fid, self.endian, "ii", [ic3_restart_codes["UGP_IO_MAGIC_NUMBER"], 2])

    def __WriteRestartConnectivity_check(self):
        # Node check
        # Header
        header = restartSectionHeader()
        header.name = "NO_CHECK"
        header.id = ic3_restart_codes["UGP_IO_NO_CHECK"]
        header.skip = header.hsize + type2nbytes["int32"] * self.params["no_count"]
        header.idata[0] = self.params["no_count"]
        header.write(self.fid, self.endian)
        # Write the nodes global ids
        BinaryWrite(
            self.fid,
            self.endian,
            "i" * self.params["no_count"],
            range(self.params["no_count"]),
        )
        # Face check
        # Header
        header = restartSectionHeader()
        header.name = "FA_CHECK"
        header.id = ic3_restart_codes["UGP_IO_FA_CHECK"]
        header.skip = header.hsize + type2nbytes["int32"] * self.params["fa_count"]
        header.idata[0] = self.params["fa_count"]
        header.write(self.fid, self.endian)
        # Write the nodes global ids
        BinaryWrite(
            self.fid,
            self.endian,
            "i" * self.params["fa_count"],
            range(self.params["fa_count"]),
        )
        # Cell check
        # Header
        header = restartSectionHeader()
        header.name = "CV_CHECK"
        header.id = ic3_restart_codes["UGP_IO_CV_CHECK"]
        header.skip = header.hsize + type2nbytes["int32"] * self.params["cv_count"]
        header.idata[0] = self.params["cv_count"]
        header.write(self.fid, self.endian)
        # Write the nodes global ids
        BinaryWrite(
            self.fid,
            self.endian,
            "i" * self.params["cv_count"],
            range(self.params["cv_count"]),
        )

    def __WriteRestartConnectivity(self):
        """
        Method writing the connectivities of a restart file.
        Also the number of nodes, faces and volumes for later checks.
        """
        # First the header for counts
        nnode, nface, ncell = (self.params[key] for key in ('no_count', 'fa_count', 'cv_count'))
        log.info(
            f"  sizes: {nnode} nodes, {nface} faces and reference to {ncell} cells",
        )
        header = restartSectionHeader()
        header.name = "NO_FA_CV_NOOFA_COUNTS"
        header.id = ic3_restart_codes["UGP_IO_NO_FA_CV_NOOFA_COUNTS"]
        header.skip = header.hsize
        header.idata[0] = nnode
        header.idata[1] = nface
        header.idata[2] = ncell
        header.idata[3] = self.params["noofa_count"]
        header.write(self.fid, self.endian)
        #
        # don't really know how useful it was; suppressed in V3
        self.__WriteRestartConnectivity_check()
        #
        # Connectivities
        # Faces-to-nodes connectivities
        # Header
        log.info("  face to node connectivity")
        header = restartSectionHeader()
        header.name = "NOOFA_I_AND_V"
        header.id = ic3_restart_codes["UGP_IO_NOOFA_I_AND_V"]
        header.skip = (
            header.hsize + type2nbytes["int32"] * nface + type2nbytes["int32"] * self.params["noofa_count"]
        )
        header.idata[0] = nface
        header.idata[1] = self.params["noofa_count"]
        header.write(self.fid, self.endian)
        # Node count per face
        BinaryWrite(self.fid, self.endian, "i" * nface, self.f2v["noofa"].tolist())  # remove first 0
        # Flattened face-to-node connectivity
        BinaryWrite(
            self.fid,
            self.endian,
            "i" * self.params["noofa_count"],
            self.f2v["noofa_v"].tolist(),
        )
        # print('noofa', self.f2v["noofa"] )
        # print('noofa_v', self.f2v["noofa_v"] )
        # Faces-to-cells connectivities
        # Header
        log.info("  face to cell connectivity")
        header = restartSectionHeader()
        header.name = "CVOFA"
        header.id = ic3_restart_codes["UGP_IO_CVOFA"]
        header.skip = header.hsize + type2nbytes["int32"] * nface * 2
        header.idata[0] = nface
        header.idata[1] = 2
        header.write(self.fid, self.endian)
        # Flattened face-to-cell connectivity
        # print("W", self.f2e)
        # print(self.f2e.ravel())
        BinaryWrite(self.fid, self.endian, "i" * nface * 2, self.f2e.ravel().tolist())
        # print('cellofface', self.f2e.ravel().tolist() )
        # Face zones, a.k.a boundary condition patches
        log.info("  marks (face based)")
        last_boco = 0
        for key, boco in self.bocos.items():
            assert key == boco.name
            assert boco.geodim in ('face', 'bdface'), "boco marks must be faces index"
            # Header
            header = restartSectionHeader()
            header.name = key
            header.id = ic3_restart_codes["UGP_IO_FA_ZONE"]
            header.skip = header.hsize
            # diff# print(self.bocos[key]["type"], type2zonekind)
            assert boco.type in type2zonekind.keys(), f"unsupported type of boco for IC3 output: {boco.type}"
            ifmin, ifmax = boco.index.range()
            log.info(f"  . ({boco.type}) {boco.name}: {ifmin}-{ifmax}")
            header.idata[0] = type2zonekind[boco.type]
            assert boco.index.type == 'range', "indexing must be a range and may need reordering"
            header.idata[1] = ifmin
            header.idata[2] = ifmax
            last_boco = max(last_boco, ifmax)
            if "periodic_transform" in boco.properties.keys():
                if boco.properties["periodic_transform"] is not None:
                    for idx, val in enumerate(boco.properties["periodic_transform"]):
                        header.rdata[idx] = val
            header.write(self.fid, self.endian)
        # Header
        header = restartSectionHeader()
        header.name = 'internal-domain'
        header.id = ic3_restart_codes["UGP_IO_FA_ZONE"]
        header.skip = header.hsize
        ifmin, ifmax = last_boco + 1, self.params['fa_count'] - 1  # last face
        log.info(f"  additional mark (FA_ZONE) for internal faces: {ifmin}-{ifmax}")
        header.idata[0] = type2zonekind['internal']
        header.idata[1] = ifmin
        header.idata[2] = ifmax
        header.rdata[:] = 0.0
        header.write(self.fid, self.endian)  # Partition information
        # Header
        log.info("  cell based partition info")
        header = restartSectionHeader()
        header.name = "CV_PART"
        header.id = ic3_restart_codes["UGP_IO_CV_PART"]
        header.skip = header.hsize + type2nbytes["int32"] * ncell
        header.idata[0] = ncell
        header.idata[1] = self.params['partition'].get('npart', 1)
        header.write(self.fid, self.endian)
        # The ranks of the processors, default to everybody 0
        BinaryWrite(
            self.fid,
            self.endian,
            "i" * ncell,
            self.params['partition'].get('icvpart', np.zeros((ncell,), dtype=np.int32)),
        )
        # Coordinates
        # Header
        log.info("  nodes coordinates")
        header = restartSectionHeader()
        header.name = "X_NO"
        header.id = ic3_restart_codes["UGP_IO_X_NO"]
        header.skip = header.hsize + type2nbytes["float64"] * nnode * 3
        header.idata[0] = nnode
        header.idata[1] = 3
        header.write(self.fid, self.endian)
        # X, Y and Z
        BinaryWrite(self.fid, self.endian, "d" * nnode * 3, self.coordinates.ravel(order='C'))

    def __WriteInformativeValues(self):
        """
        Method writing the solution information to the restart file
        """

        # Start of data
        # An extra header to announce what is coming
        header = restartSectionHeader()
        header.name = "DATA"
        header.id = ic3_restart_codes["UGP_IO_DATA"]
        header.skip = header.hsize
        header.write(self.fid, self.endian)
        # Write current number of iteration
        # Header
        key = "STEP"
        header = restartSectionHeader()
        header.name = key
        value = self.simstate[key.lower()]
        log.info(f"  write section for {key} parameter: {value}")
        header.id = ic3_restart_codes["UGP_IO_I0"]
        header.skip = header.hsize
        header.idata[0] = value
        header.write(self.fid, self.endian)
        # Write the rest, all double values
        # The monomials first
        for key in ["DT", "TIME"]:
            value = self.simstate[key.lower()]
            log.info(f"  write section for {key} parameter: {value}")
            header = restartSectionHeader()
            header.name = key
            header.id = ic3_restart_codes["UGP_IO_D0"]
            header.skip = header.hsize
            header.rdata[0] = value
            header.write(self.fid, self.endian)
        # The second level now
        for key in self.simstate["wgt"].keys():
            header = restartSectionHeader()
            header.name = key.upper() + '_WGT'
            header.id = ic3_restart_codes["UGP_IO_D0"]
            header.skip = header.hsize
            header.rdata[0] = self.simstate["wgt"][key]
            header.write(self.fid, self.endian)

    def __WriteRestartVar(self):
        """
        Method to write all the variables into a restart file.
        Scalars, vectors and tensors all together.
        """
        # Start with the node based variables
        #
        nno = self.params["no_count"]
        #
        for ndname, nddata in self.vars["nodes"].items():
            # Scalar
            if nddata.size == nddata.shape[0]:
                # Header
                log.info(f"  write node scalar data {ndname}")
                header = restartSectionHeader()
                header.name = ndname
                header.id = ic3_restart_codes["UGP_IO_NO_D1"]
                header.skip = header.hsize + type2nbytes["float64"] * nno
                header.idata[0] = nno
                header.write(self.fid, self.endian)
                # Field
                BinaryWrite(
                    self.fid,
                    self.endian,
                    "d" * nno,
                    nddata,
                )
        #
        for ndname, nddata in self.vars["nodes"].items():
            # Vector
            if len(nddata.shape) == 2:
                log.info(f"  write node vector data {ndname}")
                # Header
                header = restartSectionHeader()
                header.name = ndname
                header.id = ic3_restart_codes["UGP_IO_NO_D3"]
                header.skip = header.hsize + type2nbytes["float64"] * nno * 3
                header.idata[0] = nno
                header.idata[1] = 3
                header.write(self.fid, self.endian)
                # Field
                BinaryWrite(
                    self.fid,
                    self.endian,
                    "d" * nno * 3,
                    nddata.ravel(order='C'),
                )
        #
        for ndname, nddata in self.vars["nodes"].items():
            # Tensor
            if len(nddata.shape) == 3:
                pass

        # Then the cell based variables
        #
        if self.vars["cells"]:  # if defined
            ncv = self.params["cv_count"]
        #
        for cvname, cvdata in self.vars["cells"].items():
            # Scalar
            if cvdata.size == cvdata.shape[0]:
                # Header
                log.info(f"  write cell scalar data {cvname}")
                header = restartSectionHeader()
                header.name = cvname
                header.id = ic3_restart_codes["UGP_IO_CV_D1"]
                header.skip = header.hsize + type2nbytes["float64"] * ncv
                header.idata[0] = ncv
                header.write(self.fid, self.endian)
                # Field
                BinaryWrite(
                    self.fid,
                    self.endian,
                    "d" * ncv,
                    cvdata,
                )
        #
        for cvname, cvdata in self.vars["cells"].items():
            # Vector
            if len(cvdata.shape) == 2:
                # Header
                log.info(f"  write cell vector data {cvname}")
                header = restartSectionHeader()
                header.name = cvname
                header.id = ic3_restart_codes["UGP_IO_CV_D3"]
                header.skip = header.hsize + type2nbytes["float64"] * ncv * 3
                header.idata[0] = ncv
                header.idata[1] = 3
                header.write(self.fid, self.endian)
                # Field
                BinaryWrite(
                    self.fid,
                    self.endian,
                    "d" * ncv * 3,
                    cvdata.ravel(order='C'),
                )
        #
        for cvname, cvdata in self.vars["cells"].items():
            # Tensor
            if len(cvdata.shape) == 3:
                # Header
                log.info(f"  write cell tensor data {cvname}")
                header = restartSectionHeader()
                header.name = cvname
                header.id = ic3_restart_codes["UGP_IO_CV_D33"]
                header.skip = header.hsize + type2nbytes["float64"] * ncv * 9
                header.idata[0] = ncv
                header.idata[1] = 3
                header.idata[2] = 3
                header.write(self.fid, self.endian)
                # Field
                BinaryWrite(
                    self.fid,
                    self.endian,
                    "d" * ncv * 9,
                    cvdata.ravel(order='C'),
                )

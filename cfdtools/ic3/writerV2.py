
# Import modules
import numpy as np
import cfdtools.api as api
#import cfdtools.meshbase._mesh as _mesh
from cfdtools.ic3._ic3 import *

###################################################################################################

@api.fileformat_writer('IC3', '.ic3')
class writer():
    ''' Implementation of the writer to write ic3 restart files '''
    __version__="2"

    def __init__(self, mesh, endian='native'):
        """
        Initialization of a ic3 restart file writer.
        """
        api.io.print('std',"Initialization of IC3 writer V"+self.__version__)
        if not endian in struct_endian.keys():
            raise ValueError("unknown endian key")
        else:
            self.endian = endian
        self.vars = {"nodes":{}, "cells":{}}
        self._mesh = mesh
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
        api.io.print('std',"Setting coordinates and connectivity arrays..")

        # Keep the coordinates of the points
        self.coordinates = np.stack((self._mesh._nodes[c] 
            for c in ['x', 'y', 'z']), axis=1)

        # Compute the number of nodes and elements
        self.params = {}
        # Nodes
        #print(self.coordinates.shape)
        assert self.coordinates.shape[0] == self._mesh.nnode
        self.params["no_count"] = self.coordinates.shape[0]
        # Elements, count throughout all connectivities
        self.params["cv_count"] = self._mesh.ncell

        ### TODO : here most params are directly copied from params after IC3 reading
        ### will need an actual extract of data
        #
        # Store the information properly
        self.params["no_count"] = self._mesh._params["no_count"]
        self.params["fa_count"] = self._mesh._params["fa_count"]
        self.params["cv_count"] = self._mesh._params["cv_count"]
        self.params["noofa_count"] = self._mesh._params["noofa_count"]

        self.f2v = {}
        self.f2v["noofa"] = self._mesh._face2node['listofStarts_f2v'][1:]-self._mesh._face2node['listofStarts_f2v'][:-1]
        self.f2v["noofa_v"] = self._mesh._face2node['face2vertex']
        self.f2e = self._mesh._face2cell['cvofa']
        #print(self.f2e.shape)
        api.io.print('std',"ok.")

    def set_bocos(self, nboco=None):
        """
        Set the boundary conditions dictionary for the restart.
        It will be later written to the file using the
        __WriteRestartConnectivity method.
        The bocos slicing is expressed in terms of faces.
        """
        api.io.print('std',"Setting boundary conditions..")
        self.bocos = self._mesh._bocos
        self.bocos.pop("nfa_b")
        self.bocos.pop("nfa_bp")
        api.io.print('std',"ok.")

    def set_simstate(self, state={}):
        """
        Set the simulation state for the restart.
        It will be later written to the file using the
        __WriterInformativeValues method.
        It gives informations to CX about the current iteration number,
        the timestep, and the current time.
        """

        # Initialize the current simulation state
        self.simstate = {"wgt":{}, "step":0, "dt":0, "time":0}
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
        api.io.print('std',"Setting variables..",)
        self.vars = {"nodes":{}, "cells":{}}
        # Start with the variables stored at the vertices
        for key, item in self._mesh._nodedata.items():
            api.io.print('std',"  node data: "+key)
            self.vars["nodes"][key] = item

        # Then the variables stored at the cells:
        for key, item in self._mesh._celldata.items():
            api.io.print('std',"  cell data: "+key)
            self.vars["cells"][key] = item
        api.io.print('std',"ok.")
        
    def write_data(self, filename):
        """
        Main method of the ic3 restart file writer
        """
        self.filename = filename
        # Open the file for binary reading
        self.fid = open(self.filename, "wb")

        api.io.print('std',"> check consistency before writing")
        if not self.check():
            raise RuntimeError("Inconsistent data to write")

        api.io.print('std',"> Writing header")
        self.__WriteRestartHeader()
        #
        api.io.print('std',"> writing connectivity")
        self.__WriteRestartConnectivity()
        #
        api.io.print('std',"> writing informative values")
        self.__WriteInformativeValues()
        #
        api.io.print('std',"> writing variables")
        self.__WriteRestartVar()
        #
        api.io.print('std',"> end of file")
        header = restartSectionHeader()
        header.name = "EOF"
        header.id = ic3_restart_codes["UGP_IO_EOF"]
        header.skip = header.hsize
        header.write(self.fid, self.endian)

        # Before returning, close the file
        self.fid.close()
        del self.fid

    def check(self):
        check_error = True
        # bad field size
        keylist = []
        for key, cellitem in self.vars["cells"].items():
            if cellitem.ndof() != 1:
                keylist.append(key)
        if len(keylist) >= 1:
            check_error = False
            api.io.print('error', 'wrong size (ndof) of cell data: '+keylist.__str__())
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
        BinaryWrite(self.fid, self.endian, "i"*self.params["no_count"], range(self.params["no_count"]))
        # Face check
        # Header
        header = restartSectionHeader()
        header.name = "FA_CHECK"
        header.id = ic3_restart_codes["UGP_IO_FA_CHECK"]
        header.skip = header.hsize + type2nbytes["int32"] * self.params["fa_count"]
        header.idata[0] = self.params["fa_count"]
        header.write(self.fid, self.endian)
        # Write the nodes global ids
        BinaryWrite(self.fid, self.endian, "i"*self.params["fa_count"], range(self.params["fa_count"]))
        # Cell check
        # Header
        header = restartSectionHeader()
        header.name = "CV_CHECK"
        header.id = ic3_restart_codes["UGP_IO_CV_CHECK"]
        header.skip = header.hsize + type2nbytes["int32"] * self.params["cv_count"]
        header.idata[0] = self.params["cv_count"]
        header.write(self.fid, self.endian)
        # Write the nodes global ids
        BinaryWrite(self.fid, self.endian, "i"*self.params["cv_count"], range(self.params["cv_count"]))

    def __WriteRestartConnectivity(self):
        """
        Method writing the connectivities of a restart file.
        Also the number of nodes, faces and volumes for later checks.
        """
        # First the header for counts
        header = restartSectionHeader()
        header.name = "NO_FA_CV_NOOFA_COUNTS"
        header.id = ic3_restart_codes["UGP_IO_NO_FA_CV_NOOFA_COUNTS"]
        header.skip = header.hsize
        header.idata[0] = self.params["no_count"]
        header.idata[1] = self.params["fa_count"]
        header.idata[2] = self.params["cv_count"]
        header.idata[3] = self.params["noofa_count"]
        header.write(self.fid, self.endian)
        #
        # don't really know how useful it was; suppressed in V3
        self.__WriteRestartConnectivity_check()
        #
        # Connectivities
        # Faces-to-nodes connectivities
        # Header
        header = restartSectionHeader()
        header.name = "NOOFA_I_AND_V"
        header.id = ic3_restart_codes["UGP_IO_NOOFA_I_AND_V"]
        header.skip = header.hsize + type2nbytes["int32"] * self.params["fa_count"] + \
                                     type2nbytes["int32"] * self.params["noofa_count"]
        header.idata[0] = self.params["fa_count"]
        header.idata[1] = self.params["noofa_count"]
        header.write(self.fid, self.endian)
        # Node count per face
        BinaryWrite(self.fid, self.endian, "i"*self.params["fa_count"], self.f2v["noofa"].tolist()) # remove first 0
        # Flattened face-to-node connectivity
        BinaryWrite(self.fid, self.endian, "i"*self.params["noofa_count"], self.f2v["noofa_v"].tolist())
        # Faces-to-cells connectivities
        # Header
        header = restartSectionHeader()
        header.name = "CVOFA"
        header.id = ic3_restart_codes["UGP_IO_CVOFA"]
        header.skip = header.hsize + type2nbytes["int32"] * self.params["fa_count"] * 2
        header.idata[0] = self.params["fa_count"]
        header.idata[1] = 2
        header.write(self.fid, self.endian)
        # Flattened face-to-cell connectivity
        #print("W",self.f2e)
        # print(self.f2e.ravel())
        BinaryWrite(self.fid, self.endian, "i"*self.params["fa_count"]*2, self.f2e.ravel().tolist())
        # Face zones, a.k.a boundary condition patches
        for key in self.bocos.keys():
            # Header
            header = restartSectionHeader()
            header.name = key
            header.id = ic3_restart_codes["UGP_IO_FA_ZONE"]
            header.skip = header.hsize
            # diff# print(self.bocos[key]["type"], type2zonekind)
            header.idata[0] = type2zonekind[self.bocos[key]["type"]]
            header.idata[1] = self.bocos[key]["fa_range"][0]
            header.idata[2] = self.bocos[key]["fa_range"][1]
            if "periodic_transform" in self.bocos[key].keys():
                for idx, val in enumerate(self.bocos[key]["periodic_transform"]):
                    header.rdata[idx] = val
            header.write(self.fid, self.endian)
        # Partition information
        # Header
        header = restartSectionHeader()
        header.name = "CV_PART"
        header.id = ic3_restart_codes["UGP_IO_CV_PART"]
        header.skip = header.hsize + type2nbytes["int32"] * self.params["cv_count"]
        header.idata[0] = self.params["cv_count"]
        header.idata[1] = self._mesh._cellprop['partition'].get('npart', 1)
        header.write(self.fid, self.endian)
        # The ranks of the processors, default to everybody 0
        BinaryWrite(self.fid, self.endian, "i"*self.params["cv_count"], 
            self._mesh._cellprop['partition'].get('icvpart', np.zeros((self.params["cv_count"],), dtype=np.int32)))
        # Coordinates
        # Header
        header = restartSectionHeader()
        header.name = "X_NO"
        header.id = ic3_restart_codes["UGP_IO_X_NO"]
        header.skip = header.hsize + type2nbytes["float64"] * self.params["no_count"] * 3
        header.idata[0] = self.params["no_count"]
        header.idata[1] = 3
        header.write(self.fid, self.endian)
        # X, Y and Z
        BinaryWrite(self.fid, self.endian, "d"*self.params["no_count"]*3, self.coordinates.ravel(order='C'))

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
        api.io.print("std", "  write section for {} parameter: {}".format(key, value))
        header.id = ic3_restart_codes["UGP_IO_I0"]
        header.skip = header.hsize
        header.idata[0] = value
        header.write(self.fid, self.endian)
        # Write the rest, all double values
        # The monomials first
        for key in ["DT", "TIME"]:
            value = self.simstate[key.lower()]
            api.io.print("std", "  write section for {} parameter: {}".format(key, value))
            header = restartSectionHeader()
            header.name = key
            header.id = ic3_restart_codes["UGP_IO_D0"]
            header.skip = header.hsize
            header.rdata[0] = value
            header.write(self.fid, self.endian)
        # The second level now
        for key in self.simstate["wgt"].keys():
            header = restartSectionHeader()
            header.name = key.upper()+'_WGT'
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
        for key, item in self.vars["nodes"].items():
            # Scalar
            if item.size == item.shape[0]:
                # Header
                header = restartSectionHeader()
                header.name = key
                header.id = ic3_restart_codes["UGP_IO_NO_D1"]
                header.skip = header.hsize + type2nbytes["float64"] * self.params["no_count"]
                header.idata[0] = self.params["no_count"]
                header.write(self.fid, self.endian)
                # Field
                BinaryWrite(self.fid, self.endian, "d"*self.params["no_count"], item)
        for key, item in self.vars["nodes"].items():
            # Vector
            if len(item.shape) == 2:
                # Header
                header = restartSectionHeader()
                header.name = key
                header.id = ic3_restart_codes["UGP_IO_NO_D3"]
                header.skip = header.hsize + type2nbytes["float64"] * self.params["no_count"] * 3
                header.idata[0] = self.params["no_count"]
                header.idata[1] = 3
                header.write(self.fid, self.endian)
                # Field
                BinaryWrite(self.fid, self.endian, "d"*self.params["no_count"]*3, item.ravel(order='C'))
        for key, item in self.vars["nodes"].items():
            # Tensor
            if len(item.shape) == 3:
                pass

        # Then the cell based variables
        for key, cellitem in self.vars["cells"].items():
            item = cellitem.data()
            # Scalar
            if item.size == item.shape[0]:
                # Header
                header = restartSectionHeader()
                header.name = key
                header.id = ic3_restart_codes["UGP_IO_CV_D1"]
                header.skip = header.hsize + type2nbytes["float64"] * self.params["cv_count"]
                header.idata[0] = self.params["cv_count"]
                header.write(self.fid, self.endian)
                # Field
                BinaryWrite(self.fid, self.endian, "d"*self.params["cv_count"], item)
        for key, cellitem in self.vars["cells"].items():
            item = cellitem.data()
            # Vector
            if len(item.shape) == 2:
                # Header
                header = restartSectionHeader()
                header.name = key
                header.id = ic3_restart_codes["UGP_IO_CV_D3"]
                header.skip = header.hsize + type2nbytes["float64"] * self.params["cv_count"] * 3
                header.idata[0] = self.params["cv_count"]
                header.idata[1] = 3
                header.write(self.fid, self.endian)
                # Field
                BinaryWrite(self.fid, self.endian, "d"*self.params["cv_count"]*3, item.ravel(order='C'))
        for key, cellitem in self.vars["cells"].items():
            item = cellitem.data()
            # Tensor
            if len(item.shape) == 3:
                # Header
                header = restartSectionHeader()
                header.name = key
                header.id = ic3_restart_codes["UGP_IO_CV_D33"]
                header.skip = header.hsize + type2nbytes["float64"] * self.params["cv_count"] * 9
                header.idata[0] = self.params["cv_count"]
                header.idata[1] = 3
                header.idata[2] = 3
                header.write(self.fid, self.endian)
                # Field
                BinaryWrite(self.fid, self.endian, "d"*self.params["cv_count"]*9, item.ravel(order='C'))


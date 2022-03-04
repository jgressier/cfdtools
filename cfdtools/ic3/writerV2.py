
# Import modules
import struct
import numpy as np
import numpy.ma as npma
import collections
import sys
import cfdtools.api as api
import cfdtools.meshbase._mesh as _mesh
from cfdtools.ic3._ic3 import *

struct_endian = { 'native':'@', 'little':'<', 'big':'>'}

###################################################################################################
def BinaryWrite(bfile, endian, form, varargs):
    '''
    Method to encapsulate the few lines necessary for the translation
    of formatted components to a piece of formatted binary.
    input   : handle on an open file [type file identifier]
            format of the binary token [type string]
            list of variables to make the bytearray
    '''

    # Assume big-endian to write the file - anyways charlesx can swap
    form = struct_endian[endian] + form

    # Create a packer
    s = struct.Struct(form)
    #print(len(form), form,len(varargs),':')
    # Actually pack the string now
    packed_ba = s.pack(*varargs)

    # Try to write the bytes to the file otherwise crash
    try:
        bfile.write(packed_ba)
    except IOError:
        api.io.print('error',"Fatal error. Could not write to %s. Exiting."%(bfile.name))
        exit()

class restartSectionHeader():
    '''
    This class is designed to handle the header that is present
    before every variable/section/category saved into the restart file.
    '''

    def __init__(self):
        """
        Initialize a section header class, as it is always structured
        in the same way, i.e. with a name, an id, a skip bytes count,
        an information array, and an auxiliary array needed in specific
        cases like periodicity.
        """

        self.name   = " "*ic3_restart_codes["UGP_IO_HEADER_NAME_LEN"]
        self.id     = np.zeros((1,), dtype=np.int32)
        self.skip   = np.zeros((1,), dtype=np.int64)
        self.idata  = np.zeros((8,), dtype=np.int64)
        self.rdata  = np.zeros((16,), dtype=np.float64)

        # Initialize the format string
        self.binstr = ""
        # Fill it with prior knowledge of its content
        for i in range(ic3_restart_codes["UGP_IO_HEADER_NAME_LEN"]) :
            self.binstr += "c"               # name
        self.binstr += "i"                   # id
        self.binstr += "q"                   # skip
        self.binstr += "qqqqqqqq"            # idata (long long)
        self.binstr += "dddddddddddddddd"    # rdata (double)

        # Keep in mind the total byte size of the header
        self.hsize = ic3_restart_codes["UGP_IO_HEADER_NAME_LEN"] * type2nbytes["char"] +\
                type2nbytes["int32"] +\
                type2nbytes["int64"] +\
                self.idata.size * type2nbytes["int64"] +\
                self.rdata.size * type2nbytes["float64"]

    def write(self, bfile, endian):
        """
        Once initialization is done, this method actually writes
        the header data from the packed binary formatted string.
        """

        # Make a list from all the arguments
        varargs = []
        for i in range(ic3_restart_codes["UGP_IO_HEADER_NAME_LEN"]):
            if i < len(self.name):
                varargs.append(bytes(self.name[i], 'utf-8'))
            else:
                varargs.append(b'\0')
        # for i in varargs:
        #     print(">",i, type(i))
        varargs.append(self.id)
        varargs.append(self.skip)
        for kk in range(8):
            varargs.append(self.idata[kk])
        for kk in range(16):
            varargs.append(self.rdata[kk])
        #print(varargs)
        # Now write everything at once
        BinaryWrite(bfile, endian, self.binstr, varargs)

    def __str__(self):
        mystring="Header:\n"
        mystring +="Name:%s\n"%self.name

        mystring +="Id:%i\n"%self.id
        mystring +="Skip:%i\n"%self.skip
        mystring +="idata:("
        for i in range(8):
            mystring +="%i,"%(self.idata[i])
        mystring +=")\n"
        mystring +="rdata:("
        for i in range(16):
            mystring +="%f,"%(self.rdata[i])
        mystring +=")"
        return mystring

###################################################################################################

@api.fileformat_writer('IC3', '.ic3')
class writer():
    ''' Implementation of the writer to write ic3 restart files '''

    def __init__(self, mesh, endian='native'):
        """
        Initialization of a ic3 restart file writer.
        """
        api.io.print('std',"Initialization of IC3 writer")
        if not endian in struct_endian.keys():
            raise ValueError("unknown endian key")
        else:
            self.endian = endian
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

    def set_simstate(self, state=None):
        """
        Set the simulation state for the restart.
        It will be later written to the file using the
        __WriterInformativeValues method.
        It gives informations to CX about the current iteration number,
        the timestep, and the current time.
        """

        # Initialize the current simulation state
        self.simstate = {"wgt":{}, "step":0, "dt":0, "time":0}

        # If the state dictionary is not None then can replace
        if state is not None:
            for key, item in state.items():
                self.simstate[key] = item

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

        api.io.print('std',"Writing header ..")
        self.__WriteRestartHeader()
        #
        api.io.print('std',"Writing connectivity ..")
        self.__WriteRestartConnectivity()
        #
        api.io.print('std',"Writing informative values ..")
        self.__WriteInformativeValues()
        #
        api.io.print('std',"Writing variables ..")
        self.__WriteRestartVar()
        #
        api.io.print('std',"End of file ..")
        header = restartSectionHeader()
        header.name = "EOF"
        header.id = ic3_restart_codes["UGP_IO_EOF"]
        header.skip = header.hsize
        header.write(self.fid, self.endian)

        # Before returning, close the file
        self.fid.close()
        del self.fid

    def __WriteRestartHeader(self):
        """
        Method writing the header of a restart file.
        It is composed of two integers, the "magic number" used as a flag for endianness
        and the CharlesX version number.
        input   : handle on an open restart file, [type file identifier]
        """
        # Write the two integers
        BinaryWrite(self.fid, self.endian, "ii", 
                                    [ic3_restart_codes["UGP_IO_MAGIC_NUMBER"],
                                     ic3_restart_codes["UGP_IO_VERSION"]])

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
        # print(self.f2e)
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
        header.idata[1] = 1
        header.write(self.fid, self.endian)
        # The ranks of the processors, default to everybody 0
        BinaryWrite(self.fid, self.endian, "i"*self.params["cv_count"], np.zeros((self.params["cv_count"],), dtype=np.int32))
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
        header = restartSectionHeader()
        header.name = "STEP"
        header.id = ic3_restart_codes["UGP_IO_I0"]
        header.skip = header.hsize
        header.idata[0] = self.simstate["step"]
        header.write(self.fid, self.endian)
        # Write the rest, all double values
        # The monomials first
        for key in ["DT", "TIME"]:
            header = restartSectionHeader()
            header.name = key
            header.id = ic3_restart_codes["UGP_IO_D0"]
            header.skip = header.hsize
            header.rdata[0] = self.simstate[key.lower()]
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

        # Check the variable dictionary has been set
        try:
            self.vars
        except:
            self.vars = {"nodes":{}, "cells":{}}

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
        for key, item in self.vars["cells"].items():
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
        for key, item in self.vars["cells"].items():
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
        for key, item in self.vars["cells"].items():
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


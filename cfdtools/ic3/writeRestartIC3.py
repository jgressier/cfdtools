#!/usr/bin/env python

# Import modules
import struct
import numpy as np
import numpy.ma as npma
import collections
import sys

# Copy the restart codes from CharlesX to match its restart routine
ic3_restart_codes = {"UGP_IO_MAGIC_NUMBER":123581321,
                     "UGP_IO_VERSION":2,
                     "UGP_IO_HEADER_NAME_LEN":52, # select this to make header size 256 bytes
                     "UGP_IO_NO_FA_CV_NOOFA_COUNTS":10,
                     "UGP_IO_I0":11,
                     "UGP_IO_D0":12,
                     "UGP_IO_FA_CHECK":20,
                     "UGP_IO_NOOFA_I_AND_V":21,
                     "UGP_IO_CVOFA":22,
                     "UGP_IO_FA_ZONE":23,
                     "UGP_IO_FA_D1":24,
                     "UGP_IO_FA_D3":25,
                     "UGP_IO_FA_II1":26,
                     "UGP_IO_NO_CHECK":30,
                     "UGP_IO_X_NO":31,
                     "UGP_IO_NO_D1":32,
                     "UGP_IO_NO_D3":33,
                     "UGP_IO_NO_II1":34,
                     "UGP_IO_CV_CHECK":40,
                     "UGP_IO_CV_PART":41,
                     "UGP_IO_CV_D1":42,
                     "UGP_IO_CV_D3":43,
                     "UGP_IO_CV_D33":44,
                     "UGP_IO_CV_II1":45,
                     "UGP_IO_DATA":50,
                     "UGP_IO_EOF":51,}
# Dictionary to convert type of data [string] to a number of bytes for clean binary parsing
type2nbytes = {"char":1,
               "int32":4,
               "int64":8,
               "float32":4,
               "float64":8,
               "float128":16,}
# Actual number of vertices for a given cell type
nodes_per_cell = {
    'bi': 2,
    'tri': 3,
    'qua': 4,
    'tet': 4,
    'hex': 8,
    'pri': 6,
    'pyr': 5,}
cell_from_nodes = {
    2:'bi',
    3:'tri',
    4:'qua',
    4:'tet',
    8:'hex',
    6:'pri',
    5:'pyr',}
faces_per_cell = {
    'hex': 6,
    'pri': 5,
    'pyr': 5,
    'tet': 4,}
faces_of_cell = {
    'hex':['qua']*6,
    'pri':['tri', 'qua', 'qua', 'qua', 'tri'],
    'pyr':['qua'] + ['tri']*4,
    'tet':['tri']*4,}
ifaces_of_cell = {
    'hex':[[0,3,2,1], [0,1,5,4], [1,2,6,5], [2,3,7,6], [3,0,4,7], [4,5,6,7]],
    'pri':[[0,2,1], [0,1,4,3], [1,2,5,4], [2,0,3,5], [3,4,5]],
    'pyr':[[0,3,2,1], [0,1,4], [1,2,4], [2,3,4], [3,0,4]],
    'tet':[[0,2,1], [0,1,3], [1,2,3], [2,0,3]],}
type2zonekind = {"boundary":1,
                 "periodic_cart":2,
                 "periodic_cylx":3,
                 "periodic_cyly":4,
                 "periodic_cylz":5,
                 "internal":6}

###################################################################################################

class restartSectionHeader():
    '''
    This class is designed to handle the header that is present
    before every variable/section/category saved into the restart file.
    '''

    def __init__(self):
        '''
        Initialize a section header class, as it is always structured
        in the same way, i.e. with a name, an id, a skip bytes count,
        an information array, and an auxiliray array needed in specific
        cases like periodicity.
        '''

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
        self.binstr += "qqqqqqqq"            # idata
        self.binstr += "dddddddddddddddd"    # rdata

        # Keep in mind the total byte size of the header
        self.hsize = ic3_restart_codes["UGP_IO_HEADER_NAME_LEN"] * type2nbytes["char"] +\
                type2nbytes["int32"] +\
                type2nbytes["int64"] +\
                self.idata.size * type2nbytes["int64"] +\
                self.rdata.size * type2nbytes["float64"]

    def write(self, bfile):
        '''
        Once initialization is done, this method actually writes
        the header data from the packed binary formatted string.
        '''

        # Make a list from all the arguments
        varargs = []
        for i in range(ic3_restart_codes["UGP_IO_HEADER_NAME_LEN"]):
            if i < len(self.name):
                varargs.append(self.name[i])
            else:
                varargs.append('\0')
        varargs.append(self.id)
        varargs.append(self.skip)
        for kk in xrange(8):
            varargs.append(self.idata[kk])
        for kk in xrange(16):
            varargs.append(self.rdata[kk])

        # Now write everything at once
        BinaryWrite(bfile, self.binstr, varargs)

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

def BinaryWrite(bfile, form, varargs):
    '''
    Method to encapsulate the few lines necessary for the translation
    of formatted components to a piece of formatted binary.
    input   : handle on an open file [type file identifier]
              format of the binary token [type string]
              list of variables to make the bytearray
    '''

    # Assume big-endian to write the file - anyways charlesx can swap
    form = '>' + form

    # Create a packer
    s = struct.Struct(form)

    # Actually pack the string now
    packed_ba = s.pack(*varargs)

    # Try to write the bytes to the file otherwise crash
    try:
        bfile.write(packed_ba)
    except IOError:
        print "Fatal error. Could not write to %s. Exiting."%(bfile.name)
        exit()

###################################################################################################

class WriterRestartIC3():
    ''' Implementation of the writer to write ic3 restart files '''

    def __init__(self, filename):
        '''
        Initialization of a ic3 restart file writer.
        '''

        print
        print ":: WRITER RESTART IC3 ::"
        print


        # Keep the filename
        self.filename = filename

        # Initialize the simulation state
        self.set_simstate()

    def set_mesh(self, xyz, co):
        '''
        Store the coordinates of the nodes and the
        element-to-vertex connectivity. Also triggers the creation
        of a face-to-element and face-to-vertex connectivity
        as per CharlesX implementation.
        '''

        print "Setting coordinates and connectivity arrays.."

        # Keep the coordinates of the points
        self.coordinates = xyz

        # Keep the element-to-vertex connectivity
        self.e2v = collections.OrderedDict({})
        for key, cco in co.items():
            self.e2v[key] = cco
        must_have_keys = ["pri", "hex", "tet", "pyr"]
        for key in must_have_keys:
            if key not in self.e2v.keys():
                self.e2v[key] = None

        # Compute the number of nodes and elements
        self.params = {}
        # Nodes
        self.params["no_count"] = self.coordinates.shape[0]
        # Elements, count throughout all connectivities
        self.params["cv_count"] = 0
        for key in self.e2v:
            if self.e2v[key] is not None:
                self.params["cv_count"] += self.e2v[key].shape[0]
            else:
                self.e2v[key] = np.empty((0,), dtype=int)

        # Call upon the C function to do the heavy lifting
        from getFaceConnectivities import getfaceco
        reverse = int(0)
        trico, quaco, fa2cv, vtx2fa, ntri, nqua, hasNegVols = getfaceco(self.coordinates,
                                                                        self.e2v["hex"], self.e2v["pri"], self.e2v["pyr"], self.e2v["tet"],
                                                                        self.params["no_count"], self.params["cv_count"],
                                                                        np.array([nodes_per_cell[key] for key in self.e2v.keys()]),
                                                                        reverse)
        if hasNegVols:
            print "\tWarning. Negative volumes detected. Inverting all the connectivities."
            reverse = int(1)
            trico, quaco, fa2cv, vtx2fa, ntri, nqua, hasNegVols = getfaceco(self.coordinates,
                                                                            self.e2v["hex"], self.e2v["pri"], self.e2v["pyr"], self.e2v["tet"],
                                                                            self.params["no_count"], self.params["cv_count"],
                                                                            np.array([nodes_per_cell[key] for key in self.e2v.keys()]),
                                                                            reverse)

        # Store the information properly
        self.params["fa_count"] = ntri + nqua
        self.params["noofa_count"] = ntri*3 + nqua*4
        self.f2v = {"tri":trico, "qua":quaco}
        self.f2e = fa2cv
        self.v2f = vtx2fa

        # Concatenate the face-to-cell connectivities
        self.f2v["triAqua"] = -np.ones((self.params["fa_count"], 4), dtype=self.f2v["qua"].dtype)
        self.f2v["triAqua"][:ntri, :3] = self.f2v["tri"]
        self.f2v["triAqua"][ntri:, :] = self.f2v["qua"]

        # Flatten the face-to-cell connectivity
        # First the tris then the quads
        self.f2v["noofa_i"] = np.zeros((self.params["fa_count"],), dtype=np.int32)
        self.f2v["noofa_i"][:ntri] = 3
        self.f2v["noofa_i"][ntri:] = 4

        print "ok."

    def set_bocos(self, bocos, nboco=None):
        '''
        Set the boundary conditions dictionary for the restart.
        It will be later written to the file using the
        __WriteRestartConnectivity method.
        The bocos slicing is expressed in terms of faces.
        '''

        print "Setting boundary conditions.."

        # Initialize the dictionary
        self.bocos = collections.OrderedDict({})

        # Get the number of boundary conditions
        gotInt = False
        c_nboco = 0
        for key in bocos.keys():
            if 'nfa' in key.lower():
                continue
            if np.any([string in key.lower() for string in ["int", "interior"]]):
                gotInt = True
            c_nboco += 1
        if not gotInt:
            c_nboco += 1
        #
        if nboco is None:
            self.nboco = c_nboco
        else:
            self.nboco = c_nboco
            if nboco != c_nboco:
                if gotInt and c_nboco == nboco+1:
                    self.nboco = c_nboco
                else:
                    if gotInt:
                        c_nboco -= 1
                    print "Fatal error. Incoherence between the bocos dictionary (%d bocos) and the number (%d) of bocos declared. Exiting."%(c_nboco, nboco)
                    exit()

        # Browse the bocos and prepare them to be written
        passed_key = None
        for key, boco in bocos.items():
            # Skip the eventual extra info introduced by a CX reader
            if 'nfa' in key.lower():
                continue
            # Store the relevant information
            self.bocos[key] = collections.OrderedDict({})
            # The type of boco, defaulting to "boundary"
            if "type" in boco.keys():
                self.bocos[key]["kind"] = boco["type"]
            else:
                self.bocos[key]["kind"] = 1 # "boundary"
            # The face slicing
            if not "slicing" in boco.keys():
                print "Fatal error. No slicing found for boundary condition %s. Exiting."%(key)
                exit()
            else:
                if "int" in key.lower() or "interior" in key.lower() or "fluid" in key.lower():
                    passed_key = key
                else:
                    print "\tHandling boundary condition %s.."%key
                    # Call upon a C function for the heavy lifting
                    from getFaceSlicings import getfacesli
                    self.bocos[key]["fa_slicing"] = getfacesli(boco["slicing"], self.v2f, self.f2v["triAqua"])
                    print "\tok."
            # Periodicity information (if any)
            if "periodic" in key.lower() or \
               "periodic" in self.bocos[key]["kind"]:
                if "periodic_transform" in boco.keys():
                    self.bocos[key]["periodic_transform"] = boco["periodic_transform"]
                else:
                    print "Fatal error. Periodic (a priori) boco (%s) found but no information on the transform can be found. Exiting."%(key)
                    exit()

        # Get the face slicing for the eventual passed interior key
        if passed_key is None:
            passed_key = "int_FLUID"
            self.bocos[passed_key] = collections.OrderedDict({})
            self.bocos[passed_key]["kind"] = "internal"
        if passed_key is not None:
            print "\tHandling boundary condition %s.."%passed_key
            all_faces = np.arange(self.params["fa_count"])
            to_delete = []
            for key in self.bocos:
                if "fa_slicing" in self.bocos[key].keys():
                    to_delete += self.bocos[key]["fa_slicing"].tolist()
            all_faces = np.delete(all_faces, to_delete)
            self.bocos[passed_key]["fa_slicing"] = all_faces
            print "\tok."

        # Reorder the bocos to put the interior faces last
        newbocos = collections.OrderedDict({})
        for key, item in self.bocos.items():
            if "int" in key.lower() or "interior" in key.lower():
                continue
            newbocos[key] = item
        newbocos[passed_key] = self.bocos[passed_key]
        self.bocos = newbocos

        # Gotta reorder all the faces according to bocos now
        noofa_i = np.zeros((self.params["fa_count"],), dtype=np.int32)
        noofa_v = []
        f2e = []
        last = 0
        for key in self.bocos.keys():
            sz = self.bocos[key]["fa_slicing"].size
            noofa_v += self.f2v["triAqua"][self.bocos[key]["fa_slicing"],:].ravel().tolist()
            noofa_i[last:last+sz] = self.f2v["noofa_i"][self.bocos[key]["fa_slicing"]]
            self.bocos[key]["fa_range"] = [last, last+sz]
            f2e += self.f2e[self.bocos[key]["fa_slicing"],:].ravel().tolist()
            last += sz
        to_delete = np.where(np.asarray(noofa_v) == -1)
        noofa_v = np.delete(noofa_v, to_delete)

        # Re-assign them now
        self.f2v["noofa_i"] = np.asarray(noofa_i)
        self.f2v["noofa_v"] = np.asarray(noofa_v)
        self.f2e = np.asarray(f2e)

        print "ok."

    def set_simstate(self, state=None):
        '''
        Set the simulation state for the restart.
        It will be later written to the file using the
        __WriterInformativeValues method.
        It gives informations to CX about the current iteration number,
        the timestep, and the current time.
        '''

        # Initialize the current simulation state
        self.simstate = {"wgt":{}, "step":0, "dt":0, "time":0}

        # If the state dictionary is not None then can replace
        if state is not None:
            for key, item in state.items():
                self.simstate[key] = item

    def set_vars(self, nodesvar=None, cellsvar=None):
        '''
        Set the variables for the restart.
        They will be later written to the file using the
        __WriteRestartVar method.
        '''

        print "Setting variables..",

        self.vars = {"nodes":{}, "cells":{}}

        # Start with the variables stored at the vertices
        if nodesvar is not None:
            print "at nodes..",
            for key, item in nodesvar.items():
                self.vars["nodes"][key] = item

        # Then the variables stored at the cells:
        if cellsvar is not None:
            print "at cell centers..",
            for key, item in cellsvar.items():
                self.vars["cells"][key] = item

        print "ok."

    def write_data(self):
        '''
        Main method of the ic3 restart file writer
        '''

        # Open the file for binary reading
        self.fid = open(self.filename, "wb")

        print "Writing header .."
        self.__WriteRestartHeader()
        #
        print "Writing connectivity .."
        self.__WriteRestartConnectivity()
        #
        print "Writing informative values .."
        self.__WriteInformativeValues()
        #
        print "Writing variables .."
        self.__WriteRestartVar()
        #
        print "End of file .."
        header = restartSectionHeader()
        header.name = "EOF"
        header.id = ic3_restart_codes["UGP_IO_EOF"]
        header.skip = header.hsize
        header.write(self.fid)

        # Before returning, close the file
        self.fid.close()
        del self.fid

    def __WriteRestartHeader(self):
        '''
        Method writing the header of a restart file.
        It is composed of two integers, the "magic number" used as a flag for endianness
        and the CharlesX version number.
        input   : handle on an open restart file, [type file identifier]
        '''

        # Write the two integers
        BinaryWrite(self.fid, "ii", [ic3_restart_codes["UGP_IO_MAGIC_NUMBER"],
                                     ic3_restart_codes["UGP_IO_VERSION"]])

    def __WriteRestartConnectivity(self):
        '''
        Method writing the connectivities of a restart file.
        Also the number of nodes, faces and volumes for later checks.
        '''

        # First the header for counts
        header = restartSectionHeader()
        header.name = "NO_FA_CV_NOOFA_COUNTS"
        header.id = ic3_restart_codes["UGP_IO_NO_FA_CV_NOOFA_COUNTS"]
        header.skip = header.hsize
        header.idata[0] = self.params["no_count"]
        header.idata[1] = self.params["fa_count"]
        header.idata[2] = self.params["cv_count"]
        header.idata[3] = self.params["noofa_count"]
        header.write(self.fid)

        # Node check
        # Header
        header = restartSectionHeader()
        header.name = "NO_CHECK"
        header.id = ic3_restart_codes["UGP_IO_NO_CHECK"]
        header.skip = header.hsize + type2nbytes["int32"] * self.params["no_count"]
        header.idata[0] = self.params["no_count"]
        header.write(self.fid)
        # Write the nodes global ids
        BinaryWrite(self.fid, "i"*self.params["no_count"], range(self.params["no_count"]))

        # Face check
        # Header
        header = restartSectionHeader()
        header.name = "FA_CHECK"
        header.id = ic3_restart_codes["UGP_IO_FA_CHECK"]
        header.skip = header.hsize + type2nbytes["int32"] * self.params["fa_count"]
        header.idata[0] = self.params["fa_count"]
        header.write(self.fid)
        # Write the nodes global ids
        BinaryWrite(self.fid, "i"*self.params["fa_count"], range(self.params["fa_count"]))

        # Cell check
        # Header
        header = restartSectionHeader()
        header.name = "CV_CHECK"
        header.id = ic3_restart_codes["UGP_IO_CV_CHECK"]
        header.skip = header.hsize + type2nbytes["int32"] * self.params["cv_count"]
        header.idata[0] = self.params["cv_count"]
        header.write(self.fid)
        # Write the nodes global ids
        BinaryWrite(self.fid, "i"*self.params["cv_count"], range(self.params["cv_count"]))

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
        header.write(self.fid)
        # Node count
        BinaryWrite(self.fid, "i"*self.params["fa_count"], self.f2v["noofa_i"].tolist())
        # Flattened face-to-node connectivity
        BinaryWrite(self.fid, "i"*self.params["noofa_count"], self.f2v["noofa_v"].tolist())
        # Faces-to-cells connectivities
        # Header
        header = restartSectionHeader()
        header.name = "CVOFA"
        header.id = ic3_restart_codes["UGP_IO_CVOFA"]
        header.skip = header.hsize + type2nbytes["int32"] * self.params["fa_count"] * 2
        header.idata[0] = self.params["fa_count"]
        header.idata[1] = 2
        header.write(self.fid)
        # Flattened face-to-cell connectivity
        BinaryWrite(self.fid, "i"*self.params["fa_count"]*2, self.f2e.tolist())

        # Face zones, a.k.a boundary condition patches
        for key in self.bocos.keys():
            # Header
            header = restartSectionHeader()
            header.name = key
            header.id = ic3_restart_codes["UGP_IO_FA_ZONE"]
            header.skip = header.hsize
            header.idata[0] = type2zonekind[self.bocos[key]["kind"]]
            header.idata[1] = self.bocos[key]["fa_range"][0]
            header.idata[2] = self.bocos[key]["fa_range"][1] - 1
            if "periodic_transform" in self.bocos[key].keys():
                for idx, val in enumerate(self.bocos[key]["periodic_transform"]):
                    header.rdata[idx] = val
            header.write(self.fid)

        # Partition information
        # Header
        header = restartSectionHeader()
        header.name = "CV_PART"
        header.id = ic3_restart_codes["UGP_IO_CV_PART"]
        header.skip = header.hsize + type2nbytes["int32"] * self.params["cv_count"]
        header.idata[0] = self.params["cv_count"]
        header.idata[1] = 1
        header.write(self.fid)
        # The ranks of the processors, default to everybody 0
        BinaryWrite(self.fid, "i"*self.params["cv_count"], np.zeros((self.params["cv_count"],), dtype=np.int32))

        # Coordinates
        # Header
        header = restartSectionHeader()
        header.name = "X_NO"
        header.id = ic3_restart_codes["UGP_IO_X_NO"]
        header.skip = header.hsize + type2nbytes["float64"] * self.params["no_count"] * 3
        header.idata[0] = self.params["no_count"]
        header.idata[1] = 3
        header.write(self.fid)
        # X, Y and Z
        BinaryWrite(self.fid, "d"*self.params["no_count"]*3, self.coordinates.ravel(order='C'))

    def __WriteInformativeValues(self):
        '''
        Method writing the solution information to the restart file
        '''

        # Start of data
        # An extra header to announce what is coming
        header = restartSectionHeader()
        header.name = "DATA"
        header.id = ic3_restart_codes["UGP_IO_DATA"]
        header.skip = header.hsize
        header.write(self.fid)

        # Write current number of iteration
        # Header
        header = restartSectionHeader()
        header.name = "STEP"
        header.id = ic3_restart_codes["UGP_IO_I0"]
        header.skip = header.hsize
        header.idata[0] = self.simstate["step"]
        header.write(self.fid)

        # Write the rest, all double values
        # The monomials first
        for key in ["DT", "TIME"]:
            header = restartSectionHeader()
            header.name = key
            header.id = ic3_restart_codes["UGP_IO_D0"]
            header.skip = header.hsize
            header.rdata[0] = self.simstate[key.lower()]
            header.write(self.fid)
        # The second level now
        for key in self.simstate["wgt"].keys():
            header = restartSectionHeader()
            header.name = key.upper()+'_WGT'
            header.id = ic3_restart_codes["UGP_IO_D0"]
            header.skip = header.hsize
            header.rdata[0] = self.simstate["wgt"][key]
            header.write(self.fid)

    def __WriteRestartVar(self):
        '''
        Method to write all the variables into a restart file.
        Scalars, vectors and tensors all together.
        '''

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
                header.write(self.fid)
                # Field
                BinaryWrite(self.fid, "d"*self.params["no_count"], item)
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
                header.write(self.fid)
                # Field
                BinaryWrite(self.fid, "d"*self.params["no_count"]*3, item.ravel(order='C'))
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
                header.write(self.fid)
                # Field
                BinaryWrite(self.fid, "d"*self.params["cv_count"], item)
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
                header.write(self.fid)
                # Field
                BinaryWrite(self.fid, "d"*self.params["cv_count"]*3, item.ravel(order='C'))
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
                header.write(self.fid)
                # Field
                BinaryWrite(self.fid, "d"*self.params["cv_count"]*9, item.ravel(order='C'))

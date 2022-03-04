# coding: utf8

# Import modules
import struct
import numpy as np
import numpy.ma as npma
import collections
import sys,os
import cfdtools.api as api
import cfdtools.meshbase._mesh as _mesh
from cfdtools.ic3._ic3 import *

def _printreadable(string, value):
    if isinstance(value, (int, float, str, np.int32, np.int64)):
        print(string+':',value)
    elif isinstance(value, np.ndarray):
        if value.size <= 10:
            print(string+': ndarray',value.shape, value)
        else:
            print(string+': ndarray',value.shape)
    else:
        print(string+': '+str(type(value)))

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

        self.name   = ""
        self.id     = np.zeros((1,), dtype=np.int32)
        self.skip   = np.zeros((1,), dtype=np.int64)
        self.idata  = np.zeros((8,), dtype=np.int64)
        self.rdata  = np.zeros((16,), dtype=np.float64)

    def readVar(self, bfile, byte_swap,nametypes, reset_offset=True):
        '''
        Once initialization is done, this method actually reads
        the header data from the packed binary formatted string.
        '''
        id_list=[]
        for nametype in nametypes:
            id_list.append(ic3_restart_codes[nametype])
        if(reset_offset is True): bfile.seek(8, os.SEEK_SET)
        #print self.skip[0],self.skip
        while (self.id[0] not in id_list) and (self.id[0] !=ic3_restart_codes["UGP_IO_EOF"]):
            
            self.name= ""
            self.id     = np.zeros((1,), dtype=np.int32)
            self.idata  = np.zeros((8,), dtype=np.int64)
            self.rdata  = np.zeros((16,), dtype=np.float64)
            
            api.io.print('debug', "skipping data %d"%self.skip[0])
            if (self.skip[0]>0):
                bfile.seek(self.skip[0]-256, os.SEEK_CUR)
            
            self.skip   = np.zeros((1,), dtype=np.int64)
            # Initialize the format string
            binstr = ""
            # Fill it with prior knowledge of its content
            for i in range(ic3_restart_codes["UGP_IO_HEADER_NAME_LEN"]) :
                binstr += "c"               # name
            binstr += "i"                   # id
            binstr += "q"                   # skip
            binstr += "qqqqqqqq"            # idata
            binstr += "dddddddddddddddd"    # rdata
            # Keep in mind the total byte size of the header
            self.hsize = ic3_restart_codes["UGP_IO_HEADER_NAME_LEN"] * type2nbytes["char"] +\
                    type2nbytes["int32"] +\
                    type2nbytes["int64"] +\
                    self.idata.size * type2nbytes["int64"] +\
                    self.rdata.size * type2nbytes["float64"]


            # Split the packed binary string into its tokens
            s = list(BinaryRead(bfile, binstr, byte_swap, self.hsize))

            # Remove trailing space from header name (h.name)
            i = 0
            while (s[i] != b'\x00') and  ( i < ic3_restart_codes["UGP_IO_HEADER_NAME_LEN"] ) : #remove trailing spaces
                self.name += s[i].decode() # .decode() for python3 portage
                i += 1
            if(self.name.startswith("DEOF")): break 
            # Store the rest of the tokens in the right namespace
            nlen = ic3_restart_codes["UGP_IO_HEADER_NAME_LEN"]
            self.id[0] = s[nlen]
            self.skip[0] = s[nlen+self.id.size]
            for i in range(self.idata.size):
                self.idata[i] = s[nlen+self.id.size+self.skip.size+i]
            for i in range(self.rdata.size):
                self.rdata[i] = s[nlen+self.id.size+self.skip.size+self.idata.size+i]
            #print(self)
        return (self.id[0] in id_list)
        
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

def BinaryRead(bfile, form, byte_swap, size):
    '''
    Method to encapsulate the few lines necessary for the translation
    of a piece of formatted binary onto its components given its format.
    input   : handle on an open restart file [type file identifier]
              format of the binary token [type string]
              endianness flag [type boolean]
              size of the binary token to be read [type int]
    '''

    # Choose the right prefix to the format string based on endianness
    if byte_swap is True: # little-endian
        form = '<' + form
    else: # big-endian
        form = '>' + form

    # Create a "packed-binary reader structure"
    s = struct.Struct(form)

    # Try to read 'size' bytes from the file otherwise crash
    try:
        while True:
            record = bfile.read(size)
            if len(record) != size:
                break;
            return s.unpack(record)
    except IOError:
        api.io.print('error', "Fatal error. Could not read %d bytes from %s. Exiting."%(size, bfile.name))
        exit()

###################################################################################################

@api.fileformat_reader('IC3', '.ic3')
class reader(api._files):
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

    def __str__(self):
        s = '    filename: '+self.filename
        s+= '\n   simulation: '+str(list(self.simulation_state.keys()))
        s+= '\n    mesh keys: '+str(list(self.mesh.keys()))
        s+= '\nvariable keys: '+str(list(self.variables.keys()))
        return s

    def printinfo(self):
        print(self)
        print("- mesh properties")
        for key,item in self.mesh.items():
            if isinstance(item,dict):
                for key2,item2 in item.items():
                    if isinstance(item2,dict):
                        for key3,item3 in item2.items():
                            _printreadable('  mesh.'+key+'.'+key2+'.'+key3, item3)
                    else:
                        _printreadable('  mesh.'+key+'.'+key2, item2)
            else:
                _printreadable('  mesh.'+key, item)
        print("- variable properties")
        for key,item in self.variables.items():
            for key2,value in item.items():
                _printreadable('variables.'+key+'.'+key2, value)
                # for key3,value3 in value.items():
                #     print('variables.'+key+'.'+key2+'.'+key3+':'+str(type(value3)))

    def read_data(self):
        '''
        Main method of the IC3 restart reader.
        Parses in order the file using sub-methods described below.
        output  : the mesh itself
                  the lot of variables stored in the restart file
                  information on the state of the simulation
        '''
        api.io.print('std',":: READER RESTART IC3 ::")

        if not self._exists:
            print("Fatal error. File %s cannot be found."%(self.filename))
            exit()

        # Open the file for binary reading
        api.io.print('debug','reading ',self.filename)
        self.fid = open(self.filename, "rb")

        api.io.print('std', "Reading header ..")
        self.__ReadRestartHeader()
        #
        api.io.print('std', "Reading connectivity ..")
        self.__ReadRestartConnectivity()
        #
        api.io.print('std', "Reading informative values ..")
        self.__ReadInformativeValues()
        #
        api.io.print('std', "Reading variables ..")
        self.__ReadRestartVar()
        #

        # Before returning, close the file
        self.fid.close()
        del self.fid

        #return self.mesh["coordinates"], self.mesh["connectivity"]["e2v"], self.mesh["bocos"], self.variables["nodes"], self.variables["cells"], (self.simulation_state, self.mesh["params"])
        meshdata = _mesh.mesh(self.mesh['params']['cv_count'], self.mesh['params']['no_count'])
        meshdata.set_nodescoord_nd(self.mesh['coordinates'])
        meshdata.set_face2cell(self.mesh['connectivity']['cvofa'])
        meshdata.set_face2node(self.mesh['connectivity']['noofa'])
        meshdata.set_bocos(self.mesh['bocos'])
        meshdata.set_celldata(self.variables['cells'])
        meshdata.set_nodedata(self.variables['nodes'])
        meshdata.set_facedata(self.variables['faces'])
        meshdata.set_params(self.mesh['params'])
        meshdata.update_params(self.simulation_state)

        return meshdata

    def __ReadRestartHeader(self):
        '''
        Method reading the header of a restart file.
        It is composed of two integers, the "magic number" used as a flag for endianness
        and the CharlesX version number.
        input   : handle on an open restart file, [type file identifier]
        output  : the endianness of the open restart file [type boolean]
        '''

        # By default suppose big-endian format
        self.byte_swap = False

        # Read the first integer (int64)
        s = list(BinaryRead(self.fid, "ii", False, 2*type2nbytes["int32"]))
        # If, with big-endian assumption, the first integer comes out wrong, swap to little-endian
        if s[0] != ic3_restart_codes["UGP_IO_MAGIC_NUMBER"]:
            # Change the flag
            self.byte_swap=True
            # Transform the second integer of the list to match the version number
            aux_struct = struct.Struct(">i")
            packed_version = aux_struct.pack(s[1])
            del aux_struct
            aux_struct = struct.Struct("<i")
            s[1] = aux_struct.unpack(packed_version)[0]

        # Some info for the user
        api.io.print('std', "\t Restart file is version %s, current IC3 version is %d."%(str(s[1]).strip(), ic3_restart_codes["UGP_IO_VERSION"]))

    def __ReadRestartConnectivity(self):
        '''
        Method reading the first blocks passed the header, containing informations
        on the nodes, the faces, the cells and the connectivity in between those
        input   : handle on an open restart file [type file identifier]
                  endianness flag [boolean]
        output  : mesh structure containing all the necessary information to build the grid
        '''

        # Initialize the mesh
        self.mesh = {"params":{"no_count":0, "fa_count":0, "cv_count":0, "noofa_count":0, "nboco":0},
                     "connectivity":{"noofa":{}, "cvofa":{}, "nkeys":0},
                     "coordinates":None,
                     "bocos":{"nfa_b":0, "nfa_bp":0},
                     "partition":None}

        # Get the header
        h = restartSectionHeader()
        if (not h.readVar(self.fid, self.byte_swap,["UGP_IO_NO_FA_CV_NOOFA_COUNTS"])): exit()

        # Store the size informations at the right place
        self.mesh["params"]["no_count"] = h.idata[0]
        self.mesh["params"]["fa_count"] = h.idata[1]
        self.mesh["params"]["cv_count"] = h.idata[2]
        self.mesh["params"]["noofa_count"] = h.idata[3]
        del h

        # Integrity check
        if self.check_integrity:
            # Check the restart is whole by counting the global ids of no, fa and cv
            # For nodes
            api.io.print('std', "\t Checking nodes integrity .."); sys.stdout.flush()
            h = restartSectionHeader()
            if(not h.readVar(self.fid, self.byte_swap,["UGP_IO_NO_CHECK"])): exit()

            assert h.idata[0] == self.mesh["params"]["no_count"]
            assert h.id[0] == ic3_restart_codes["UGP_IO_NO_CHECK"]
            nodes_id = np.empty((self.mesh["params"]["no_count"],), dtype=np.int32)
            for loopi in range(self.mesh["params"]["no_count"]):  # xrange to range (python3 portage)
                s = BinaryRead(self.fid, "i", self.byte_swap, type2nbytes["int32"])
                nodes_id[loopi] = s[0]
            assert np.allclose(nodes_id, np.arange(self.mesh["params"]["no_count"]))
            api.io.print('std', "ok.") ; sys.stdout.flush()
            del nodes_id, h
            # For faces
            api.io.print('std', "\t Checking faces integrity .."); sys.stdout.flush()
            h = restartSectionHeader()
            if(not h.readVar(self.fid, self.byte_swap,["UGP_IO_FA_CHECK"])): exit()
            faces_id = np.empty((self.mesh["params"]["fa_count"],), dtype=np.int32)
            for loopi in range(self.mesh["params"]["fa_count"]):  # xrange to range (python3 portage)
                s = BinaryRead(self.fid, "i", self.byte_swap, type2nbytes["int32"])
                faces_id[loopi] = s[0]
            assert np.allclose(faces_id, np.arange(self.mesh["params"]["fa_count"]))
            api.io.print('std', "ok."); sys.stdout.flush()
            del faces_id, h
            # For cells
            api.io.print('std', "\t Checking cells integrity .."); sys.stdout.flush()
            h = restartSectionHeader()
            if(not h.readVar(self.fid, self.byte_swap,["UGP_IO_CV_CHECK"])): exit()
            assert h.idata[0] == self.mesh["params"]["cv_count"]
            assert h.id[0] == ic3_restart_codes["UGP_IO_CV_CHECK"]
            cells_id = np.empty((self.mesh["params"]["cv_count"],), dtype=np.int32)
            for loopi in range(self.mesh["params"]["cv_count"]):  # xrange to range (python3 portage)
                s = BinaryRead(self.fid, "i", self.byte_swap, type2nbytes["int32"])
                cells_id[loopi] = s[0]
            assert np.allclose(cells_id, np.arange(self.mesh["params"]["cv_count"]))
            api.io.print('std', "ok."); sys.stdout.flush()
            del cells_id, h

        # The two connectivities now
        #
        #- First, NOOFA
        #
        api.io.print('std', "\t Parsing face to node connectivity .."); sys.stdout.flush()
        h = restartSectionHeader()
        if (not h.readVar(self.fid, self.byte_swap,["UGP_IO_NOOFA_I_AND_V"])): exit()
        #
        assert h.idata[0] == self.mesh["params"]["fa_count"]
        assert h.idata[1] == self.mesh["params"]["noofa_count"]
        # Get the node count per face
        nno_per_face = np.empty((self.mesh["params"]["fa_count"],), dtype=np.int32)
        for loopi in range(self.mesh["params"]["fa_count"]): # xrange to range (python3 portage)
            s = BinaryRead(self.fid, "i", self.byte_swap, type2nbytes["int32"])
            nno_per_face[loopi] = s[0]
        print("nooface:",nno_per_face)
        uniq, counts = np.unique(nno_per_face, return_counts=True)
        print("uniq:",uniq)
        uniq = [nno2fatype[val] for val in uniq]
        api.io.print('std', "found %s faces .."%(" and ".join(uniq))); sys.stdout.flush()
        # Initialize the proper connectivity arrays in self.mesh
        self.mesh["connectivity"]["noofa"]["listofStarts_f2v"] = np.concatenate(([0,], np.cumsum(nno_per_face)))
        self.mesh["connectivity"]["noofa"]["face2vertex"] = np.zeros((np.sum(nno_per_face),), dtype=np.int64)
        assert self.mesh["connectivity"]["noofa"]["listofStarts_f2v"][0] == 0
        assert self.mesh["connectivity"]["noofa"]["listofStarts_f2v"][-2] == self.mesh["connectivity"]["noofa"]["face2vertex"].size - nno_per_face[-1]
        # Now loop on the restart file to fill the connectivities
        for loopi in range(self.mesh["params"]["fa_count"]):  # xrange to range (python3 portage)
            sta, sto = self.mesh["connectivity"]["noofa"]["listofStarts_f2v"][loopi], self.mesh["connectivity"]["noofa"]["listofStarts_f2v"][loopi+1]
            s = BinaryRead(self.fid, "i"*nno_per_face[loopi], self.byte_swap, type2nbytes["int32"]*nno_per_face[loopi])
            self.mesh["connectivity"]["noofa"]["face2vertex"][sta:sto] = np.asarray(s).astype(np.int64)
        api.io.print('std', "ok."); sys.stdout.flush()
        del nno_per_face, h, uniq, counts
        #
        #- Second, CVOFA
        #
        api.io.print('std', "\t Parsing face to cell connectivity .."); sys.stdout.flush()
        h = restartSectionHeader()
        if (not h.readVar(self.fid, self.byte_swap,["UGP_IO_CVOFA"])): exit()

        assert h.idata[0] == self.mesh["params"]["fa_count"]
        assert h.idata[1] == 2
        # Initialize the proper connectivity arrays in self.mesh
        self.mesh["connectivity"]["cvofa"]["cvofa"] = np.zeros((self.mesh["params"]["fa_count"], 2), dtype=np.int64)
        # Now loop on the restart file to fill the connectivities
        for loopi in range(self.mesh["params"]["fa_count"]):  # xrange to range (python3 portage)
            s = BinaryRead(self.fid, "ii", self.byte_swap, type2nbytes["int32"]*self.mesh["connectivity"]["cvofa"]["cvofa"].shape[1])
            self.mesh["connectivity"]["cvofa"]["cvofa"][loopi, :] = np.asarray(s).astype(np.int64)
        api.io.print('std', "ok."); sys.stdout.flush()
        del h
        # Checks and a few associations
        assert self.mesh["connectivity"]["cvofa"]["cvofa"].max() < self.mesh["params"]["cv_count"]
        assert self.mesh["connectivity"]["cvofa"]["cvofa"].max() == self.mesh["params"]["cv_count"]-1
        uniq, counts = np.unique(self.mesh["connectivity"]["cvofa"]["cvofa"][:,1], return_counts=True)
        # Number of assigned boundary faces
        try:
            iwhere = np.where(uniq==-1)[0][0]
        except IndexError:
            self.mesh["bocos"]["nfa_b"] = 0
        else:
            self.mesh["bocos"]["nfa_b"] = counts[iwhere]
        del uniq, counts
        # Number of periodic boundary faces
        self.mesh["bocos"]["nfa_bp"] = np.count_nonzero(self.mesh["connectivity"]["cvofa"]["cvofa"][:,1] < -1)
        # All periodic boundary faces to -1
        mask = self.mesh["connectivity"]["cvofa"]["cvofa"][:,1] < -1
        self.mesh["connectivity"]["cvofa"]["cvofa"][:,1][mask] = -1
        #
        #- Convert CVOFA to element2face
        #
        # api.io.print('std', "\t Converting CVOFA to an element2face connectivity .."); sys.stdout.flush()
        # self.mesh["connectivity"]["cvofa"]["cvofa"] += 1
        # uniq, unique_count = np.unique(self.mesh["connectivity"]["cvofa"]["cvofa"], return_counts=True)
        # mask = uniq > 0 # >= 0
        # element2face_size = np.sum(unique_count[mask]) + unique_count[mask].size
        # from getFACEEL import getfaceel
        # e2f = getfaceel(self.mesh["params"]["fa_count"], self.mesh["params"]["cv_count"],
        #                 element2face_size,
        #                 unique_count[mask],
        #                 self.mesh["connectivity"]["cvofa"]["cvofa"][:,0], self.mesh["connectivity"]["cvofa"]["cvofa"][:,1])
        # e2f = np.hstack( ( unique_count[mask][:, None], e2f ) )
        # e2f_ma = npma.masked_where(e2f==-1, e2f)
        # if e2f_ma.mask!=False:
        #     api.io.print('std', "Masking",e2f_ma.mask)
        #     element2face = e2f[~e2f_ma.mask].squeeze().ravel()
        # else:
        #     element2face = e2f.squeeze().ravel()

        # listofStarts_e2f = np.insert(np.cumsum(unique_count[mask]+1), 0, 0)[:-1]
        # e2f = np.delete(element2face, listofStarts_e2f)
        # self.mesh["connectivity"]["cvofa"]["element2face"] = e2f
        # self.mesh["connectivity"]["cvofa"]["element2face"] -= 1
        # self.mesh["connectivity"]["cvofa"]["listofStarts_e2f"] = listofStarts_e2f - np.arange(len(listofStarts_e2f))
        # self.mesh["connectivity"]["cvofa"]["listofStarts_e2f"] = np.concatenate((self.mesh["connectivity"]["cvofa"]["listofStarts_e2f"], [e2f.shape[0],]))
        # self.mesh["connectivity"]["cvofa"]["cvofa"] -= 1
        # api.io.print('std', "ok."); sys.stdout.flush()
        # del uniq, element2face_size, e2f, e2f_ma, element2face, listofStarts_e2f

        # The boundary conditions now
        api.io.print('std', "\t Parsing boundary conditions .."); sys.stdout.flush()
        while True:
            h = restartSectionHeader()
            if(not h.readVar(self.fid, self.byte_swap,["UGP_IO_FA_ZONE"],reset_offset=False)): break

            self.mesh["params"]["nboco"] += 1
            self.mesh["bocos"][h.name] = {}
            self.mesh["bocos"][h.name]["type"] = zonekind2type[h.idata[0]]
            self.mesh["bocos"][h.name]["fa_range"] = np.array([h.idata[1], h.idata[2]])
            self.mesh["bocos"][h.name]["periodic_transform"] = h.rdata
            #
            famin, famax = self.mesh["bocos"][h.name]["fa_range"]
            sta = self.mesh["connectivity"]["noofa"]["listofStarts_f2v"][famin]
            try:
                sto = self.mesh["connectivity"]["noofa"]["listofStarts_f2v"][famax+1]
            except IndexError:
                sto = self.mesh["connectivity"]["noofa"]["face2vertex"].size
            #
            self.mesh["bocos"][h.name]["slicing"] = np.unique(self.mesh["connectivity"]["noofa"]["face2vertex"][sta:sto])
            if h.idata[0] == 6:
                break
        api.io.print('standard', "ok.")
        sys.stdout.flush()

        # Parse the header of the partition information
        api.io.print('std', "\t Parsing partitioning information .."); sys.stdout.flush()
        h = restartSectionHeader()
        if(not h.readVar(self.fid, self.byte_swap,["UGP_IO_CV_PART"])): exit()

        self.mesh["partition"] = np.zeros((self.mesh["params"]["cv_count"],), dtype=np.int32)
        for loopi in range(self.mesh["params"]["cv_count"]):  # xrange to range (python3 portage)
            s = BinaryRead(self.fid, "i", self.byte_swap, type2nbytes["int32"])
            self.mesh["partition"][loopi] = s[0]
        api.io.print('std', "ok.")
        sys.stdout.flush()

        # The coordinates of the vertices finally
        api.io.print('std', "\t Parsing vertices coordinates .."); sys.stdout.flush()
        h = restartSectionHeader()
        if(not h.readVar(self.fid, self.byte_swap,["UGP_IO_X_NO"])): exit()

        assert h.idata[0] == self.mesh["params"]["no_count"]
        assert h.idata[1] == 3
        self.mesh["coordinates"] = np.zeros((self.mesh["params"]["no_count"], 3), dtype=np.float64)
        for loopi in range(self.mesh["params"]["no_count"]):  # xrange to range (python3 portage)
            s = BinaryRead(self.fid, "ddd", self.byte_swap, type2nbytes["float64"]*self.mesh["coordinates"].shape[1])
            self.mesh["coordinates"][loopi, :] = np.asarray(s)
        self.mesh["coordinates"] = np.ascontiguousarray(self.mesh["coordinates"])
        api.io.print('std', "ok."); sys.stdout.flush()

        # Convert the two connectivity tables to a classical element2vertex connectivity
        # api.io.print('std', "\t Converting element2face and face2vertex to an element2vertex connectivity .."); sys.stdout.flush()
        # from getNOoCV import getnoocv
        # los_e2v, e2v = getnoocv(self.mesh["params"]["cv_count"], self.mesh["params"]["no_count"], self.mesh["params"]["fa_count"],
        #                         self.mesh["connectivity"]["noofa"]["listofStarts_f2v"],
        #                         self.mesh["connectivity"]["cvofa"]["listofStarts_e2f"],
        #                         self.mesh["connectivity"]["noofa"]["face2vertex"],
        #                         self.mesh["connectivity"]["cvofa"]["element2face"],
        #                         self.mesh["connectivity"]["cvofa"]["cvofa"][:,0], self.mesh["connectivity"]["cvofa"]["cvofa"][:,1],
        #                         self.mesh["coordinates"][:,0], self.mesh["coordinates"][:,1], self.mesh["coordinates"][:,2])
        # _nvTotalperElt = np.diff(los_e2v)
        # nvTotalperElt = np.repeat(_nvTotalperElt, _nvTotalperElt)
        # uniqueNoV = np.unique(nvTotalperElt)
        # connectivities = [None]*len(uniqueNoV)
        # self.mesh["connectivity"]["e2v"] = collections.OrderedDict()
        # self.mesh["connectivity"]["cell_indices"] = collections.OrderedDict()
        # for idx, nov in enumerate(uniqueNoV):
        #     self.mesh["connectivity"]["nkeys"] += 1
        #     _maskNoV = _nvTotalperElt==nov
        #     maskNoV = nvTotalperElt==nov
        #     connectivities[idx] = e2v[maskNoV].reshape((-1, nov))
        #     el_type = cell_from_nodes[nov]
        #     self.mesh["connectivity"]["e2v"][el_type] = connectivities[idx]
        #     self.mesh["connectivity"]["cell_indices"][el_type] = np.arange(self.mesh["params"]["cv_count"])[_maskNoV]
        # api.io.print('std', "ok."); sys.stdout.flush()

    def __ReadInformativeValues(self):
        '''
        Method reading all the values also stored in a restart file,
        i.e. the step number, the time, the timestep.
        input   : handle on an open restart file, [type file identifier]
                  endianness flag [boolean]
        output  : simulation state structure containing informations about the current state
                  of the simulation
        '''

        # Initialize the state dictionary
        self.simulation_state = {"step":0, "dt":0, "time":0, "wgt":{}}

        # First, a header saying data to introduce to this block we are parsing now


        #removed as not used at the moment
        #h = restartSectionHeader()
        #if(not readVar(self.fid, self.byte_swap,"UGP_IO_DAT)A")
        #if (not varfound): exit()

        reset_offset=True
        while True:
            h = restartSectionHeader()
            if (not h.readVar(self.fid, self.byte_swap,["UGP_IO_I0"],reset_offset=reset_offset)): break
            else: reset_offset=False
            print("I0 var "+h.name.lower())
            self.simulation_state[h.name.lower()] = h.idata[0]
        reset_offset=True
        while True:
            h = restartSectionHeader()
            if (not h.readVar(self.fid, self.byte_swap,["UGP_IO_D0"],reset_offset=reset_offset)): break
            else: reset_offset=False
            print("D0 var "+h.name.lower())
            self.simulation_state[h.name.lower()] = h.rdata[0]

    def __ReadRestartVar(self):
        '''
        Method reading the variables from the restart file
        input   : handle on an open restart file, [type file identifier]
                  endianness flag [boolean]
                  mesh structure
        output  : structure containing all the variables
        '''

        # Some extra modules
        import copy

        # Initialize the variable dictionary
        self.variables = {"nodes":{}, "cells":{}, "faces":{}}

        # First come the scalars
        api.io.print('std', "\t First the scalars ..")
        reset_offset=True
        while True:
            h = restartSectionHeader()
            if(not h.readVar(self.fid, self.byte_swap,["UGP_IO_NO_D1","UGP_IO_FA_D1","UGP_IO_CV_D1"],reset_offset=reset_offset)): break
            reset_offset=False

            if h.idata[0] == self.mesh["params"]["no_count"]:
                self.variables["nodes"][h.name] = np.zeros((self.mesh["params"]["no_count"],), dtype=np.float64)
                s = BinaryRead(self.fid, "d"*self.mesh["params"]["no_count"], self.byte_swap, type2nbytes["float64"]*self.mesh["params"]["no_count"])
                self.variables["nodes"][h.name] = np.asarray(s)
                api.io.print('std', "\t %s%s:\t %+.5e / %+.5e / %+.5e (min/mean/max)."%(h.name, ' '*(20-len(h.name)), np.asarray(s).min(), np.mean(np.asarray(s)), np.asarray(s).max()))
            elif h.idata[0] == self.mesh["params"]["fa_count"]:
                self.variables["faces"][h.name] = np.zeros((self.mesh["params"]["fa_count"],), dtype=np.float64)
                s = BinaryRead(self.fid, "d"*self.mesh["params"]["fa_count"], self.byte_swap, type2nbytes["float64"]*self.mesh["params"]["fa_count"])
                self.variables["faces"][h.name] = np.asarray(s)
                api.io.print('std', "\t %s%s:\t %+.5e / %+.5e / %+.5e (min/mean/max)."%(h.name, ' '*(20-len(h.name)), np.asarray(s).min(), np.mean(np.asarray(s)), np.asarray(s).max()))
            elif h.idata[0] == self.mesh["params"]["cv_count"]:
                self.variables["cells"][h.name] = np.zeros((self.mesh["params"]["cv_count"],), dtype=np.float64)
                s = BinaryRead(self.fid, "d"*self.mesh["params"]["cv_count"], self.byte_swap, type2nbytes["float64"]*self.mesh["params"]["cv_count"])
                self.variables["cells"][h.name] = np.asarray(s)
                # If multiple connectivities, gotta order the tables correctly
                if self.mesh["connectivity"]["nkeys"] > 1:
                    aux = copy.deepcopy(self.variables["cells"][h.name])
                    self.variables["cells"][h.name] = np.empty((0,), dtype=aux.dtype)
                    for uns_type, indices in self.mesh["connectivity"]["cell_indices"].items():
                        self.variables["cells"][h.name] = np.concatenate((self.variables["cells"][h.name],
                                                                          aux[indices]))
                api.io.print('std', "\t %s%s:\t %+.5e / %+.5e / %+.5e (min/mean/max)."%(h.name, ' '*(20-len(h.name)), np.asarray(s).min(), np.mean(np.asarray(s)), np.asarray(s).max()))
            else:
                api.io.print('std', "Fatal error. Incoherence in dataset %s. Exiting."%(h.name))
                exit()
        api.io.print('std', "\t ok.")
        

        # Then the vectors
        api.io.print('std', "\t Then the vectors ..")
        reset_offset=True
        while True:
            h = restartSectionHeader()
            if(not h.readVar(self.fid, self.byte_swap,["UGP_IO_NO_D3","UGP_IO_FA_D3","UGP_IO_CV_D3"],reset_offset=reset_offset)): break
            reset_offset=False

            if h.idata[0] == self.mesh["params"]["no_count"]:
                self.variables["nodes"][h.name] = np.zeros((self.mesh["params"]["no_count"], 3), dtype=np.float64)
                s = BinaryRead(self.fid, "ddd"*self.mesh["params"]["no_count"], self.byte_swap, type2nbytes["float64"]*self.mesh["params"]["no_count"]*3)
                self.variables["nodes"][h.name] = np.asarray(s).reshape((-1, 3))
                api.io.print('std', "\t %s%s:\t %+.5e / %+.5e / %+.5e (min/mean/max)."%(h.name, ' '*(20-len(h.name)), np.asarray(s).min(), np.mean(np.asarray(s)), np.asarray(s).max()))
            elif h.idata[0] == self.mesh["params"]["fa_count"]:
                self.variables["faces"][h.name] = np.zeros((self.mesh["params"]["fa_count"], 3), dtype=np.float64)
                s = BinaryRead(self.fid, "ddd"*self.mesh["params"]["fa_count"], self.byte_swap, type2nbytes["float64"]*self.mesh["params"]["fa_count"]*3)
                self.variables["faces"][h.name] = np.asarray(s).reshape((-1, 3))
                api.io.print('std', "\t %s%s:\t %+.5e / %+.5e / %+.5e (min/mean/max)."%(h.name, ' '*(20-len(h.name)), np.asarray(s).min(), np.mean(np.asarray(s)), np.asarray(s).max()))
            elif h.idata[0] == self.mesh["params"]["cv_count"]:
                self.variables["cells"][h.name] = np.zeros((self.mesh["params"]["cv_count"], 3), dtype=np.float64)
                s = BinaryRead(self.fid, "ddd"*self.mesh["params"]["cv_count"], self.byte_swap, type2nbytes["float64"]*self.mesh["params"]["cv_count"]*3)
                self.variables["cells"][h.name] = np.asarray(s).reshape((-1, 3))
                # If multiple connectivities, gotta order the tables correctly
                if self.mesh["connectivity"]["nkeys"] > 1:
                    aux = copy.deepcopy(self.variables["cells"][h.name])
                    self.variables["cells"][h.name] = np.empty((0, 3), dtype=aux.dtype)
                    for uns_type, indices in self.mesh["connectivity"]["cell_indices"].items():
                        self.variables["cells"][h.name] = np.concatenate((self.variables["cells"][h.name],
                                                                          aux[indices, :]), axis=0)
                api.io.print('std', "\t %s%s:\t %+.5e / %+.5e / %+.5e (min/mean/max)."%(h.name, ' '*(20-len(h.name)), np.asarray(s).min(), np.mean(np.asarray(s)), np.asarray(s).max()))
            else:
                api.io.print('std', "Fatal error. Incoherence in dataset %s. Exiting."%(h.name))
                exit()
        api.io.print('std', "\t ok.")
        

        # Then the tensors
        api.io.print('std', "\t Then the tensors ..")
        reset_offset=True
        while True:
            h = restartSectionHeader()
            if(not h.readVar(self.fid, self.byte_swap,["UGP_IO_CV_D33"],reset_offset=reset_offset)): break
            reset_offset=False

            if h.idata[0] == self.mesh["params"]["cv_count"]:
                self.variables["cells"][h.name] = np.zeros((self.mesh["params"]["cv_count"], 3, 3), dtype=np.float64)
                s = BinaryRead(self.fid, "d"*self.mesh["params"]["cv_count"]*3*3, self.byte_swap, type2nbytes["float64"]*self.mesh["params"]["cv_count"]*3*3)
                self.variables["cells"][h.name] = np.asarray(s).reshape((-1, 3, 3))
                # If multiple connectivities, gotta order the tables correctly
                if self.mesh["connectivity"]["nkeys"] > 1:
                    aux = copy.deepcopy(self.variables["cells"][h.name])
                    self.variables["cells"][h.name] = np.empty((0, 3, 3), dtype=aux.dtype)
                    for uns_type, indices in self.mesh["connectivity"]["cell_indices"].items():
                        self.variables["cells"][h.name] = np.concatenate((self.variables["cells"][h.name],
                                                                          aux[indices, :, :]), axis=0)
                api.io.print('std', "\t %s%s:\t %+.5e / %+.5e / %+.5e (min/mean/max)."%(h.name, ' '*(20-len(h.name)), np.asarray(s).min(), np.mean(np.asarray(s)), np.asarray(s).max()))
            else:
                api.io.print('std', "Fatal error. Incoherence in dataset %s. Exiting."%(h.name))
                exit()
        api.io.print('std', "\t ok.")
        

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

#     # Module import for script use
#     import copy
#     import os
#     import sys
#     import vtk as _vtk
#     from vtk.util import numpy_support

#     # Parse command line arguments seeking a restart file name
#     defaultFName = "restart_test.out"
#     if len(sys.argv) == 1:
#         fName = defaultFName
#         cIntegrity = False
#         api.io.print('std', "Use requires the name of a restart file to be input in argument.")
#         api.io.print('std', "Using a default name, %s."%(defaultFName))
#     else:
#         fName = defaultFName
#         cIntegrity = False
#         for arg in sys.argv[1:]:
#             if '--restartname' in arg:
#                 fName = arg.split('=')[-1]
#             elif '--checkIntegrity' in arg:
#                 _cIntegrity = arg.split('=')[-1]
#                 if _cIntegrity == "False":
#                     cIntegrity = False
#                 elif _cIntegrity == "True":
#                     cIntegrity = True
#     if fName == defaultFName:
#         api.io.print('std', "Use requires the name of a restart file to be input in argument.")
#         api.io.print('std', "Using a default name, %s."%(defaultFName))

#     # Check given file name
#     if not os.path.isfile(fName):
#         api.io.print('std', "Fatal error: %s cannot be found. Exiting."%(fName))
#         exit()

#     # Else proceed
#     readr = ReaderRestartIC3(fName, cIntegrity)
#     xyz, co, bocos, simulation_state, nodesvar, cellsvar, params = readr.read_data()

#     ###################################################################################################

#     # Now, we can export the data under a more human-friendly format
#     #
#     #- IN VTK
#     #
#     # Build point coordinates for VTK
#     nb_point = xyz.shape[0]
#     if hasattr(_vtk, 'vtkSOADataArrayTemplate'):
#         # VTK 8.1.0+: use SOA coordinate array data structure (zero-copy)
#         points_vtk = _vtk.vtkSOADataArrayTemplate[np.float64]()  # shall avoid copy
#         points_vtk.SetNumberOfComponents(3)
#         points_vtk.charlesx_arrays = []
#         for index in range(3):
#             np_array = copy.deepcopy(xyz[:, index])
#             points_vtk.SetArray(index, np_array, nb_point, True, True)
#             points_vtk.charlesx_arrays.append(np_array)  # record array ref so it won't be gc'ed
#     else:
#         # VTK 8.0.1-: build a new AOS coordinate array
#         points_vtk = numpy_support.numpy_to_vtk(mesh["coordinates"], deep=True)
#     vtkpoints = _vtk.vtkPoints()
#     vtkpoints.SetData(points_vtk)
#     # Build connectivity for VTK
#     vtkprimitive = {
#         'bi': _vtk.vtkLine(),
#         'tri': _vtk.vtkTriangle(),
#         'qua': _vtk.vtkQuad(),
#         'tet': _vtk.vtkTetra(),
#         'hex': _vtk.vtkHexahedron(),
#         'pri': _vtk.vtkWedge(),
#         'pyr': _vtk.vtkPyramid(),
#     }
#     cells = np.empty((0, ), dtype=np.int64)
#     cell_types = np.empty((0, ), dtype=np.int64)
#     offsets = np.empty((0, ), dtype=np.int64)
#     offset_start = 0
#     total_nb_cells = 0

#     for uns_type, connect in co.items():

#         nb_cells = connect.shape[0]

#         # record cell-type for each individual cell
#         if uns_type in vtkprimitive:
#             cell_type = vtkprimitive[uns_type].GetCellType()
#         else:
#             raise ValueError(str_error('Unknown cell type'))
#         cell_types = np.append(cell_types, np.tile(cell_type, (nb_cells, 1)))

#         # put number of vertices before each cell
#         cells = np.append(cells, np.concatenate((np.tile(nodes_per_cell[uns_type], (nb_cells, 1)), connect),
#                                                 axis=1).flat)

#         # start offset of each cell in connectivity array
#         offset_stop = offset_start + nb_cells * (nodes_per_cell[uns_type] + 1)
#         offsets = np.append(offsets, np.arange(offset_start, offset_stop, nodes_per_cell[uns_type] + 1))
#         offset_start = offset_stop
#         total_nb_cells += nb_cells
#     idtype_vtk = _vtk.vtkIdTypeArray().GetDataType()
#     uchartype_vtk = _vtk.vtkUnsignedCharArray().GetDataType()

#     cells_vtk = numpy_support.numpy_to_vtk(cells, deep=True, array_type=idtype_vtk)
#     cell_array = _vtk.vtkCellArray()
#     cell_array.SetCells(total_nb_cells, cells_vtk)
#     # Build VTK unstructured Grid
#     if len(cell_types) == 1:
#         vtk_obj = _vtk.vtkUnstructuredGrid()
#         vtk_obj.SetPoints(vtkpoints)
#         vtk_obj.SetCells(cell_types[0], cell_array)
#     else:
#         cell_types_vtk = numpy_support.numpy_to_vtk(cell_types, deep=True, array_type=uchartype_vtk)
#         offsets_vtk = numpy_support.numpy_to_vtk(offsets, deep=True, array_type=idtype_vtk)
#         vtk_obj = _vtk.vtkUnstructuredGrid()
#         vtk_obj.SetPoints(vtkpoints)
#         vtk_obj.SetCells(cell_types_vtk, offsets_vtk, cell_array)

#     elemdata = {
#         "nodes": vtk_obj.GetPointData(),
#         "cells": vtk_obj.GetCellData(),
#     }
#     # And add the variables
#     for var in cellsvar:
#         np_array = cellsvar[var]
#         # cast without copy if possible
#         if np_array.flags.contiguous:
#             np_array = np_array.astype(np.float64, copy=False)
#         else:
#             np_array = np_array.astype(np.float64)
#         # zero-copy if contiguous array
#         vtkarray = numpy_support.numpy_to_vtk(np_array, deep=False)
#         vtkarray.charlesx_array = np_array  # avoid garbage collection
#         vtkarray.SetName(var)
#         #
#         elemdata["cells"].AddArray(vtkarray)
#     for var in nodesvar:
#         np_array = nodesvar[var]
#         # cast without copy if possible
#         if np_array.flags.contiguous:
#             np_array = np_array.astype(np.float64, copy=False)
#         else:
#             np_array = np_array.astype(np.float64)
#         # zero-copy if contiguous array
#         vtkarray = numpy_support.numpy_to_vtk(np_array, deep=False)
#         vtkarray.charlesx_array = np_array  # avoid garbage collection
#         vtkarray.SetName(var)
#         #
#         elemdata["nodes"].AddArray(vtkarray)
#     # Finally, write the file
#     writer = _vtk.vtkXMLUnstructuredGridWriter()
#     if _vtk.vtkVersion.GetVTKMajorVersion() >= 6:
#         writer.SetInputData(vtk_obj)
#     else:
#         writer.SetInput(vtk_obj)
#     filename = "restart.%d"%(simulation_state["step"])
#     writer.SetFileName(filename + '.vtu')
#     writer.SetDataModeToBinary()
#     writer.Write()



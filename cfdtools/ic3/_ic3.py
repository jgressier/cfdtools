import struct
import os
import numpy as np
import cfdtools.api as api

# Copy the restart codes from CharlesX to match its restart routine
ic3_restart_codes = {"UGP_IO_MAGIC_NUMBER":123581321,
                     "UGP_IO_VERSION":-1, # not specific here
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

properties_ugpcode = {
    ic3_restart_codes["UGP_IO_FA_D1"]:  { 'structcode': "d", 'size': 8, 'numpytype': np.float64 },
    ic3_restart_codes["UGP_IO_FA_D3"]:  { 'structcode': "d", 'size': 8, 'numpytype': np.float64 },
    ic3_restart_codes["UGP_IO_NO_D1"]:  { 'structcode': "d", 'size': 8, 'numpytype': np.float64 },
    ic3_restart_codes["UGP_IO_NO_D3"]:  { 'structcode': "d", 'size': 8, 'numpytype': np.float64 },
    ic3_restart_codes["UGP_IO_CV_D1"]:  { 'structcode': "d", 'size': 8, 'numpytype': np.float64 },
    ic3_restart_codes["UGP_IO_CV_D3"]:  { 'structcode': "d", 'size': 8, 'numpytype': np.float64 },
    ic3_restart_codes["UGP_IO_CV_D33"]: { 'structcode': "d", 'size': 8, 'numpytype': np.float64 },
    ic3_restart_codes["UGP_IO_FA_II1"]:  { 'structcode': "q", 'size': 8, 'numpytype': np.int64 },
    ic3_restart_codes["UGP_IO_NO_II1"]:  { 'structcode': "q", 'size': 8, 'numpytype': np.int64 },
    ic3_restart_codes["UGP_IO_CV_II1"]:  { 'structcode': "q", 'size': 8, 'numpytype': np.int64 },
}

# Dictionary to convert type of data [string] to a number of bytes for clean binary parsing
type2nbytes = {"char":1,
               "int32":4,
               "int64":8,
               "float32":4,
               "float64":8,
               "float128":16,}
# Dictionary to convert number of nodes per face to type of face
nno2fatype = {2:"line",
              3:"tri",
              4:"qua",}

fatype2nno = {"line":2,
              "tri":3,
              "qua":4,}
# Actual number of vertices for a given cell type
nodes_per_cell = {
    'bi': 2,
    'tri': 3,
    'qua': 4,
    'tet': 4,
    'hex': 8,
    'pri': 6,
    'pyr': 5,
}

# cell_from_nodes = {
#     2:'bi',
#     3:'tri',
#     4:'qua',
#     4:'tet',
#     8:'hex',
#     6:'pri',
#     5:'pyr',
# }

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
# Dictionary to convert zone kind to zone type (as a string)
zonekind2type = {1:"boundary",
                 2:"periodic_cart",
                 3:"periodic_cylx",
                 4:"periodic_cyly",
                 5:"periodic_cylz",
                 6:"internal"}

struct_endian = { 'native':'@', 'little':'<', 'big':'>'}

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

    def readVar(self, bfile, byte_swap, nametypes, reset_offset=True):
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
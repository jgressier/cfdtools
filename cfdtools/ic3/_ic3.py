import logging
import struct
import os

import numpy as np

import cfdtools.api as api

log = logging.getLogger(__name__)

# Copy the restart codes from CharlesX to match its restart routine
ic3_restart_codes = {
    "UGP_IO_MAGIC_NUMBER": 123581321,
    "UGP_IO_VERSION": -1,  # not specific here
    "UGP_IO_HEADER_NAME_LEN": 52,  # select this to make header size 256 bytes
    "UGP_IO_NO_FA_CV_NOOFA_COUNTS": 10,
    "UGP_IO_I0": 11,
    "UGP_IO_D0": 12,
    "UGP_IO_FA_CHECK": 20,
    "UGP_IO_NOOFA_I_AND_V": 21,
    "UGP_IO_CVOFA": 22,
    "UGP_IO_FA_ZONE": 23,
    "UGP_IO_FA_D1": 24,
    "UGP_IO_FA_D3": 25,
    "UGP_IO_FA_II1": 26,
    "UGP_IO_NO_CHECK": 30,
    "UGP_IO_X_NO": 31,
    "UGP_IO_NO_D1": 32,
    "UGP_IO_NO_D3": 33,
    "UGP_IO_NO_II1": 34,
    "UGP_IO_CV_CHECK": 40,
    "UGP_IO_CV_PART": 41,
    "UGP_IO_CV_D1": 42,
    "UGP_IO_CV_D3": 43,
    "UGP_IO_CV_D33": 44,
    "UGP_IO_CV_II1": 45,
    "UGP_IO_DATA": 50,
    "UGP_IO_EOF": 51,
}

scode, sz, ntype = 'structcode', 'size', 'numpytype'

properties_ugpcode = {
    ic3_restart_codes["UGP_IO_FA_D1"]: {scode: "d", sz: 8, ntype: np.float64},
    ic3_restart_codes["UGP_IO_FA_D3"]: {scode: "d", sz: 8, ntype: np.float64},
    ic3_restart_codes["UGP_IO_NO_D1"]: {scode: "d", sz: 8, ntype: np.float64},
    ic3_restart_codes["UGP_IO_NO_D3"]: {scode: "d", sz: 8, ntype: np.float64},
    ic3_restart_codes["UGP_IO_CV_D1"]: {scode: "d", sz: 8, ntype: np.float64},
    ic3_restart_codes["UGP_IO_CV_D3"]: {scode: "d", sz: 8, ntype: np.float64},
    ic3_restart_codes["UGP_IO_CV_D33"]: {scode: "d", sz: 8, ntype: np.float64},
    ic3_restart_codes["UGP_IO_FA_II1"]: {scode: "q", sz: 8, ntype: np.int64},
    ic3_restart_codes["UGP_IO_NO_II1"]: {scode: "q", sz: 8, ntype: np.int64},
    ic3_restart_codes["UGP_IO_CV_II1"]: {scode: "q", sz: 8, ntype: np.int64},
}

del scode, sz, ntype

# Dictionary to convert type of data [string] to a number of bytes for clean binary parsing
type2nbytes = {
    "char": 1,
    "int32": 4,
    "int64": 8,
    "float32": 4,
    "float64": 8,
    "float128": 16,
}


# fatype2nno = {"line":2,
#               "tri":3,
#               "qua":4,}
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

# faces_per_cell = {
#     'hex': 6,
#     'pri': 5,
#     'pyr': 5,
#     'tet': 4,}

# faces_of_cell = {
#     'hex':['qua']*6,
#     'pri':['tri', 'qua', 'qua', 'qua', 'tri'],
#     'pyr':['qua'] + ['tri']*4,
#     'tet':['tri']*4,}

# ifaces_of_cell = {
#     'hex':[[0,3,2,1], [0,1,5,4], [1,2,6,5], [2,3,7,6], [3,0,4,7], [4,5,6,7]],
#     'pri':[[0,2,1], [0,1,4,3], [1,2,5,4], [2,0,3,5], [3,4,5]],
#     'pyr':[[0,3,2,1], [0,1,4], [1,2,4], [2,3,4], [3,0,4]],
#     'tet':[[0,2,1], [0,1,3], [1,2,3], [2,0,3]],}

# #define FA_ZONE_UNKNOWN          -1
# #define FA_ZONE_PERIODIC_UNKNOWN -2

# #define FA_ZONE_BOUNDARY          1

# #define FA_ZONE_PERIODIC_CART     2
# #define FA_ZONE_PERIODIC_CYL_X    3
# #define FA_ZONE_PERIODIC_CYL_Y    4
# #define FA_ZONE_PERIODIC_CYL_Z    5

# #define FA_ZONE_INTERNAL          6

# // range for all and periodic zones (tommie needs this?)...
# #define FA_ZONE_FIRST             1
# #define FA_ZONE_LAST              6
# #define FA_ZONE_PERIODIC_FIRST    2
# #define FA_ZONE_PERIODIC_LAST     5
type2zonekind = {
    "boundary": 1,
    "perio_cart": 2,
    "perio_cylx": 3,
    "perio_cyly": 4,
    "perio_cylz": 5,
    "internal": 6,
}
# Dictionary to convert zone kind to zone type (as a string)
zonekind2type = {itype: type for type, itype in type2zonekind.items()}

struct_endian = {'native': '@', 'little': '<', 'big': '>'}


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
    if byte_swap:  # little-endian
        form = '<' + form
    else:  # big-endian
        form = '>' + form

    # Create a "packed-binary reader structure"
    s = struct.Struct(form)

    # Try to read 'size' bytes from the file otherwise crash
    try:
        while True:
            record = bfile.read(size)
            if len(record) != size:
                log.error(
                    "mismatched record ({}) and expected ({}) sizes".format(
                        len(record),
                        size,
                    ),
                )
                break
            return s.unpack(record)
    except IOError:
        api.error_stop(f"Fatal error. Could not read {size} bytes from {bfile.name!r}. Exiting.")


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
    # Actually pack the string now
    packed_ba = s.pack(*varargs)

    # Try to write the bytes to the file otherwise crash
    try:
        bfile.write(packed_ba)
    except IOError:
        api.error_stop("Fatal error. Could not write to {bfile.name!r}. Exiting.")


class restartSectionHeader:
    '''
    This class is designed to handle the header that is present
    before every variable/section/category saved into the restart file.
    '''

    def __init__(self, skip=0):
        """
        Initialize a section header class, as it is always structured
        in the same way, i.e. with a name, an id, a skip bytes count,
        an information array, and an auxiliary array needed in specific
        cases like periodicity.
        """
        header_name_length = ic3_restart_codes["UGP_IO_HEADER_NAME_LEN"]
        self.name = " " * header_name_length
        self.id = np.zeros((1,), dtype=np.int32)
        self._skip = np.array([skip], dtype=np.int64)
        self.idata = np.zeros((8,), dtype=np.int64)
        self.rdata = np.zeros((16,), dtype=np.float64)

        # Initialize the format string
        self.binstr = ""
        # Fill it with prior knowledge of its content
        self.binstr += "c" * header_name_length  # name
        self.binstr += "i"  # id
        self.binstr += "q"  # skip
        self.binstr += "qqqqqqqq"  # idata (long long)
        self.binstr += "dddddddddddddddd"  # rdata (double)

        # Keep in mind the total byte size of the header
        self.hsize = (
            header_name_length * type2nbytes["char"]
            + type2nbytes["int32"]
            + type2nbytes["int64"]
            + self.idata.size * type2nbytes["int64"]
            + self.rdata.size * type2nbytes["float64"]
        )

    def skip(self):
        return self._skip[0]  # numpy array size 1 to handle int type

    def readVar(self, bfile, byte_swap, nametypes, reset_offset=True, required=False):
        """
        Once initialization is done, this method actually reads
        the header data from the packed binary formatted string.
        """
        id_list = [ic3_restart_codes[nametype] for nametype in nametypes]
        if reset_offset:
            bfile.seek(8, os.SEEK_SET)
        while True:
            self.name = ""
            self.id = np.zeros((1,), dtype=np.int32)
            self.idata = np.zeros((8,), dtype=np.int64)
            self.rdata = np.zeros((16,), dtype=np.float64)

            log.debug("skipping data %d", self.skip())
            if self.skip() > 0:
                bfile.seek(self.skip() - 256, os.SEEK_CUR)

            self._skip = np.zeros((1,), dtype=np.int64)

            # Split the packed binary string into its tokens
            s = list(BinaryRead(bfile, self.binstr, byte_swap, self.hsize))

            # Remove trailing space from header name (h.name)
            i = 0
            while (s[i] != b'\x00') and (
                i < ic3_restart_codes["UGP_IO_HEADER_NAME_LEN"]
            ):  # remove trailing spaces
                self.name += s[i].decode()  # .decode() for python3 portage
                i += 1
            if self.name.startswith("DEOF"):
                break
            # Store the rest of the tokens in the right namespace
            nlen = ic3_restart_codes["UGP_IO_HEADER_NAME_LEN"]
            self.id[0] = s[nlen]
            self._skip[0] = s[nlen + self.id.size]  # store size in skip for next read
            ibeg = nlen + self.id.size + self._skip.size
            iend = ibeg + self.idata.size
            self.idata = s[ibeg:iend]
            ibeg = iend
            iend = ibeg + self.rdata.size
            self.rdata = s[ibeg:iend]
            if self.id[0] in id_list:
                return True
            if self.id[0] == ic3_restart_codes["UGP_IO_EOF"]:
                break
        if required:
            api.error_stop(f"Fatal error. Section(s) {nametypes} not found in dataset {self.name}. Exiting.")
        return False

    def readReqVar(self, *args, **kwargs):
        """
        error_stop if required data is not found.
        """
        return self.readVar(*args, **kwargs, required=True)

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
        varargs.append(self.id)
        varargs.append(self.skip)
        for kk in range(8):
            varargs.append(self.idata[kk])
        for kk in range(16):
            varargs.append(self.rdata[kk])
        # Now write everything at once
        BinaryWrite(bfile, endian, self.binstr, varargs)

    def __str__(self):
        mystring = "> Header:"
        mystring += "\n"
        mystring += "Name : %s" % self.name
        mystring += "\n"
        mystring += f"Id   : {self.id} {list(dict(filter(lambda items: items[1] == self.id[0], ic3_restart_codes.items())).keys())}" 
        mystring += "\n"
        mystring += "hsize: %i" % self.hsize
        mystring += "\n"
        mystring += "skip : %i" % self.skip()
        mystring += "\n"
        mystring += "idata: (" + ", ".join(str(i) for i in self.idata) + ")"
        mystring += "\n"
        mystring += "rdata: (" + ", ".join(str(r) for r in self.rdata) + ")"
        return mystring


class binreader(api._files):
    '''Implementation of the reader to read IC3 restart files.'''

    # def __init__(self, filename):
    #     '''
    #     Initialization of an IC3 restart reader.
    #     Just save the filename and a boolean for an integrity check.
    #     input   : IC3 restart file name [type string]
    #               whether to check the integrity of the file beforehand [type boolean]
    #     '''
    #     super().__init__(filename)
    #     self.check_integrity = cIntegrity
    #     self.ic3_version = -1

    def read_headers(self):
        """
        Main method of the IC3 restart reader.
        Parses in order the file using sub-methods described below.
        """
        log.info("READER RESTART IC3 - only headers")

        if not self.exists():
            raise FileNotFoundError("Fatal error. File %s cannot be found." % (self.filename))

        # Open the file for binary reading
        log.debug('opening %s', self.filename)
        with open(self.filename, "rb") as self.fid:
            log.info("reading header (first section)")
            self._ReadRestartHeader()
            #
            reset_offset = True
            skip = 0
            while True:
                h = restartSectionHeader(skip=skip)
                if not h.readVar(
                    self.fid,
                    self.byte_swap,
                    ic3_restart_codes.keys(),
                    reset_offset=reset_offset,
                ):
                    break
                reset_offset = False  # continue
                skip = h.skip()
                print(h)
                if h.id[0] == ic3_restart_codes['UGP_IO_EOF']:
                    log.info('UGP EOF reached')
                    break
        log.debug('%s closed', self.filename)
        del self.fid

        return

    def _ReadRestartHeader(self):
        '''
        Method reading the header of a restart file.
        It is composed of two integers, the "magic number" used as a flag for endianness
        and the IC3 version number.
        input:  handle on an open restart file, [type file identifier]
        output: the endianness of the open restart file [type boolean]
        '''
        # By default suppose big-endian format
        self.byte_swap = False

        # Read the first integer (int64)
        s = list(BinaryRead(self.fid, "ii", False, 2 * type2nbytes["int32"]))
        # If, with big-endian assumption, the first integer comes out wrong, swap to little-endian
        if s[0] != ic3_restart_codes["UGP_IO_MAGIC_NUMBER"]:
            # Change the flag
            self.byte_swap = True
            # Transform the second integer of the list to match the version number
            aux_struct = struct.Struct(">i")
            packed_version = aux_struct.pack(s[1])
            del aux_struct
            aux_struct = struct.Struct("<i")
            s[1] = aux_struct.unpack(packed_version)[0]

        # Some info for the user
        self.ic3_version = s[1]

        log.info(
            f"  version: {self.ic3_version} " + ("little-endian" if self.byte_swap else "big-endian"),
        )

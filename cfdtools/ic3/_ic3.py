
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
cell_from_nodes = {
    2:'bi',
    3:'tri',
    4:'qua',
    4:'tet',
    8:'hex',
    6:'pri',
    5:'pyr',
}
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


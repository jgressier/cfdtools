import cfdtools.api as api
import cfdtools.meshbase._elements as ele
from collections import defaultdict, OrderedDict
from itertools import groupby
import numpy as np

class singleindex():
    """class of different implementation of list of index

    Returns:
        _type_: _description_
    """
    __available_types = ['list', 'range']

    def __init__(self):
        self._type = None

    def _delete(self):
        self._list = []
        self._range = []

    def set_list(self, ilist):
        self._delete()
        self._type = 'list'
        self._list = ilist

    def set_range(self, imin, imax):
        self._delete()
        self._type = 'range'
        self._range = [imin, imax]

    def range(self):
        if self._type == 'range':
            return self._range
        else:
            api.error_stop("unable to get range from list connectivity")

    def list(self):
        if self._type == 'range':
            return list(range(self._range[0], self._range[1]+1))
        elif self._type == 'list':
            return self._list
        else:
            api.error_stop("unknown type")

class doubleindex():
    """class for a genuine connectivity, regular which size = nelem x ndim 
    """
    def __init__(self, nelem=0, dim=0):
        self._conn = None
        self._nelem = nelem
        self._dim   = dim

    @property
    def conn(self):
        return self._conn

    @conn.setter
    def conn(self, array):
        self._conn = array
        self._elem = array.shape[0]
        self._dim  = array.shape[1]

    def append(self, array):
        if self._conn is None:
            self.conn = array
        else:
            assert(array.shape[1]==self._dim)
            self._conn.append(array, axis=0)
            self._nelem += array.shape[0]

class compressed_doubleindex():
    """class for a compressed double index (CSR like)

    Returns:
        _type_: _description_
    """
    def __init__(self, index, value) -> None:
        self._index = index
        self._value = value

class elem_connectivity():
    def __init__(self):
        self._nelem = 0
        self._elem2node = OrderedDict()

    def add_elems(self, etype: str, elem2node: np.ndarray):
        #print("elem2node", elem2node)
        self._elem2node[etype] = { 'starts' : self._nelem, 'elem2node' : elem2node }
        self._nelem += elem2node.shape[0]

    def index_cell2node(self, etype):
        ist = self._elem2node[etype]['starts']
        cell2node = self._elem2node[etype]['cell2node']
        return range(ist, ist+cell2node.shape[0]), cell2node

    def nelem(self):
        return self._nelem

    def print(self):
        for celltype, elemco in self._elem2node.items():
            api.io.print("std", f"  {celltype}: {elemco['elem2node'].shape}")

    def create_faces_from_elems(self):
        """build a dict of face type to a list of tuples of each (oriented) face and its neighbor

        Args:
            elems (dict): dict of elements with node definition

        Returns:
            _type_: dict of face type 
        """
        faces_neighbor = defaultdict(list)
        for elemtype, elemsdict in self._elem2node.items(): # elemtype: 'hexa8', elemsarray: ndarray[nelem,8]
            istart = elemsdict['starts']
            elemsarray = elemsdict['elem2node']
            for ielem in range(elemsarray.shape[0]):
                for ftype, listfaces in ele.elem2faces[elemtype].items():
                    for face in listfaces:
                        faces_neighbor[ftype].append( (tuple(elemsarray[ielem, face]), istart+ielem) ) 
        return faces_neighbor

def find_duplicates(faces_neighbor: dict):
    """find duplicated faces and build unique face/node face/cell connectivity

    Args:
        faces_neighbor (dict): _description_
    """
    internalfaces = elem_connectivity()
    boundaryfaces = elem_connectivity()
    iface2cell = doubleindex(0,2) # 
    bface2cell = doubleindex(0,2) # 

    def face_from_ufacedict(uface_dict):
        return np.array(list(map(lambda flist: flist[0][0], uface_dict.values())))

    for ftype, listfaces in faces_neighbor.items(): # find pairs for a given face type
        nf_all = len(listfaces)
        # build a dict of "sorted node" face with list of tuple (face, elem)
        # face_pairs = dict()
        # for uface, facepair in groupby(listfaces, lambda tup: tuple(sorted(tup[0]))):
        #     face_pairs[uface] = list(facepair)
        face_pairs = defaultdict(list)
        for tface in listfaces:
            face_pairs[tuple(sorted(tface[0]))].append(tface)            
        nf_unique = len(face_pairs)
        assert (nf_unique < nf_all)
        assert (2*nf_unique >= nf_all)
        # extract all unique face
        uniqueface_dict = dict(filter(lambda tup: len(tup[1])==1, face_pairs.items())) # tup[1] is the value of key
        #for k,t in face_pairs.items():
        #    print(k,t,len(t))
        #print(dict(map(lambda tup: (tup[0],len(tup[1])), uniqueface_dict.items())))
        boundaryfaces.add_elems(ftype, face_from_ufacedict(uniqueface_dict))
        f2c = np.full((len(uniqueface_dict),2), -1)
        # get index of connected cells
        f2c[:,0] = list(map(lambda flist: flist[0][1], uniqueface_dict.values()))
        bface2cell.append(f2c)
        # remove these faces
        print(len(uniqueface_dict),'/', len(face_pairs))
        for uface in uniqueface_dict.keys():
            face_pairs.pop(uface)
        print(len(uniqueface_dict), len(face_pairs))
        # get all first face of each pair of tuple (face,ielem)
        intfaces = face_from_ufacedict(face_pairs)
        internalfaces.add_elems(ftype, intfaces)
        # get all elements connections via faces
        # get index of connected cells
        f2c = np.array(list(map(lambda flist: [flist[0][1], flist[1][1]], face_pairs.values())))
        print(f2c.shape)
        iface2cell.append(f2c)

    return internalfaces, iface2cell, boundaryfaces, bface2cell
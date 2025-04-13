from collections import defaultdict, OrderedDict
from itertools import chain
import logging

import numpy as np

import cfdtools.api as api
import cfdtools.meshbase._elements as _elem

log = logging.getLogger(__name__)


class indexlist:
    """class of different implementation of list of index

    Returns:
        _type_: _description_
    """

    __available_types = ['list', 'range']

    def __init__(self, irange=None, ilist=None):
        self._type = None
        assert (irange is None) or (ilist is None)
        if irange is not None:
            self.set_range(irange)
        if ilist is not None:
            self.set_list(ilist)

    @property
    def type(self):
        return self._type

    @property
    def size(self):
        return len(self._list) if self.type == 'list' else self._range[1] - self._range[0] + 1

    def _delete(self):
        self._type = None
        self._list = None
        self._range = None

    # @property
    def range(self):
        if self._type == 'range':
            return self._range
        else:
            api.error_stop("unable to get range from list connectivity")

    def set_range(self, irange):
        """define range of index with first and last included"""
        self._delete()
        self._type = 'range'
        self._range = [*irange]

    def list(self):
        if self._type == 'range':
            return list(range(self._range[0], self._range[1] + 1))
        elif self._type == 'list':
            return list(self._list)  # ensure list if may be
        else:
            api.error_stop(f"unknown type: {self._type}")

    def _listof(self, ilist):
        if isinstance(ilist, list):
            rlist = ilist
        elif isinstance(ilist, np.ndarray):
            rlist = self._list = ilist.tolist()
        else:
            log.error(f"{type(ilist)}")
            api.error_stop("unknow type in indexlist class")
        return rlist

    def set_list(self, ilist):
        self._delete()
        self._type = 'list'
        self._list = self._listof(ilist)

    def append(self, ilist: list):
        # assert self._type == 'list'
        # self._list.extend(self._listof(ilist))
        self.set_list(self.list() + self._listof(ilist))

    def shift(self, i):
        if self._type == 'range':
            return indexlist(irange=[self._range[0] + i, self._range[1] + i])
        elif self._type == 'list':
            return indexlist(ilist=[j + i for j in self._list])

    def compress(self):
        """try to make it a range"""
        if self._type == 'list':
            if np.all(self._list == np.arange(self._list[0], self._list[0] + len(self._list))):
                self.set_range([self._list[0], self._list[-1]])
        # else no error, keep list

    def __getitem__(self, indices):  # caution: may be highly costly
        return self.list()[indices]

    def __add__(self, other):
        # TODO: may be optimized to keep ranges
        return indexlist(ilist=self.list() + other.list())

    def __str__(self):
        if self._type == 'range':
            return f"range {self._range[0]} to {self._range[1]}"
        elif self._type == 'list':
            n = np.array(self._list)
            return "list min:max:size = {}:{}:{}".format(n.min(), n.max(), n.size)


class indexindirection:
    """class for a genuine connectivity, regular which size = nelem x 2"""

    def __init__(self, array: np.ndarray = None):
        if array is None:
            self._nelem = 0
            self._conn = None
        else:
            self.conn = array  # setter

    @property
    def conn(self):
        return self._conn

    @conn.setter
    def conn(self, array: np.ndarray):
        assert array.shape[1] == 2
        self._conn = array
        self._nelem = array.shape[0]

    @property
    def nelem(self):
        return self._nelem

    def __getitem__(self, indices):
        # direct access to index without check
        return self._conn.__getitem__(indices)

    def append(self, array: np.ndarray):
        if self._conn is None:
            self.conn = array
        else:
            assert array.shape[1] == 2
            self.conn = np.concatenate((self.conn, array), axis=0)


class compressed_listofindex:
    """class for a compressed list of index (CSR like)

    Returns:
        _type_: _description_
    """

    def __init__(self, index, value) -> None:
        self.set(index, value)

    def set(self, index, value):
        self._index = index
        self._value = value
        self._nelem = self._index.size
        self._size = self._value.size

    @property
    def size(self):
        return self._size

    def check(self):
        assert self._index[0] == 0
        assert self._index.max() == self._size
        return True


class elem_connectivity:
    def __init__(self):
        self._nelem = 0
        self._elem2node = OrderedDict()

    def add_elems(self, etype: str, elem2node: np.ndarray, index: indexlist = None):
        dim = elem2node.shape[0]
        ind = indexlist(irange=[self._nelem, self._nelem + dim - 1]) if index is None else index
        self._nelem += dim
        if etype in self.elems():
            self._elem2node[etype]['index'].append(ind.list())
            self._elem2node[etype]['elem2node'] = np.concatenate(
                (self._elem2node[etype]['elem2node'], elem2node), axis=0
            )
        else:
            self._elem2node[etype] = {'index': ind, 'elem2node': elem2node}

    # def index_elem2node(self, etype):
    #     ist = self._elem2node[etype]['starts']
    #     cell2node = self._elem2node[etype]['cell2node']
    #     return range(ist, ist+cell2node.shape[0]), cell2node

    @property
    def nelem(self):
        return self._nelem

    def items(self):
        return self._elem2node.items()

    def elems(self):
        return self._elem2node.keys()

    def keys(self):  # duplicate but needed to iterate as a dict
        return self._elem2node.keys()

    def __getitem__(self, key):
        return self._elem2node[key]['elem2node']

    def check(self):
        # check uniqueness of all index
        index = np.concatenate(tuple(e2n['index'].list() for _, e2n in self._elem2node.items()))
        uniq = np.unique(index)
        assert index.min() == 0
        assert index.max() == index.size - 1
        assert np.all(uniq == index)
        return True

    # def reindexed_elems(self, etype):
    #     ind = self._elem2node[etype]['index']
    #     e2n[ind,:] = self._elem2node[etype]['elem2node'] # dangerous, index may be too large
    #     return e2n

    def print(self, prefix="", detailed=False):
        for elemtype, elemco in self._elem2node.items():
            log.info(
                prefix + f"{elemtype}: {elemco['elem2node'].shape} with index {elemco['index']}",
            )
            if detailed:
                log.info(prefix + f"  index: {elemco['index'].list()}")
                log.info(prefix + f"  faces: {elemco['elem2node']}")

    def all_index(self):
        return list(sum([econ['index'].list() for _, econ in self.items()], []))

    # def index_elem_tuples(self):
    #     # optim: here, .list() is not mandatory but avoid massively calling .list().getitem()
    #     return list( (i, face.ravel().tolist())
    #                 for _, e2n in self.items()
    #                     for i, face in zip(e2n['index'].list(), np.vsplit(e2n['elem2node'], e2n['elem2node'].shape[0]))
    #             )
    def index_elem_tuples(self):
        """creates a list of tuple (index of face, nodes of faces)"""
        list_of_tuples = []
        for _, e2n in self.items():
            # optim: here, .list() is not mandatory but avoid massively calling .list().getitem()
            ind = e2n['index'].list()
            f2n = e2n['elem2node']
            list_of_tuples.extend([(ind[i], f2n[i, :].ravel().tolist()) for i in range(f2n.shape[0])])
        return list_of_tuples

    def nodes_of_indexlist(self, elemlist):
        """get list of nodes given list of index of elements"""
        ind = []
        f2n = []
        for _, e2n in self.items():
            ind.extend(e2n['index'].list())
            f2n.extend(e2n['elem2node'].tolist())
        # joinlist = list(
        #     chain.from_iterable([tnod[1] for tnod in filter(lambda tup: tup[0] in elemlist, zip(ind, f2n))])
        # )
        joinlist = list(chain.from_iterable(map(dict(zip(ind, f2n)).get, elemlist)))
        return list(set(joinlist))  # make unique

    def importfrom_compressedindex(self, zconn: compressed_listofindex):
        # there is no test but must only applied to faces
        nodeperface = zconn._index[1:] - zconn._index[:-1]
        uniq, counts = np.unique(nodeperface, return_counts=True)
        for facesize, nface in zip(uniq, counts):
            typef = _elem.face_from_nnode[facesize]
            index = np.argwhere(nodeperface == facesize)
            zind = zconn._index[index]
            nodes = np.hstack(tuple(zconn._value[zind + i] for i in range(facesize)))
            self.add_elems(typef, nodes, indexlist(ilist=index))

    # @profile
    def exportto_compressedindex(self) -> compressed_listofindex:
        concat = sorted(self.index_elem_tuples())
        nnode_perface = [len(face) for _, face in concat]
        index = np.concatenate(([0], np.cumsum(nnode_perface)))
        value = np.array(list(chain(*[face for _, face in concat])))
        zconn = compressed_listofindex(index, value)
        return zconn

    def importfrom_merge(self, list_elem):
        # computes size of all elem_connectivity
        sizes = [
            np.sum([np.array(e2n['index'].list()).size for _, e2n in elemcon.items()])
            for elemcon in list_elem
        ]
        sizes = [0] + sizes[:-1]  # start with 0, last is useless
        mergedict = defaultdict(dict)
        for shift, elemcon in zip(sizes, list_elem):
            for key, elemtype in elemcon.items():
                if key in mergedict.keys():
                    mergedict[key]['index'] = mergedict[key]['index'] + elemtype['index'].shift(shift)
                    mergedict[key]['elem2node'] = np.concatenate(
                        (mergedict[key]['elem2node'], elemtype['elem2node']), axis=0
                    )
                else:
                    mergedict[key]['index'] = elemtype['index'].shift(shift)
                    mergedict[key]['elem2node'] = elemtype['elem2node']

        for elem, elemtype in mergedict.items():
            self.add_elems(elem, elemtype['elem2node'], elemtype['index'])

    # @profile
    def create_faces_from_elems(self):
        # @profile
        def __build_face_and_neighbour():
            """build a dict of face type to a list of tuples of each (oriented) face and its neighbor

            Args:
                elems (dict): dict of elements with node definition

            Returns:
                _type_: dict of face type
            """
            faces_neighbour = defaultdict(list)
            for (
                elemtype,
                elemsdict,
            ) in self._elem2node.items():  # elemtype: 'hexa8', elemsarray: ndarray[nelem,8]
                index = elemsdict['index'].list()  # call export to list now
                elemsarray = elemsdict['elem2node']
                # V0
                # for ielem in range(elemsarray.shape[0]):
                #     for ftype, listfaces in _elem.elem2faces[elemtype].items():
                #         # for face in listfaces:
                #         #     faces_neighbour[ftype].append( (tuple(elemsarray[ielem, face]), index[ielem]) )
                #         faces_neighbour[ftype].extend(
                #             [ (tuple(elemsarray[ielem, face]), index[ielem]) for face in listfaces ] )
                # V1 (-3%)
                # for ftype, listfaces in _elem.elem2faces[elemtype].items():
                #     for face in listfaces:
                #         faces_neighbour[ftype].extend(
                #             [ (tuple(elemsarray[ielem, face]), index[ielem]) for ielem in range(elemsarray.shape[0]) ] )
                # V2 (-30%)
                # NODE ORDER of face IS REVERSED
                for ftype, face_of_elem in _elem.elem2faces[elemtype].items():
                    for eface in face_of_elem:
                        reindex_f = elemsarray[:, list(reversed(eface))].tolist()
                        faces_neighbour[ftype].extend(
                            [(tuple(fnodes), ind) for fnodes, ind in zip(reindex_f, index)]
                        )
            return faces_neighbour

        # @profile
        def __find_duplicates(faces_neighbour: dict):
            """find duplicated faces and build unique face/node face/cell connectivity

            Args:
                faces_neighbour (dict): _description_
            """
            internalfaces = elem_connectivity()
            boundaryfaces = elem_connectivity()
            iface2cell = indexindirection()  #
            bface2cell = indexindirection()  #

            def face_from_ufacedict(uface_dict):
                return np.array(list(map(lambda flist: flist[0][0], uface_dict.values())))

            # find pairs for a given face type
            for ftype, listfaces in faces_neighbour.items():
                nf_all = len(listfaces)
                # build a dict of "sorted node" face with list of tuple (face, elem)
                # face_pairs = dict()
                # for uface, facepair in groupby(listfaces, lambda tup: tuple(sorted(tup[0]))):
                #     face_pairs[uface] = list(facepair)
                face_pairs = defaultdict(list)
                for tface in listfaces:
                    face_pairs[tuple(sorted(tface[0]))].append(tface)  # 50% COST
                nf_unique = len(face_pairs)
                assert nf_unique < nf_all
                assert 2 * nf_unique >= nf_all
                # extract all unique face
                #   since reversed when created, boco faces are pointing outward
                uniqueface_dict = dict(
                    filter(lambda tup: len(tup[1]) == 1, face_pairs.items())
                )  # tup[1] is the value of key
                boundaryfaces.add_elems(ftype, face_from_ufacedict(uniqueface_dict))
                f2c = np.full((len(uniqueface_dict), 2), -1)
                # get index of connected cells
                f2c[:, 0] = list(map(lambda flist: flist[0][1], uniqueface_dict.values()))
                bface2cell.append(f2c)
                # remove these faces
                for uface in uniqueface_dict.keys():
                    face_pairs.pop(uface)
                # get all first face of each pair of tuple (face,ielem)
                intfaces = face_from_ufacedict(face_pairs)  # 10% COST
                internalfaces.add_elems(ftype, intfaces)
                # get all elements connections via faces
                # get index of connected cells # 25% COST
                f2c = np.array(
                    list(
                        map(
                            lambda flist: [flist[0][1], flist[1][1]],
                            face_pairs.values(),
                        )
                    )
                )
                iface2cell.append(f2c)

            return internalfaces, iface2cell, boundaryfaces, bface2cell

        faces_neighbour = __build_face_and_neighbour()
        return __find_duplicates(faces_neighbour)

    def nodelist(self):
        """create list of (unique) nodes from element connectivity""" 
        nodeset = set()
        for _, econ in self.items():
            nodeset.update(econ['elem2node'].ravel().tolist())
        return list(nodeset)

    def extrude(self, nplanes: int, inodeshift: int):
        """_summary_

        Args:
            n (int): number of planes
            inodeshift (int): number of nodes of each plane
        """
        assert nplanes > 1
        ncell = nplanes - 1
        newcon = elem_connectivity()
        for etype, econ in self.items():
            nelem = econ['elem2node'].shape[0]
            elemcon = np.tile(econ['elem2node'], (ncell, 2))
            fnnode = _elem.nnode_elem[etype]
            elemcon[:, fnnode : 2 * fnnode] += inodeshift
            for i in range(ncell):
                elemcon[i * nelem : (i + 1) * nelem, :] += i * inodeshift
            # print(etype, econ['elem2node'], elemcon)
            index = np.tile(econ['index'].list(), (ncell))
            for i in range(ncell):
                index[i * nelem : (i + 1) * nelem] += i * inodeshift
            newcon.add_elems(_elem.extruded_face[etype], elemcon, indexlist(ilist=index.tolist()))
        return newcon

import builtins
import logging
from collections import OrderedDict, defaultdict
from itertools import chain

import numpy as np

import cfdtools.api as api
import cfdtools.meshbase._elements as _elem

log = logging.getLogger(__name__)
_mprof_profile = getattr(builtins, "profile", lambda function: function)


@_mprof_profile
def _group_face_candidates(faces, ncorner, return_boundary_keys=False):
    """Group ordered face candidates without changing face numbering.

    This is the common topological kernel used by the monoblock and
    multiblock builders.  Candidate order defines face orientation and final
    numbering; canonical corner-node keys are used only to find duplicates.
    """
    canonical_nodes = np.sort(faces[:, :ncorner], axis=1)
    _, first, inverse, counts = np.unique(
        canonical_nodes,
        axis=0,
        return_index=True,
        return_inverse=True,
        return_counts=True,
    )
    if np.any(counts > 2):
        raise ValueError("non-manifold face found while building connectivity")

    # NumPy sorts canonical keys. Restore first-occurrence order to preserve
    # orientation and numbering from the historical dictionary algorithm.
    group_order = np.argsort(first, kind="mergesort")
    ordered_counts = counts[group_order]
    boundary_groups = group_order[ordered_counts == 1]
    internal_groups = group_order[ordered_counts == 2]
    boundary = first[boundary_groups]
    internal_first = first[internal_groups]

    candidates_by_group = np.argsort(inverse, kind="mergesort")
    group_starts = np.cumsum(counts) - counts
    internal_second = candidates_by_group[group_starts[internal_groups] + 1]
    boundary_keys = (
        canonical_nodes[boundary].copy() if return_boundary_keys else None
    )
    return boundary, internal_first, internal_second, boundary_keys


def _concatenate_or_single(arrays):
    """Return one candidate array directly, concatenating only when needed."""
    return arrays[0] if len(arrays) == 1 else np.concatenate(arrays)


def _concatenate_and_release(arrays):
    """Concatenate axis-zero fragments and release each source after copying."""
    if len(arrays) == 1:
        return arrays.pop()

    dtype = arrays[0].dtype
    size = 0
    trailing_shape = arrays[0].shape[1:]
    for array in arrays:
        dtype = np.promote_types(dtype, array.dtype)
        size += array.shape[0]
    result = np.empty((size, *trailing_shape), dtype=dtype)

    offset = 0
    for index in range(len(arrays)):
        array = arrays[index]
        last = offset + array.shape[0]
        result[offset:last] = array
        arrays[index] = None
        offset = last
    arrays.clear()
    return result


def _assemble_ordered_fragments(arrays, destination_by_source):
    """Copy fragments directly into final order, releasing them one by one."""
    dtype = arrays[0].dtype
    trailing_shape = arrays[0].shape[1:]
    size = 0
    for array in arrays:
        dtype = np.promote_types(dtype, array.dtype)
        size += array.shape[0]
    if size != destination_by_source.size:
        raise ValueError("inconsistent multiblock assembly sizes")

    result = np.empty((size, *trailing_shape), dtype=dtype)
    offset = 0
    for index in range(len(arrays)):
        array = arrays[index]
        last = offset + array.shape[0]
        result[destination_by_source[offset:last]] = array
        arrays[index] = None
        offset = last
    arrays.clear()
    return result


class indexlist:
    """class of different implementation of list of index

    Returns:
        _type_: _description_
    """

    __available_types = ["list", "range"]

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
        return (
            len(self._list)
            if self.type == "list"
            else self._range[1] - self._range[0] + 1
        )

    def _delete(self):
        self._type = None
        self._list = None
        self._range = None

    # @property
    def range(self):
        if self._type == "range":
            return self._range
        else:
            api.error_stop("unable to get range from list connectivity")

    def set_range(self, irange):
        """define range of index with first and last included"""
        self._delete()
        self._type = "range"
        self._range = [*irange]

    def list(self):
        if self._type == "range":
            return list(range(self._range[0], self._range[1] + 1))
        elif self._type == "list":
            return list(self._list)  # ensure list if may be
        else:
            api.error_stop(f"unknown type: {self._type}")

    def to_numpy(self, dtype=None):
        """Return indices as a NumPy array without an intermediate Python list.

        Range-backed indices are generated directly with :func:`numpy.arange`.
        List-backed indices necessarily allocate an array but avoid making an
        additional copy of the Python list first.
        """
        if self._type == "range":
            return np.arange(self._range[0], self._range[1] + 1, dtype=dtype)
        if self._type == "list":
            return np.asarray(self._list, dtype=dtype)
        api.error_stop(f"unknown type: {self._type}")

    def _listof(self, ilist):
        if isinstance(ilist, list):
            rlist = ilist
        elif isinstance(ilist, np.ndarray):
            rlist = self._list = ilist.tolist()
        else:
            log.error(f"{type(ilist)}")
            api.error_stop("unknown type in indexlist class")
        return rlist

    def set_list(self, ilist):
        self._delete()
        self._type = "list"
        self._list = self._listof(ilist)

    def append(self, ilist: list):
        # assert self._type == 'list'
        # self._list.extend(self._listof(ilist))
        self.set_list(self.list() + self._listof(ilist))

    def shift(self, i):
        if self._type == "range":
            return indexlist(irange=[self._range[0] + i, self._range[1] + i])
        elif self._type == "list":
            return indexlist(ilist=[j + i for j in self._list])

    def compress(self):
        """try to make it a range"""
        if self._type == "list":
            if np.all(
                self._list == np.arange(self._list[0], self._list[0] + len(self._list))
            ):
                self.set_range([self._list[0], self._list[-1]])
        # else no error, keep list

    def __getitem__(self, indices):  # caution: may be highly costly
        return self.list()[indices]

    def __add__(self, other):
        # TODO: may be optimized to keep ranges
        return indexlist(ilist=self.list() + other.list())

    def __str__(self):
        if self._type == "range":
            return f"range {self._range[0]} to {self._range[1]}"
        elif self._type == "list":
            n = np.array(self._list)
            return "list min:max:size = {}:{}:{}".format(n.min(), n.max(), n.size)


class indexindirection:
    """class for a genuine connectivity, regular which size = nelem x 2"""

    def __init__(self, array: np.ndarray = None):
        self.conn = array  # setter

    def __str__(self):
        return str(self.conn)

    @property
    def conn(self):
        return self._conn

    @conn.setter
    def conn(self, array: np.ndarray):
        if array is None:
            self._nelem = 0
            self._conn = None
        else:
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
        if self.conn is None:
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

    def __str__(self) -> str:
        return str(self._elem2node)

    def add_elems(self, etype: str, elem2node: np.ndarray, index: indexlist = None):
        dim = elem2node.shape[0]
        ind = (
            indexlist(irange=[self._nelem, self._nelem + dim - 1])
            if index is None
            else index
        )
        self._nelem += dim
        if etype in self.elems():
            self._elem2node[etype]["index"].append(ind.list())
            self._elem2node[etype]["elem2node"] = np.concatenate(
                (self._elem2node[etype]["elem2node"], elem2node), axis=0
            )
        else:
            self._elem2node[etype] = {"index": ind, "elem2node": elem2node}

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
        return self._elem2node[key]["elem2node"]

    def check(self):
        # check uniqueness of all index
        index = np.concatenate(
            tuple(e2n["index"].list() for _, e2n in self._elem2node.items())
        )
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
                prefix
                + f"{elemtype}: {elemco['elem2node'].shape} with index {elemco['index']}",
            )
            if detailed:
                log.info(prefix + f"  index: {elemco['index'].list()}")
                log.info(prefix + f"  faces: {elemco['elem2node']}")

    def all_index(self):
        return list(sum([econ["index"].list() for _, econ in self.items()], []))

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
            ind = e2n["index"].list()
            f2n = e2n["elem2node"]
            list_of_tuples.extend(
                [(ind[i], f2n[i, :].ravel().tolist()) for i in range(f2n.shape[0])]
            )

        return list_of_tuples

    def nodes_of_indexlist(self, elemlist):
        """get list of nodes given list of index of elements"""
        ind = []
        f2n = []
        for _, e2n in self.items():
            ind.extend(e2n["index"].list())
            f2n.extend(e2n["elem2node"].tolist())
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
        concat = self.index_elem_tuples()
        concat.sort()
        nnode_perface = [len(face) for _, face in concat]
        index = np.concatenate(([0], np.cumsum(nnode_perface)))
        value = np.array(list(chain(*[face for _, face in concat])))
        zconn = compressed_listofindex(index, value)
        return zconn

    def importfrom_merge(self, list_elem):
        # computes size of all elem_connectivity
        sizes = [
            np.sum([np.array(e2n["index"].list()).size for _, e2n in elemcon.items()])
            for elemcon in list_elem
        ]
        sizes = [0] + sizes[:-1]  # start with 0, last is useless
        mergedict = defaultdict(dict)
        for shift, elemcon in zip(sizes, list_elem):
            for key, elemtype in elemcon.items():
                if key in mergedict.keys():
                    mergedict[key]["index"] = mergedict[key]["index"] + elemtype[
                        "index"
                    ].shift(shift)
                    mergedict[key]["elem2node"] = np.concatenate(
                        (mergedict[key]["elem2node"], elemtype["elem2node"]), axis=0
                    )
                else:
                    mergedict[key]["index"] = elemtype["index"].shift(shift)
                    mergedict[key]["elem2node"] = elemtype["elem2node"]

        for elem, elemtype in mergedict.items():
            self.add_elems(elem, elemtype["elem2node"], elemtype["index"])

    # @profile
    def _create_faces_from_elems_legacy(self):
        # @profile
        def __build_face_and_neighbour():
            """Build a dict of face type to a list of tuples of each (oriented) face and its neighbors.

            Args:
                elems (dict): dict of elements with node definition

            Returns:
                _type_: dict of face type
            """
            faces_neighbour = defaultdict(list)
            # elemtype: 'hexa8', elemsarray: ndarray[nelem,8]
            for elemtype, elemsdict in self._elem2node.items():
                # cell/element indices
                index = elemsdict["index"].list()  # call export to list now
                elemsarray = elemsdict["elem2node"]
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
                for ftype, face_of_elem in _elem.elem2faces[elemtype].items():
                    # Loop over the definition of each face of the element type
                    for eface in face_of_elem:
                        # Get all face definition from the mesh for each face.
                        index_f = elemsarray[:, eface]
                        # Store the face indices with the element index
                        faces_neighbour[ftype].extend(
                            [
                                (tuple(fnodes), ind)
                                for fnodes, ind in zip(index_f, index)
                            ]
                        )

            return faces_neighbour

        # @profile
        def __find_duplicates(faces_neighbour: dict):
            """Find duplicated faces and build unique face/node face/cell connectivity.

            Parameters
            ----------
            faces_neighbour: dict
                 description

            Returns
            -------
            internalfaces: elem_connectivity
                Connectivity of internal faces.
            iface2cell: indexindirection
                Cells of each internal face.
            boundaryfaces: elem_connectivity
                Connectivity of boundary faces.
            bface2cell: indexindirection
                Cells of each boundary face.
            """
            internalfaces = elem_connectivity()
            boundaryfaces = elem_connectivity()
            iface2cell = indexindirection()  #
            # give the two cells attached to a face
            bface2cell = indexindirection()

            def face_from_ufacedict(uface_dict):
                """Get the tuple of node indices for all first face."""
                return np.array([flist[0][0] for flist in uface_dict.values()])

            # find pairs for a given face type
            for ftype, listfaces in faces_neighbour.items():
                nf_all = len(listfaces)
                # quad interface only, could be adapted to other ftype
                nb_linear_face_pts = 4
                # build a dict of "sorted node" face with list of tuple (face, elem)
                # face_pairs = dict()
                # for uface, facepair in groupby(listfaces, lambda tup: tuple(sorted(tup[0]))):
                #     face_pairs[uface] = list(facepair)
                face_pairs = defaultdict(list)
                for tface in listfaces:
                    # nodes of faces are sorted to get a dict key
                    face_pairs[tuple(sorted(tface[0][:nb_linear_face_pts]))].append(
                        tface
                    )  # 50% COST
                nf_unique = len(face_pairs)
                assert nf_unique < nf_all
                assert 2 * nf_unique >= nf_all

                # extract all unique faces: boundary faces
                #   since reversed when created, boco faces are pointing outward
                unique_faces = {
                    key: pair for key, pair in face_pairs.items() if len(pair) == 1
                }
                boundaryfaces.add_elems(ftype, face_from_ufacedict(unique_faces))
                # get cells of each boundary face. Left cell is always defined.
                # Right cell is -1.
                # get index of connected cells
                f2c = np.array([[flist[0][1], -1] for flist in unique_faces.values()])
                bface2cell.append(f2c)

                # extract all double faces: internal faces
                face_pairs = {
                    key: pair for key, pair in face_pairs.items() if len(pair) == 2
                }
                # get all first face of each pair of tuple (face,ielem)
                intfaces = face_from_ufacedict(face_pairs)  # 10% COST
                internalfaces.add_elems(ftype, intfaces)
                # get cells of each internal face.
                # get all elements connections via faces
                # get index of connected cells # 25% COST
                f2c = np.array(
                    [[flist[0][1], flist[1][1]] for flist in face_pairs.values()]
                )
                iface2cell.append(f2c)

            return internalfaces, iface2cell, boundaryfaces, bface2cell

        faces_neighbour = __build_face_and_neighbour()
        return __find_duplicates(faces_neighbour)

    def create_faces_from_elems(self):
        """Create face connectivity with NumPy grouping by canonical face nodes.

        The first oriented occurrence of a face is kept, matching the former
        insertion-order behaviour while avoiding millions of Python objects.
        """
        # Preallocate the complete candidate arrays.  The previous list-plus-
        # concatenate approach temporarily retained both every source block
        # and the concatenated copy, nearly doubling candidate storage.
        candidate_counts = OrderedDict()
        node_dtypes = {}
        cell_dtypes = {}
        cell_indices_by_type = {}
        for elemtype, elemsdict in self._elem2node.items():
            cell_nodes = elemsdict["elem2node"]
            cell_indices = elemsdict["index"].to_numpy()
            cell_indices_by_type[elemtype] = cell_indices
            for ftype, face_of_elem in _elem.elem2faces[elemtype].items():
                if not face_of_elem:
                    continue
                candidate_counts[ftype] = candidate_counts.get(ftype, 0) + (
                    cell_nodes.shape[0] * len(face_of_elem)
                )
                node_dtypes.setdefault(ftype, []).append(cell_nodes.dtype)
                cell_dtypes.setdefault(ftype, []).append(cell_indices.dtype)

        faces_by_type = {
            ftype: np.empty(
                (count, _elem.nnode_elem[ftype]),
                dtype=np.result_type(*node_dtypes[ftype]),
            )
            for ftype, count in candidate_counts.items()
        }
        cells_by_type = {
            ftype: np.empty(
                count, dtype=np.result_type(*cell_dtypes[ftype])
            )
            for ftype, count in candidate_counts.items()
        }
        offsets = {ftype: 0 for ftype in candidate_counts}
        for elemtype, elemsdict in self._elem2node.items():
            cell_nodes = elemsdict["elem2node"]
            cell_indices = cell_indices_by_type[elemtype]
            nelemtype = cell_nodes.shape[0]
            for ftype, face_of_elem in _elem.elem2faces[elemtype].items():
                for local_face in face_of_elem:
                    first = offsets[ftype]
                    last = first + nelemtype
                    faces_by_type[ftype][first:last] = cell_nodes[:, local_face]
                    cells_by_type[ftype][first:last] = cell_indices
                    offsets[ftype] = last
        del cell_indices_by_type

        internalfaces = elem_connectivity()
        boundaryfaces = elem_connectivity()
        iface2cell = indexindirection()
        bface2cell = indexindirection()

        for ftype, faces in faces_by_type.items():
            cells = cells_by_type[ftype]
            boundary, internal_first, internal_second, _ = (
                _group_face_candidates(
                    faces, _elem.ncorner_linear[ftype]
                )
            )
            boundaryfaces.add_elems(ftype, faces[boundary])
            bface2cell.append(
                np.column_stack(
                    (
                        cells[boundary],
                        -np.ones(boundary.size, dtype=cells.dtype),
                    )
                )
            )
            internalfaces.add_elems(ftype, faces[internal_first])
            iface2cell.append(
                np.column_stack(
                    (cells[internal_first], cells[internal_second])
                )
            )

        return internalfaces, iface2cell, boundaryfaces, bface2cell

    @staticmethod
    def _append_internal_face_records(records, ranks, faces, neighbours):
        """Append one independently releasable batch of internal faces."""
        records["ranks"].append(ranks)
        records["faces"].append(faces)
        records["neighbours"].append(neighbours)

    @staticmethod
    @_mprof_profile
    def _assemble_multiblock_connectivity(face_order, internal_records, pending):
        """Assemble final arrays while releasing block fragments progressively."""
        internalfaces = elem_connectivity()
        boundaryfaces = elem_connectivity()
        iface2cell = indexindirection()
        bface2cell = indexindirection()

        for ftype in face_order:
            records = internal_records[ftype]
            if records["ranks"]:
                ranks = _concatenate_and_release(records["ranks"])
                order = np.argsort(ranks, kind="mergesort")
                del ranks
                destination_by_source = np.empty_like(order)
                destination_by_source[order] = np.arange(
                    order.size, dtype=order.dtype
                )
                del order

                ordered_faces = _assemble_ordered_fragments(
                    records["faces"], destination_by_source
                )
                internalfaces.add_elems(ftype, ordered_faces)
                del ordered_faces

                ordered_neighbours = _assemble_ordered_fragments(
                    records["neighbours"], destination_by_source
                )
                del destination_by_source
                iface2cell.append(ordered_neighbours)
                del ordered_neighbours
            else:
                # Preserve the monoblock representation: a face type which
                # occurs in elements is present even when it has no internal
                # faces (for example prism triangular faces on a boundary).
                internalfaces.add_elems(
                    ftype, np.empty((0, _elem.nnode_elem[ftype]), dtype=int)
                )
                iface2cell.append(np.empty((0, 2), dtype=int))

            open_faces = pending[ftype]
            if open_faces:
                nface = len(open_faces)
                sample = next(iter(open_faces.values()))
                rank_dtype = np.asarray(sample[0]).dtype
                face_dtype = sample[1].dtype
                cell_dtype = np.asarray(sample[2]).dtype
                for record in open_faces.values():
                    rank_dtype = np.promote_types(
                        rank_dtype, np.asarray(record[0]).dtype
                    )
                    face_dtype = np.promote_types(
                        face_dtype, record[1].dtype
                    )
                    cell_dtype = np.promote_types(
                        cell_dtype, np.asarray(record[2]).dtype
                    )
                ranks = np.empty(nface, dtype=rank_dtype)
                faces = np.empty(
                    (nface, _elem.nnode_elem[ftype]),
                    dtype=face_dtype,
                )
                cells = np.empty(nface, dtype=cell_dtype)
                for index, record in enumerate(open_faces.values()):
                    ranks[index], faces[index], cells[index] = record
                open_faces.clear()

                order = np.argsort(ranks, kind="mergesort")
                del ranks
                ordered_faces = faces[order]
                del faces
                ordered_cells = cells[order]
                del cells, order
                boundaryfaces.add_elems(ftype, ordered_faces)
                del ordered_faces
                bface2cell.append(
                    np.column_stack(
                        (
                            ordered_cells,
                            -np.ones(
                                ordered_cells.size,
                                dtype=ordered_cells.dtype,
                            ),
                        )
                    )
                )
                del ordered_cells
            else:
                boundaryfaces.add_elems(
                    ftype, np.empty((0, _elem.nnode_elem[ftype]), dtype=int)
                )
                bface2cell.append(np.empty((0, 2), dtype=int))

        return internalfaces, iface2cell, boundaryfaces, bface2cell

    @staticmethod
    @_mprof_profile
    def _validate_unique_final_faces(internalfaces, boundaryfaces):
        """Reject non-manifold faces with bounded comparison temporaries."""
        comparison_chunk = 1_000_000
        for ftype, face_data in internalfaces.items():
            face_arrays = [face_data["elem2node"]]
            if ftype in boundaryfaces.elems():
                face_arrays.append(boundaryfaces[ftype])
            ncorner = _elem.ncorner_linear[ftype]
            nface = sum(array.shape[0] for array in face_arrays)
            canonical_nodes = np.empty(
                (nface, ncorner), dtype=face_arrays[0].dtype
            )
            offset = 0
            for array in face_arrays:
                last = offset + array.shape[0]
                canonical_nodes[offset:last] = array[:, :ncorner]
                canonical_nodes[offset:last].sort(axis=1)
                offset = last

            order = np.lexsort(
                tuple(
                    canonical_nodes[:, column]
                    for column in range(ncorner - 1, -1, -1)
                )
            )
            for first in range(0, max(nface - 1, 0), comparison_chunk):
                last = min(first + comparison_chunk, nface - 1)
                left = canonical_nodes[order[first:last]]
                right = canonical_nodes[order[first + 1 : last + 1]]
                if np.any(np.all(left == right, axis=1)):
                    raise ValueError(
                        "non-manifold face found while building connectivity"
                    )
            del canonical_nodes, order

    @_mprof_profile
    def create_faces_from_elems_multiblock(self, nblocks, verbose=False):
        """Create face connectivity by processing consecutive cell blocks.

        Unlike :meth:`create_faces_from_elems`, this method only retains face
        candidates for one block at a time.  A compact dictionary records the
        faces which are still open at a block boundary; a later block either
        reconnects one of them or adds another open face.  Final face and
        neighbour arrays are, of course, retained as they are the result of
        the operation.

        Cells are partitioned in increasing cell-index order. Consequently,
        locality of these indices matters for memory use: a spatially ordered
        mesh keeps a small interface between blocks, whereas a scrambled
        numbering can make the open-face registry large.

        Args:
            nblocks: Number of non-empty consecutive cell-index blocks.
            verbose: Log per-block candidate, creation and reconnection
                counts at ``INFO`` level.
        """
        if not isinstance(nblocks, (int, np.integer)) or nblocks < 1:
            raise ValueError("nblocks must be a positive integer")
        if self.nelem == 0:
            raise ValueError("cannot build faces from an empty connectivity")
        if nblocks > self.nelem:
            raise ValueError("nblocks cannot exceed the number of cells")

        if verbose:
            log.info("multiblock face connectivity (%d blocks)", nblocks)
            log.info(
                "  %9s %15s %10s %12s %10s %12s %10s",
                "block",
                "cell ranks",
                "cells",
                "candidates",
                "opened",
                "reconnected",
                "pending",
            )

        # The monoblock implementation is the established, fully vectorized
        # reference.  Keep it on this path both for speed and exact legacy
        # ordering.
        if nblocks == 1:
            result = self.create_faces_from_elems()
            if verbose:
                ncandidate = sum(
                    len(local_faces) * elemsdict["elem2node"].shape[0]
                    for elemtype, elemsdict in self._elem2node.items()
                    for local_faces in _elem.elem2faces[elemtype].values()
                )
                log.info(
                    "  %9s %15s %10d %12d %10d %12d %10d",
                    "1/1",
                    "0:{}".format(self.nelem),
                    self.nelem,
                    ncandidate,
                    result[2].nelem,
                    result[0].nelem,
                    result[2].nelem,
                )
            return result

        # ``rank`` reproduces the first-occurrence ordering of the monoblock
        # implementation even though candidates are visited block by block.
        # It is used only for the final result, not to retain candidates.
        face_order = OrderedDict()
        rank_bases = {}
        for elemtype, elemsdict in self._elem2node.items():
            nelemtype = elemsdict["elem2node"].shape[0]
            for ftype, local_faces in _elem.elem2faces[elemtype].items():
                if not local_faces:
                    continue
                rank_bases[(elemtype, ftype)] = face_order.get(ftype, 0)
                face_order[ftype] = face_order.get(ftype, 0) + (
                    len(local_faces) * nelemtype
                )

        ordered_cell_indices = np.sort(
            np.concatenate(
                [
                    elemsdict["index"].to_numpy()
                    for elemsdict in self._elem2node.values()
                ]
            )
        )
        block_edges = np.linspace(0, self.nelem, nblocks + 1, dtype=int)
        pending = {ftype: {} for ftype in face_order}
        internal_records = {
            ftype: {"ranks": [], "faces": [], "neighbours": []}
            for ftype in face_order
        }

        for iblock, (first_rank, last_rank) in enumerate(
            zip(block_edges[:-1], block_edges[1:]), start=1
        ):
            block_cell_indices = ordered_cell_indices[first_rank:last_rank]
            face_blocks = {ftype: [] for ftype in face_order}
            cell_blocks = {ftype: [] for ftype in face_order}
            rank_blocks = {ftype: [] for ftype in face_order}

            for elemtype, elemsdict in self._elem2node.items():
                all_cells = elemsdict["index"].to_numpy()
                selected = np.isin(all_cells, block_cell_indices)
                if not np.any(selected):
                    continue
                row_indices = np.flatnonzero(selected)
                cell_nodes = elemsdict["elem2node"][selected]
                cell_indices = all_cells[selected]
                for ftype, local_faces in _elem.elem2faces[elemtype].items():
                    if not local_faces:
                        continue
                    face_blocks[ftype].append(
                        np.concatenate(
                            [
                                cell_nodes[:, local_face]
                                for local_face in local_faces
                            ],
                            axis=0,
                        )
                    )
                    cell_blocks[ftype].append(
                        np.tile(cell_indices, len(local_faces))
                    )
                    rank_base = rank_bases[(elemtype, ftype)]
                    rank_blocks[ftype].append(
                        np.concatenate(
                            [
                                rank_base
                                + ilocal * all_cells.size
                                + row_indices
                                for ilocal in range(len(local_faces))
                            ]
                        )
                    )

            block_candidates = block_opened = block_reconnected = 0
            for ftype in face_order:
                if not face_blocks[ftype]:
                    continue
                faces = _concatenate_or_single(face_blocks[ftype])
                cells = _concatenate_or_single(cell_blocks[ftype])
                ranks = _concatenate_or_single(rank_blocks[ftype])
                # Release source arrays immediately when concatenation made a
                # copy. Candidate ranks are already in monoblock order.
                face_blocks[ftype] = None
                cell_blocks[ftype] = None
                rank_blocks[ftype] = None
                (
                    boundary,
                    internal_first,
                    internal_second,
                    boundary_keys,
                ) = _group_face_candidates(
                    faces,
                    _elem.ncorner_linear[ftype],
                    return_boundary_keys=True,
                )
                if internal_first.size:
                    self._append_internal_face_records(
                        internal_records[ftype],
                        ranks[internal_first],
                        faces[internal_first],
                        np.column_stack(
                            (
                                cells[internal_first],
                                cells[internal_second],
                            )
                        ),
                    )
                    block_reconnected += internal_first.size

                matched_ranks = []
                matched_faces = []
                matched_neighbours = []
                for candidate, canonical_key in zip(boundary, boundary_keys):
                    key = tuple(canonical_key)
                    current = (
                        ranks[candidate], faces[candidate], cells[candidate]
                    )
                    previous = pending[ftype].pop(key, None)
                    if previous is None:
                        # A row slice is a view of the full candidate block.
                        # Copying this small face is essential: retaining the
                        # view would retain the entire block in memory.
                        pending[ftype][key] = (
                            current[0], current[1].copy(), current[2]
                        )
                        block_opened += 1
                    else:
                        if previous[0] <= current[0]:
                            first_face, second_face = previous, current
                        else:
                            first_face, second_face = current, previous
                        matched_ranks.append(first_face[0])
                        matched_faces.append(first_face[1])
                        matched_neighbours.append(
                            (first_face[2], second_face[2])
                        )
                        block_reconnected += 1
                if matched_ranks:
                    self._append_internal_face_records(
                        internal_records[ftype],
                        np.asarray(matched_ranks),
                        np.asarray(matched_faces),
                        np.asarray(matched_neighbours),
                    )
                block_candidates += faces.shape[0]
                # Do not retain views or temporary arrays from the last face
                # type until final assembly.
                current = previous = first_face = second_face = None
                canonical_key = key = None
                del (
                    faces,
                    cells,
                    ranks,
                    boundary,
                    internal_first,
                    internal_second,
                    boundary_keys,
                    matched_ranks,
                    matched_faces,
                    matched_neighbours,
                )

            if verbose:
                log.info(
                    "  %9s %15s %10d %12d %10d %12d %10d",
                    "{}/{}".format(iblock, nblocks),
                    "{}:{}".format(first_rank, last_rank),
                    last_rank - first_rank,
                    block_candidates,
                    block_opened,
                    block_reconnected,
                    sum(len(open_faces) for open_faces in pending.values()),
                )

        # Python loop variables remain alive after a loop. Explicitly release
        # the last block selections and candidate containers before allocating
        # final arrays.
        block_cell_indices = face_blocks = cell_blocks = rank_blocks = None
        all_cells = selected = row_indices = None
        cell_nodes = cell_indices = None
        ordered_cell_indices = block_edges = None

        result = self._assemble_multiblock_connectivity(
            face_order, internal_records, pending
        )
        internal_records = pending = None

        # A face completed in an earlier block may be encountered again only
        # later. Detect that non-manifold case after final assembly.
        self._validate_unique_final_faces(result[0], result[2])
        return result

    def nodelist(self):
        """Create list of (unique) nodes from element connectivity."""
        nodeset = set()
        for _, econ in self.items():
            nodeset.update(econ["elem2node"].ravel().tolist())
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
            nelem = econ["elem2node"].shape[0]
            elemcon = np.tile(econ["elem2node"], (ncell, 2))
            fnnode = _elem.nnode_elem[etype]
            elemcon[:, fnnode : 2 * fnnode] += inodeshift
            for i in range(ncell):
                elemcon[i * nelem : (i + 1) * nelem, :] += i * inodeshift
            # print(etype, econ['elem2node'], elemcon)
            index = np.tile(econ["index"].list(), (ncell))
            for i in range(ncell):
                index[i * nelem : (i + 1) * nelem] += i * inodeshift
            newcon.add_elems(
                _elem.extruded_face[etype], elemcon, indexlist(ilist=index.tolist())
            )
        return newcon

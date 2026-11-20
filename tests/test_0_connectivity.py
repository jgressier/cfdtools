import cfdtools.meshbase._connectivity as _conn
import cfdtools.meshbase._elements as _elem
import cfdtools.meshbase.simple as _simple
import numpy as np
import pytest


def test_indexlist_initlist0():
    ilist = [10, 12, 15, 20]
    iconn = _conn.indexlist()
    iconn.set_list(ilist)
    assert iconn.list() == ilist


def test_indexlist_initlist1():
    ilist = [10, 12, 15, 20]
    iconn = _conn.indexlist(ilist=ilist)
    assert iconn.list() == ilist  # test property


def test_indexlist_initrange0():
    imin, imax = 10, 20
    iconn = _conn.indexlist()
    iconn.set_range([imin, imax])
    assert iconn.range() == [imin, imax]
    assert iconn.list() == list(range(imin, imax + 1))


@pytest.mark.parametrize(
    "iconn, expected",
    [
        (_conn.indexlist(irange=[10, 20]), np.arange(10, 21)),
        (_conn.indexlist(ilist=[10, 12, 15, 20]), np.array([10, 12, 15, 20])),
    ],
)
def test_indexlist_to_numpy(iconn, expected):
    assert np.array_equal(iconn.to_numpy(), expected)


def test_indexlist_range_to_numpy_does_not_materialize_list(monkeypatch):
    iconn = _conn.indexlist(irange=[10, 20])
    monkeypatch.setattr(iconn, "list", lambda: pytest.fail("list() was called"))
    assert np.array_equal(iconn.to_numpy(), np.arange(10, 21))


def test_element_corner_metadata_is_complete_and_consistent():
    assert _elem.ncorner_linear.keys() == _elem.nnode_elem.keys()
    assert all(
        0 < _elem.ncorner_linear[elemtype] <= nnode
        for elemtype, nnode in _elem.nnode_elem.items()
    )
    assert _elem.ncorner_linear["quad9"] == 4


def test_indexindirection_void():
    dconn = _conn.indexindirection()
    assert dconn.nelem == 0


def test_indexindirection_init():
    i = np.arange(10)
    j = np.arange(10) * 2
    dconn = _conn.indexindirection(array=np.array([i, j]).T)
    assert np.all(dconn.conn[:, 1] == j)
    assert np.all(dconn[:, 1] == j)


def test_indexindirection_append():
    i = np.arange(10)
    j = np.arange(10) * 2
    dconn = _conn.indexindirection(array=np.array([i, j]).T)
    i = np.arange(10) + 10
    j = np.arange(10) * 3
    dconn.append(np.array([i, j]).T)
    assert dconn.nelem == 20


def test_compressed_listofindex():
    i = np.arange(11) * 3  # 10 elements of 3 values, add last index 30
    v = 100 + np.arange(30)
    zcon = _conn.compressed_listofindex(i, v)
    assert zcon.check()


class TestElem:
    dict_basiccon = {
        '2quads': ('quad4', np.array([[0, 1, 2, 3], [2, 1, 4, 5]])),
        '3quadsintri': ('quad4', np.array([[0, 5, 6, 4], [5, 1, 3, 6], [2, 4, 6, 3]])),
        '2hexas': (
            'hexa8',
            np.array([[0, 1, 2, 3, 4, 5, 6, 7], [4, 5, 6, 7, 8, 9, 10, 11]]),
        ),
    }

    def set_elemcon(self, econ):
        elemc = _conn.elem_connectivity()
        elemc.add_elems(*self.dict_basiccon[econ])  # type and connectivy
        return elemc

    @staticmethod
    def two_cells_sharing_first_face(elemtype, corners_only=False):
        nnode = _elem.nnode_elem[elemtype]
        cell0 = np.arange(nnode)
        cell1 = np.arange(nnode, 2 * nnode)
        local_face = next(
            face
            for faces in _elem.elem2faces[elemtype].values()
            for face in faces
        )
        shared = local_face[:4] if corners_only else local_face
        cell1[shared] = cell0[shared]
        elemc = _conn.elem_connectivity()
        elemc.add_elems(elemtype, np.vstack((cell0, cell1)))
        return elemc

    @staticmethod
    def assert_results_equal(reference, transformed):
        for reference_connectivity, transformed_connectivity in zip(reference, transformed):
            if isinstance(reference_connectivity, _conn.elem_connectivity):
                assert list(reference_connectivity.elems()) == list(
                    transformed_connectivity.elems()
                )
                for elemtype, reference_data in reference_connectivity.items():
                    assert np.array_equal(
                        reference_data['elem2node'],
                        transformed_connectivity[elemtype],
                    )
                    assert reference_data['index'].list() == transformed_connectivity._elem2node[
                        elemtype
                    ]['index'].list()
            else:
                assert np.array_equal(reference_connectivity.conn, transformed_connectivity.conn)

    @pytest.mark.parametrize('econ', ['2quads', '2hexas'])
    def test_elem_init(self, econ):
        elemc = self.set_elemcon(econ)
        assert elemc.check()

    def test_destructive_fragment_assembly_releases_sources(self):
        fragments = [
            np.array([[10], [20]], dtype=np.int32),
            np.array([[30]], dtype=np.int64),
        ]
        concatenated = _conn._concatenate_and_release(fragments)
        assert fragments == []
        assert concatenated.dtype == np.int64
        assert np.array_equal(concatenated[:, 0], [10, 20, 30])

        fragments = [np.array([[10], [20]]), np.array([[30]])]
        destination_by_source = np.array([2, 0, 1])
        ordered = _conn._assemble_ordered_fragments(
            fragments, destination_by_source
        )
        assert fragments == []
        assert np.array_equal(ordered[:, 0], [20, 30, 10])

    @pytest.mark.parametrize('econ', ['2quads', '2hexas'])
    def test_elem_createface(self, econ):
        elemc = self.set_elemcon(econ)
        (
            intfaces,
            iface2cell,
            boundaryfaces,
            bface2cell,
        ) = elemc.create_faces_from_elems()
        assert iface2cell.nelem == intfaces.nelem
        assert bface2cell.nelem == boundaryfaces.nelem

    @pytest.mark.parametrize('econ', ['2quads', '3quadsintri', '2hexas'])
    def test_elem_createface_matches_legacy_ordering(self, econ):
        elemc = self.set_elemcon(econ)
        self.assert_results_equal(
            elemc._create_faces_from_elems_legacy(),
            elemc.create_faces_from_elems(),
        )

    @pytest.mark.parametrize(
        'elemtype',
        ['quad4', 'tetra4', 'pyra5', 'prism6', 'hexa8', 'hexa27'],
    )
    def test_current_matches_legacy_for_all_supported_elements(self, elemtype):
        """Compare every face type within the narrower legacy input contract."""
        nnode = _elem.nnode_elem[elemtype]
        cells = []
        face_groups = [
            faces for faces in _elem.elem2faces[elemtype].values() if faces
        ]
        for pair, faces in enumerate(face_groups):
            cell0 = np.arange(2 * pair * nnode, (2 * pair + 1) * nnode)
            cell1 = np.arange((2 * pair + 1) * nnode, (2 * pair + 2) * nnode)
            cell1[faces[0]] = cell0[faces[0]]
            cells.extend((cell0, cell1))
        cells = np.asarray(cells)
        cell_indices = np.arange(10, 10 + 7 * len(cells), 7)
        elemc = _conn.elem_connectivity()
        elemc.add_elems(
            elemtype,
            cells,
            _conn.indexlist(ilist=cell_indices),
        )
        self.assert_results_equal(
            elemc._create_faces_from_elems_legacy(),
            elemc.create_faces_from_elems(),
        )

    def test_current_matches_legacy_on_structured_cube(self):
        """Compare realistic internal and boundary faces with range-backed indices."""
        cube_mesh = _simple.Cube(3, 2, 2).export_mesh()
        elemc = _conn.elem_connectivity()
        elemc.add_elems('hexa8', cube_mesh._cell2node['hexa8'])
        assert elemc._elem2node['hexa8']['index'].type == 'range'
        self.assert_results_equal(
            elemc._create_faces_from_elems_legacy(),
            elemc.create_faces_from_elems(),
        )

    @pytest.mark.parametrize('nblocks', [2, 3, 10])
    def test_multiblock_matches_monoblock_on_structured_cube(self, nblocks):
        """Block boundaries reconnect to the exact monoblock connectivity."""
        cube_mesh = _simple.Cube(3, 2, 2).export_mesh()
        elemc = _conn.elem_connectivity()
        elemc.add_elems('hexa8', cube_mesh._cell2node['hexa8'])
        self.assert_results_equal(
            elemc.create_faces_from_elems(),
            elemc.create_faces_from_elems_multiblock(nblocks),
        )

    def test_multiblock_matches_monoblock_on_mixed_elements(self):
        """A face shared by different element types may cross a block edge."""
        # Different node-index dtypes also exercise promotion while pending
        # faces from separate blocks are assembled.
        hexa = np.arange(_elem.nnode_elem['hexa8'], dtype=np.int32)
        prism = np.arange(20, 20 + _elem.nnode_elem['prism6'])
        hexa_face = _elem.elem2faces['hexa8']['quad4'][0]
        prism_face = _elem.elem2faces['prism6']['quad4'][0]
        prism[prism_face] = hexa[hexa_face]
        elemc = _conn.elem_connectivity()
        elemc.add_elems('hexa8', hexa[np.newaxis, :])
        elemc.add_elems('prism6', prism[np.newaxis, :])
        self.assert_results_equal(
            elemc.create_faces_from_elems(),
            elemc.create_faces_from_elems_multiblock(2),
        )

    def test_multiblock_reports_block_progress(self, caplog):
        elemc = self.set_elemcon('2hexas')
        with caplog.at_level('INFO', logger='cfdtools.meshbase._connectivity'):
            elemc.create_faces_from_elems_multiblock(2, verbose=True)
        assert 'block' in caplog.text
        assert 'candidates' in caplog.text
        assert 'opened' in caplog.text
        assert 'reconnected' in caplog.text
        assert 'pending' in caplog.text
        assert '1/2' in caplog.text

    @pytest.mark.parametrize('nblocks', [0, 3])
    def test_multiblock_rejects_invalid_block_count(self, nblocks):
        elemc = self.set_elemcon('2hexas')
        with pytest.raises(ValueError, match='nblocks'):
            elemc.create_faces_from_elems_multiblock(nblocks)

    @pytest.mark.parametrize(
        'elemtype',
        ['quad4', 'tetra4', 'pyra5', 'prism6', 'hexa8', 'hexa27'],
    )
    def test_all_supported_elements_transform_cells_to_faces(self, elemtype):
        elemc = self.two_cells_sharing_first_face(elemtype)
        transformed = elemc.create_faces_from_elems()
        internalfaces, iface2cell, boundaryfaces, bface2cell = transformed
        assert internalfaces.nelem == 1
        assert np.array_equal(iface2cell.conn, [[0, 1]])
        assert boundaryfaces.nelem == sum(
            2 * len(faces) for faces in _elem.elem2faces[elemtype].values()
        ) - 2
        assert np.all(bface2cell.conn[:, 1] == -1)
        face_type, local_face = next(
            (face_type, face)
            for face_type, faces in _elem.elem2faces[elemtype].items()
            for face in faces
        )
        assert np.array_equal(
            internalfaces[face_type][0],
            np.arange(_elem.nnode_elem[elemtype])[local_face],
        )

    def test_high_order_faces_keep_linear_corner_matching_rule(self):
        elemc = self.two_cells_sharing_first_face('hexa27', corners_only=True)
        legacy = elemc._create_faces_from_elems_legacy()
        transformed = elemc.create_faces_from_elems()
        self.assert_results_equal(legacy, transformed)
        assert transformed[0].nelem == 1

    def test_mixed_elements_can_share_a_face(self):
        hexa = np.arange(_elem.nnode_elem['hexa8'])
        prism = np.arange(20, 20 + _elem.nnode_elem['prism6'])
        hexa_face = _elem.elem2faces['hexa8']['quad4'][0]
        prism_face = _elem.elem2faces['prism6']['quad4'][0]
        prism[prism_face] = hexa[hexa_face]
        elemc = _conn.elem_connectivity()
        elemc.add_elems('hexa8', hexa[np.newaxis, :])
        elemc.add_elems('prism6', prism[np.newaxis, :])
        transformed = elemc.create_faces_from_elems()
        assert transformed[0].nelem == 1
        assert np.array_equal(transformed[1].conn, [[0, 1]])

    @pytest.mark.parametrize(
        'elemtype',
        ['quad4', 'tetra4', 'pyra5', 'prism6', 'hexa8', 'hexa27'],
    )
    def test_single_cell_has_boundary_faces_only(self, elemtype):
        elemc = _conn.elem_connectivity()
        elemc.add_elems(
            elemtype,
            np.arange(_elem.nnode_elem[elemtype])[np.newaxis, :],
        )
        internalfaces, iface2cell, boundaryfaces, bface2cell = (
            elemc.create_faces_from_elems()
        )
        assert internalfaces.nelem == 0
        assert iface2cell.nelem == 0
        expected_faces = sum(
            len(faces) for faces in _elem.elem2faces[elemtype].values()
        )
        assert boundaryfaces.nelem == expected_faces
        assert np.array_equal(bface2cell.conn, [[0, -1]] * expected_faces)

    def test_non_manifold_face_is_rejected(self):
        elemc = _conn.elem_connectivity()
        elemc.add_elems('quad4', np.tile(np.arange(4), (3, 1)))
        with pytest.raises(ValueError, match='non-manifold face'):
            elemc.create_faces_from_elems()
        with pytest.raises(ValueError, match='non-manifold face'):
            elemc.create_faces_from_elems_multiblock(2)

    def test_cell_indices_are_preserved_in_face_neighbours(self):
        elemc = _conn.elem_connectivity()
        cells = self.dict_basiccon['2hexas'][1]
        elemc.add_elems(
            'hexa8',
            cells,
            _conn.indexlist(ilist=[10, 42]),
        )
        original = cells.copy()
        internalfaces, iface2cell, _, bface2cell = elemc.create_faces_from_elems()
        assert internalfaces.nelem == 1
        assert set(iface2cell.conn[0]) == {10, 42}
        assert set(bface2cell.conn[:, 0]) == {10, 42}
        assert np.array_equal(elemc['hexa8'], original)
        self.assert_results_equal(
            elemc.create_faces_from_elems(),
            elemc.create_faces_from_elems_multiblock(2),
        )

    def test_face_transformation_is_deterministic(self):
        elemc = self.set_elemcon('3quadsintri')
        first = elemc.create_faces_from_elems()
        second = elemc.create_faces_from_elems()
        self.assert_results_equal(first, second)

    @pytest.mark.parametrize('econ', ['3quadsintri', '2hexas'])
    def test_elem_merge(self, econ):
        elemc = self.set_elemcon(econ)
        intfaces, _, boundaryfaces, _ = elemc.create_faces_from_elems()
        mergedfaces = _conn.elem_connectivity()
        mergedfaces.importfrom_merge((boundaryfaces, intfaces))
        for ec in [boundaryfaces, intfaces, mergedfaces]:
            ec.print()
        mergedfaces.check()

    def test_elemtocompress(self):
        elemc = _conn.elem_connectivity()
        elemc.add_elems(*self.dict_basiccon['2quads'])  # type and connectivy
        zconn = elemc.exportto_compressedindex()
        assert np.all(zconn._index == [0, 4, 8])

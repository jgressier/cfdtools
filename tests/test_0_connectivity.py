import cfdtools.meshbase._connectivity as _conn
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


def test_iter_and_list_conversion():
    g = _conn.GeneralizedIndex([range(3), [10, 11], 20])
    assert list(g) == [0, 1, 2, 10, 11, 20]

def test_indexing():
    g = _conn.GeneralizedIndex([range(3), [10, 11], 20])
    assert g[0] == 0
    assert g[3] == 10
    assert g[5] == 20

def test_negative_indexing():
    g = _conn.GeneralizedIndex([1, 2, 3, 4])
    assert g[-1] == 4
    assert g[-2] == 3
    assert g[-4] == 1
    with pytest.raises(IndexError):
        _ = g[-5]

def test_length():
    g = _conn.GeneralizedIndex([range(3), [7], [8, 9]])
    assert len(g) == 6

def test_list_method():
    g = _conn.GeneralizedIndex([range(2), [99]])
    assert g.list() == [0, 1, 99]

def test_append():
    g = _conn.GeneralizedIndex([1, 2])
    g.append(3)
    assert list(g) == [1, 2, 3]

def test_extend_with_list():
    g = _conn.GeneralizedIndex([1])
    g.extend([2, 3])
    assert list(g) == [1, 2, 3]

def test_extend_with_range():
    g = _conn.GeneralizedIndex([5])
    g.extend(range(6, 8))
    assert list(g) == [5, 6, 7]

def test_invalid_input_raises():
    with pytest.raises(ValueError):
        _conn.GeneralizedIndex(123)  # not list or range

    with pytest.raises(ValueError):
        _conn.GeneralizedIndex([object()])  # not int, list or range

def test_out_of_bounds_index():
    g = _conn.GeneralizedIndex([1, 2, 3])
    with pytest.raises(IndexError):
        _ = g[10]
    with pytest.raises(IndexError):
        _ = g[-4]

def test_repr():
    g = _conn.GeneralizedIndex([1, 2, 3])
    assert "GeneralizedIndex" in repr(g)
    assert repr(g).endswith("[1, 2, 3])")

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

    @pytest.mark.parametrize('econ', ['2quads', '2hexas'])
    def test_elem_init(self, econ):
        elemc = self.set_elemcon(econ)
        assert elemc.check()

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

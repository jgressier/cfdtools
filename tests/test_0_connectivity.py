import cfdtools.meshbase._connectivity as conn
import numpy as np
import pytest

def test_indexlist_list():
    ilist = [ 10, 12, 15, 20]
    iconn = conn.indexlist()
    iconn.set_list(ilist)
    assert iconn.list() == ilist

def test_indexlist_range():
    imin, imax = 10, 20
    iconn = conn.indexlist()
    iconn.set_range(imin, imax)
    assert iconn.range() == [imin, imax]
    assert iconn.list() == list(range(imin, imax+1))

def test_indexindirection_void():
    dconn = conn.indexindirection()
    assert dconn.nelem == 0

def test_indexindirection_init():
    i = np.arange(10)
    j = np.arange(10)*2
    dconn = conn.indexindirection(array=np.array([i, j]).T)
    assert np.all(dconn.conn[:,1] == j)
    assert np.all(dconn[:,1] == j)

def test_indexindirection_append():
    i = np.arange(10)
    j = np.arange(10)*2
    dconn = conn.indexindirection(array=np.array([i, j]).T)
    i = np.arange(10)+10
    j = np.arange(10)*3
    dconn.append(np.array([i, j]).T)
    assert dconn.nelem == 20

def test_compressed_listofindex():
    i = np.arange(11)*3 # 10 elements of 3 values, add last index 30
    v = 100+np.arange(30)
    zcon = conn.compressed_listofindex(i, v)
    assert zcon.check()

class TestElem():
    dict_basiccon = { 
        '2quads': ('quad4', np.array([[0, 1, 2, 3], [2, 1, 4, 5]])),
    }

    @pytest.mark.parametrize('econ', ['2quads'])
    def test_elem_init(self, econ):
        elemc = conn.elem_connectivity()
        elemc.add_elems(*self.dict_basiccon[econ]) # type and connectivy
        assert elemc.check()
        return elemc

    @pytest.mark.parametrize('econ', ['2quads'])
    def test_elem_createface(self, econ):
        elemc = self.test_elem_init(econ)
        faces = elemc.create_faces_from_elems()
        intfaces, iface2cell, boundaryfaces, bface2cell = conn.find_duplicates(faces)
        assert iface2cell.nelem == intfaces.nelem
        assert bface2cell.nelem == boundaryfaces.nelem
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

#def test_compressed_listofindex():

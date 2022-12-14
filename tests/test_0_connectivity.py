import cfdtools.meshbase._connectivity as conn
import pytest

def test_singleindex_list():
    ilist = [ 10, 12, 15, 20]
    iconn = conn.singleindex()
    iconn.set_list(ilist)
    assert iconn.list() == ilist

def test_singleindex_range():
    imin, imax = 10, 20
    iconn = conn.singleindex()
    iconn.set_range(imin, imax)
    assert iconn.range() == [imin, imax]
    assert iconn.list() == list(range(imin, imax+1))
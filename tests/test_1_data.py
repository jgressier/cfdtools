import cfdtools.meshbase._data as _data
import numpy as np
import pytest

def test_DataSetBase():
    uv = _data.DataSetBase()
    assert uv.Xrep == 'cellaverage'
    assert uv.Trep == 'instant'
    assert uv.ndof == 1

def test_DataSet_default():
    uv = _data.DataSet()
    assert uv.Xrep == 'cellaverage'
    assert uv.Trep == 'instant'
    assert uv.ndof == 1
    uv.add_data('U', 10.)
    uv.add_data('V', 15.)
    for var in ('U', 'V'):
        assert var in uv.keys()
    assert uv['U'] == 10.
    sum = np.sum([val for _, val in uv.items()])
    assert sum == pytest.approx(25.)

def test_DataSetList_default():
    uv = _data.DataSetList(10)
    assert uv.Xrep == 'cellaverage'
    assert uv.Trep == 'timeevol'
    assert uv.ndof == 1

# @pytest.mark.parametrize("filename", ["cavity-degen.hdf"])
# def test_reader(filename):
#     input = cgns.cgnsMesh(_datadir.joinpath(filename))
#     input.read_data()
#     rmesh = input.export_mesh()
#     assert rmesh.check()


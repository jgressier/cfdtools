import cfdtools.data as data
import numpy as np
import pytest

def test_DataSetBase():
    uv = data.DataSetBase()
    assert uv.Xrep == 'cellaverage'
    assert uv.Trep == 'instant'
    assert uv.ndof == 1

def test_DataSet_default():
    uv = data.DataSet()
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

def test_DataSet_ndof():
    uv = data.DataSet(Xrep='spectralcell', ndof=4)
    assert uv.Xrep == 'spectralcell'
    assert uv.Trep == 'instant'
    assert uv.ndof == 4
    uv.add_data('U', 10.) # not consistent data, not tested
    uv.add_data('V', 15.)
    for var in ('U', 'V'):
        assert var in uv.keys()
    assert uv['U'] == 10.
    sum = np.sum([val for _, val in uv.items()])
    assert sum == pytest.approx(25.)

def test_DataSet_3Dline():
    line = data.DataSet(Xrep='nodal')
    assert line.Xrep == 'nodal'
    assert line.Trep == 'instant'
    assert line.ndof == 1
    x = np.linspace(100., 200., 501)
    line.set_nodes(x)
    line.add_data('U', 10+np.sin(2*np.pi*x))
    assert line.geoprop['x'].min() == pytest.approx(100.)
    assert line.geoprop['x'].max() == pytest.approx(200.)
    assert np.average(line.data['U']) == pytest.approx(10.)

class TestDataSet_timeevol():

    def uv_dataset(self, n)-> data.DataSet:
        t = np.linspace(0., 1., n)
        uv = data.DataSet(Trep='timeevol')
        uv.add_data('time', t)
        uv.add_data('U', 10+np.sin(2*np.pi*100*t))
        uv.add_data('V', np.sin(2*np.pi*(50+20*t)*t))
        return uv

    def test_set(self):
        n = 10001
        uv = self.uv_dataset(n)
        dtavg, dtdev, nn = uv.dtstats()
        assert n == nn
        assert dtavg == pytest.approx(1./(n-1))
        assert dtdev == pytest.approx(0.)

    def test_spectrum(self):
        n = 10001
        uv = self.uv_dataset(n)
        sp = uv.dataSet_spectrum()
        assert sp.Trep == 'spectrum'

    def test_spectrogram(self):
        n = 10001
        uv = self.uv_dataset(n)
        sp = uv.dataSet_spectrogram()
        assert sp.Trep == 'spectrogram'
        assert 'time' in sp.keys()

def test_DataSetList_default():
    uv = data.DataSetList(10)
    assert uv.Xrep == 'cellaverage'
    assert uv.Trep == 'timeevol'
    assert uv.ndof == 1

# @pytest.mark.parametrize("filename", ["cavity-degen.hdf"])
# def test_reader(filename):
#     input = cgns.cgnsMesh(_datadir.joinpath(filename))
#     input.read_data()
#     rmesh = input.export_mesh()
#     assert rmesh.check()


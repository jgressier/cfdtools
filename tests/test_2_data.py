import cfdtools.data as data
import numpy as np
#import pytest


# class TestDataSet_import():

#     def __init__(self) -> None:
#         pass
    
#     def uv_dataset(self, n)-> data.DataSet:
#         t = np.linspace(0., 1., n)
#         uv = data.DataSet(Trep='timeevol')
#         uv.add_data('time', t)
#         uv.add_data('U', 10+np.sin(2*np.pi*100*t))
#         uv.add_data('V', np.sin(2*np.pi*(50+20*t)*t))
#         return uv

#     def test_set(self):
#         n = 10001
#         uv = self.uv_dataset(n)
#         dtavg, dtdev, nn = uv.dtstats()
#         assert n == nn
#         assert dtavg == pytest.approx(1./(n-1))
#         assert dtdev == pytest.approx(0.)

#     def test_spectrum(self):
#         n = 10001
#         uv = self.uv_dataset(n)
#         sp = uv.dataSet_spectrum()
#         assert sp.Trep == 'spectrum'

#     def test_spectrogram(self):
#         n = 10001
#         uv = self.uv_dataset(n)
#         sp = uv.dataSet_spectrogram()
#         assert sp.Trep == 'spectrogram'
#         assert 'time' in sp.keys()

# class TestDataSetList_import():

#     def __init__(self) -> None:
#         pass
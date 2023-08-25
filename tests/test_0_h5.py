from cfdtools import __version__
import cfdtools.hdf5._hdf5 as _hdf5
import cfdtools.hdf5 as hdf5 # check user availability
import numpy as np
import pytest

def test_WriteOpen(builddir):
    for i, dtype in enumerate(_hdf5._available_types):
        pname = builddir / f"test{i}.h5"
        f = hdf5.h5File(pname)
        f.open(mode="w", datatype=dtype, version=0) # version is dummy but needed
        f.close()
        g = hdf5.h5File(pname)
        g.open() # default is 'r'
        a = g._h5file.attrs
        assert a['cfdtools_version'] == __version__
        assert a['cfd_datatype'] == dtype
        assert a['data_version'] == 0
        assert g.datatype == dtype
        g.printinfo()
        g.close()

def test_readaccess(datadir):
    f = hdf5.h5File(datadir / "cavity-degen.hdf")
    f.open()
    f.printinfo()
    x = f["/BASE#1/cavity-degen/GridCoordinates/CoordinateX/ data"]
    assert np.mean(x) == pytest.approx(1.166317)
    f.close()

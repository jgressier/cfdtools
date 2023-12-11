import pathlib

import h5py
import numpy as np
import pytest

import cfdtools.cgns as cgns
import cfdtools.ic3.writerV4 as ic3writer
import cfdtools.data as _data


@pytest.mark.parametrize(
    'filename',
    [
        'cavity-degen.hdf',
        'cavity-degen-facebc.hdf',
    ],
)
def test_convert_ic3(datadir: pathlib.Path, builddir: pathlib.Path, filename: str):
    cgmesh = cgns.cgnsMesh(datadir / filename)
    cgmesh.read_data()
    rmesh = cgmesh.export_mesh()

    celldata = _data.DataSet('cellaverage')
    celldata.add_data('RHO', np.ones(rmesh.ncell))
    celldata.add_data('RHOU', np.empty((3 * rmesh.ncell)).reshape(-1, 3))
    celldata.add_data('RHOE', 5 * np.ones(rmesh.ncell))
    rmesh.set_celldata(celldata)
    assert rmesh.check()

    ic3write = ic3writer.writer(rmesh)
    outmesh = builddir / filename
    outmesh = outmesh.with_suffix('.h5')
    outsol = pathlib.Path(str(outmesh.with_suffix('')) + '_solution.h5')
    ic3write.write_data(str(outmesh), str(outsol))

    with h5py.File(outmesh, 'r') as fid:
        assert fid.attrs['cv_count'] == 16
        assert fid.attrs['no_count'] == 27
        assert fid.attrs['fa_count'] == 42

        assert fid['coordinates'].shape[0] == 81
        assert fid['coordinates'].dtype == 'f8'
        assert fid['Connectivities']['Face']['cvofa'].shape[0] == 84
        assert fid['Connectivities']['Face']['cvofa'].dtype == 'i4'
        assert fid['Connectivities']['Face']['noofa_i'].shape[0] == 42

    with h5py.File(outsol, 'r') as fid:
        assert fid.attrs['step'] == 0
        assert fid.attrs['time'] == 0
        assert fid['rho'].shape[0] == 16
        assert fid['rhou'].shape[0] == 48
        assert all(fid['rhoe'][()] == [5.0] * 16)

    outmesh.unlink()
    outsol.unlink()

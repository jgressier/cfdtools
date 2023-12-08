import pathlib

import h5py
import numpy as np
import pytest

import cfdtools.cgns as cgns
import cfdtools.ic3.readerV4 as ic3reader
import cfdtools.ic3.writerV4 as ic3writer
import cfdtools.data as _data


def create_data(datadir: pathlib.Path, builddir: pathlib.Path, filename: str):
    """Create data for tests."""
    cgmesh = cgns.cgnsMesh(datadir / filename)
    cgmesh.read_data()
    rmesh = cgmesh.export_mesh()

    celldata = _data.DataSet('cellaverage')
    celldata.add_data('RHO', np.ones(rmesh.ncell))
    celldata.add_data('RHOU', np.empty((3 * rmesh.ncell)).reshape(-1, 3))
    celldata.add_data('RHOE', 5 * np.ones(rmesh.ncell))
    celldata.add_data('H_AVG', np.ones(rmesh.ncell))
    rmesh.set_celldata(celldata)
    assert rmesh.check()

    return rmesh


@pytest.mark.parametrize(
    'filename',
    [
        'cavity-degen.hdf',
        'cavity-degen-facebc.hdf',
    ],
)
def test_convert_ic3(datadir: pathlib.Path, builddir: pathlib.Path, filename: str):
    rmesh = create_data(datadir, builddir, filename)

    outmesh = builddir / filename
    outmesh = outmesh.with_suffix('.h5')
    outsol = pathlib.Path(str(outmesh.with_suffix('')) + '_solution.h5')
    ic3write = ic3writer.writer(rmesh)
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


@pytest.mark.parametrize(
    'filename',
    [
        'cavity-degen.hdf',
        'cavity-degen-facebc.hdf',
    ],
)
def test_read_ic3v4_mesh(datadir: pathlib.Path, builddir: pathlib.Path, filename: str):
    rmesh = create_data(datadir, builddir, filename)

    outmesh = builddir / filename
    outmesh = outmesh.with_suffix('.h5')
    outsol = pathlib.Path(str(outmesh.with_suffix('')) + '_solution.h5')
    ic3write = ic3writer.writer(rmesh)
    ic3write.write_data(str(outmesh), str(outsol))

    ic3read = ic3reader.reader(str(outmesh), str(outsol))
    ic3read.read_data()
    ic3read.printinfo()
    rmesh = ic3read.export_mesh()
    assert rmesh.check()

    outmesh.unlink()
    outsol.unlink()

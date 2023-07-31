import cfdtools.cgns as cgns
import cfdtools.ic3.writerV3 as ic3writer
import cfdtools.api as api
from pathlib import Path
import pytest

_datadir = Path("./tests/data")
_builddir = Path("./tests/build")


@pytest.mark.parametrize(
    "filename",
    [
        "cavity-degen.hdf",
        "cavity-degen-facebc.hdf",
    ],
)
def test_reader(filename):
    cgmesh = cgns.cgnsMesh(_datadir / filename)
    cgmesh.read_data()
    rmesh = cgmesh.export_mesh()
    assert rmesh.check()


@pytest.mark.parametrize(
    "filename",
    [
        "cavity-degen.hdf",
        "cavity-degen-facebc.hdf",
    ],
)
def test_convert_ic3(filename):
    cgmesh = cgns.cgnsMesh(_datadir / filename)
    cgmesh.read_data()
    rmesh = cgmesh.export_mesh()
    assert rmesh.check()
    ic3write = ic3writer.writer(rmesh)
    _builddir.mkdir(exist_ok=True)
    outfile = api._files(_builddir / filename)
    outfile.change_suffix('.ic3')
    ic3write.write_data(outfile.filename)
    Path(outfile.filename).unlink()

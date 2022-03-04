import pytest
import cfdtools.ic3.readerV2 as ic3readv2
import cfdtools.ic3.writerV2 as ic3writv2
from pathlib import Path
import filecmp

_datadir=Path("./tests/data")
_builddir=Path("./tests/build")

@pytest.mark.parametrize("filename", ["Box3x3x2v2.ic3", "Box3x3x2v3.ic3"])
def test_reader(filename):
    ic3mesh = ic3readv2.reader(_datadir.joinpath(filename))
    rmesh = ic3mesh.read_data()
    assert rmesh.check()

@pytest.mark.parametrize("filename", ["Box3x3x2v2.ic3"])
def test_writer_v2_litend(filename):
    _builddir.mkdir(exist_ok = True)
    basefile = _datadir / filename
    outfile = _builddir / filename
    ic3read = ic3readv2.reader(basefile)
    rmesh = ic3read.read_data()
    assert rmesh.check()
    ic3write = ic3writv2.writer(rmesh, endian='little')
    ic3write.write_data(outfile)
    assert filecmp.cmp(basefile, outfile)

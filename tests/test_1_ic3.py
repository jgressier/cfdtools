import pytest
import cfdtools.ic3.reader_legacy as ic3reader
import cfdtools.ic3.writerV2 as ic3wv2
# import cfdtools.ic3.writerV3 as ic3wv3
from pathlib import Path
import filecmp


@pytest.mark.parametrize(
    "filename", ["Box3x3x2v2.ic3", "Box3x3x2v3.ic3", "nrg-tinycube-v2.ic3"]
)
def test_reader(datadir, filename):
    ic3mesh = ic3reader.reader(datadir / filename, cIntegrity=True)
    ic3mesh.read_data()
    ic3mesh.printinfo()
    rmesh = ic3mesh.export_mesh()
    assert rmesh.check()


@pytest.mark.parametrize("filename", ["Box3x3x2v2.ic3", "nrg-tinycube-v2.ic3"])
def test_writer_v2_litend(datadir, builddir, filename):
    builddir.mkdir(exist_ok=True)
    basefile = datadir / filename
    outfile = builddir / filename
    ic3read = ic3reader.reader(basefile)
    ic3read.read_data()
    rmesh = ic3read.export_mesh()
    assert rmesh.check()
    ic3write = ic3wv2.writer(rmesh, endian='little')
    ic3write.write_data(outfile)
    assert filecmp.cmp(basefile, outfile)
    outfile.unlink()


# @pytest.mark.parametrize("filename", ["sam_sd3.ic3"])
# def test_writer_v3_litend(filename):
#     _builddir.mkdir(exist_ok = True)
#     basefile = _datadir / filename
#     outfile = _builddir / filename
#     ic3read = ic3reader.reader(basefile)
#     rmesh = ic3read.read_data()
#     assert rmesh.check()
#     ic3write = ic3wv3.writer(rmesh, endian='little')
#     ic3write.write_data(outfile)
#     assert 1 # safe run
#     #assert filecmp.cmp(basefile, outfile)

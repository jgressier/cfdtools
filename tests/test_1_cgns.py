import cfdtools.cgns as cgns
import cfdtools.ic3.writerV3 as ic3writer
import cfdtools.api as api
from pathlib import Path
import pytest

_datadir=Path("./tests/data")
_builddir=Path("./tests/build")

@pytest.mark.parametrize("filename", ["cavity-degen.hdf"])
def test_reader(filename):
    input = cgns.cgnsMesh(_datadir.joinpath(filename))
    input.read_data()
    rmesh = input.export_mesh()
    assert rmesh.check()

# @pytest.mark.parametrize("filename", ["cavity-degen.hdf"
#     input = cgns.reader(_datadir.joinpath(filename))
#     input.read_data()
#     rmesh = input.export_mesh()
#     ic3write = ic3writer.writer(rmesh)
#     _builddir.mkdir(exist_ok=True)
#     outfile = api._files(_builddir / Path(filename))
#     outfile.change_suffix('.ic3')
#     ic3write.write_data(outfile.filename)
# # def test_reader3dv22():
#     rmesh = gmsh.reader(_datadir+'box3d-v22.msh')
#     rmesh = gmshmesh.read_data()
#     assert rmesh.check()

# def test_reader3dv41():
#     gmshmesh = gmsh.reader(_datadir+'box3d-v41.msh')
#     rmesh = gmshmesh.read_data()
#     assert rmesh.check()
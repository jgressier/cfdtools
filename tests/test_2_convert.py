import cfdtools.gmsh as gmsh
import cfdtools.ic3.writerV3 as ic3writer
import cfdtools.api as api
from pathlib import Path
import pytest

_datadir=Path("./tests/data")
_builddir=Path("./tests/build")

@pytest.mark.parametrize("filename", ["box3d-v22.msh", "box3d-v41.msh", "test_3d.msh"])
def test_gmsh_to_ic3(filename):
    input = gmsh.reader(_datadir.joinpath(filename))
    input.read_data()
    rmesh = input.export_mesh()
    ic3write = ic3writer.writer(rmesh)
    outfile = api._files(_builddir / Path(filename))
    outfile.change_suffix('.ic3')
    ic3write.write_data(outfile.filename)

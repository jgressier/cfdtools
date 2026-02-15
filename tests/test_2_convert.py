import cfdtools.gmsh as gmsh
import cfdtools.ic3.writerV3 as ic3writer
import cfdtools.api as api
from pathlib import Path
import pytest


@pytest.mark.parametrize("filename", ["box3d-v22.msh", "box3d-v41.msh", "test_3d.msh", "multi.msh"])
def test_convert_ic3(datadir, builddir, filename):
    gmesh = gmsh.reader(datadir / filename)
    gmesh.read_data()
    rmesh = gmesh.export_mesh()
    assert rmesh.check()
    ic3write = ic3writer.writer(rmesh)
    builddir.mkdir(exist_ok=True)
    outfile = api._files(builddir / filename)
    outfile.change_suffix('.ic3')
    ic3write.write_data(outfile.filename)
    Path(outfile.filename).unlink()

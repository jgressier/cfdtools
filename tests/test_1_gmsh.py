import cfdtools.gmsh as gmsh
from pathlib import Path
import pytest

_datadir = Path("./tests/data")
_builddir = Path("./tests/build")


@pytest.mark.parametrize(
    "filename",
    [
        "box3d-v22.msh",
        "box3d-v41.msh",
        "test_2d.msh",
        "test_2d_small.msh",
        "test_3d.msh",
    ],
)
def test_reader(filename):
    gmesh = gmsh.reader(_datadir / filename)
    gmesh.read_data()
    rmesh = gmesh.export_mesh()
    assert rmesh.check()

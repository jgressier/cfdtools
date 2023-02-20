import cfdtools.gmsh as gmsh
from pathlib import Path
import pytest

_datadir=Path("./tests/data")

@pytest.mark.parametrize("filename", ["box3d-v22.msh", "box3d-v41.msh", "test_2d.msh", "test_2d_small.msh", "test_3d.msh"])
def test_reader(filename):
    input = gmsh.reader(_datadir.joinpath(filename))
    input.read_data()
    rmesh = input.export_mesh()
    assert rmesh.check()

# def test_reader3dv22():
#     rmesh = gmsh.reader(_datadir+'box3d-v22.msh')
#     rmesh = gmshmesh.read_data()
#     assert rmesh.check()

# def test_reader3dv41():
#     gmshmesh = gmsh.reader(_datadir+'box3d-v41.msh')
#     rmesh = gmshmesh.read_data()
#     assert rmesh.check()
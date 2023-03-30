import cfdtools.probes as prb
from pathlib import Path
import pytest

_datadir = Path("./tests/data")

# @pytest.mark.parametrize("filename", ["box3d-v22.msh", "box3d-v41.msh", "test_2d.msh", "test_2d_small.msh", "test_3d.msh"])
# def test_reader(filename):
#     input = gmsh.reader(_datadir.joinpath(filename))
#     input.read_data()
#     rmesh = input.export_mesh()
#     assert rmesh.check()

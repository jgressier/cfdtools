# import pytest
import cfdtools.meshbase.simple as simplemesh
from pathlib import Path
# import filecmp

# _datadir = Path("./tests/data")
# _builddir = Path("./tests/build")


def test_cube_mini():
    cube = simplemesh.Cube(2, 2, 2)
    rmesh = cube.export_mesh()
    rmesh.printinfo()
    assert rmesh.check()


def test_cube_large():
    cube = simplemesh.Cube(100, 100, 100)
    rmesh = cube.export_mesh()
    rmesh.printinfo()
    assert rmesh.check()

import cfdtools._cli as cli
from pathlib import Path
import pytest

_datadir = Path("./tests/data")
_builddir = Path("./tests/build")

def outdirlist(nameaslist):
    return nameaslist + ["--outpath", str(_builddir)]


@pytest.mark.parametrize(
    "filename",
    [
        "./tests/data/Box3x3x2v2.ic3",
        "./tests/data/nrg-o1hllc-vf-v3.ic3",
    ],
)
def test_info(filename):
    assert cli.info([filename])


def test_cube_ic3():
    filename = "cube.ic3"
    _builddir.mkdir(exist_ok=True)
    file1 = cli.writecube(outdirlist([filename]))
    assert file1
    Path(file1).unlink()


@pytest.mark.parametrize(
    "filename",
    [
        "./tests/data/Box3x3x2v2.ic3",
        "./tests/data/nrg-tinycube-v2.ic3",
        "./tests/data/nrg-o1hllc-vf-v3.ic3",
    ],
)
def test_ic3brief(filename):
    assert cli.ic3brief([filename])


@pytest.mark.parametrize(
    "filename",
    [
        "./tests/data/Box3x3x2v2.ic3",
    ],
)
def test_ic3writev2(filename):
    _builddir.mkdir(exist_ok=True)
    file1 = cli.write_ic3v2(outdirlist([filename]))
    file2 = cli.write_ic3v2(outdirlist([filename])) # write twice to test safe new name
    assert file1 != filename
    assert file2 != filename
    assert file1 != file2
    Path(file1).unlink()
    Path(file2).unlink()


@pytest.mark.parametrize(
    "filename",
    [
        "./tests/data/Box3x3x2v2.ic3",
        "./tests/data/nrg-o1hllc-vf-v3.ic3",
    ],
)
def test_ic3writev3(filename):
    _builddir.mkdir(exist_ok=True)
    file1 = cli.write_ic3v3(outdirlist([filename]))
    assert file1 != filename
    Path(file1).unlink()


@pytest.mark.parametrize("args", [["--fmt", "CGNS", "./tests/data/cavity-degen.hdf"]])
def test_vtkwrite(args):
    _builddir.mkdir(exist_ok=True)
    file1 = cli.write_vtk(outdirlist(args))
    assert file1
    Path(file1).unlink()

def test_vtkwrite_extrude():
    args = ["--fmt", "CGNS", "--extrude", "5", "./tests/data/cavity-degen.hdf"]
    _builddir.mkdir(exist_ok=True)
    file1 = cli.write_vtk(outdirlist(args))
    assert file1
    Path(file1).unlink()

def test_vtkwrite_scale():
    args = ["--fmt", "CGNS", "--scale", "2.", ".5", "1.", "./tests/data/cavity-degen.hdf"]
    _builddir.mkdir(exist_ok=True)
    file1 = cli.write_vtk(outdirlist(args))
    assert file1
    Path(file1).unlink()

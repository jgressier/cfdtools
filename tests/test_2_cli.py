import cfdtools._cli as cli
from pathlib import Path
import pytest


def outdirlist(nameaslist):
    _builddir = "./tests/build/"
    return nameaslist + ["--outpath", str(_builddir)]


@pytest.mark.parametrize("filename", ["Box3x3x2v2.ic3", "nrg-o1hllc-vf-v3.ic3"])
def test_info(datadir, filename):
    filepath = str(datadir / filename)
    assert cli.cfdinfo([filepath])


def test_cube_ic3(builddir):
    filename = "cube.ic3"
    builddir.mkdir(exist_ok=True)
    file1 = cli.cfdwritecube(outdirlist([filename]))
    assert file1
    Path(file1).unlink()


@pytest.mark.parametrize(
    "filename",
    [
        "Box3x3x2v2.ic3",
        "nrg-tinycube-v2.ic3",
        "nrg-o1hllc-vf-v3.ic3",
    ],
)
def test_ic3brief(datadir, filename):
    filepath = str(datadir / filename)
    assert cli.ic3brief([filepath])


@pytest.mark.parametrize("filename", ["Box3x3x2v2.ic3"])
def test_ic3writev2(datadir, builddir, filename):
    builddir.mkdir(exist_ok=True)
    filepath = str(datadir / filename)
    file1 = cli.cfdwrite_ic3v2(outdirlist([filepath]))
    file2 = cli.cfdwrite_ic3v2(outdirlist([filepath]))  # write twice to test safe new name
    assert file1 != filepath
    assert file2 != filepath
    assert file1 != file2
    Path(file1).unlink()
    Path(file2).unlink()


@pytest.mark.parametrize( "filename", [ "Box3x3x2v2.ic3", "nrg-o1hllc-vf-v3.ic3"])
def test_ic3writev3(datadir, builddir, filename):
    builddir.mkdir(exist_ok=True)
    filepath = str(datadir / filename)
    file1 = cli.cfdwrite_ic3v3(outdirlist([filepath]))
    assert file1 != filepath
    Path(file1).unlink()


def test_vtkbrief(datadir: Path):
    file = "cubemixed0000.vtu"
    assert cli.vtkbrief([str(datadir / file)])


def test_vtkpack(datadir: Path):
    filelist = map(str, sorted(list(datadir.glob("cubemixed*.vtu"))))
    outfile = cli.vtkpack(filelist)
    assert Path(outfile).exists()
    Path(outfile).unlink()



@pytest.mark.parametrize("args", [["--fmt", "CGNS", "./tests/data/cavity-degen.hdf"]])
def test_vtkwrite(builddir, args):
    builddir.mkdir(exist_ok=True)
    file1 = cli.cfdwrite_vtk(outdirlist(args))
    assert file1
    Path(file1).unlink()


def test_vtkwrite_extrude(builddir):
    args = ["--fmt", "CGNS", "--extrude", "5", "./tests/data/cavity-degen.hdf"]
    builddir.mkdir(exist_ok=True)
    file1 = cli.cfdwrite_vtk(outdirlist(args))
    assert file1
    Path(file1).unlink()


def test_vtkwrite_scale(builddir):
    args = ["--fmt", "CGNS", "--scale", "2.", ".5", "1.", "./tests/data/cavity-degen.hdf"]
    builddir.mkdir(exist_ok=True)
    file1 = cli.cfdwrite_vtk(outdirlist(args))
    assert file1
    Path(file1).unlink()

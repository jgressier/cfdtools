from pathlib import Path

import pytest

import cfdtools.api as api
import cfdtools.gmsh as gmsh
import cfdtools.ic3.writerV3 as ic3writer
import cfdtools.ic3.writerV4 as ic3writer_v4


@pytest.mark.parametrize(
    "filename",
    ["box3d-v22.msh", "box3d-v41.msh", "test_3d.msh", "multi.msh"],
)
def test_convert_ic3(datadir: Path, builddir: Path, filename):
    """Writer IC3 V3 For all meshes."""
    gmesh = gmsh.reader(datadir / filename)
    gmesh.read_data()
    rmesh = gmesh.export_mesh()
    assert rmesh.check()
    ic3write = ic3writer.writer(rmesh)
    outfile = api._files(builddir / filename)
    outfile.change_suffix(".ic3")
    ic3write.write_data(outfile.filename)
    Path(outfile.filename).unlink()


@pytest.mark.parametrize("filename", ["mesh3_o2.msh"])
def test_convert_ic3_quadratic(datadir: Path, builddir: Path, filename):
    """Writer IC3 V4 for quadratic meshes."""
    gmesh = gmsh.reader(datadir / filename)
    gmesh.read_data(exclude_center_points=True)
    rmesh = gmesh.export_mesh()
    assert rmesh.check()
    ic3write = ic3writer_v4.writer(rmesh)
    outfile = api._files(builddir / filename)
    outfile.change_suffix(".ic3")
    ic3write.write_data(outfile.filename)
    Path(outfile.filename).unlink()

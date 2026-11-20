from pathlib import Path

import pytest

import cfdtools.api as api
import cfdtools.cgns as cgns
import cfdtools.ic3.writerV3 as ic3writer
import cfdtools.ic3.writerV4 as ic3writer_v4


@pytest.mark.parametrize(
    "filename",
    [
        "cavity-degen.hdf",
        "cavity-degen-facebc.hdf",
        "mesh3_o2.cgns",  # quadratic mesh: HEXA_27 volumes, QUAD_9 boundary faces
    ],
)
def test_reader(datadir, filename):
    cgmesh = cgns.cgnsMesh(datadir / filename)
    cgmesh.read_data()
    rmesh = cgmesh.export_mesh()
    assert rmesh.check()


@pytest.mark.parametrize(
    "filename",
    [
        "cavity-degen.hdf",
        "cavity-degen-facebc.hdf",
        "mesh3_o2.cgns",  # quadratic mesh: HEXA_27 volumes, QUAD_9 boundary faces
    ],
)
def test_convert_ic3(datadir, builddir, filename):
    cgmesh = cgns.cgnsMesh(datadir / filename)
    cgmesh.read_data()
    rmesh = cgmesh.export_mesh()
    assert rmesh.check()
    ic3write = ic3writer.writer(rmesh)
    builddir.mkdir(exist_ok=True)
    outfile = api._files(builddir / filename)
    outfile.change_suffix(".ic3")
    ic3write.write_data(outfile.filename)
    Path(outfile.filename).unlink()


def test_convert_ic3_quadratic(datadir, builddir):
    """CGNS quadratic mesh, center nodes removed, written with the IC3 V4 writer."""
    cgmesh = cgns.cgnsMesh(datadir / "mesh3_o2.cgns")
    cgmesh.read_data(exclude_center_points=True)
    rmesh = cgmesh.export_mesh()
    assert rmesh.check()
    ic3write = ic3writer_v4.writer(rmesh)
    builddir.mkdir(exist_ok=True)
    outfile = api._files(builddir / "mesh3_o2.ic3")
    ic3write.write_data(outfile.filename)
    Path(outfile.filename).unlink()

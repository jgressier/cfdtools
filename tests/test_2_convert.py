from pathlib import Path

import pytest

import cfdtools.api as api
import cfdtools.gmsh as gmsh
import cfdtools.ic3.writerV3 as ic3writer
import cfdtools.ic3.writerV4 as ic3writer_v4


def _prepare_mesh_for_ic3(datadir: Path, filename: str):
    """Read a gmsh mesh ready for IC3 reindexing.

    Runs the face-connectivity / boco steps the IC3 writer performs before
    writing.
    """
    gmesh = gmsh.reader(datadir / filename)
    gmesh.read_data()
    rmesh = gmesh.export_mesh()
    assert rmesh.check()
    rmesh.make_face_connectivity(method="legacy")
    if any(boco.nodebased() for boco in rmesh._bocos.values()):
        rmesh.bocomarks_set_node_to_face()
    return rmesh


@pytest.mark.parametrize("filename", ["box3d-v22.msh", "box3d-v41.msh"])
def test_reindex_boundaryfaces_yields_range_bocos(datadir: Path, filename):
    """Reindexing yields a contiguous range for every boco.

    After node->face conversion and reindexing, every boco must be a
    contiguous range. The IC3 writers rely on this: they assert
    ``boco.index.type == "range"`` and store only the [first, last] bounds.
    """
    rmesh = _prepare_mesh_for_ic3(datadir, filename)
    rmesh.reindex_boundaryfaces()
    assert rmesh._bocos, "expected at least one boco mark"
    for name, boco in rmesh._bocos.items():
        assert boco.index.type == "range", f"boco {name} is not a range"
        lo, hi = boco.index.range()
        assert lo <= hi, f"boco {name} has inverted bounds {lo}..{hi}"


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

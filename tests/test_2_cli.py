import cfdtools._cli as cli
import cfdtools.meshbase.simple as simplemesh
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


@pytest.mark.parametrize(
    'command, writer_class',
    [
        (cli.cfdwrite_ic3v2, cli.ic3.writerV2.writer),
        (cli.cfdwrite_ic3v3, cli.ic3.writerV3.writer),
    ],
)
def test_cfdwrite_ic3_selects_legacy_connectivity(
    monkeypatch, command, writer_class
):
    captured = {}

    def fake_write_generic(
        argv, ext, writer, fname=None, writer_options=None
    ):
        captured.update(
            argv=argv,
            ext=ext,
            writer=writer,
            writer_options=writer_options,
        )
        return 'output.ic3'

    monkeypatch.setattr(cli, 'write_generic', fake_write_generic)
    assert command(['input.mesh']) == 'output.ic3'
    assert captured['ext'] == '.ic3'
    assert captured['writer'] is writer_class
    assert captured['writer_options'] == {
        'face_connectivity_method': 'legacy'
    }


@pytest.mark.parametrize(
    'writer_class',
    [cli.ic3.writerV2.writer, cli.ic3.writerV3.writer],
)
def test_ic3_writer_propagates_legacy_connectivity(monkeypatch, writer_class):
    mesh = simplemesh.Cube(2, 2, 2).export_mesh()
    selected_methods = []
    make_face_connectivity = mesh.make_face_connectivity

    def record_method(method='current', nblocks=None):
        selected_methods.append((method, nblocks))
        return make_face_connectivity(method=method, nblocks=nblocks)

    monkeypatch.setattr(mesh, 'make_face_connectivity', record_method)
    writer_class(mesh, face_connectivity_method='legacy')
    assert selected_methods == [('legacy', None)]


@pytest.mark.parametrize(
    'writer_class',
    [cli.ic3.writerV2.writer, cli.ic3.writerV3.writer],
)
def test_ic3_writer_propagates_multiblock_connectivity(monkeypatch, writer_class):
    mesh = simplemesh.Cube(2, 2, 2).export_mesh()
    selected = []
    make_face_connectivity = mesh.make_face_connectivity

    def record_method(method='current', nblocks=None):
        selected.append((method, nblocks))
        return make_face_connectivity(method=method, nblocks=nblocks)

    monkeypatch.setattr(mesh, 'make_face_connectivity', record_method)
    writer_class(
        mesh, face_connectivity_method='multiblock', face_connectivity_nblocks=2
    )
    assert selected == [('multiblock', 2)]


def test_cfdwrite_ic3v2_multiblock_cli(datadir, builddir):
    """CLI --face-method multiblock --nblocks writes a valid IC3 file."""
    builddir.mkdir(exist_ok=True)
    filepath = str(datadir / "Box3x3x2v2.ic3")
    outpath = str(builddir)
    file1 = cli.cfdwrite_ic3v2(
        [filepath, "--outpath", outpath, "--face-method", "multiblock", "--nblocks", "2"]
    )
    assert file1
    assert Path(file1).exists()
    Path(file1).unlink()


def test_cfdwrite_ic3v2_multiblock_requires_nblocks(datadir, builddir):
    """CLI rejects --face-method multiblock without --nblocks."""
    builddir.mkdir(exist_ok=True)
    filepath = str(datadir / "Box3x3x2v2.ic3")
    with pytest.raises(RuntimeError, match="--nblocks is required"):
        cli.cfdwrite_ic3v2(
            [filepath, "--outpath", str(builddir), "--face-method", "multiblock"]
        )


def test_vtkbrief(datadir: Path):
    file = "cubemixed0000.vtu"
    assert cli.vtkbrief([str(datadir / file)])


def test_vtkdiag(datadir: Path):
    file = "cubemixed0000.vtu"
    assert cli.vtkdiag([str(datadir / file)])


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

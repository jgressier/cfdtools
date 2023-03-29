import cfdtools.api as api
import cfdtools.cgns as cgns
import cfdtools.meshbase._mesh as mesh
from pathlib import Path
import pytest

_datadir=Path("./tests/data")
_builddir=Path("./tests/build")

def test_extrude():
    input = cgns.cgnsMesh(_datadir.joinpath("cavity-degen.hdf"))
    input.read_data()
    mesh2d = input.export_mesh()
    assert mesh2d.check()
    mesh3d = mesh2d.export_extruded(direction=[0., 0., 1.], extrude=[0., .5,  1.])
    assert mesh3d.check()
    mesh3d.printinfo()
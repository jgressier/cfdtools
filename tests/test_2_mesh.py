import cfdtools.api as api
import cfdtools.cgns as cgns
from cfdtools.meshbase._mesh import _default_domain_name
#import pytest


def test_extrude(datadir):
    input = cgns.cgnsMesh(datadir.joinpath("cavity-degen.hdf"))
    input.read_data()
    mesh2d = input.export_mesh()
    assert mesh2d.check()
    mesh3d = mesh2d.export_extruded(direction=[0.0, 0.0, 1.0], extrude=[0.0, 0.5, 1.0])
    assert mesh3d.check()
    mesh3d.printinfo()
    assert mesh3d.get_mark(_default_domain_name+"0") is not None
    assert mesh3d.get_mark(_default_domain_name+"1") is not None

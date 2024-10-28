import numpy as np
import cfdtools.meshbase._mesh as mesh
import cfdtools.meshbase.simple as simplemesh
import pytest


def test_meshconnection_translate():
    vec = [0.0, 0.0, 5.0]
    con = mesh.meshconnection()
    con.set_translation(vec)
    assert np.all(con['translation vector'] == np.array(vec))


def test_meshconnection_rotate():
    con = mesh.meshconnection()
    con.set_rotation('rotz', angle=30.0)
    assert np.all(con['axis'] == np.array([0.0, 0.0, 1.0]))
    assert con['angle'] == 30.0


def test_submeshmark():
    mark = mesh.submeshmark("mymark")
    assert mark.name == "mymark"
    assert mark.geodim is None
    mark.geodim = 'node'
    assert mark.nodebased()
    mark.geodim = 'face'
    assert mark.facebased()


def test_perio_translation_auto():
    cube = simplemesh.Cube(5, 5, 5)
    rmesh = cube.export_mesh()
    meshco = rmesh.build_perio(mark1="imin", mark2='imax')
    assert meshco.transform == 'translate'
    assert meshco.contype == None
    assert np.allclose(meshco['translation vector'], [1., 0, 0.])

def test_perio_rotation():
    cube = simplemesh.Cube(5, 5, 5)
    rmesh = cube.export_mesh()
    # map y, z to r, theta ; rotation is along x
    rmesh.morph(lambda x, y, z: (x, (1+y)*np.cos(np.pi/2*z), (1+y)*np.sin(np.pi/2*z)))
    con = mesh.meshconnection()
    con.set_rotation('rotx', angle=90.0)
    # build_perio uses con:meshconnection and returns a meshconnection filled with index
    meshco = rmesh.build_perio(mark1="kmin", mark2='kmax', connection=con)
    assert meshco.is_rotation()
    assert meshco.contype == None
    assert np.allclose(meshco['axis'], [1., 0, 0.])
    assert meshco['angle'] == pytest.approx(90.)

# def test_unmarked():
#     cube = simplemesh.Cube(2, 2, 2)
#     rmesh = cube.export_mesh()
#     imin_marks = rmesh.pop_mark("imin") # check it exists and pop
#     assert imin_marks is not None
#     print(imin_marks)
#     assert rmesh.make_unmarked_BC(name = "new_faces") # not compatible with node marks
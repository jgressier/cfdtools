import numpy as np
import cfdtools.meshbase._mesh as mesh
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

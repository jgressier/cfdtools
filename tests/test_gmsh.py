import cfdtools.gmsh as gmsh

_datadir="./tests/data/"

def test_reader3dv22():
    rmesh = gmsh.reader(_datadir+'box3d-v22.msh')
    rmesh.read_data()
    assert True

def test_reader3dv41():
    rmesh = gmsh.reader(_datadir+'box3d-v41.msh')
    rmesh.read_data()
    assert True
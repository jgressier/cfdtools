import cfdtools.meshbase.simple as sm
from cfdtools.vtk import vtkMesh, vtkList
from cfdtools.hdf5 import h5File

# import cfdtools.vtk as vtk
# import cfdtools.api as api
from pathlib import Path
import pytest


def test_cube_vtk(builddir):
    filename = "cube.vtu"
    filepath = builddir / filename
    cube = sm.Cube(10, 10, 10)
    mesh = cube.export_mesh()
    vtkmesh = vtkMesh(mesh)
    vtkmesh.write_data(filepath)
    Path(filepath).unlink()


def test_vtkread(datadir):
    name = datadir / "cubemixed0000.vtu"
    vtkfile = vtkMesh()
    vtkfile.read(name)
    assert vtkfile.pyvista_grid.n_cells == 1000


def test_vtkList(datadir):
    namelist = list(datadir.glob("cubemixed00*.vtu"))
    vtklist = vtkList(namelist)
    assert vtklist.nfile == len(namelist)
    assert vtklist.allexist()


def test_vtkList_check(datadir):
    namelist = list(datadir.glob("cubemixed00*.vtu"))
    vtklist = vtkList(namelist, verbose=True)
    assert not vtklist.check_order('cellcenter')  # known to be mixed


def test_vtkList_read(datadir):
    namelist = list(datadir.glob("cubemixed00*.vtu"))
    vtklist = vtkList(namelist, verbose=True)
    vtklist.read()


def test_vtkList_dump(datadir, builddir):
    namelist = list(datadir.glob("cubemixed00*.vtu"))
    vtklist = vtkList(namelist, verbose=True)
    vtklist.read()
    h5filename = builddir / "vtklist.hdf"
    vtklist.dumphdf(h5filename)
    h5file = h5File(h5filename)
    h5file.open()
    assert h5file.datatype == 'datalist'
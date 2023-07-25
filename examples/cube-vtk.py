import cfdtools.meshbase.simple as sm
from cfdtools.vtk import vtkMesh

cube = sm.Cube(10, 10, 10)
mesh = cube.export_mesh()

w = vtkMesh(mesh)
w.write_data("pipo.vtu")
print('done')

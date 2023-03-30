import cfdtools.meshbase.simple as sm
import cfdtools.vtk.writer as vtkw

cube = sm.Cube(10, 10, 10)
mesh = cube.export_mesh()

w = vtkw.writer(mesh)
w.set_mesh()
w.write_mesh("pipo.vtu")
print('done')

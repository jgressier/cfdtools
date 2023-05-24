
import cfdtools.meshbase.simple as sm
from cfdtools.vtk import vtkMesh
import pyvista as pv
import numpy as np

cube = sm.Cube(10, 10, 10)
mesh = cube.export_mesh()

w = vtkMesh(mesh)

def f(xyz, t):
    return np.sin(xyz[:,0]-t)+xyz[:,2]+xyz[:,2]**2

for i in range(10):
    newcon = dict()
    for itype, icon in w.pyvista_grid().cells_dict.items():
        newcon[itype] = np.roll(icon, i*7, axis=0)
    new = pv.UnstructuredGrid(newcon, w._coords)
    ctr = new.cell_centers().points
    new.cell_data['Q'] = f(ctr, i)
    new.save(f"cubemixed{i:04d}.vtu")
print('done')

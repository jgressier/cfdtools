#from cfdtools.meshbase._mesh import meshconnection
import cfdtools.meshbase.simple as sm
from cfdtools import log

import logging
log.setLevel(logging.DEBUG)

# transx = meshconnection()
# transx.set_translation([0.0, 1.0, 0.0])

cube = sm.Cube(5, 5, 5)
mesh = cube.export_mesh()
meshco = mesh.build_perio(mark1="jmin", mark2='jmax')

print(meshco)
mesh.printinfo(detailed=True)

print('done')

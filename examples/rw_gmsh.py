import cfdtools._cli as cli

import cfdtools.gmsh as gmsh
import cfdtools.ic3.writerV3 as ic3writer
import cfdtools.api as api

_datadir="../tests/data/"
_builddir="../tests/build/"

#filename = "small_cube.msh" # 4.1
filename = "test_3d.msh"

api.io.set_modes(api.io._available) # all outputs
api.io.set_modes(api.io._available.remove("debug")) # all outputs but debug

# ic3read = ic3reader.reader("./tests/data/Box3x3x2v3.ic3")
# ic3read = ic3reader.reader("./examples/restart_SD_21.out")
reader = gmsh.reader(_datadir+filename)
reader.read_data()
rmesh = reader.export_mesh()
#reader.printinfo()
rmesh._make_face_connectivity()
rmesh.printinfo()
rmesh.check()

# ic3write = ic3writer.writer(rmesh)
# ic3write.write_data(_builddir+filename)
# print("done")
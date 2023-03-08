import cfdtools._cli as cli

import cfdtools.gmsh as gmsh
import cfdtools.ic3.writerV3 as ic3writer
import cfdtools.api as api

_datadir="./tests/data/"
_builddir="./tests/build/"

api.io.set_modes(api.io._available) # all outputs
api.io.set_modes(api.io._available.remove("debug")) # all outputs but debug

filename = "box3d-v41.msh" # 4.1
filename = "box2x2.msh" # 4.1
filename = "circle-6quads.msh" # 4.1
#filename = "test_2d.msh"
#filename = "naca4412-1M.msh"
#filename = "SingleJet_VF.msh"
# ic3read = ic3reader.reader("./tests/data/Box3x3x2v3.ic3")
# ic3read = ic3reader.reader("./examples/restart_SD_21.out")
reader = gmsh.reader(_datadir+filename)
#reader = gmsh.reader("./examples/"+filename)
reader.read_data()
rmesh = reader.export_mesh()
#reader.printinfo()
#api.io.print('std','PRINT INFO')
rmesh.check()
rmesh.printinfo()
api.io.print('std','CHECK')
rmesh.check()

api.io.print('std','EXPORT TO IC3')
ic3write = ic3writer.writer(rmesh)
outfile = api._files(_builddir+filename)
outfile.change_suffix('.ic3')
ic3write.write_data(outfile.filename)
print("done")
import cfdtools.ic3.reader_legacy as ic3reader
import cfdtools.ic3.writerV3 as ic3writer

_datadir = "./tests/data/"
_builddir = "./tests/build/"

filename = "maxime_sd2.ic3"

# ic3read = ic3reader.reader("./tests/data/Box3x3x2v3.ic3")
# ic3read = ic3reader.reader("./examples/restart_SD_21.out")
ic3read = ic3reader.reader(_datadir + filename)
ic3read.read_data()
ic3read.printinfo()
rmesh = ic3read.export_mesh()
rmesh.printinfo()
ic3write = ic3writer.writer(rmesh)
ic3write.write_data(_builddir + filename)
print("done")

from cfdtools import log
import cfdtools._cli as cli
from cfdtools.ic3 import reader
# cli.info(["examples/restart_gene_21.out"])
# cli.info(["examples/restart_hades_21.out"])
cli.cfdinfo(["tests/data/restart-perio.ic3"])
# x/y/z dimensions = 3/5/2
# periodic faces are y (2x6 faces) and z (2x15 faces) 

r = reader("tests/data/restart-perio.ic3")
log.info("> Reading IC3 file")
r.read_data()
log.info("> export IC3 data to cfdtools mesh")
mesh = r.export_mesh()
log.info("> info")
mesh.printinfo()
for name, bc in mesh._bocos.items():
    log.info(f"> boco {bc.name} - {name}")
    log.info(f">   type: {bc.type}")
    log.info(f">   facebased: {bc.facebased()}")
    log.info(f">   nodebased: {bc.nodebased()}")
    log.info(f">   index: {bc.index}")
    log.info(f">   slice: {bc._properties['slicing']}")

import cfdtools._cli as cli
from cfdtools.ic3 import reader
# cli.info(["examples/restart_gene_21.out"])
# cli.info(["examples/restart_hades_21.out"])
cli.info(["tests/data/restart-perio.ic3"])

r = reader("tests/data/restart-perio.ic3")
r.read_data()

mesh = r.export_mesh()
mesh.printinfo()
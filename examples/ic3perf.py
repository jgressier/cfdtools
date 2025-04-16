from cfdtools import log
import logging
from cfdtools.utils._dev import TraceMemory
import cfdtools.api as api
from cfdtools.ic3.writerV3 import writer as ic3writer
import cfdtools.meshbase.simple as sm

log.setLevel(logging.DEBUG)

for size in [100, 200, 400]:
    with TraceMemory():
        log.info(f"> size {size}")
        with api.Timer(f"create cube {size}", nelem=size**3):
            cube = sm.Cube(size, size, size)
        with api.Timer(f"export cube", nelem=size**3):
            mesh = cube.export_mesh()
            #meshco = mesh.build_perio(mark1="jmin", mark2='jmax')
        with api.Timer(f"write cube", nelem=size**3):
            writer = ic3writer(mesh)
            writer.write_data(f"cube-{size}.ic3")
        log.info(f"done size {size}")
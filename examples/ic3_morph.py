import cfdtools._cli as cli

import cfdtools.ic3.reader_legacy as ic3reader
import cfdtools.ic3.writerV3 as ic3writer
import cfdtools.api as api

import aerokit.aero.degree as deg
import numpy as np
import psutil

api.io.set_modes(api.io._available)  # all outputs
api.io.set_modes(api.io._available.remove("debug"))  # all outputs but debug

mem = psutil.virtual_memory().percent
print(f"memory: {mem}")

ic3read = ic3reader.reader("./tan25D2-185k+data.restart")
ic3read.read_data()
rmesh = ic3read.export_mesh()

mem = psutil.virtual_memory().percent
print(f"memory: {mem}")

# ic3read.printinfo()

rmesh.printinfo()


def fmorph(x, y, z):
    yd, yD, delta = 1.0, 2.0, -0.4
    angref, angnew = 25.0, 15.0
    nx = x * (deg.tan(angref) / deg.tan(angnew))
    ny = y + (y > yd) * (y < 2 * yD - yd) * 0.5 * delta * (
        1.0 - np.cos(np.pi * (y - yd) / (yD - yd))
    )
    nz = z
    return nx, ny, nz


rmesh.morph(fmorph)

ic3write = ic3writer.writer(rmesh)
ic3write.write_data("./new_dwedge.restart")

mem = psutil.virtual_memory().percent
print(f"memory: {mem}")

del ic3read
mem = psutil.virtual_memory().percent
print(f"del ic3read  / memory: {mem}")

del ic3write
mem = psutil.virtual_memory().percent
print(f"del ic3write / memory: {mem}")

print("done")

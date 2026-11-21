# simple subpackage
from itertools import product

import numpy as np

import cfdtools.meshbase._connectivity as _conn
from cfdtools.meshbase import _mesh


class Cube:
    def __init__(self, nx, ny, nz=1) -> None:
        self._nx = nx
        self._ny = ny
        self._nz = nz
        self.ncell = nx * ny * nz
        self.nnode = (nx + 1) * (ny + 1) * (nz + 1)

    def elem_ijk(self):
        return list(product(range(self._nx), range(self._ny), range(self._nz)))

    def nodes_ijk(self):
        return list(
            product(range(self._nx + 1), range(self._ny + 1), range(self._nz + 1))
        )

    def nodes_xyz(self, i, j, k):
        return (i / float(self._nx), j / float(self._ny), k / float(self._nz))

    def nodeglobindex_ijk(self, i, j, k):
        # k lines, j-k planes
        return k + (self._nz + 1) * (j + (self._ny + 1) * i)

    def elemglobindex_ijk(self, i, j, k):
        # k lines, j-k planes
        return k + self._nz * (j + self._ny * i)

    def _elems(self):
        ei, ej, ek = np.indices((self._nx, self._ny, self._nz), dtype=np.int64)
        ei, ej, ek = (array.ravel() for array in (ei, ej, ek))
        hexa2node = np.column_stack(
            [
                self.nodeglobindex_ijk(ei, ej, ek),
                self.nodeglobindex_ijk(ei + 1, ej, ek),
                self.nodeglobindex_ijk(ei + 1, ej + 1, ek),
                self.nodeglobindex_ijk(ei, ej + 1, ek),
                self.nodeglobindex_ijk(ei, ej, ek + 1),
                self.nodeglobindex_ijk(ei + 1, ej, ek + 1),
                self.nodeglobindex_ijk(ei + 1, ej + 1, ek + 1),
                self.nodeglobindex_ijk(ei, ej + 1, ek + 1),
            ]
        )
        index = self.elemglobindex_ijk(ei, ej, ek)
        return hexa2node, index

    def export_mesh(self):
        meshdata = _mesh.Mesh(ncell=self.ncell, nnode=self.nnode)
        # set cell connectivity
        cell2node = _conn.elem_connectivity()
        hexanode, ielem = self._elems()
        cell2node.add_elems("hexa8", hexanode, _conn.indexlist(ilist=ielem))
        meshdata.set_cell2node(cell2node)
        # set node coordinates
        ni, nj, nk = np.indices(
            (self._nx + 1, self._ny + 1, self._nz + 1), dtype=np.int64
        )
        ni, nj, nk = (array.ravel() for array in (ni, nj, nk))
        x, y, z = self.nodes_xyz(ni, nj, nk)
        meshdata.set_nodescoord_xyz(x, y, z)
        # boco
        boundary_values = {
            "imin": (0, 0),
            "imax": (0, self._nx),
            "jmin": (1, 0),
            "jmax": (1, self._ny),
            "kmin": (2, 0),
            "kmax": (2, self._nz),
        }
        for name, (axis, value) in boundary_values.items():
            bcmark = _mesh.submeshmark(name)
            bcmark.geodim = "bdnode"
            bcmark.type = "boundary"
            ranges = [
                np.arange(self._nx + 1),
                np.arange(self._ny + 1),
                np.arange(self._nz + 1),
            ]
            ranges[axis] = np.array([value])
            ni, nj, nk = (
                array.ravel() for array in np.meshgrid(*ranges, indexing="ij")
            )
            nodes = self.nodeglobindex_ijk(ni, nj, nk)
            bcmark.index = _conn.indexlist(ilist=nodes)
            meshdata.add_boco(bcmark)
        return meshdata


if __name__ == "__main__":
    box = Cube(2, 1, 1)
    mesh = box.export_mesh()
    mesh.printinfo()

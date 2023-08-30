# simple subpackage
import cfdtools.api as api
import cfdtools.meshbase._mesh as _mesh
import cfdtools.meshbase._connectivity as _conn
from itertools import product
import numpy as np


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
        return list(product(range(self._nx + 1), range(self._ny + 1), range(self._nz + 1)))

    def nodes_xyz(self, i, j, k):
        return (i / float(self._nx), j / float(self._ny), k / float(self._nz))

    def nodeglobindex_ijk(self, i, j, k):
        # k lines, j-k planes
        return k + (self._nz + 1) * (j + (self._ny + 1) * i)

    def elemglobindex_ijk(self, i, j, k):
        # k lines, j-k planes
        return k + self._nz * (j + self._ny * i)

    def _elems(self):
        ei, ej, ek = map(np.array, zip(*self.elem_ijk()))
        hexa2node = np.array(
            list(
                zip(
                    *tuple(
                        c.tolist()
                        for c in [
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
                )
            )
        )
        index = self.elemglobindex_ijk(ei, ej, ek)
        return hexa2node, index.tolist()

    def export_mesh(self):
        meshdata = _mesh.Mesh(ncell=self.ncell, nnode=self.nnode)
        # set cell connectivity
        cell2node = _conn.elem_connectivity()
        hexanode, ielem = self._elems()
        cell2node.add_elems('hexa8', np.array(hexanode), _conn.indexlist(ilist=ielem))
        meshdata.set_cell2node(cell2node)
        # set node coordinates
        ni, nj, nk = np.array(self.nodes_ijk()).T
        x, y, z = self.nodes_xyz(ni, nj, nk)
        meshdata.set_nodescoord_xyz(x, y, z)
        # boco
        tasks = {
            'imin': {'ijknodes': product([0], range(self._ny + 1), range(self._nz + 1))},
            'imax': {'ijknodes': product([self._nx], range(self._ny + 1), range(self._nz + 1))},
            'jmin': {'ijknodes': product(range(self._nx + 1), [0], range(self._nz + 1))},
            'jmax': {'ijknodes': product(range(self._nx + 1), [self._ny], range(self._nz + 1))},
            'kmin': {'ijknodes': product(range(self._nx + 1), range(self._ny + 1), [0])},
            'kmax': {'ijknodes': product(range(self._nx + 1), range(self._ny + 1), [self._nz])},
        }
        for name, task in tasks.items():
            bcmark = _mesh.submeshmark(name)
            bcmark.geodim = 'bdnode'
            bcmark.type = 'boundary'
            ni, nj, nk = np.array(list(task['ijknodes'])).T
            nodes = self.nodeglobindex_ijk(ni, nj, nk)
            bcmark.index = _conn.indexlist(ilist=nodes.tolist())
            meshdata.add_boco(bcmark)
        return meshdata


if __name__ == "__main__":
    box = Cube(2, 1, 1)
    mesh = box.export_mesh()
    mesh.printinfo()

from pathlib import Path

import numpy as np
import pytest

import cfdtools.api as api
import cfdtools.cgns as cgns
import cfdtools.meshbase._connectivity as _conn
import cfdtools.meshbase._mesh as _mesh

_datadir = Path("./tests/data")
_builddir = Path("./tests/build")


def test_extrude():
    input = cgns.cgnsMesh(_datadir.joinpath("cavity-degen.hdf"))
    input.read_data()
    mesh2d = input.export_mesh()
    assert mesh2d.check()
    mesh3d = mesh2d.export_extruded(direction=[0.0, 0.0, 1.0], extrude=[0.0, 0.5, 1.0])
    assert mesh3d.check()
    mesh3d.printinfo()


def _make_quad9_2d():
    """Build a single quad9 element on the unit square [0,1]x[0,1] at z=0."""
    # CGNS quad9 node order: corners(0-3), edge mids(4-7), face center(8)
    coords = np.array(
        [
            [0.0, 0.0, 0.0],  # 0
            [1.0, 0.0, 0.0],  # 1
            [1.0, 1.0, 0.0],  # 2
            [0.0, 1.0, 0.0],  # 3
            [0.5, 0.0, 0.0],  # 4 bottom edge mid
            [1.0, 0.5, 0.0],  # 5 right edge mid
            [0.5, 1.0, 0.0],  # 6 top edge mid
            [0.0, 0.5, 0.0],  # 7 left edge mid
            [0.5, 0.5, 0.0],  # 8 face center
        ]
    )
    m = _mesh.Mesh(ncell=1, nnode=9)
    m.set_nodescoord_nd(coords)
    c2n = _conn.elem_connectivity()
    c2n.add_elems("quad9", np.array([[0, 1, 2, 3, 4, 5, 6, 7, 8]]))
    m.set_cell2node(c2n)
    return m


def test_extrude_quadratic():
    mesh2d = _make_quad9_2d()
    mesh3d = mesh2d.export_extruded(direction=[0.0, 0.0, 1.0], extrude=[0.0, 0.5, 1.0])
    # nrange=3 planes -> 2 hexa27 stacked, 45 nodes
    assert mesh3d.ncell == 2
    assert mesh3d.nnode == 45
    # element type must be hexa27
    assert "hexa27" in mesh3d._cell2node.elems()
    con = mesh3d._cell2node["hexa27"]
    assert con.shape == (2, 27)
    # element 0 first node is the bottom corner 0
    assert con[0, 0] == 0
    coords = mesh3d.nodescoord(ndarray=True)
    # volume center of element 0 (CGNS node 26) -> node 17 = mid plane center
    # at (0.5, 0.5, 0.25)
    vol_center = con[0, 26]
    np.testing.assert_allclose(coords[vol_center], [0.5, 0.5, 0.25])
    # top corner of element 0 (CGNS node 4) -> global node 18 = plane1 corner 0
    top_corner = con[0, 4]
    np.testing.assert_allclose(coords[top_corner], [0.0, 0.0, 0.5])
    mesh3d.printinfo()

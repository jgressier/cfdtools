"""Tests for quadratic (P2) mesh reading, shared by the CGNS and gmsh readers.

``mesh3_o2.cgns`` and ``mesh3_o2.msh`` describe the *same* mesh: a unit cube
meshed with 8 HEXA27 cells (125 P2 nodes) and 6 named boundary patches of 25
(5x5) nodes each. Because the two readers are fully independent code paths,
reading both and asserting they agree is a strong correctness check on top of
the per-format content assertions below.

Comparisons are made on order-invariant quantities (sorted point clouds, cell
centroids, per-patch node sets) so they survive the node/cell renumbering that
legitimately differs between the two formats.
"""

import numpy as np
import pytest

import cfdtools.cgns as cgns
import cfdtools.gmsh as gmsh

EXPECTED_BOCOS = {"front", "back", "left", "right", "top", "bottom"}
NCELL = 8
NNODE = 125
NBOCO_NODE = 25  # 5x5 P2 nodes per cube face


def _read(datadir, fmt, exclude_center_points=False):
    """Read mesh3_o2 in the requested format and return the exported Mesh."""
    if fmt == "cgns":
        reader = cgns.cgnsMesh(datadir / "mesh3_o2.cgns")
    elif fmt == "gmsh":
        reader = gmsh.reader(datadir / "mesh3_o2.msh")
    else:
        raise ValueError(f"unknown format {fmt!r}")
    reader.read_data(exclude_center_points=exclude_center_points)
    return reader.export_mesh()


def _points(mesh):
    x, y, z = mesh.nodescoord(ndarray=False)
    return np.column_stack((x, y, z))


def _sorted_points(pts, decimals=9):
    """Lexicographically sort a point cloud so it can be compared up to ordering."""
    pts = np.round(np.asarray(pts), decimals)
    return pts[np.lexsort((pts[:, 2], pts[:, 1], pts[:, 0]))]


def _centroids(mesh):
    """Cell centroids, taken as the mean of the 8 corner nodes of each HEXA27."""
    corners = mesh._cell2node["hexa27"][:, :8]
    return _points(mesh)[corners].mean(axis=1)


def _boco_point_sets(mesh):
    """Map each boundary patch to its (order-invariant) set of node coordinates."""
    pts = _points(mesh)
    return {
        name: _sorted_points(pts[np.asarray(boco.index.list())])
        for name, boco in mesh._bocos.items()
    }


@pytest.mark.parametrize("fmt", ["cgns", "gmsh"])
def test_reader_content(datadir, fmt):
    """The reader returns the expected cells, nodes, element type and patches."""
    mesh = _read(datadir, fmt)
    assert mesh.check()
    assert mesh.ncell == NCELL
    assert mesh.nnode == NNODE
    # a single quadratic volume element type
    assert set(mesh._cell2node.elems()) == {"hexa27"}
    assert mesh._cell2node["hexa27"].shape == (NCELL, 27)
    # six node-based boundary patches
    assert set(mesh._bocos) == EXPECTED_BOCOS
    for boco in mesh._bocos.values():
        assert boco.nodebased()
        assert boco.index.size == NBOCO_NODE
    # geometry: unit cube, with every patch lying in a coordinate plane
    pts = _points(mesh)
    assert np.allclose(pts.min(axis=0), [0.0, 0.0, 0.0])
    assert np.allclose(pts.max(axis=0), [1.0, 1.0, 1.0])
    for patch in _boco_point_sets(mesh).values():
        constant_axes = np.isclose(patch.min(axis=0), patch.max(axis=0))
        assert constant_axes.sum() >= 1, "boundary patch is not planar"


def test_cross_format_equivalence(datadir):
    """CGNS and gmsh describe the same mesh; the two readers must agree."""
    cg = _read(datadir, "cgns")
    gm = _read(datadir, "gmsh")
    assert cg.ncell == gm.ncell
    assert cg.nnode == gm.nnode
    assert set(cg._cell2node.elems()) == set(gm._cell2node.elems())
    # node coordinates, invariant to node renumbering
    assert np.allclose(_sorted_points(_points(cg)), _sorted_points(_points(gm)))
    # cell centroids, invariant to cell ordering
    assert np.allclose(_sorted_points(_centroids(cg)), _sorted_points(_centroids(gm)))
    # boundary patches address the same physical nodes
    cg_bocos, gm_bocos = _boco_point_sets(cg), _boco_point_sets(gm)
    assert set(cg_bocos) == set(gm_bocos)
    for name in cg_bocos:
        assert np.allclose(cg_bocos[name], gm_bocos[name])


@pytest.mark.parametrize("fmt", ["cgns", "gmsh"])
def test_exclude_center_points(datadir, fmt):
    """Dropping HEXA27 centre nodes removes one node per cell and renumbers safely."""
    full = _read(datadir, fmt)
    reduced = _read(datadir, fmt, exclude_center_points=True)
    assert reduced.check()
    # exactly one (interior) node removed per cell, connectivity reduced to 26
    assert reduced.nnode == full.nnode - full.ncell
    assert reduced._cell2node["hexa27"].shape == (NCELL, 26)
    # the removed nodes are the cell centres, preserved as cell data
    point27 = np.asarray(reduced.pop_celldata("point27"))
    assert point27.shape == (NCELL, 3)
    assert np.allclose(_sorted_points(point27), _sorted_points(_centroids(full)))
    # renumbering preserves the boundary patches (centres are interior nodes)
    full_bocos, reduced_bocos = _boco_point_sets(full), _boco_point_sets(reduced)
    assert set(full_bocos) == set(reduced_bocos)
    for name in full_bocos:
        assert np.allclose(full_bocos[name], reduced_bocos[name])

# coding: utf8

# Import modules
import collections
import logging

import numpy as np

import cfdtools.api as api
import cfdtools.data as _data
import cfdtools.meshbase._connectivity as _conn
import cfdtools.meshbase._elements as _elem
import cfdtools.meshbase._mesh as _mesh
from cfdtools.gmsh._gmsh import gmshtype_elem, nodes_per_gmshtype  # , nodes_per_cell

log = logging.getLogger(__name__)


@api.fileformat_reader("GMSH", ".msh")
class reader(api._files):
    """Implementation of the reader to read Gmsh meshes."""

    def __init__(self, filename, cIntegrity=False):
        """Initialize of a GMSH reader.

        param: filename: file name [type string]
        """
        super().__init__(filename)
        self.check_integrity = cIntegrity
        self.celldata = _data.DataSet("cellaverage")

    @property
    def ncell(self):
        return self._ncell

    def read_data(self, exclude_center_points=False):
        log.info(f"> GMSH reader: starts reading {self.filename}")
        # Check file exists
        if not self.exists():
            api.error_stop(f"File not found: {self.filename!r}")

        # Set the node ordering convention
        _elem.elem2faces = _elem.gmsh_elem2faces

        # Read file (streaming for v4, legacy whole-file parse for v2).
        version = self.__peek_version(self.filename)
        if version >= 4:
            (
                families,
                bctype,
                x,
                y,
                z,
                connectivity,
                connect_bc,
                bnd,
                bnd2fam,
                self._maxdim,
                self._ncell,
            ) = self.__read_v4_stream(self.filename)
        elif version <= 2:
            families, bctype, x, y, z, elts = self.__read_v2_stream(self.filename)
            (connectivity, connect_bc, bnd, bnd2fam, self._maxdim, self._ncell) = (
                self.__build_connectivity(elts)
            )
        else:
            api.error_stop(
                f"unexpected mesh version: {version} (expected <= 2 or >= 4)"
            )

        # Check for 3D
        if self._maxdim < 2 or self._maxdim > 3:
            api.error_stop(f"unexpected max element dimension: {self._maxdim}")
        log.info("%dD mesh", self._maxdim)

        # Choose the integer dtype for connectivity (node and cell indices).
        # Signed int32 halves the memory and is safe as long as every node and
        # cell index fits (indices are 0-based, so the largest is count - 1);
        # it must stay *signed* because -1 is used as a sentinel (boundary face
        # neighbour, 27th-point removal). Fall back to int64 above the limit.
        nnode = len(x)
        idtype = (
            np.int32
            if max(nnode, self._ncell) - 1 <= np.iinfo(np.int32).max
            else np.int64
        )
        if idtype != np.int64:
            for elt_type in connectivity:
                connectivity[elt_type] = connectivity[elt_type].astype(
                    idtype, copy=False
                )
            for bnd_tag in connect_bc:
                for elt_type in connect_bc[bnd_tag]:
                    connect_bc[bnd_tag][elt_type] = connect_bc[bnd_tag][
                        elt_type
                    ].astype(idtype, copy=False)

        # Initialize returned variables
        boundaries = {}
        self._nboco = 0

        if exclude_center_points:
            log.info("  exclude center points of volume connectivity")
            x = np.asarray(x)
            y = np.asarray(y)
            z = np.asarray(z)
            x, y, z, point27 = self.__remove_27th_point_hexa27(
                connect_bc, connectivity, x, y, z
            )
            self.celldata.add_data("point27", point27)

        # Fill the dictionary of boundary conditions
        self._famm = []
        if families is not None:  # B.C. are defined
            for bnd_tag in bnd:
                family = families[str(bnd2fam[bnd_tag])]
                if family != "fluid":
                    self.__create_bnd(
                        boundaries,
                        np.array(bnd[bnd_tag]),
                        family,
                        bctype,
                        connect_bc[bnd_tag],
                    )
        else:  # define B.C.
            for bnd_tag in bnd.keys():
                family = "fluid"
                self.__create_bnd(
                    boundaries,
                    np.array(bnd[bnd_tag]),
                    family,
                    bctype,
                    connect_bc[bnd_tag],
                )

        # Add the fluid boundary condition with the internal type now
        boundaries["int_fluid"] = {
            "slicing": np.array([], dtype=np.int64),
            "type": "internal",
            "periodic_transform": np.zeros((16,), dtype=np.float64),
        }
        boundaries["int_fluid"]["slicing"] = np.unique(
            np.concatenate(
                [connectivity[elt_type].ravel() for elt_type in connectivity]
            )
        )

        self._coords = (x, y, z)
        self._cellconnectivity = connectivity
        self._boundaries = boundaries

    def export_mesh(self):
        log.info("Export gmsh mesh to cfdtools mesh data")
        meshdata = _mesh.Mesh(nnode=len(self._coords[0]))
        meshdata.set_nodescoord_xyz(*self._coords)
        # meshdata.set_face2node(self.mesh['connectivity']['noofa'])
        cellconn = _conn.elem_connectivity()
        # extract cell connectivity only but not face connectivity
        for etype, econn in self._cellconnectivity.items():
            if _elem.dim_elem[etype] == self._maxdim:
                cellconn.add_elems(etype, econn)
        meshdata.set_cell2node(cellconn)

        for name, bc_dict in self._boundaries.items():
            if bc_dict["type"] == "boundary":
                boco = _mesh.submeshmark(name)
                boco.geodim = "bdnode"
                boco.properties["type"] = bc_dict["type"]
                boco.properties["periodic_transform"] = bc_dict["periodic_transform"]
                boco.index = _conn.indexlist(ilist=bc_dict["slicing"])
                meshdata.add_boco(boco)
        meshdata.set_celldata(self.celldata)
        # meshdata.set_nodedata(self.variables['nodes'])
        # meshdata.set_facedata(self.variables['faces'])
        # meshdata.set_params(self.mesh['params'])
        return meshdata

    def __create_bnd(self, boundaries, window, family, bctype, connectivity):
        fff = family
        self._nboco += 1

        # Loop to concatenate multiple faces/entities into the same boundary (For complex geometries)

        if fff in self._famm:
            bc = boundaries[family]
            bc["type"] = "None"
            bc["periodic_transform"] = "None"
            # merge new node ids into the existing slicing array
            new_nodes = np.concatenate(
                [connectivity[elt_type].ravel() for elt_type in connectivity]
            )
            bc["slicing"] = np.unique(np.concatenate([bc["slicing"], new_nodes]))
            bc["type"] = "boundary"
            bc["periodic_transform"] = np.zeros((16,), dtype=np.float64)

        else:
            boundaries[family] = {
                "slicing": None,
                "type": None,
                "periodic_transform": None,
            }

            bc = boundaries[family]
            bc["slicing"] = np.unique(
                np.concatenate(
                    [connectivity[elt_type].ravel() for elt_type in connectivity]
                )
            )
            bc["type"] = "boundary"
            bc["periodic_transform"] = np.zeros((16,), dtype=np.float64)
            self._famm.append(fff)

    def __peek_version(self, filename):
        """Read the MSH format version from the header without parsing the file."""
        with open(filename) as f:
            while True:
                line = f.readline()
                if line == "":
                    api.error_stop(f"no $MeshFormat in {filename!r}")
                if line.rstrip() == "$MeshFormat":
                    break
            version = f.readline().split()
        self.version = int(float(version[0]))
        assert int(version[1]) == 0, "only ASCII version is supported"
        assert int(version[2]) == 8, "size of float must be 8 (64bits)"
        return self.version

    def __build_connectivity(self, elts):
        """Build volume/boundary connectivity from a legacy ``elts`` list.

        elt = [index, gmsh_type, ntag, phys_group, entity_tag, *node_ids].
        Returns (connectivity, connect_bc, bnd, bnd2fam, maxdim, ncell), the
        connectivity arrays in 0-based numbering. Used by the v2 path (the v4
        path fills the same structures directly while streaming).
        """
        log.info("Analyze...")
        maxdim = np.amax([_elem.dim_elem[gmshtype_elem[e[1]]] for e in elts])
        mesh_elt = [etype for etype, dim in _elem.dim_elem.items() if dim == maxdim]
        bc_elt = [etype for etype, dim in _elem.dim_elem.items() if dim == maxdim - 1]

        connectivity = collections.OrderedDict({})
        bnd, bnd2fam = collections.defaultdict(list), {}

        ncell = 0
        log.info("  extract volume connectivity")
        bnd_connect = 1
        for elt in elts:
            if elt[1] not in gmshtype_elem:
                log.error("%d not in available GMSH types %r.", elt[1], gmshtype_elem)
            elt_type = gmshtype_elem[elt[1]]
            if elt_type not in mesh_elt:
                continue
            ncell += 1
            if elt_type not in connectivity:
                connectivity[elt_type] = []
            ntag = elt[2]  # expected to be 2
            if ntag != 2:
                api.error_stop(f"Expected number of tags to be 2, got {ntag}")
            connectivity[elt_type].append(elt[3 + ntag :])
            bnd_connect += 1

        log.info("  extract boundary connectivity")
        connect_bc = {}
        for elt in elts:
            if elt[1] not in gmshtype_elem:
                continue
            elt_type = gmshtype_elem[elt[1]]
            if elt_type not in bc_elt:
                continue
            ntag = elt[2]  # expected to be 2
            if ntag != 2:
                api.error_stop(f"Expected number of tags to be 2, got {ntag}")
            bnd_tag = elt[4]
            bnd_fam = elt[3]
            if bnd_tag not in connect_bc:
                connect_bc[bnd_tag] = collections.defaultdict(list)
            connect_bc[bnd_tag][elt_type].append(elt[3 + ntag :])
            bnd[bnd_tag].append(bnd_connect)
            bnd2fam[bnd_tag] = bnd_fam
            bnd_connect += 1
        log.info(f"    tags are {bnd2fam}")

        for elt_type in connectivity:
            connectivity[elt_type] = np.array(connectivity[elt_type]) - 1
        for bnd_tag in connect_bc:
            for elt_type in connect_bc[bnd_tag]:
                connect_bc[bnd_tag][elt_type] = (
                    np.array(connect_bc[bnd_tag][elt_type]) - 1
                )

        return connectivity, connect_bc, bnd, bnd2fam, maxdim, ncell

    def __read_v2_stream(self, filename):
        """Streaming reader for GMSH 2.x ASCII meshes.

        Builds the same (fam, bctype, x, y, z, elts) as the legacy v2 parser but
        reads the file line by line instead of loading it whole into memory,
        then hands ``elts`` to __build_connectivity (unchanged logic). The v2
        flat element section carries no per-type counts, so ``elts`` is still
        materialised; the point here is to drop the whole-file token list.
        """
        log.info("Running version 2.0 reader (streaming)")

        def read_section(f):
            out = []
            while True:
                line = f.readline()
                if line == "":
                    api.error_stop("unexpected end of file in a GMSH section")
                s = line.rstrip()
                if s.startswith("$End"):
                    return out
                out.append(s.split())

        fam = bctype = None
        x, y, z, elts = [], [], [], []
        with open(filename) as f:
            while True:
                line = f.readline()
                if line == "":
                    break
                tag = line.rstrip()

                if tag == "$PhysicalNames":
                    pn = read_section(f)
                    fam, bctype = {}, {}
                    for entry in pn[1:]:
                        fam[entry[1]] = entry[2][1:-1]
                        bctype[entry[2][1:-1]] = entry[0]

                elif tag == "$Nodes":
                    nnodes = int(f.readline())
                    for _ in range(nnodes):
                        c = f.readline().split()
                        x.append(float(c[1]))
                        y.append(float(c[2]))
                        z.append(float(c[3]))
                    f.readline()  # consume $EndNodes

                elif tag == "$Elements":
                    nelems = int(f.readline())
                    for _ in range(nelems):
                        elts.append([int(j) for j in f.readline().split()])
                    f.readline()  # consume $EndElements
                    break  # nothing we need after $Elements

        return fam, bctype, x, y, z, elts

    def __find_ent_v4(self, entities, j):
        """Physical group for entity ``j`` -- same logic as the legacy find_ent."""
        addd = 0
        bot = int(entities[0][0]) + int(entities[0][1])
        top = int(entities[0][2]) + bot
        for entity in entities[bot : top + 1]:
            if j == int(entity[0]) and len(entity) >= 8:
                addd = int(entity[8])
                break
        return addd

    def __read_v4_stream(self, filename):
        """Streaming reader for GMSH 4.x ASCII meshes.

        Single pass with readline() only: fills preallocated numpy arrays and
        never builds the whole-file token list nor the per-element ``elts``
        list. Reproduces the legacy v4 output exactly (same node ordering, same
        find_ent boundary labelling) at a fraction of the peak memory.
        """
        log.info("Running version 4.x reader (streaming)")

        def read_section(f):
            out = []
            while True:
                line = f.readline()
                if line == "":
                    api.error_stop("unexpected end of file in a GMSH section")
                s = line.rstrip()
                if s.startswith("$End"):
                    return out
                out.append(s.split())

        fam = bctype = None
        entities = None
        x = y = z = None
        blocks = []  # (elt_type_name, entity_tag, conn[n, K] int64 1-based)

        with open(filename) as f:
            while True:
                line = f.readline()
                if line == "":
                    break
                tag = line.rstrip()

                if tag == "$PhysicalNames":
                    pn = read_section(f)
                    fam, bctype = {}, {}
                    for family in pn[1:]:
                        fam[family[1]] = family[2][1:-1]
                        bctype[family[2][1:-1]] = family[0]
                    log.info(f"    found Physical names {fam}")

                elif tag == "$Entities":
                    entities = read_section(f)

                elif tag == "$Nodes":
                    nnodes = int(f.readline().split()[1])
                    x = np.empty(nnodes, dtype=np.float64)
                    y = np.empty(nnodes, dtype=np.float64)
                    z = np.empty(nnodes, dtype=np.float64)
                    done = 0
                    while done < nnodes:
                        cnt = int(f.readline().split()[3])
                        tags = [int(f.readline()) for _ in range(cnt)]
                        for t in tags:
                            c = f.readline().split()
                            pos = t - 1
                            x[pos] = float(c[0])
                            y[pos] = float(c[1])
                            z[pos] = float(c[2])
                        done += cnt
                    f.readline()  # consume $EndNodes

                elif tag == "$Elements":
                    nblk = int(f.readline().split()[0])
                    for _ in range(nblk):
                        h = f.readline().split()
                        ent_tag, etype_code, n = int(h[1]), int(h[2]), int(h[3])
                        name = gmshtype_elem[etype_code]
                        K = nodes_per_gmshtype[etype_code]
                        conn = np.empty((n, K), dtype=np.int64)
                        for r in range(n):
                            p = f.readline().split()
                            conn[r] = [int(v) for v in p[1 : K + 1]]
                        blocks.append((name, ent_tag, conn))
                    f.readline()  # consume $EndElements
                    break  # nothing we need after $Elements

        log.info("Analyze...")
        maxdim = max(_elem.dim_elem[name] for name, _, _ in blocks)

        connectivity = collections.OrderedDict()
        connect_bc = {}
        bnd, bnd2fam = collections.defaultdict(list), {}
        ncell = 0
        for name, ent_tag, conn in blocks:
            dim = _elem.dim_elem[name]
            if dim == maxdim:  # volume element
                ncell += conn.shape[0]
                if name in connectivity:
                    connectivity[name] = np.concatenate((connectivity[name], conn - 1))
                else:
                    connectivity[name] = conn - 1
            elif dim == maxdim - 1:  # boundary element
                if ent_tag not in connect_bc:
                    connect_bc[ent_tag] = collections.OrderedDict()
                if name in connect_bc[ent_tag]:
                    connect_bc[ent_tag][name] = np.concatenate(
                        (connect_bc[ent_tag][name], conn - 1)
                    )
                else:
                    connect_bc[ent_tag][name] = conn - 1
                bnd[ent_tag].extend([0] * conn.shape[0])
                bnd2fam[ent_tag] = self.__find_ent_v4(entities, ent_tag)
        log.info(f"    tags are {bnd2fam}")

        return (
            fam,
            bctype,
            x,
            y,
            z,
            connectivity,
            connect_bc,
            bnd,
            bnd2fam,
            maxdim,
            ncell,
        )

    def __remove_27th_point_hexa27(self, conn_bc, conn_vol, x, y, z):
        """Remove the 27th local node of each HEXA27 element.

        Delete coordinates only for nodes used exclusively as 27th nodes.

        Parameters
        ----------
        conn_vol : (nel,27) array_like
            HEXA27 connectivity
        x, y, z : (nnode,ndim) array_like
            Coordinates
        Returns
        -------
        conn : (nel,26) ndarray
            Connectivity without node #27
        x : ndarray
            Corrected coordinates
        removed_nodes : ndarray
            Deleted global node ids
        """
        assert "hexa27" in conn_vol, f"HEXA27 connectivity missing in {conn_vol.keys()}"
        conn = conn_vol["hexa27"]

        # Keep first 26 nodes of the initial connectivity
        conn26_old = conn[:, :26]

        # Nodes in 27th column
        removed_nodes = conn[:, 26]

        # Store the actual coordinates of those centre points (for celldata)
        removed_coords = np.column_stack(
            (x[removed_nodes], y[removed_nodes], z[removed_nodes])
        )

        # Keep mask
        initial_nb_nodes = len(x)
        keep = np.ones(initial_nb_nodes, dtype=bool)
        keep[removed_nodes] = False

        # New coordinates
        x = x[keep]
        y = y[keep]
        z = z[keep]

        # Renumber map old -> new (keep the connectivity dtype; -1 sentinel
        # requires a signed type)
        new_id = -np.ones(initial_nb_nodes, dtype=conn.dtype)
        new_id[keep] = np.arange(np.sum(keep))

        # Renumber connectivity
        conn_vol["hexa27"] = new_id[conn26_old]

        for bnd_tag in conn_bc:
            for elt_type in conn_bc[bnd_tag]:
                assert elt_type == "quad9", f"Unexpected BC element type {elt_type}"
                conn_bc[bnd_tag][elt_type] = new_id[conn_bc[bnd_tag][elt_type]]

        return x, y, z, removed_coords

#!/usr/bin/env python3
"""Compare face-connectivity implementations on a structured cube.

The default 100 x 100 x 100 cube has one million hexahedra.  Building its
face connectivity creates six million face candidates before duplicates are
removed, so the legacy Python-dictionary implementation can require several
GiB of process memory.

Run from an installed checkout:

    python examples/perf/cube_connectivity.py

Use a smaller cube while developing:

    python examples/perf/cube_connectivity.py --size 20
"""

import argparse
import gc
import hashlib
import json
import resource
import sys
import time

import numpy as np

import cfdtools.meshbase.simple as simplemesh


def peak_rss_bytes():
    """Return the process peak resident set size in bytes on Unix platforms."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS reports bytes, while Linux and most other Unix systems report KiB.
    return peak if sys.platform == "darwin" else peak * 1024


def mib(value):
    return value / 1024.0 ** 2


class Measurements:
    def __init__(self):
        self._reference_rss = peak_rss_bytes()
        self.phases = []

    def measure(self, name, action):
        gc.collect()
        start = time.perf_counter()
        result = action()
        elapsed = time.perf_counter() - start
        rss = peak_rss_bytes()
        self.phases.append(
            {
                "name": name,
                "seconds": elapsed,
                "peak_rss_mib": mib(rss),
                "peak_rss_delta_mib": mib(rss - self._reference_rss),
            }
        )
        return result


def expected_face_counts(nx, ny, nz):
    """Return expected (internal, boundary) faces for a Cartesian hexa grid."""
    internal = (nx - 1) * ny * nz + nx * (ny - 1) * nz + nx * ny * (nz - 1)
    boundary = 2 * (ny * nz + nx * nz + nx * ny)
    return internal, boundary


def connectivity_signature(result):
    """Return face counts and a digest of face and neighbour connectivity."""
    internalfaces, iface2cell, boundaryfaces, bface2cell = result
    digest = hashlib.sha256()
    for label, connectivity in (
        ("internal", internalfaces),
        ("boundary", boundaryfaces),
    ):
        digest.update(label.encode())
        for elemtype, elemdata in connectivity.items():
            digest.update(elemtype.encode())
            digest.update(np.ascontiguousarray(elemdata["elem2node"]).tobytes())
            digest.update(np.ascontiguousarray(elemdata["index"].to_numpy()).tobytes())
    for neighbours in (iface2cell.conn, bface2cell.conn):
        digest.update(np.ascontiguousarray(neighbours).tobytes())
    return (internalfaces.nelem, boundaryfaces.nelem), digest.hexdigest()


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--size",
        type=int,
        default=128,
        help="cells in each direction; default: %(default)s",
    )
    parser.add_argument(
        "--implementation",
        choices=("new", "legacy", "both"),
        default="both",
        help="face-connectivity implementation to benchmark; default: %(default)s",
    )
    parser.add_argument(
        "--skip-faces",
        action="store_true",
        help="measure cube allocation only",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="emit measurements as JSON after the human-readable summary",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.size < 1:
        raise SystemExit("--size must be positive")

    nx = ny = nz = args.size
    ncell = nx * ny * nz
    nnode = (nx + 1) * (ny + 1) * (nz + 1)
    ninternal, nboundary = expected_face_counts(nx, ny, nz)
    measurements = Measurements()

    print("Cube connectivity benchmark")
    print("  cells: {:,}; nodes: {:,}".format(ncell, nnode))
    print(
        "  expected faces: {:,} internal + {:,} boundary = {:,}".format(
            ninternal, nboundary, ninternal + nboundary
        )
    )

    mesh = measurements.measure(
        "build structured cube",
        lambda: simplemesh.Cube(nx, ny, nz).export_mesh(),
    )

    signatures = {}
    if not args.skip_faces:
        methods = {
            "new": mesh._cell2node.create_faces_from_elems,
            "legacy": mesh._cell2node._create_faces_from_elems_legacy,
        }
        selected = ("new", "legacy") if args.implementation == "both" else (args.implementation,)
        for implementation in selected:
            result = measurements.measure(
                "{} face connectivity".format(implementation),
                methods[implementation],
            )
            counts, signature = connectivity_signature(result)
            assert counts == (ninternal, nboundary), (
                "unexpected {} face counts: got {}, expected {}".format(
                    implementation, counts, (ninternal, nboundary)
                )
            )
            signatures[implementation] = signature
            del result
            gc.collect()

        if args.implementation == "both":
            assert signatures["new"] == signatures["legacy"], (
                "new and legacy face connectivities differ"
            )

    print("\nMeasurements (process peak RSS includes NumPy and Python allocations):")
    for phase in measurements.phases:
        print(
            "  {name:<25} {seconds:8.3f} s  peak RSS {peak_rss_mib:9.1f} MiB"
            "  (delta {peak_rss_delta_mib:+9.1f} MiB)".format(**phase)
        )

    if args.implementation == "both" and not args.skip_faces:
        elapsed = {phase["name"]: phase["seconds"] for phase in measurements.phases}
        speedup = elapsed["legacy face connectivity"] / elapsed["new face connectivity"]
        print("\nComparison:")
        print("  results: new and legacy connectivities are identical")
        print("  new implementation speed-up: {:.2f}x".format(speedup))

    if args.json:
        print("\n" + json.dumps(measurements.phases, indent=2))


if __name__ == "__main__":
    main()

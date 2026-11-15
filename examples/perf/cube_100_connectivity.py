#!/usr/bin/env python3
"""Benchmark face-connectivity construction on a structured cube.

The default 100 x 100 x 100 cube has one million hexahedra.  Building its
face connectivity creates six million face candidates before duplicates are
removed, so it can require several GiB of process memory with the current
Python-dictionary implementation.

Run from an installed checkout:

    python examples/perf/cube_100_connectivity.py

Use a smaller cube while developing:

    python examples/perf/cube_100_connectivity.py --size 20
"""

import argparse
import gc
import json
import resource
import sys
import time

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


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--size",
        type=int,
        default=128,
        help="cells in each direction; default: %(default)s",
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

    if not args.skip_faces:
        measurements.measure("build face connectivity", mesh.make_face_connectivity)
        internal = mesh._faces["internal"]["face2node"].nelem
        boundary = mesh._faces["boundary"]["face2node"].nelem
        assert (internal, boundary) == (ninternal, nboundary), (
            "unexpected face counts: got ({}, {}), expected ({}, {})".format(
                internal, boundary, ninternal, nboundary
            )
        )

    print("\nMeasurements (process peak RSS includes NumPy and Python allocations):")
    for phase in measurements.phases:
        print(
            "  {name:<25} {seconds:8.3f} s  peak RSS {peak_rss_mib:9.1f} MiB"
            "  (delta {peak_rss_delta_mib:+9.1f} MiB)".format(**phase)
        )

    if args.json:
        print("\n" + json.dumps(measurements.phases, indent=2))


if __name__ == "__main__":
    main()

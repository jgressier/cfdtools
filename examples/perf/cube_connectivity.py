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

Select implementations and explicit multiblock counts:

    python examples/perf/cube_connectivity.py \
        --implementation mono multi --blocks 1 2 5 10
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
            digest.update(np.ascontiguousarray(elemdata["index"].to_list()).tobytes())
    for neighbours in (iface2cell.conn, bface2cell.conn):
        digest.update(np.ascontiguousarray(neighbours).tobytes())
    return (internalfaces.nelem, boundaryfaces.nelem), digest.hexdigest()


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--size",
        type=int,
        default=100,
        help="cells in each direction; default: %(default)s",
    )
    parser.add_argument(
        "--implementation",
        choices=("legacy", "mono", "multi"),
        nargs="+",
        default=("mono", "multi"),
        help=(
            "one or more face-connectivity implementations to benchmark; "
            "default: mono multi"
        ),
    )
    parser.add_argument(
        "--skip-faces",
        action="store_true",
        help="measure cube allocation only",
    )
    parser.add_argument(
        "--blocks",
        type=int,
        nargs="+",
        default=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
        help=(
            "block counts used by the multi implementation; "
            "default: 1 2 3 4 5 6 7 8 9 10"
        ),
    )
    parser.add_argument(
        "--quiet-blocks",
        action="store_true",
        help="do not show candidate/open/reconnected counts for every block",
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
    if any(nblocks < 1 for nblocks in args.blocks):
        raise SystemExit("--blocks values must be positive")

    nx = ny = nz = args.size
    ncell = nx * ny * nz
    nnode = (nx + 1) * (ny + 1) * (nz + 1)
    ninternal, nboundary = expected_face_counts(nx, ny, nz)
    selected = tuple(dict.fromkeys(args.implementation))
    if "multi" in selected and any(nblocks > ncell for nblocks in args.blocks):
        raise SystemExit("--blocks values cannot exceed the number of cells")
    measurements = Measurements()

    print("Cube connectivity benchmark")
    print("  cells: {:,}; nodes: {:,}".format(ncell, nnode))
    print(
        "  expected faces: {:,} internal + {:,} boundary = {:,}".format(
            ninternal, nboundary, ninternal + nboundary
        )
    )
    sys.stdout.flush()

    mesh = measurements.measure(
        "build structured cube",
        lambda: simplemesh.Cube(nx, ny, nz).export_mesh(),
    )

    signatures = {}
    if not args.skip_faces:
        methods = {
            "mono": mesh._cell2node.create_faces_from_elems,
            "legacy": mesh._cell2node._create_faces_from_elems_legacy,
        }
        for implementation in (name for name in selected if name != "multi"):
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

        if "multi" in selected:
            for nblocks in args.blocks:
                label = "multi:{}".format(nblocks)
                result = measurements.measure(
                    "multiblock ({} {})".format(
                        nblocks, "block" if nblocks == 1 else "blocks"
                    ),
                    lambda nblocks=nblocks: (
                        mesh._cell2node.create_faces_from_elems_multiblock(
                            nblocks, verbose=not args.quiet_blocks
                        )
                    ),
                )
                counts, signature = connectivity_signature(result)
                assert counts == (ninternal, nboundary), (
                    "unexpected multiblock face counts with {} blocks: got {}, "
                    "expected {}".format(nblocks, counts, (ninternal, nboundary))
                )
                signatures[label] = signature
                del result
                gc.collect()

        reference_name, reference_signature = next(iter(signatures.items()))
        for name, signature in signatures.items():
            assert signature == reference_signature, (
                "{} connectivity differs from {}".format(name, reference_name)
            )

    print("\nMeasurements (process peak RSS includes NumPy and Python allocations):")
    for phase in measurements.phases:
        print(
            "  {name:<25} {seconds:8.3f} s  peak RSS {peak_rss_mib:9.1f} MiB"
            "  (delta {peak_rss_delta_mib:+9.1f} MiB)".format(**phase)
        )

    if {"legacy", "mono"}.issubset(selected) and not args.skip_faces:
        elapsed = {phase["name"]: phase["seconds"] for phase in measurements.phases}
        speedup = elapsed["legacy face connectivity"] / elapsed["mono face connectivity"]
        print("\nComparison:")
        print("  results: mono and legacy connectivities are identical")
        print("  mono implementation speed-up: {:.2f}x".format(speedup))

    if "multi" in selected and not args.skip_faces:
        print(
            "  multiblock results: identical for block counts {}".format(
                ", ".join(str(nblocks) for nblocks in args.blocks)
            )
        )
        print(
            "  note: RSS is a process high-water mark; run one configuration "
            "per process when comparing memory peaks strictly."
        )

    if args.json:
        print("\n" + json.dumps(measurements.phases, indent=2))


if __name__ == "__main__":
    main()

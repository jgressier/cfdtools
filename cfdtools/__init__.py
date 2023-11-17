# coding: utf8
__version__ = "0.5.2"

import argparse
import logging

parser = argparse.ArgumentParser()
parser.add_argument(
    "-log",
    "--log",
    default="info",
    help=("provide logging level. Example: --log debug', default='info'"),
)

try:
    args = parser.parse_args()
except SystemExit:
    # silently forget unknown arguments
    args = None

if args is not None:
    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(
            f"Invalid log level: {args.log}. Must be one of: {' | '.join(logging._nameToLevel.keys())}"
        )

    log = logging.getLogger(__name__)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    log.addHandler(handler)
    log.setLevel(args.log.upper())

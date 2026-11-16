# coding: utf8
import logging
import sys

__version__ = "0.6.0"

level = 'info'
if '--log' in sys.argv:
    index = sys.argv.index('--log')
    if not sys.argv[index + 1].startswith('-'):
        level = sys.argv[index + 1]

if level.upper() not in logging._nameToLevel:
    raise ValueError(f"Invalid log level: {level!r}."
                     "\n            "
                     f"Must be in: ( {' | '.join(logging._nameToLevel.keys())} )")

log = logging.getLogger(__name__)

handler = logging.StreamHandler()

#handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
handler.setFormatter(logging.Formatter("%(levelname)s %(message)s"))
log.addHandler(handler)
log.setLevel(level.upper())

import cfdtools._cli as cli
from pathlib import Path
import pytest

_datadir=Path("./tests/data")
_builddir=Path("./tests/build")

@pytest.mark.parametrize("args", [["./tests/data/Box3x3x2v2.ic3"]])
def test_info(args):
    assert cli.info(args)
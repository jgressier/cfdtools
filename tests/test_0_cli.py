import cfdtools._cli as cli
from pathlib import Path
import pytest

_datadir=Path("./tests/data")
_builddir=Path("./tests/build")

@pytest.mark.parametrize("args", [["./tests/data/Box3x3x2v2.ic3"], ["./tests/data/nrg-o1hllc-vf-v3.ic3"]])
def test_info(args):
    assert cli.info(args)

@pytest.mark.parametrize("args", [["./tests/data/Box3x3x2v2.ic3"], 
    ["./tests/data/nrg-tinycube-v2.ic3"], ["./tests/data/nrg-o1hllc-vf-v3.ic3"]])
def test_ic3brief(args):
    assert cli.ic3brief(args)

@pytest.mark.parametrize("args", [["./tests/data/Box3x3x2v2.ic3"]])
def test_ic3writev2(args):
    cli.write_ic3v2(args) # write 2 times to test safe new name
    assert cli.write_ic3v2(args)

@pytest.mark.parametrize("args", [["./tests/data/Box3x3x2v2.ic3"], ["./tests/data/nrg-o1hllc-vf-v3.ic3"]])
def test_ic3writev3(args):
    assert cli.write_ic3v3(args)

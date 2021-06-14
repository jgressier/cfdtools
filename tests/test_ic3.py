import cfdtools.ic3.readerV2 as ic3readv2
import pytest

_datadir="./tests/data/"

@pytest.mark.parametrize("filename", ["Box3x3x2v2.ic3", "Box3x3x2v3.ic3"])
def test_reader(filename):
    ic3mesh = ic3readv2.reader(_datadir+filename)
    rmesh = ic3mesh.read_data()
    assert rmesh.check()

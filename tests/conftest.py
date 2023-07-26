# coding: utf8
import pytest
from pathlib import Path

@pytest.fixture(scope='session')
def datadir():
    return Path("./tests/data")

@pytest.fixture(scope='session')
def builddir():
    p = Path("./tests/build")
    p.mkdir(exist_ok=True)
    return p

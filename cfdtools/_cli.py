import argparse
from pathlib import Path

import cfdtools.api as api
# readers
import cfdtools.ic3 as ic3 #.reader_legacy # needed to map readers
#import cfdtools.gmsh as gmsh

# To add a command line tool, just add the function pyproject.toml in section
# [tool.poetry.scripts]
# cfdinfo = 'cfdtools._cli:info'

#print(api._fileformat_map)

class cli_argparser(argparse.ArgumentParser):
    def __init__(self):
        super().__init__(self)
        self._available_readers = list(filter(lambda fname: 'reader' in api._fileformat_map[fname].keys(), api._fileformat_map))

    def addarg_filenameformat(self):
        self.add_argument('filename', help="file")
        self.add_argument('--fmt', help="input format", choices=self._available_readers)

    def parse_args(self, argv):
        self._args = super().parse_args(argv)

    def args(self):
        return self._args

    def parse_filenameformat(self):
        """parse args to get filename, automatic or specified format
        """
        ext = Path(self._args.filename).suffix
        thisfmt = list(filter(lambda n: api._fileformat_map[n]['ext']==ext, api._fileformat_map))
        if len(thisfmt)==0:
            api.io.print('error','no extension found')
            exit()
        elif len(thisfmt)>1:
            api.io.print('error','too many extensions found, must specify format')
            exit()
        self._fileformat = thisfmt[0]
        self._reader = api._fileformat_map[self._fileformat]['reader']

def info(argv=None):
    """call specific printinfo function from reader

    Args:
        argv (_type_, optional): _description_. Defaults to None.
    """
    #api.io.set_modes(api.io._available)
    parser = cli_argparser()
    parser.addarg_filenameformat()
    parser.parse_args(argv)
    parser.parse_filenameformat()
    #
    inputfile = Path(parser.args().filename)
    r = parser._reader(str(inputfile))
    r.read_data()
    r.printinfo()
    return True # needed for pytest

def write_ic3v2(argv=None):
    parser = cli_argparser()
    parser.addarg_filenameformat()
    parser.parse_args(argv)
    parser.parse_filenameformat()
    #
    inputfile = Path(parser.args().filename)
    r = parser._reader(str(inputfile))
    cfdmesh = r.read_data()
    outputfile = Path(inputfile.stem).with_suffix('.ic3')
    output = ic3.writerV2.writer(cfdmesh)
    output.write_data(str(outputfile))
    return True # needed for pytest

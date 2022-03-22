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

class cli_argparser():
    def __init__(self, **kwargs):
        self._parser = argparse.ArgumentParser(**kwargs)
        self._available_readers = list(filter(lambda fname: 'reader' in api._fileformat_map[fname].keys(), api._fileformat_map))

    def addarg_filenameformat(self):
        self._parser.add_argument('filename', help="file")
        self._parser.add_argument('--fmt', help="input format", choices=self._available_readers)

    def addarg_data(self):
        self._parser.add_argument('--remove-node-data', nargs='+', help="list of data to remove",)
        self._parser.add_argument('--remove-face-data', nargs='+', help="list of data to remove",)
        self._parser.add_argument('--remove-cell-data', nargs='+', help="list of data to remove",)

    def parse_cli_args(self, argv):
        self._args = self._parser.parse_args(argv)

    def args(self):
        return self._args

    def parse_filenameformat(self):
        """parse args to get filename, automatic or specified format
        """
        if self.args().fmt is None:
            ext = Path(self._args.filename).suffix
            thisfmt = list(filter(lambda n: api._fileformat_map[n]['ext']==ext, api._fileformat_map))
        else:
            thisfmt = [self.args().fmt]
        if len(thisfmt)==0:
            api.io.print('error','no extension found')
            exit()
        elif len(thisfmt)>1:
            api.io.print('error','too many extensions found, must specify format with --fmt')
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
    parser.parse_cli_args(argv)
    parser.parse_filenameformat()
    #
    inputfile = Path(parser.args().filename)
    r = parser._reader(str(inputfile))
    r.read_data()
    r.printinfo()
    return True # needed for pytest

def ic3brief(argv=None):
    parser = cli_argparser()
    parser.addarg_filenameformat()
    parser.parse_cli_args(argv)
    parser.parse_filenameformat()
    #
    r = ic3.binreader(parser.args().filename)
    r.read_headers()
    return True

def write_generic(argv, ext, writer):
    parser = cli_argparser()
    parser.addarg_filenameformat()
    parser.addarg_data()
    parser.parse_cli_args(argv)
    parser.parse_filenameformat()
    #
    file = api._files(parser.args().filename)
    r = parser._reader(file.filename)
    cfdmesh = r.read_data()
    #
    if parser.args().remove_cell_data:
        for var in parser.args().remove_cell_data:
            if cfdmesh.pop_celldata(var) is None:
                api.io.print('std', f'  cannot find cell data {var}')
            else:
                api.io.print('std', f'  pop cell data {var}')
    #
    file.remove_dir()
    file.change_suffix(ext)
    if file.find_safe_newfile() > 0:
        api.io.print("std","change output to safe new name "+file.filename)
    output = writer(cfdmesh)
    output.write_data(file.filename)
    return True # needed for pytest

def write_ic3v2(argv=None):
    return write_generic(argv, '.ic3', ic3.writerV2.writer)

def write_ic3v3(argv=None):
    return write_generic(argv, '.ic3', ic3.writerV3.writer)

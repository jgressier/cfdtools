import argparse
from pathlib import Path

import cfdtools.api as api
# readers
import cfdtools.ic3 as ic3 #.reader_legacy # needed to map readers
import cfdtools.gmsh as gmsh
import cfdtools.meshbase.simple as simplemesh
import cfdtools.probes.plot as probeplot
import cfdtools.probes.data as probedata

# To add a command line tool, just add the function pyproject.toml in section
# [tool.poetry.scripts]
# cfdinfo = 'cfdtools._cli:info'

#print(api._fileformat_map)

def cli_header(name):
    api.io.print('std',"CFDTOOLS - "+name)

class cli_argparser():
    def __init__(self, **kwargs):
        self._parser = argparse.ArgumentParser(**kwargs)
        self._available_readers = list(filter(lambda fname: 'reader' in api._fileformat_map[fname].keys(), api._fileformat_map))

    def add_argument(self, option, **kwargs):
        return self._parser.add_argument(option, **kwargs)

    def addarg_filenameformat(self):
        self.add_argument('filename', help="file")
        self.add_argument('--fmt', help="input format", choices=self._available_readers)
        self.add_argument('--outpath', help="output folder")

    def addarg_prefix(self):
        self.add_argument('prefix', help="prefix of files")

    def addarg_data(self):
        self.add_argument('--remove-node-data', nargs='+', help="list of data to remove",)
        self.add_argument('--remove-face-data', nargs='+', help="list of data to remove",)
        self.add_argument('--remove-cell-data', nargs='+', help="list of data to remove",)

    def parse_cli_args(self, argv):
        self._args = self._parser.parse_args(argv)

    def args(self, key=None):
        return self._args if key is None else vars(self._args)[key]

    def argsdict(self):
        return vars(self._args)

    def parse_filenameformat(self):
        """parse args to get filename, automatic or specified format
        """
        if self.args().fmt is None:
            ext = Path(self._args.filename).suffix
            thisfmt = list(filter(lambda n: api._fileformat_map[n]['ext']==ext, api._fileformat_map))
        else:
            thisfmt = [self.args().fmt]
        if len(thisfmt)==0:
            api.error_stop('no extension found')
        elif len(thisfmt)>1:
            api.error_stop('too many extensions found, must specify format with --fmt')
        self._fileformat = thisfmt[0]
        self._reader = api._fileformat_map[self._fileformat].get('reader', None)
        self._writer = api._fileformat_map[self._fileformat].get('writer', None)

def info(argv=None):
    """call specific printinfo function from reader

    Args:
        argv (_type_, optional): _description_. Defaults to None.
    """
    cli_header("cfdinfo")
    #api.io.set_modes(api.io._available)
    parser = cli_argparser()
    parser.addarg_filenameformat()
    parser.parse_cli_args(argv)
    parser.parse_filenameformat()
    #
    inputfile = Path(parser.args('filename'))
    r = parser._reader(str(inputfile))
    r.read_data()
    mesh = r.export_mesh()
    mesh.printinfo()
    return True # needed for pytest

def ic3brief(argv=None):
    cli_header("ic3brief")
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
    r.read_data()
    cfdmesh = r.export_mesh()
    #
    if parser.args().remove_cell_data:
        for var in parser.args().remove_cell_data:
            if cfdmesh.pop_celldata(var) is None:
                api.io.print('std', f'  cannot find cell data {var}')
            else:
                api.io.print('std', f'  pop cell data {var}')
    #
    if parser.args().outpath is None:
        file.remove_dir()
    else:
        file.change_dir(parser.args().outpath)
    file.change_suffix(ext)
    if file.find_safe_newfile() > 0:
        api.io.print("std","change output to safe new name "+file.filename)
    output = writer(cfdmesh)
    output.write_data(file.filename)
    return True # needed for pytest

def write_ic3v2(argv=None):
    cli_header("cfdwrite_ic3v2")
    return write_generic(argv, '.ic3', ic3.writerV2.writer)

def write_ic3v3(argv=None):
    cli_header("cfdwrite_ic3v3")
    return write_generic(argv, '.ic3', ic3.writerV3.writer)

def writecube(argv=None):
    cli_header("cfdwritecube")
    """call specific printinfo function from reader

    Args:
        argv (_type_, optional): _description_. Defaults to None.
    """
    #api.io.set_modes(api.io._available)
    parser = cli_argparser()
    parser.addarg_filenameformat()
    parser.add_argument("--nx", action="store", dest="nx", default=10, type=int, help="number of cells in x direction")
    parser.add_argument("--ny", action="store", dest="ny", default=10, type=int, help="number of cells in y direction")
    parser.add_argument("--nz", action="store", dest="nz", default=10, type=int, help="number of cells in z direction")
    parser.parse_cli_args(argv)
    parser.parse_filenameformat()
    nx, ny, nz = parser.args().nx, parser.args().ny, parser.args().nz
    #
    api.io.print('std', f"> create Cube {nx}x{ny}x{nz}")
    cube = simplemesh.Cube(nx, ny, nz)
    mesh = cube.export_mesh()
    #
    file = api._files(parser.args().filename)
    w = parser._writer(mesh)
    w.write_data(file.filename)
    return True # needed for pytest

def ic3probe_plotline(argv=None):
    cli_header("ic3probe_plotline")
    parser = cli_argparser(description="Process line probes from IC3")
    parser.addarg_prefix()
    #parser.add_argument("filenames", nargs="*", help="list of files")
    parser.add_argument("-v", action="store_true", dest="verbose", help="verbose output")
    parser.add_argument("--data", action="store", dest="datalist", default="P", help="quantity to plot")
    parser.add_argument("--axis", action="store", dest="axis", default="X", help="axis to follow")
    parser.add_argument("--map", action="store", dest="map", default="time", choices=["time", "freq"], help="type of map")
    parser.add_argument("--check", action="store_true", dest="check", help="process some checks")
    parser.add_argument("--cmap", action="store", dest="cmap", default="turbo", help="colormap")
    parser.add_argument("--cmaplevels", action="store", dest="nlevels", default=30, type=int, help="colormap number of levels")
    parser.parse_cli_args(argv)
    #parser.parse_filenameformat()
    #basename, ext = os.path.splitext(parser.args().filenames[0])
    var = parser.args('datalist')[0] #ext[1:]
    basename = parser.args('prefix')
    expected_data = [var, parser.args().axis] # axis must be the last to get right time size
    # check files and read data
    data = probedata.phydata(basename, verbose=parser.args().verbose)

    api.io.print('std', "> read data ")
    for ivar in expected_data:
        data.check_data(ivar, prefix=basename)

    # --- read all expected data ---
    api.io.print('std', "> processing " + parser.args().map + " map of " + var)
    run_plot = { "time": probeplot.plot_timemap, "freq": probeplot.plot_freqmap}

    # run
    run_plot[parser.args().map](data, **parser.argsdict())
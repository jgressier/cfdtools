import argparse
import logging
from pathlib import Path

import numpy as np

# Command-Line Interface

import cfdtools.api as api

# readers and writers - must be imported to update api format dict
import cfdtools.ic3 as ic3  # .reader_legacy # needed to map readers
import cfdtools.gmsh as gmsh  # .reader_legacy # needed to map readers
import cfdtools.cgns as cgns  # .reader_legacy # needed to map readers
import cfdtools.vtk as vtk

#
import cfdtools.meshbase.simple as simplemesh
import cfdtools.probes.plot as probeplot
import cfdtools.probes.data as probedata

log = logging.getLogger(__name__)

# To add a command line tool, just add the function pyproject.toml in section
# [project.scripts]
# cfdinfo = 'cfdtools._cli:info'

# print(api._fileformat_map)


# W/o  argument: @cli_header()
#                 --> <funcname>
# With altfname: @cli_header(altfname="<altfn>")
#                 --> <altfname>
# With a prefix: @cli_header("<pfx>")
#                 --> <pfx><funcname>
# With altfname: @cli_header("<pfx>", altfname="<altfn>")
#                 --> <pfx><altfname>
def cli_header(prefix=None, altfname=None):
    def name_decorator(func):
        def decorator(*args, **kwargs):
            fname = func.__name__ if altfname is None else altfname
            if prefix is not None:
                fname = prefix + fname
            func.__globals__['__fname__'] = fname  # need noqa: F821 for flake8 if using __fname__
            log.info(f"CFDTOOLS - {fname}")
            return func(*args, **kwargs)

        return decorator

    return name_decorator


class cli_argparser:
    def __init__(self, **kwargs):
        self._parser = argparse.ArgumentParser(**kwargs)
        self._available_readers = list(
            filter(lambda fname: 'reader' in api._fileformat_map[fname].keys(), api._fileformat_map)
        )

    def add_argument(self, option, **kwargs):
        return self._parser.add_argument(option, **kwargs)

    def addarg_filenameformat(self, format=None):
        self.add_argument('filename', help="file")
        self.add_argument('--fmt', help="input file format", choices=self._available_readers, default=format)
        self.add_argument('--outpath', help="output folder path")
        self.add_argument('--check', action="store_true", dest="check", help="process some checks")
        self.add_argument('--info', action="store_true", dest="info", help="print information")
        self.add_argument(
            '--log', default='info', help=("provide logging level. Default: info. Example: --log debug.")
        )

    def addarg_prefix(self):
        self.add_argument('prefix', help="prefix of files")

    def addarg_filelist(self):
        self.add_argument('filelist', nargs="+", help="list of files")

    def addarg_removedata(self):
        self.add_argument(
            '--remove-node-data',
            nargs='+',
            help="list of node data to remove",
        )
        self.add_argument(
            '--remove-face-data',
            nargs='+',
            help="list of face data to remove",
        )
        self.add_argument(
            '--remove-cell-data',
            nargs='+',
            help="list of cell data to remove",
        )
        self.add_argument(
            '--',
            dest='',
            help="may be used to separate data to remove from following filename\n\n",
        )

    def addarg_transform(self):
        self.add_argument('--extrude', type=int, help="total number of planes")
        self.add_argument('--scale', nargs=3, type=float, help="x, y, z scaling coefficients")

    def parse_cli_args(self, argv):
        self._args = self._parser.parse_args(argv)

    def args(self, key=None):
        return self._args if key is None else vars(self._args)[key]

    def argsdict(self):
        return vars(self._args)

    def parse_filenameformat(self):
        """parse args to get filename, automatic or specified format"""
        if self.args().fmt is None:
            ext = Path(self._args.filename).suffix
            thisfmt = list(filter(lambda n: api._fileformat_map[n]['ext'] == ext, api._fileformat_map))
        else:
            thisfmt = [self.args().fmt]
        if len(thisfmt) == 0:
            api.error_stop("no extension found")
        elif len(thisfmt) > 1:
            api.error_stop("too many extensions found\n" "must specify format with --fmt")
        self._fileformat = thisfmt[0]
        self._reader = api._fileformat_map[self._fileformat].get('reader', None)
        self._writer = api._fileformat_map[self._fileformat].get('writer', None)


@cli_header()
def cfdinfo(argv=None):
    """ fully reads all supported formats, 
    converts to an internal mesh and data structure 
    and prints a sum up of available information.
    """
    parser = cli_argparser(prog=__fname__)  # noqa: F821
    parser.addarg_filenameformat()
    parser.parse_cli_args(argv)
    parser.parse_filenameformat()
    #
    inputfile = Path(parser.args('filename'))
    r = parser._reader(str(inputfile), cIntegrity=parser.args('check'))
    r.read_data()
    mesh = r.export_mesh()
    mesh.printinfo()
    return True  # needed for pytest


@cli_header()
def ic3brief(argv=None):
    """reads IC3 related format (v2, v3 and soon v4)
    and print information on headers (mesh and variables)
    without reading data itself (faster)
    """
    parser = cli_argparser(prog=__fname__)  # noqa: F821
    parser.addarg_filenameformat(format='IC3')
    parser.parse_cli_args(argv)
    parser.parse_filenameformat()
    #
    r = ic3.binreader(parser.args().filename)
    r.read_headers()
    return True  # needed for pytest


@cli_header()
def vtkbrief(argv=None):
    """reads all pyvista-available formats and prints information on mesh and data
    """
    parser = cli_argparser()
    parser.addarg_filenameformat(format="VTK")
    parser.parse_cli_args(argv)
    parser.parse_filenameformat()
    #
    r = vtk.vtkMesh()
    r.read(parser.args().filename)
    r.brief()
    return True  # needed for pytest


@cli_header()
def vtkpack(argv=None):
    """reads a list of VTK-like files and 
    packs it to an cfdtools hdf5 format
    (mesh consistency will be checked and data reordered if necessary)
    """
    parser = cli_argparser()
    parser.addarg_filelist()
    parser.parse_cli_args(argv)
    #
    log.info(f"> number of files: {len(parser.args().filelist)}")
    vtklist = vtk.vtkList(parser.args().filelist, verbose=True)
    if vtklist.allexist():
        log.info("  all files exist")
    else:
        api.error_stop("some files are missing")
    vtklist.read()
    outfilename = vtklist.dumphdf("dumped.h5")
    log.info(f"> mesh and data dumped to {outfilename}")
    return outfilename  # needed for pytest


def write_generic(argv, ext, writer, fname=None):
    parser = cli_argparser(prog=fname)
    parser.addarg_filenameformat()
    parser.addarg_removedata()
    parser.addarg_transform()
    parser.parse_cli_args(argv)
    parser.parse_filenameformat()
    #
    file = api._files(parser.args().filename)
    log.info(f"> read mesh file {file.filename}")
    timer = api.Timer()
    timer.start()
    r = parser._reader(file.filename)
    r.read_data()
    ncell = r.ncell
    timer.stop(nelem=ncell)
    log.info("> export mesh ")
    cfdmesh = r.export_mesh()
    #
    if parser.args().remove_cell_data:
        for var in parser.args().remove_cell_data:
            if cfdmesh.pop_celldata(var) is None:
                log.info(f"  cannot find cell data {var}")
            else:
                log.info(f"  pop cell data {var}")
    #
    if parser.args().outpath is None:
        file.remove_dir()
    else:
        file.change_dir(parser.args().outpath)
    file.change_suffix(ext)
    if file.find_safe_newfile() > 0:
        log.info("change output to safe new name " + file.filename)
    if parser.args().extrude:
        timer.start()
        nz = parser.args().extrude
        log.info(f"> extrusion along nz={nz} cells, {nz*ncell} total cells")
        cfdmesh = cfdmesh.export_extruded(extrude=np.linspace(0.0, 1.0, nz + 1, endpoint=True))
        timer.stop(nelem=nz * ncell)
    if parser.args().scale:
        cfdmesh.scale(parser.args().scale)
    if parser.args().info:
        cfdmesh.printinfo()
    output = writer(cfdmesh)
    output.write_data(file.filename)
    log.info(f"file {file.filename} written")
    return file.filename  # filename needed for pytest (for eventual rm)


@cli_header()
def cfdwrite_ic3v2(argv=None):
    """reads all available formats and transform it to IC3.v2 format
    if available
    """
    return write_generic(argv, '.ic3', ic3.writerV2.writer, fname=__fname__)  # noqa: F821


@cli_header()
def cfdwrite_ic3v3(argv=None):
    """reads all available formats and transform it to IC3.v3 format
    if available
    `cfdwrite_ic3` is a shortname for last current IC3 writer, namely `cfdwrite_ic3v3`.
    """
    return write_generic(argv, '.ic3', ic3.writerV3.writer, fname=__fname__)  # noqa: F821


@cli_header()
def cfdwrite_vtk(argv=None):
    """reads all available formats and transform it to VTU format
    if available (IC3.v3 to VTK not yet available)
    """
    return write_generic(argv, '.vtu', vtk.vtkMesh, fname=__fname__)  # noqa: F821


@cli_header()
def cfdwritecube(argv=None):
    """creates a cartesian mesh, converts it to an unstrutured hexa mesh
    and save it to specific format

        Options:

            --nx --ny --nz : mesh sizes (default: 10)
            --fmt : file format
    """
    parser = cli_argparser(prog=__fname__)  # noqa: F821
    parser.addarg_filenameformat()
    parser.add_argument(
        "--nx",
        action="store",
        dest="nx",
        default=10,
        type=int,
        help="number of cells in x direction",
    )
    parser.add_argument(
        "--ny",
        action="store",
        dest="ny",
        default=10,
        type=int,
        help="number of cells in y direction",
    )
    parser.add_argument(
        "--nz",
        action="store",
        dest="nz",
        default=10,
        type=int,
        help="number of cells in z direction",
    )
    parser.parse_cli_args(argv)
    parser.parse_filenameformat()
    nx, ny, nz = parser.args().nx, parser.args().ny, parser.args().nz
    #
    log.info(f"> create Cube {nx}x{ny}x{nz}")
    cube = simplemesh.Cube(nx, ny, nz)
    mesh = cube.export_mesh()
    #
    file = api._files(parser.args().filename)
    w = parser._writer(mesh)
    w.write_data(file.filename)
    return file.filename  # filename needed for pytest (for eventual rm)


@cli_header()
def ic3probe_plotline(argv=None):
    """parse an IC3 set of (csv files) probe line and a time or frequency map

        Options:

            --data: variable to analyze and plot
            --axis: coordinate on x-axis of the map
            --map: 'time' or 'frequency'
            --cmap: matplotlib name of colormap
            --cmaplevels: number of levels
            --check: performs some verifications

    """
    parser = cli_argparser(description="Process line probes from IC3")
    parser.addarg_prefix()
    # parser.add_argument("filenames", nargs="*", help="list of files")
    parser.add_argument(
        "-v",
        action="store_true",
        dest="verbose",
        help="verbose output",
    )
    parser.add_argument(
        "--data",
        action="store",
        dest="datalist",
        default="P",
        help="quantity to plot",
    )
    parser.add_argument(
        "--axis",
        action="store",
        dest="axis",
        default="X",
        help="axis to follow",
    )
    parser.add_argument(
        "--map",
        action="store",
        dest="map",
        default="time",
        choices=["time", "freq"],
        help="type of map",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        dest="check",
        help="process some checks",
    )
    parser.add_argument(
        "--cmap",
        action="store",
        dest="cmap",
        default="turbo",
        help="colormap",
    )
    parser.add_argument(
        "--cmaplevels",
        action="store",
        dest="nlevels",
        default=30,
        type=int,
        help="colormap number of levels",
    )
    parser.parse_cli_args(argv)
    # parser.parse_filenameformat()
    # basename, ext = os.path.splitext(parser.args().filenames[0])
    var = parser.args('datalist')  # ext[1:]
    basename = parser.args('prefix')
    # axis must be the last to get right time size
    expected_data = [
        var,
        parser.args().axis,
    ]
    # check files and read data
    data = probedata.phydata(basename, verbose=parser.args().verbose)

    log.info(f"> read data in {basename}")
    for ivar in expected_data:
        data.check_data(ivar, prefix=basename)

    # --- read all expected data ---
    log.info(f"> processing {parser.args().map} map of {var}")
    run_plot = {
        "time": probeplot.plot_timemap,
        "freq": probeplot.plot_freqmap,
    }

    # run
    run_plot[parser.args().map](data, **parser.argsdict())

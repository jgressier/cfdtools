import argparse
import cfdtools.api as api

def info(argv=None):
    #api.io.set_modes(api.io._available)
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help="file")
    # parser.add_argument('-n',        action='store',  dest='npts', type=int, default=100, help="number of points")
    # parser.add_argument('--data',    action='append', dest='datalist', choices=quantities, default=['Mach'], help="quantity to plot")
    # parser.add_argument('--auto',    action='store_true', default=True, help="choose line according to bounding box") # not necessary, only for compatibility
    # parser.add_argument('--verbose', action='store_true', help="print messages") # not necessary, only for compatibility
    args = parser.parse_args(argv)
    print(args)
    #
    import cfdtools.ic3.readRestartIC3 as ic3r
    #
    r = ic3r.ReaderRestartIC3(args.filename, False)
    #xyz, co, bocos, nodesvar, cellsvar, extras = r.read_data()
    r.read_data()
    r.printinfo()
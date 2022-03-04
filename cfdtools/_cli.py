import argparse
from pathlib import Path

import cfdtools.api as api
# readers
import cfdtools.ic3 as ic3
#import cfdtools.gmsh as gmsh

#print(api._fileformat_map)

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
    ext = Path(args.filename).suffix
    names= list(filter(lambda n: api._fileformat_map[n]['ext']==ext, api._fileformat_map))
    print('names',names)
    if len(names)==0:
        api.io.print('error','no extension found')
        exit()
    elif len(names)>1:
        api.io.print('error','no extension found')
        exit()
    else: # only one name
        r = api._fileformat_map[names[0]]['reader'](args.filename)
        r.read_data()
        r.printinfo()
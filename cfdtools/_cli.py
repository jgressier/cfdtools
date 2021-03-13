import argparse

def info(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', nargs='*', help="file")
    # parser.add_argument('-n',        action='store',  dest='npts', type=int, default=100, help="number of points")
    # parser.add_argument('--data',    action='append', dest='datalist', choices=quantities, default=['Mach'], help="quantity to plot")
    # parser.add_argument('--auto',    action='store_true', default=True, help="choose line according to bounding box") # not necessary, only for compatibility
    # parser.add_argument('--verbose', action='store_true', help="print messages") # not necessary, only for compatibility
    args = parser.parse_args()
    print(args)
    #
    import ic3.readRestartIC3
    #
    ,
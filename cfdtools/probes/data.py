import os
from shutil import Error
import numpy as np
import cfdtools.api as api

# options definition

varname_syn = {  # name for line_probe
    "X": [ "x" ],
    "Y": [ "y" ],
    "Z": [ "z" ],
    "P": [ "p", "ps", "Ps", "PS" ],
    "U-X": [ "UX", "Ux", "ux", "U_X", "U_x" ],
    "X": [ "X", "x" ]
}

def minavgmax(d):
    return (f(d) for f in [np.min, np.average, np.max])

class phydata():
    dependency_vars = {
        "Asound": ["P", "RHO"],
        "S": ["P", "RHO"],
        "Mach": ["U", "Asound"],
        "U": ["U-X", "U-Y", "U-Z"],
        "INVX+": ["U-X", "Asound" ],
        "INVX-": ["U-X", "Asound" ]
    }

    def __init__(self, basename, verbose=False):
        self.basename = basename
        self.alldata = dict()
        self.verbose = verbose
        self.compute_varname = {
            "U" : self.compute_U,
            "Mach" : self.compute_Mach,
            "Asound" : self.compute_Asound,
            "INVX+" : self.compute_invxp,
            "INVX-" : self.compute_invxm,
            "S" : self.compute_entropy,
        }

    # TODO: this part should be moved to physics module
    def compute_U(self):
        self.alldata['U'] = np.sqrt(self.alldata['U-X']**2 + self.alldata['U-Y']**2 + self.alldata['U-Z']**2 )

    def compute_Mach(self):
        self.alldata['Mach'] = self.alldata['U']/self.alldata['Asound']

    def compute_Asound(self):
        self.alldata['Asound'] = np.sqrt(1.4*self.alldata['P']/self.alldata['RHO'])

    def compute_invxp(self):
        self.alldata['INVX+'] = 5*self.alldata['Asound']+self.alldata['U-X']

    def compute_invxm(self):
        self.alldata['INVX-'] = 5*self.alldata['Asound']-self.alldata['U-X']

    def compute_entropy(self):
        self.alldata['S'] = 1./.4*np.log(self.alldata['P']/self.alldata['RHO']**1.4)

    def check_data(self, varname, prefix=""):
        if self.verbose: print("- request "+varname)
        success = varname in self.alldata
        if not success:
            # try to directly read data
            success = self.read_data(varname, prefix)
        if not success: # try to compute it
            if varname in self.dependency_vars:
                success = np.all([self.check_data(depvar) for depvar in self.dependency_vars[varname]])
            if success:
                if self.verbose: print("- compute "+varname)
                self.compute_varname[varname]()
            else:
                raise NameError(varname + " missing or unable to compute")
        if success:
            api.io.print('std', "- "+varname+" min:avg:max = {:.3f} : {:.3f} : {:.3f}".format(*minavgmax(self.alldata[varname])))
        return success

    def read_data(self, varname, prefix=""):
        fname = prefix + "." + varname
        if os.path.exists(fname):
            if self.verbose: print("- read "+varname+" in "+fname)
            rdata = np.genfromtxt(fname, delimiter=" ")
            if rdata.ndim == 1:  # supposed to be coordinate
                self.alldata[varname] = rdata[3:]  # extract only coordinate (remove time and it)
            elif rdata.ndim == 2:  # supposed to be data
                self.alldata[varname] = rdata[:, 3:]  # extract data  (remove time and it)
                if (
                    "time" not in self.alldata
                ):  # if time missing, get it from current data, no consistency test with other data
                    if self.verbose: print(" . define 'time'")
                    self.alldata["time"] = rdata[:, 1]
            else:
                raise Error("unexpected data size " + varname)
            return True # success
        else:
            return False # no file

# --------------
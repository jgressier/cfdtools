import logging
import os
from shutil import Error

import numpy as np

log = logging.getLogger(__name__)

# options definition

varname_syn = {  # name for line_probe
    "X": ["x", "X"],
    "Y": ["y"],
    "Z": ["z"],
    "P": ["p", "ps", "Ps", "PS"],
    "U-X": ["UX", "Ux", "ux", "U_X", "U_x"],
}


def minavgmax(d):
    return (f(d) for f in [np.min, np.average, np.max])


class phydata:
    dependency_vars = {
        "Asound": ["P", "RHO"],
        "S": ["P", "RHO"],
        "Mach": ["U", "Asound"],
        "U": ["U-X", "U-Y", "U-Z"],
        "INVX+": ["U-X", "Asound"],
        "INVX-": ["U-X", "Asound"],
    }

    def __init__(self, basename, verbose=False):
        self.basename = basename
        self.alldata = dict()
        self.verbose = verbose
        self.compute_varname = {
            "U": self.compute_U,
            "Mach": self.compute_Mach,
            "Asound": self.compute_Asound,
            "INVX+": self.compute_invxp,
            "INVX-": self.compute_invxm,
            "S": self.compute_entropy,
        }

    # TODO: this part should be moved to physics module
    def compute_U(self):
        self.alldata['U'] = np.sqrt(
            self.alldata['U-X'] ** 2 + self.alldata['U-Y'] ** 2 + self.alldata['U-Z'] ** 2
        )

    def compute_Mach(self):
        self.alldata['Mach'] = self.alldata['U'] / self.alldata['Asound']

    def compute_Asound(self):
        self.alldata['Asound'] = np.sqrt(1.4 * self.alldata['P'] / self.alldata['RHO'])

    def compute_invxp(self):
        self.alldata['INVX+'] = 5 * self.alldata['Asound'] + self.alldata['U-X']

    def compute_invxm(self):
        self.alldata['INVX-'] = 5 * self.alldata['Asound'] - self.alldata['U-X']

    def compute_entropy(self):
        self.alldata['S'] = 1.0 / 0.4 * np.log(self.alldata['P'] / self.alldata['RHO'] ** 1.4)

    def check_data(self, varname, prefix=""):
        coordinates_set = { 'X', 'Y', 'Z'}
        if self.verbose:
            log.info("- request " + varname)
        success = varname in self.alldata
        if not success:
            # try to directly read data
            success = self.read_data(varname, prefix, force_coordinate=(varname in coordinates_set))
        if not success:  # try to compute it
            if varname in self.dependency_vars:
                success = np.all([self.check_data(depvar, prefix) for depvar in self.dependency_vars[varname]])
            if success:
                if self.verbose:
                    log.info("- compute " + varname)
                self.compute_varname[varname]()
            else:
                raise NameError(varname + " missing or unable to compute")
        if success:
            log.info(
                "- "
                + varname
                + " min:avg:max = {:.3f} : {:.3f} : {:.3f}".format(*minavgmax(self.alldata[varname])),
            )
        return success

    def read_data(self, varname, prefix="", force_coordinate=False):
        """Get data from a CSV file

        The CSV file is supposed to have at least 4 columns (it, time, size, data[size])
        `force_coordinate` is used to get only the first line
        The file is name `prefix.varname` 

        Args:
            varname (str): variable name in self.alldata[]
            prefix (str, optional): prefix of the file 'prefix.varname'
            force_coordinate (bool, optional): _description_. Defaults to False.

        Raises:
            Error: _description_

        Returns:
            Bool: success in reading file
        """
        fname = prefix + "." + varname
        if os.path.exists(fname):
            if self.verbose:
                log.info("- read " + varname + " in " + fname)
            rdata = np.genfromtxt(fname, delimiter=" ")
            if rdata.ndim == 1:  # supposed to be coordinate
                # extract only coordinate (remove time and it)
                self.alldata[varname] = rdata[3:]
            elif rdata.ndim == 2:  # supposed to be data
                # extract data  (remove time and it)
                if force_coordinate:
                    self.alldata[varname] = rdata[0, 3:]
                    assert np.allclose(self.alldata[varname], np.average(rdata[:, 3:], axis=0))
                else:
                    self.alldata[varname] = rdata[:, 3:]
                    if "time" not in self.alldata:
                        # if time missing, get it from current data, no consistency test with other data
                        if self.verbose:
                            log.info(" . define 'time'")
                        self.alldata["time"] = rdata[:, 1]
            else:
                raise Error("unexpected data size " + varname)
            return True  # success
        else:
            return False  # no file


# --------------

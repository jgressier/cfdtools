# Import modules
import logging

# import numpy as np
import cfdtools.api as api

from cfdtools.ic3._ic3 import (
    ic3_restart_codes,
    properties_ugpcode,
    type2nbytes,
    BinaryWrite,
    restartSectionHeader,
)
from cfdtools.ic3 import writerV2

log = logging.getLogger(__name__)

###################################################################################################


@api.fileformat_writer('IC3', '.ic3')
class writer(writerV2.writer):
    '''Implementation of the writer to write ic3 restart files'''

    __version__ = "3"

    def __init__(self, mesh, endian='native'):
        """
        Initialization of a ic3 restart file writer.
        """
        super().__init__(mesh, endian)

    def check(self):
        return True

    def __WriteRestartHeader(self):
        """
        Method writing the header of a restart file.
        It is composed of two integers, the "magic number" used as a flag for endianness
        and the CharlesX version number.
        input   : handle on an open restart file, [type file identifier]
        """
        # Write the two integers
        BinaryWrite(
            self.fid,
            self.endian,
            "ii",
            [ic3_restart_codes["UGP_IO_MAGIC_NUMBER"], 3],
        )

    def __WriteRestartConnectivity_check(self):
        # nothing to do in V3
        return

    def __WriteRestartVar(self):
        """
        Method to write all the variables into a restart file.
        Scalars, vectors and tensors all together.
        """
        # Start with the node based variables
        #
        nno = self.params["no_count"]
        #
        for ndname, nddata in self.vars["nodes"].items():
            # Scalar
            if nddata.size == nddata.shape[0]:
                # Header
                log.info(f"  write node scalar data {ndname}")
                header = restartSectionHeader()
                header.name = ndname
                header.id = {
                    "float64": ic3_restart_codes["UGP_IO_NO_D1"],
                    "int64": ic3_restart_codes["UGP_IO_NO_II1"],
                }[nddata.dtype.name]
                header.skip = header.hsize + nddata.itemsize * nno
                header.idata[0] = nno
                header.write(self.fid, self.endian)
                # Field
                chartype = properties_ugpcode[header.id]['structcode']
                BinaryWrite(
                    self.fid,
                    self.endian,
                    chartype * nno,
                    nddata,
                )
        #
        for ndname, nddata in self.vars["nodes"].items():
            # Vector
            if len(nddata.shape) == 2:
                log.info(f"  write node vector data {ndname}")
                # Header
                header = restartSectionHeader()
                header.name = ndname
                header.id = ic3_restart_codes["UGP_IO_NO_D3"]
                header.skip = header.hsize + type2nbytes["float64"] * nno * 3
                header.idata[0] = nno
                header.idata[1] = 3
                header.write(self.fid, self.endian)
                # Field
                BinaryWrite(
                    self.fid,
                    self.endian,
                    "d" * nno * 3,
                    nddata.ravel(order='C'),
                )
        #
        for ndname, nddata in self.vars["nodes"].items():
            # Tensor
            if len(nddata.shape) == 3:
                pass

        # Then the cell based variables
        #
        if self.vars["cells"]:  # if defined
            ndof = self.vars["cells"].ndof
            ncv = self.params["cv_count"]
            totsize = ndof * ncv
        #
        for cvname, cvdata in self.vars["cells"].items():
            # Scalar
            if cvdata.size == cvdata.shape[0]:
                # Header
                log.info(f"  write cell scalar data {cvname}")
                header = restartSectionHeader()
                header.name = cvname
                header.id = {
                    "float64": ic3_restart_codes["UGP_IO_CV_D1"],
                    "int64": ic3_restart_codes["UGP_IO_CV_II1"],
                }[cvdata.dtype.name]
                header.skip = header.hsize + type2nbytes[cvdata.dtype.name] * totsize
                header.idata[0] = ncv
                header.idata[1] = ndof
                header.write(self.fid, self.endian)
                # Field
                chartype = properties_ugpcode[header.id]['structcode']
                BinaryWrite(
                    self.fid,
                    self.endian,
                    chartype * totsize,
                    cvdata,
                )
        #
        for cvname, cvdata in self.vars["cells"].items():
            # Vector
            if len(cvdata.shape) == 2:
                # Header
                log.info(f"  write cell vector data {cvname}")
                header = restartSectionHeader()
                header.name = cvname
                header.id = ic3_restart_codes["UGP_IO_CV_D3"]
                header.skip = header.hsize + type2nbytes["float64"] * totsize * 3
                header.idata[0] = ncv
                header.idata[1] = ndof
                # header.idata[1] = 3
                header.write(self.fid, self.endian)
                # Field
                BinaryWrite(
                    self.fid,
                    self.endian,
                    "d" * totsize * 3,
                    cvdata.ravel(order='C'),
                )
        #
        for cvname, cvdata in self.vars["cells"].items():
            # Tensor
            if len(cvdata.shape) == 3:
                # Header
                log.info(f"  write cell tensor data {cvname}")
                header = restartSectionHeader()
                header.name = cvname
                header.id = ic3_restart_codes["UGP_IO_CV_D33"]
                header.skip = header.hsize + type2nbytes["float64"] * totsize * 9
                header.idata[0] = ncv
                header.idata[1] = ndof
                # header.idata[1] = 3
                header.idata[2] = 3
                header.write(self.fid, self.endian)
                # Field
                BinaryWrite(
                    self.fid,
                    self.endian,
                    "d" * totsize * 9,
                    cvdata.ravel(order='C'),
                )

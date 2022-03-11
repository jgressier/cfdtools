
# Import modules
import numpy as np
import cfdtools.api as api
import cfdtools.meshbase._mesh as _mesh
from cfdtools.ic3._ic3 import *
from cfdtools.ic3.writerV2 import writer as writer_v2

###################################################################################################

@api.fileformat_writer('IC3', '.ic3')
class writer(writer_v2):
    ''' Implementation of the writer to write ic3 restart files '''
    __version__="3"


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
        BinaryWrite(self.fid, self.endian, "ii", [ic3_restart_codes["UGP_IO_MAGIC_NUMBER"], 3])

    def __WriteRestartConnectivity_check(self):
        # nothing to do in V3
        return

    def __WriteRestartVar(self):
        """
        Method to write all the variables into a restart file.
        Scalars, vectors and tensors all together.
        """
        # Start with the node based variables
        for key, item in self.vars["nodes"].items():
            # Scalar
            if item.size == item.shape[0]:
                # Header
                api.io.print('std', '  write node scalar data '+key)
                header = restartSectionHeader()
                header.name = key
                header.id = {"float64": ic3_restart_codes["UGP_IO_NO_D1"],
                    "int64": ic3_restart_codes["UGP_IO_NO_II1"]}[item.dtype.name]
                header.skip = header.hsize + item.itemsize * self.params["no_count"]
                header.idata[0] = self.params["no_count"]
                header.write(self.fid, self.endian)
                # Field
                chartype = properties_ugpcode[header.id]['structcode']
                BinaryWrite(self.fid, self.endian, chartype*self.params["no_count"], item)
        for key, item in self.vars["nodes"].items():
            # Vector
            if len(item.shape) == 2:
                api.io.print('std', '  write node vector data '+key)
                # Header
                header = restartSectionHeader()
                header.name = key
                header.id = ic3_restart_codes["UGP_IO_NO_D3"]
                header.skip = header.hsize + type2nbytes["float64"] * self.params["no_count"] * 3
                header.idata[0] = self.params["no_count"]
                header.idata[1] = 3
                header.write(self.fid, self.endian)
                # Field
                BinaryWrite(self.fid, self.endian, "d"*self.params["no_count"]*3, item.ravel(order='C'))
        for key, item in self.vars["nodes"].items():
            # Tensor
            if len(item.shape) == 3:
                pass

        # Then the cell based variables
        for key, item in self.vars["cells"].items():
            ndof = item.ndof()
            npdata = item.data()
            ncv = self.params["cv_count"]
            totsize = ndof*ncv
            # Scalar
            if npdata.size == npdata.shape[0]:
                # Header
                api.io.print('std', '  write cell scalar data '+key)
                header = restartSectionHeader()
                header.name = key
                header.id = {"float64": ic3_restart_codes["UGP_IO_CV_D1"],
                    "int64": ic3_restart_codes["UGP_IO_CV_II1"]}[npdata.dtype.name]
                header.skip = header.hsize + type2nbytes[npdata.dtype.name] * totsize
                header.idata[0] = ncv
                header.idata[1] = ndof
                header.write(self.fid, self.endian)
                # Field
                chartype = properties_ugpcode[header.id]['structcode']
                BinaryWrite(self.fid, self.endian, chartype*totsize, npdata)
        for key, item in self.vars["cells"].items():
            ndof = item.ndof()
            npdata = item.data()
            totsize = ndof*self.params["cv_count"]
            # Vector
            if len(npdata.shape) == 2:
                # Header
                api.io.print('std', '  write cell vector data '+key)
                header = restartSectionHeader()
                header.name = key
                header.id = ic3_restart_codes["UGP_IO_CV_D3"]
                header.skip = header.hsize + type2nbytes["float64"] * totsize * 3
                header.idata[0] = ncv
                header.idata[1] = ndof
                #header.idata[1] = 3
                header.write(self.fid, self.endian)
                # Field
                BinaryWrite(self.fid, self.endian, "d"*totsize*3, npdata.ravel(order='C'))
        for key, item in self.vars["cells"].items():
            ndof = item.ndof()
            npdata = item.data()
            totsize = ndof*self.params["cv_count"]
            # Tensor
            if len(npdata.shape) == 3:
                # Header
                api.io.print('std', '  write cell tensor data '+key)
                header = restartSectionHeader()
                header.name = key
                header.id = ic3_restart_codes["UGP_IO_CV_D33"]
                header.skip = header.hsize + type2nbytes["float64"] * totsize * 9
                header.idata[0] = ncv
                header.idata[1] = ndof
                #header.idata[1] = 3
                header.idata[2] = 3
                header.write(self.fid, self.endian)
                # Field
                BinaryWrite(self.fid, self.endian, "d"*totsize*9, npdata.ravel(order='C'))


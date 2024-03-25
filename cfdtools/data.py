import logging

import numpy as np
import numpy.fft as fftm
import matplotlib.pyplot as plt

import cfdtools.api as api
import cfdtools.hdf5 as hdf5
import cfdtools.meshbase as meshbase

log = logging.getLogger(__name__)


class DataSetBase:
    _available_Xrep = ('nodal', 'cellaverage')
    _available_Trep = ('instant',)

    def __init__(self, Xrep='cellaverage', ndof=1, Trep='instant'):
        self.Trep = Trep
        self.Xrep = Xrep
        self.ndof = ndof  # must be placed after Xrep set up
        self.xsize = 0
        self._properties = dict()
        self._geoprop = dict()

    @property
    def ndof(self):
        return self._ndof

    @ndof.setter
    def ndof(self, ndof):
        assert ndof == 1 or self._Xrep == 'spectralcell'
        self._ndof = ndof

    @property
    def Xrep(self):
        return self._Xrep

    @Xrep.setter
    def Xrep(self, Xrep):
        assert Xrep in self._available_Xrep
        self._Xrep = Xrep

    @property
    def Trep(self):
        return self._Trep

    @Trep.setter
    def Trep(self, Trep):
        # self._available_Xrep is essential to self-adapt to class
        assert Trep in self._available_Trep
        self._Trep = Trep

    @property
    def geoprop(self, key=None):
        """return geometric properties (dict or value)

        Args:
            key (_type_, optional): _description_. Defaults to None.

        Returns:
            _type_: return dict or value of dict if key provided
        """
        return self._geoprop[key] if key else self._geoprop

    def set_nodes(self, x, y=None, z=None):
        """ """
        self._geoprop['x'] = x
        if y:
            self._geoprop['y'] = y
        if z:
            self._geoprop['z'] = z

    def set_mesh(self, mesh):
        self._geoprop['mesh'] = mesh


class DataSet(DataSetBase):
    _available_Xrep = ('nodal', 'cellaverage', 'spectralcell')
    _available_Trep = ('instant', 'timeevol', 'spectrum', 'spectrogram')

    def __init__(self, Xrep='cellaverage', ndof=1, Trep='instant'):
        super().__init__(Xrep, ndof, Trep)
        self._data = dict()

    def add_data(self, name, data, time=None):
        if time:
            self._properties['time'] = time
        self._data[name] = data

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        self._data = value

    def keys(self):
        return self._data.keys()

    def items(self):
        return self._data.items()

    def __getitem__(self, name):
        return self._data[name]

    def dtstats(self, imin=None, imax=None):
        assert self.Trep == 'timeevol' and 'time' in self.keys()
        t = self._data['time'][imin:imax]
        dt = t[1:] - t[:-1]
        return np.mean(dt), np.std(dt), t.shape[0]

    def plot(self, x, y, axes=plt, **kwargs):
        axes.plot(self._data[x], self._data[y], **kwargs)

    def contourf(self, x, y, v, axes=plt, **kwargs):
        axes.contourf(self._data[x], self._data[y], self._data[v].T, **kwargs)

    def dataSet_spectrum(self, datafilter=None, imin=None, imax=None):
        dtavg, dtdev, n = self.dtstats(imin, imax)  # check 'timeevol'
        if dtdev / dtavg > 1.0e-6:
            log.warning(f"dt standard deviation is significant: {dtdev / dtavg}")
        f = fftm.fftfreq(n, dtavg)
        newdataset = DataSet(self.Xrep, self.ndof, Trep='spectrum')
        newdataset.add_data('freq', f)
        datalist = datafilter if datafilter else list(self.keys())
        datalist.remove('time')
        for name in datalist:
            newdataset.add_data(name, np.abs(fftm.fft(self[name])))
        return newdataset

    def dataSet_spectrogram(self, datafilter=None, window=None):
        lowcrop = 2
        datalist = datafilter if datafilter else list(self.keys())
        datalist.remove('time')
        dtavg, dtdev, ntot = self.dtstats()  # check 'timeevol'
        if dtdev / dtavg > 1.0e-6:
            log.warning(f"dt standard deviation is significant: {100*dtdev / dtavg:.4f}%")
        newdataset = DataSet(self.Xrep, self.ndof, Trep='spectrogram')
        nwind = window if window else int(ntot / 100)
        hwin = int(nwind / 2)
        f = fftm.fftfreq(nwind, dtavg)
        newdataset.add_data('freq', f[lowcrop : hwin + 1])
        nit = int(2 * ntot / nwind) - 1
        newdataset.add_data('time', (self['time'][hwin::hwin])[:-1])
        for name in datalist:
            v = np.zeros((nit, hwin + 1 - lowcrop))
            # print(ntot, nit, nwind, hwin, newdataset['freq'].shape, f.shape)
            for i in range(nit):
                v[i, :] = np.log(
                    np.abs(fftm.rfft(self[name][hwin * i : hwin * (i + 2) + 1] * np.hanning(2 * hwin + 1)))
                )[lowcrop:]
            # print(v.shape, newdataset['freq'].shape, newdataset['time'].shape)
            newdataset.add_data(name, v)
        return newdataset


# class manufactured_DataSet(DataSet):
#
#     def __init__(self, datafunctions, Xrep='cellaverage', ndof=1, Trep='instant'):
#         super().__init__(Xrep, ndof, Trep)


class DataSetList(DataSetBase):
    _version = 1
    _available_Xrep = ('nodal', 'cellaverage', 'spectralcell')
    _available_Trep = ('instant', 'timeevol', 'pod')
    # defines the names of data items that should be written as attributes; not hdf5 dataset
    _property_names = ('time',)

    def __init__(self, ndataset, Xrep='cellaverage', ndof=1, Trep='timeevol'):
        super().__init__(Xrep, ndof, Trep)
        self._ndataset = ndataset
        self._datalist = []

    def add_datalist(self, datalist: dict, time=None):
        if time:
            datalist['time'] = time
        self._datalist.append(datalist)

    def dataSet(self, datafilter=None):
        datalist = datafilter if datafilter else list(self.keys())
        datalist.remove('time')
        dtavg, dtdev, ntot = self.dtstats()  # check 'timeevol'
        if dtdev / dtavg > 1.0e-6:
            log.warning(f"dt standard deviation is significant: {100*dtdev / dtavg:.4f}%")
        newdataset = DataSet(self.Xrep, self.ndof, Trep='spectrogram')
        return newdataset

    def _dumphdfgroup(self, hgroup: hdf5.Group, **options):
        """dump all data to an hdf group, should be used as an internal function of cfdtools, not the user

        Args:
            hgroup (hdf5.Group): _description_
        """
        for i, datadict in enumerate(self._datalist):
            datagroup = hgroup.create_group(f"i{i:06}")
            for vname, var in datadict.items():
                if vname in self._property_names:
                    datagroup.attrs[vname] = var
                else:
                    datagroup.create_dataset(vname, data=var, **options)

    def xdmf_content(self, filename, geometry_content):
        """Create the XDMF content associated to the data set.

        :param str filename: Name of the output HDF5 file.
        :param list(str) geometry_content: XDMF content for the geometry and topology of the mesh.
          The mesh is time-invariant.
        :return: The XDMF content.
        :rtype: list(str)
        """
        lines = ['<Grid Name="IC3" GridType="Collection" CollectionType="Temporal">']

        for i, datadict in enumerate(self._datalist):
            lines += ['<Grid Name="Unstructured Mesh">']
            lines += [f'<Time Value="{i}"/>']
            lines += geometry_content
            data_wo_properties = {
                vname: datadict[vname] for vname in datadict if vname not in self._property_names
            }
            for vname, var in data_wo_properties.items():
                if len(var.shape) == 1:
                    typ = "Scalar"
                elif len(var.shape) == 2:
                    typ = "Vector"
                else:
                    raise ValueError(f"Type of data cannot be deduced with shape dimension {len(var.shape)}")
                dim = " ".join(str(j) for j in var.shape)

                lines += [f'<Attribute AttributeType="{typ}" Center="Cell" Name="{vname}">']
                lines += [f'<DataItem Dimensions="{dim}" Format="HDF">']
                lines += [f"{filename}:/datalist/i{i:06}/{vname}"]
                lines += ["</DataItem>"]
                lines += ["</Attribute>"]
            lines += ["</Grid>"]
        lines += ["</Grid>"]

        return lines

    def dumphdf(self, filename, overwrite=False, **options):
        file = hdf5.h5File(filename)
        if not overwrite:
            file.find_safe_newfile()
        file.open(mode="w", datatype='datalist', version=self._version)
        # map element types from vtk type to cfdtools type
        if isinstance(self._geoprop['mesh'], meshbase.Mesh):
            celldict = dict(self._geoprop['mesh']._cell2node)
        else:
            api.error_stop("not yet implemented to dump mesh with instance {self._geoprop['mesh']}")
        file.write_unsmesh(celldict, self._geoprop['mesh'].nodescoord(ndarray=True), **options)
        #
        hdata = file._h5file.create_group("datalist")
        self._dumphdfgroup(hdata, **options)
        file.close()
        return file.filename

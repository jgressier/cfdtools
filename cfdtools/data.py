import numpy as np
import numpy.fft as fftm
import matplotlib.pyplot as plt
import cfdtools.api as api
import cfdtools.hdf5 as hdf5

# import numpy as np


class DataSetBase:

    _available_Xrep = ('nodal', 'cellaverage')
    _available_Trep = 'instant'

    def __init__(self, Xrep='cellaverage', ndof=1, Trep='instant'):
        self.Trep = Trep
        self.Xrep = Xrep
        self.ndof = ndof # must be placed after Xrep set up
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
        self._geoprop['x'] = x
        if y: self._geoprop['y'] = y
        if z: self._geoprop['z'] = z


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
    def data(self, key=None):
        return self._data[key] if key else self._data

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
        dtavg, dtdev, n = self.dtstats(imin, imax) # check 'timeevol'
        if dtdev / dtavg > 1.e-6:
            api.io.warning(f"dt standard deviation is significant: {dtdev / dtavg}")
        f = fftm.fftfreq(n, dtavg)
        newdataset = DataSet(self.Xrep, self.ndof, Trep='spectrum')
        newdataset.add_data('freq', f)
        datalist = datafilter if datafilter else list(self.keys())
        datalist.remove('time')
        for name in datalist:
            newdataset.add_data(name, np.abs(fftm.fft(self[name]-dtavg)))
        return newdataset

    def dataSet_spectrogram(self, datafilter=None, window=None):
        lowcrop = 2
        datalist = datafilter if datafilter else list(self.keys())
        datalist.remove('time')
        dtavg, dtdev, ntot = self.dtstats() # check 'timeevol'
        #print(dtavg, dtdev)
        if dtdev / dtavg > 1.e-6:
            api.io.warning(f"dt standard deviation is significant: {100*dtdev / dtavg:.4f}%")
        newdataset = DataSet(self.Xrep, self.ndof, Trep='spectrogram')
        nwind = window if window else int(ntot/100)
        hwin = int(nwind/2)
        f = fftm.fftfreq(nwind, dtavg)
        newdataset.add_data('freq', f[lowcrop:hwin+1])
        nit = int(2*ntot/nwind)-1
        newdataset.add_data('time', (self['time'][hwin::hwin])[:-1])
        for name in datalist:
            v = np.zeros((nit, hwin+1-lowcrop))
            #print(ntot, nit, nwind, hwin, newdataset['freq'].shape, f.shape)
            for i in range(nit):
                v[i,:] = np.log(np.abs(fftm.rfft(self[name][hwin*i:hwin*(i+2)+1]*np.hanning(2*hwin+1))))[lowcrop:]  
            #print(v.shape, newdataset['freq'].shape, newdataset['time'].shape)
            newdataset.add_data(name, v)
        return newdataset

# class manufactured_DataSet(DataSet):

#     def __init__(self, datafunctions, Xrep='cellaverage', ndof=1, Trep='instant'):
#         super().__init__(Xrep, ndof, Trep)


class DataSetList(DataSetBase):

    _available_Xrep = ('nodal', 'cellaverage', 'spectralcell')
    _available_Trep = ('instant', 'timeevol', 'pod')

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
        dtavg, dtdev, ntot = self.dtstats() # check 'timeevol'
        #print(dtavg, dtdev)
        if dtdev / dtavg > 1.e-6:
            api.io.warning(f"dt standard deviation is significant: {100*dtdev / dtavg:.4f}%")
        newdataset = DataSet(self.Xrep, self.ndof, Trep='spectrogram')
        return newdataset

    def _dumphdfgroup(self, hgroup: hdf5.Group, **options):
        """dump all data to an hdf group, should be used as an internal function of cfdtools, not the user

        Args:
            hgroup (hdf5.Group): _description_
        """
        n = len(self._datalist)
        for i, datadict in enumerate(self._datalist):
            datagroup = hgroup.create_group(f"i{i:06}")
            for vname, var in datadict.items():
                datagroup.create_dataset(vname, data=var, **options)

import logging

import matplotlib.pyplot as plt
import numpy as np
import numpy.fft as fftm

import cfdtools.plot as cfdplt
from cfdtools.utils.maths import minavgmax

log = logging.getLogger(__name__)


def check_axis(axisdata):
    if axisdata.ndim > 1:
        axis = np.mean(axisdata, axis=0)
        # print(axis.shape)
        err = np.mean(axisdata**2, axis=0) - axis**2
        log.info(f"min:avg:max change of several axis {err.shape}: {list(*minavgmax(err))}")
        log.warning(
            "several lines for axis: they are merged (average change is {:.2e})".format(
                np.sqrt(np.mean(np.abs(err)))
            )
        )
    else:
        axis = axisdata
    return axis


def plot_timemap(data, **kwargs):
    basename = kwargs.get('prefix')
    axis = kwargs.get('axis')
    var = kwargs.get('datalist')[0]
    cmap, nlevels = kwargs['cmap'], kwargs['nlevels']
    figname = basename + "." + var + ".time.png"
    fig = plt.figure(1, figsize=(10, 8))
    # fig.suptitle('', fontsize=12, y=0.93)
    # labels = []
    # plt.plot(x[0], qdata[0])
    # labels.append(file)
    # plt.legend(labels, loc='upper left',prop={'size':10})
    # plt.axis([0., 50., 0., 90.])
    plt.xlabel(axis, fontsize=10)
    plt.ylabel("time", fontsize=10)
    colmap = cfdplt.normalizeCmap(cmap, nlevels)
    if kwargs['verbose']:
        log.info(
            "- fields sizes are (axis, time, data) %r %r %r",
            data.alldata[axis].shape,
            data.alldata["time"].shape,
            data.alldata[var].shape,
        )
    axis = check_axis(data.alldata[axis])
    plt.contourf(axis, data.alldata["time"], data.alldata[var], levels=nlevels, cmap=colmap)
    plt.colorbar()
    # plt.minorticks_on()
    # plt.grid(which='major', linestyle='-', alpha=0.8)
    # plt.grid(which='minor', linestyle=':', alpha=0.5)
    log.info("> saving figure " + figname)
    fig.savefig(figname, bbox_inches="tight")
    fig.clf()


def plot_freqmap(data, **kwargs):
    basename = kwargs.get('prefix')
    axis = kwargs.get('axis')
    var = kwargs.get('datalist')[0]
    cmap, nlevels = kwargs['cmap'], kwargs['nlevels']
    figname = basename + "." + var + ".freq.png"
    t = data.alldata["time"]
    dtmin, dtavg, dtmax = minavgmax(t[1:] - t[:-1])
    log.info("- dt min:avg:max = {:.3f}:{:.3f}:{:.3f}".format(dtmin, dtavg, dtmax))
    if kwargs['check']:
        log.info("    t min:max = {:.3f}:{:.3f}".format(t.min(), t.max()))
        log.info("    dt < 0    = %r", np.where(t[1:] - t[:-1] < 0.0))
    # fig.suptitle('', fontsize=12, y=0.93)
    # labels = []
    # plt.plot(x[0], qdata[0])
    # labels.append(file)
    # plt.legend(labels, loc='upper left',prop={'size':10})
    # plt.axis([0., 50., 0., 90.])
    fig = plt.figure(1, figsize=(10, 8))
    plt.xlabel(axis, fontsize=10)
    plt.ylabel("frequency", fontsize=10)
    n = data.alldata[var].shape[0]
    f = fftm.fftfreq(n, dtavg)
    psdmap = np.abs(fftm.fft(data.alldata[var] - np.average(data.alldata[var]), axis=0))
    if kwargs['verbose']:
        log.info(data.alldata[var].shape, n, psdmap.shape, f.shape)
    colmap = cfdplt.normalizeCmap(cmap, nlevels)
    axis = check_axis(data.alldata[axis])
    plt.contourf(
        axis,
        f[1 : n // 200],
        np.abs(psdmap[1 : n // 200, :]),
        levels=nlevels,
        cmap=colmap,
    )
    plt.colorbar()
    # plt.minorticks_on()
    # plt.grid(which='major', linestyle='-', alpha=0.8)
    # plt.grid(which='minor', linestyle=':', alpha=0.5)
    log.info("> saving figure " + figname)
    fig.savefig(figname, bbox_inches="tight")
    fig.clf()

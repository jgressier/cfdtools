# plot.py
import matplotlib


def normalizeCmap(cmapStr, nBins):
    colormap = matplotlib.cm.get_cmap(cmapStr, nBins)
    if isinstance(colormap, matplotlib.colors.LinearSegmentedColormap):
        return matplotlib.colors.LinearSegmentedColormap(cmapStr, segmentdata=colormap._segmentdata, N=nBins)
    else:
        return matplotlib.colors.ListedColormap(matplotlib.cm.get_cmap(cmapStr, nBins).colors)

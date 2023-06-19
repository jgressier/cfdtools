import numpy as np


def minavgmax(d):
    return (f(d) for f in [np.min, np.average, np.max])


def distance(a, b):
    """distance between two arrays of positions"""
    return np.sqrt(np.sum(np.square(b - a), axis=1))

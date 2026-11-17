import numpy as np

def moving_average(array, window=10) -> np.ndarray:
    return np.mean(rolling_window(array, window), -1)

def moving_deviation(array, window=10) -> np.ndarray:
    return np.std(rolling_window(array, window), -1)

def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

def weighted_mean(values, errors):
    return ( values/errors**2 ).sum() / ((1/errors**2).sum())

def weighted_error(errors):
    return np.sqrt( 1 / (1/errors**2).sum() )
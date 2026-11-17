# -*- coding: utf-8 -*-
"""
@author: Olivier Couture
         olivier.couture@list.lu

Description: Small processing functions
"""

import numpy as np

def moving_average(array: np.ndarray, window=10) -> np.ndarray:
    """Does a moving average on array

    Args:
        array (np.ndarray): input data that will be averaged.
        window (int, optional): rolling window. Defaults to 10.

    Returns:
        np.ndarray: rolling average of array.
    """

    return np.mean(rolling_window(array, window), -1)

def moving_deviation(array: np.ndarray, window=10) -> np.ndarray:
    """Does a moving standard deviation on array.

    Args:
        array (np.ndarray): input data on which standard deviation is realized.
        window (int, optional): rolling window. Defaults to 10.

    Returns:
        np.ndarray: rolling standard deviation of array
    """

    return np.std(rolling_window(array, window), -1)

def rolling_window(a, window):
    """Backend for performing rolling operations. Creates rolled arrays. Could probably work with other types but needs the shape and stripe methods

    Args:
        a (np.ndarray): array to be rolled.
        window (int): rolling window.

    Returns:
        np.ndarray: rolled array.
    """

    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

def weighted_mean(values: np.ndarray, errors: np.ndarray) -> float:
    """Calculates the error-weighted mean.

    Args:
        values (np.ndarray): values to average
        errors (np.ndarray): errors of the values

    Raises:
        ValueError: values and errors need to be 1 dimensional
        ValueError: values and errors need to have the same size

    Returns:
        np.float: error-weighted mean
    """
    # Formula: \frac{\sum_i \frac{x_i}{\\sigma_i^2}}{\sum_i \frac{1}{\sigma_i^2}}

    if len(values.shape) != 1 or len(values.shape) !=1:
        raise ValueError("Arrays need to be 1 dimensional; values:{}, errors:{}".format(len(values.shape),len(errors.shape)))
    elif len(values) != len(errors):
        raise ValueError("Arrays need to be of same length; {} {}".format(len(values), len(errors)))

    return ( values/errors**2 ).sum() / ((1/errors**2).sum())

def weighted_error(errors):
    """Calculates the error of an error-weighted mean

    Args:
        errors (np.ndarray): errors of the values used in the error-weighted mean

    Raises:
        ValueError: errors need to be 1 dimensional

    Returns:
        np.float: error of error-weighted mean
    """
    # Formula: \sqrt{ \frac{1}{ \sum_i \frac{1}{\sigma_i^2} }

    if len(values.shape) !=1:
        raise ValueError("Array need to be 1 dimensional; {}".format(len(errors.shape)))

    return np.sqrt( 1 / (1/errors**2).sum() )
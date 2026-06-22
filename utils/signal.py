from typing import Tuple

import numpy as np
from scipy.signal import fftconvolve, windows

EPSILON = 1e-12 # Small constant to avoid division by zero

def generate_ar1(ar1, N, sigma = 1.0):
    """
    Generate an autoregressive (AR(1)) time series.

    Parameters:
    ar1 (float): The AR(1) coefficient, which determines the influence of the 
                    previous value in the series on the current value.
    N (int): The length of the time series to generate.
    sigma (float, optional): The standard deviation of the Gaussian noise added 
                                at each step. Default is 1.0.

    Returns:
    numpy.ndarray: A 1D array containing the generated AR(1) time series.
    """
    sig = np.zeros(N)
    for i in range(1, N):
        sig[i] = sig[i - 1] * ar1 + np.random.randn() * sigma
    return sig

def estimate_ar1(sig: np.ndarray) -> Tuple[float, float]:
    """
    Estimate the lag-1 autocorrelation coefficient and sigma of a time series.
    AR-1 model: y_t = ar1 * y_(t-1) + e_t, where e_t is white noise with sigma.

    Parameters
    ----------
    sig : numpy.ndarray
        Input time series as a one-dimensional array.

    Returns
    -------
    tuple[float, float]
        The estimated lag-1 autocorrelation coefficient and sigma.
    """
    _sig = np.copy(sig)
    _sig -= np.mean(sig)
    # _sig /= (np.std(_sig) + EPSILON)  # Add epsilon to avoid division by zero
    cov0 = (_sig * _sig).mean()
    cov1 = (_sig[:-1] * _sig[1:]).mean()
    ar1 = cov1 / (cov0 + EPSILON)  # Add epsilon to avoid division by zero
    sigma = np.sqrt(max(0, (1 - ar1 ** 2) * cov0))  # Ensure non-negative value
    return ar1, sigma

def smooth_weighted(sig: np.ndarray, weight: np.ndarray, axis = 0):
    """
    Smooth a multi-dimensional array along a specified axis using a weighted convolution.

    This function applies a one-dimensional convolution to smooth the input array `sig` 
    along the specified `axis` using the provided `weight` array. The convolution is performed 
    with 'same' mode to ensure the output has the same shape as the input. The function handles 
    NaN values in the input by treating them as zeros during convolution and normalizing the result 
    by the convolved weights.

    Parameters
    ----------
    sig : numpy.ndarray
        The input multi-dimensional array to be smoothed.
    weight : numpy.ndarray
        A one-dimensional array of weights used for convolution. It should be normalized (sum to 1).
    axis : int, optional
        The axis along which to perform the smoothing. Default is 0.

    Returns
    -------
    numpy.ndarray
        A new array of the same shape as `sig`, containing the smoothed values.

    Examples
    --------
    >>> import numpy as np
    >>> sig = np.random.rand(100, 10)
    >>> weight = np.array([0.25, 0.5, 0.25])
    >>> smoothed_sig = smooth_weighted(sig, weight, axis=0)
    """

    if not (0 <= axis < sig.ndim):
        raise ValueError(f"`axis` must be between 0 and {sig.ndim - 1}.")

    # Move the specified axis to the front for easier processing
    sig = np.moveaxis(sig, axis, 0)
    
    # Prepare an output array with the same shape as sig
    smoothed = np.empty_like(sig)

    # Convolve each slice along the first axis with the weight
    for i in range(sig.shape[1]):
        slice_ = sig[:, i]
        valid_mask = ~np.isnan(slice_)
        convolved = fftconvolve(np.where(valid_mask, slice_, 0), weight, mode='same')
        weight_convolved = fftconvolve(valid_mask.astype(float), weight, mode='same')
        smoothed[:, i] = np

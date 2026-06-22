import numpy as np
from scipy.signal import fftconvolve

def weight_boxcar(w: float) -> np.ndarray:
    """
    Build a fractional (non-integer) boxcar window for smoothing wavelet spectra along the scale axis.
    This function creates a symmetric window that approximates a continuous rectangular window 
    of a specified width in octaves, sampled on a discrete grid with a given number of voices per octave. 
    The window is normalized to ensure its values sum to one, making it suitable for convolution operations.

    Parameters
    ----------
    w : float
        Desired window width expressed in octaves (e.g., `0.6` corresponds to ±0.3 octave around each scale).
        Must be positive.

    Returns
    -------
    numpy.ndarray
        A one-dimensional, symmetric array of length `L`, where `L` is the smallest odd integer 
        greater than or equal to `w * nv`. The window values sum to one, making it directly usable 
        in convolution operations such as `numpy.convolve`, `scipy.signal.convolve`, or 
        `scipy.ndimage.uniform_filter`.

    Raises
    ------
    ValueError
        If `w` is not positive or if `nv` is not a positive integer.

    Examples
    --------
    >>> import numpy as np
    >>> k = weight_boxcar(0.6)
    """

    if w <= 0:
        raise ValueError("`w` must be positive.")
    
    # Fractional part to be split between the two edges
    residual = ((w - 1) % 2) / 2.0  # 0 ≤ residual < 1

    # Ensure an odd window length so the window is centred
    length = int(2 * ((w + 1) // 2) + 1)

    # Ones in the middle, fractional weights at the ends
    window = np.ones(length, dtype=float)
    window[0] = window[-1] = residual

    # Normalise to unit area
    window /= w
    return window

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

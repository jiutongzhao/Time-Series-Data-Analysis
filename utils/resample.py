import numpy as np
import scipy.interpolate
import scipy.signal

def sinc_interp(t, sig, t_interp, dt):
    """
    Perform sinc interpolation on a signal.

    Sinc interpolation is used to resample a signal at specified interpolation points.
    It is particularly useful for reconstructing a continuous signal from discrete samples.

    Parameters:
        t (array-like): The original time points of the signal.
        sig (array-like): The signal values corresponding to the time points `t`.
        t_interp (array-like): The desired time points for interpolation.
        dt (float, optional): The time step between samples in `t`. If not provided, 
                                it is estimated as the median of the differences in `t`.

    Returns:
        numpy.ndarray: The interpolated signal values at the specified `t_interp` points.
    """
    weight = np.sinc((t_interp[:, None] - t[None, :]) / dt)
    return weight @ sig

def fourier_interp(t, sig, t_interp):
    """
    Perform Fourier interpolation on a signal.

    Fourier interpolation is used to resample a signal at specified interpolation points
    using the Fourier transform. This method is particularly useful for periodic signals.

    Parameters:
        t (array-like): The original time points of the signal.
        sig (array-like): The signal values corresponding to the time points `t`.
        t_interp (array-like): The desired time points for interpolation.

    Returns:
        numpy.ndarray: The interpolated signal values at the specified `t_interp` points.
    """
    N = len(t)
    T = t[-1] - t[0] + (t[1] - t[0])
    freqs = np.fft.fftfreq(N, d=(t[1] - t[0]))
    coef = np.fft.fft(sig)
    sig_interp = np.zeros_like(t_interp, dtype=np.complex128)

    for i, ti in enumerate(t_interp):
        exponent = np.exp(2j * np.pi * freqs * (ti - t[0]))
        sig_interp[i] = np.sum(coef * exponent) / N

    return sig_interp.real

def resample(t, sig, t_resampled, method = {'nearest', 'linear', 'akima', 'cubic', 'sinc', 'fourier'}, dt = None):
    """
    Resample a signal using various interpolation methods.

    Parameters:
    -----------
    t : array-like
        The original time points of the signal.
    sig : array-like
        The signal values corresponding to the time points `t`.
    t_resampled : array-like
        The new time points at which the signal should be resampled.
    dt : float
        The time step used for sinc interpolation (only applicable if `method` is 'sinc').
    method : {'nearest', 'linear', 'akima', 'cubic', 'sinc'}, optional
        The interpolation method to use. Options include:
        - 'nearest': Nearest-neighbor interpolation.
        - 'linear': Linear interpolation.
        - 'akima': Akima spline interpolation.
        - 'cubic': Cubic spline interpolation.
        - 'sinc': Sinc interpolation.
        - 'fourier': Fourier transform-based interpolation.

    Returns:
    --------
    array-like
        The resampled signal values corresponding to `t_resampled`.

    Raises:
    -------
    ValueError
        If an unknown interpolation method is specified.

    Notes:
    ------
    - The 'sinc' method requires the `sinc_interp` function to be defined elsewhere.
    - The other methods rely on `scipy.interpolate` for interpolation.
    """

    if method is 'sinc':
        if dt is None:
            dt = np.median(np.diff(t))
        return sinc_interp(t, sig, t_resampled, dt)
    elif method is 'nearest':
        _interp = scipy.interpolate.interp1d(t, sig, kind = 'nearest', bounds_error = False)
    elif method is 'linear':
        _interp = scipy.interpolate.interp1d(t, sig, bounds_error = False)
    elif method is 'akima':
        _interp = scipy.interpolate.Akima1DInterpolator(t, sig)
    elif method is 'cubic':
        _interp = scipy.interpolate.CubicSpline(t, sig)
    else:
        raise ValueError(f"Unknown interpolation method: {method}")
    return _interp(t_resampled)
import numpy as np
import scipy.ndimage
import ssqueezepy
from scipy.signal import fftconvolve, windows
from scipy.optimize import brentq


def estimate_ar1(sig):
    _sig = np.copy(sig)
    _sig -= np.mean(sig)
    _sig /= np.std(_sig)
    cov0 = (_sig * _sig).mean()
    cov1 = (_sig[:-1] * _sig[1:]).mean()
    ar1 = cov1 / cov0
    # g = np.sqrt((1 - ar1 ** 2) * cov0)
    return ar1


def ar1_significance_level(scale, ar1=0, percentile=0.95):
    return (1 - ar1**2) / (1 + ar1**2 - 2 * ar1 * np.cos(2 * np.pi * scale)
                           )[:] * scipy.stats.chi2(2).ppf(percentile) / 2


def fractional_boxcar_kernel(w: float, nv: int = 16) -> np.ndarray:
    """
    Build a *fractional* (non-integer) boxcar kernel suitable for smoothing
    wavelet spectra along the scale axis.

    The routine is mathematically equivalent to convolving the data with a
    continuous rectangular window of width ``w`` *octaves*, then sampling the
    result on a discrete grid that has ``nv`` voices (scales) per octave.  
    If that width does not fall on an integer number of voices, the two edge
    weights are reduced proportionally so that the discrete kernel integrates
    to exactly one.

    Parameters
    ----------
    w : float
        Desired window width expressed **in octaves** (e.g., ``0.6`` gives
        ±0.3 octave around each scale).
    nv : int, optional
        Number of voices (wavelet scales) per octave.  Higher values give a
        finer scale grid; values ≥ 16 are typical.  Must be positive.

    Returns
    -------
    kernel : numpy.ndarray
        One-dimensional, symmetric array of length ``L``, where ``L`` is the
        smallest odd integer ≥ ``w * nv``.  The values sum to unity, so the
        kernel can be used directly in `numpy.convolve`, `scipy.signal.convolve`,
        or `scipy.ndimage.uniform_filter`-style operations.

    Examples
    --------
    >>> k = fractional_boxcar_kernel(0.6, nv=12)
    >>> k
    array([0.01388889, 0.13888889, 0.13888889, 0.13888889, 0.13888889,
           0.13888889, 0.13888889, 0.13888889, 0.01388889])
    >>> k.sum()
    1.0
    """
    if w <= 0:
        raise ValueError("`w` must be positive.")
    if nv <= 0:
        raise ValueError("`nv` must be a positive integer.")

    # Continuous width expressed on the discrete voice grid
    w_voice = w * nv  # e.g. 0.6 oct × 12 voices/oct = 7.2 voices

    # Fractional part to be split between the two edges
    residual = ((w_voice - 1) % 2) / 2.0  # 0 ≤ residual < 1

    # Ensure an odd kernel length so the window is centred
    length = int(2 * ((w_voice + 1) // 2) + 1)

    # Ones in the middle, fractional weights at the ends
    kernel = np.ones(length, dtype=float)
    kernel[0] = kernel[-1] = residual

    # Normalise to unit area
    kernel /= w_voice
    return kernel


def smooth_spec_gaussian(spec, sigmas):
    _spec = np.copy(spec)
    for i, s in enumerate(sigmas):
        L = int(2 * 3.0 * s) + 1
        if L < 50:
            _spec[i] = scipy.ndimage.gaussian_filter1d(spec[i],
                                                       sigma=s,
                                                       axis=0,
                                                       mode='constant',
                                                       cval=0.0)
        else:
            _kernel = windows.gaussian(L, std=s)
            _kernel /= _kernel.sum()
            _kernel = _kernel.reshape((L, ) + (1, ) * (spec[i].ndim - 1))
            _spec[i] = fftconvolve(spec[i], _kernel,
                                   mode='same')  # 注意边界条件与 ndimage 的 reflect 不同

    return _spec


def smooth_spec_boxcar(spec: np.ndarray, scales: np.ndarray, w: float, axis=0):
    dj = np.abs(np.diff(np.log2(scales))[0])
    _kernel = fractional_boxcar_kernel(w, nv=1 / dj)
    kshape = [1] * spec.ndim
    kshape[axis] = _kernel.size
    _kernel = _kernel.reshape(kshape)
    _spec = scipy.signal.convolve(
        spec * scales[:, np.newaxis, np.newaxis, np.newaxis],
        _kernel,
        mode='same')
    return _spec


def wavelet_spec(sig, scales, bandwidth=6.0, fs=1.0):
    # Compute the wavelet power spectral density

    coef, _ = ssqueezepy.cwt(sig, ('morlet', {
        'mu': bandwidth
    }),
                             scales=scales,
                             fs=fs,
                             l1_norm=False)
    coef = coef.transpose(1, 2, 0)
    spec = np.einsum('fti,ftj->ftij', coef, coef.conj())

    return spec

def transcendental_eq(x, omega0):
    # 2 x^2 (1 - e^{-omega0 x}) - 2 omega0 x - 1 + e^{-omega0 x} = 0
    return 2*x**2*(1 - np.exp(-omega0*x)) \
           - 2*omega0*x \
           - 1 \
           + np.exp(-omega0*x)

def morlet_scale_to_wavelength(mu):
    omega0 = np.geomspace(0.01, 100, 1000)
    x0 = np.zeros_like(omega0)
    for idx in range(len(omega0)):
        x0[idx] = brentq(transcendental_eq, np.sqrt(1.5), (np.sqrt(omega0[idx] ** 2 + 2) + omega0[idx]) / 1, args=(omega0[idx],), xtol = 1e-12)

def vector_angle(u, v):
    cos_theta = np.nansum(u * v, axis = -1) / (np.linalg.norm(u, axis=-1) * np.linalg.norm(v, axis=-1))
    cos_theta = np.clip(cos_theta, -1.0, 1.0)  # Ensure the value is within the valid range for arccos
    angle = np.arccos(cos_theta)
    return angle

def fractional_boxcar_kernel(w: float, nv: int = 16) -> np.ndarray:
    """
    Build a *fractional* (non-integer) boxcar kernel suitable for smoothing
    wavelet spectra along the scale axis.

    The routine is mathematically equivalent to convolving the data with a
    continuous rectangular window of width ``w`` *octaves*, then sampling the
    result on a discrete grid that has ``nv`` voices (scales) per octave.  
    If that width does not fall on an integer number of voices, the two edge
    weights are reduced proportionally so that the discrete kernel integrates
    to exactly one.

    Parameters
    ----------
    w : float
        Desired window width expressed **in octaves** (e.g., ``0.6`` gives
        ±0.3 octave around each scale).
    nv : int, optional
        Number of voices (wavelet scales) per octave.  Higher values give a
        finer scale grid; values ≥ 16 are typical.  Must be positive.

    Returns
    -------
    kernel : numpy.ndarray
        One-dimensional, symmetric array of length ``L``, where ``L`` is the
        smallest odd integer ≥ ``w * nv``.  The values sum to unity, so the
        kernel can be used directly in `numpy.convolve`, `scipy.signal.convolve`,
        or `scipy.ndimage.uniform_filter`-style operations.

    Examples
    --------
    >>> k = fractional_boxcar_kernel(0.6, nv=12)
    >>> k
    array([0.01388889, 0.13888889, 0.13888889, 0.13888889, 0.13888889,
           0.13888889, 0.13888889, 0.13888889, 0.01388889])
    >>> k.sum()
    1.0
    """
    if w <= 0:
        raise ValueError("`w` must be positive.")
    if nv <= 0:
        raise ValueError("`nv` must be a positive integer.")

    # Continuous width expressed on the discrete voice grid
    w_voice = w * nv                       # e.g. 0.6 oct × 12 voices/oct = 7.2 voices

    # Fractional part to be split between the two edges
    residual = ((w_voice - 1) % 2) / 2.0   # 0 ≤ residual < 1

    # Ensure an odd kernel length so the window is centred
    length = int(2 * ((w_voice + 1) // 2) + 1)

    # Ones in the middle, fractional weights at the ends
    kernel = np.ones(length, dtype=float)
    kernel[0] = kernel[-1] = residual

    # Normalise to unit area
    kernel /= w_voice
    return kernel

def estimate_ar1(sig):
    _sig = np.copy(sig)
    _sig -= np.mean(sig)
    _sig /= np.std(_sig)
    cov0 = (_sig * _sig).mean()
    cov1 = (_sig[:-1] * _sig[1:]).mean()
    ar1 = cov1 / cov0
    # g = np.sqrt((1 - ar1 ** 2) * cov0)
    return ar1

def smooth_spec_gaussian(spec, scales, a = 1, axis = 0):
    _spec = np.copy(spec)
    for i, s in enumerate(scales):
        _spec[i] = scipy.ndimage.gaussian_filter1d(spec[i], sigma = s * a, axis = axis, mode = 'constant', cval = 0.0)

    return _spec

def smooth_spec_box(spec: np.ndarray, scales: np.ndarray, w: float, axis = 0):    
    dj = np.abs(np.diff(np.log2(scales))[0])
    w = w / dj
    _kernel = fractional_boxcar_kernel(w)
    kshape = [1] * spec.ndim
    kshape[axis] = _kernel.size
    _kernel = _kernel.reshape(kshape)
    _spec = scipy.signal.convolve(spec * scales[:, np.newaxis, np.newaxis, np.newaxis], _kernel, mode = 'same')
    return _spec


def zero_padding(sig, N = None):
    if N is None:
        N = 2 ** np.ceil(np.log2(len(sig)) + 1)
    if len(sig) < N:
        return np.pad(sig, (0, int(N - len(sig))), mode='constant', constant_values=0)
    else:
        return sig

def ar1_significance_level(scale, ar1 = 0, percentile = 0.95):
    return (1 - ar1 ** 2) / (1 + ar1 ** 2 - 2 * ar1 * np.cos(2 * np.pi * scale))[:] * scipy.stats.chi2(2).ppf(percentile) / 2


"""
Rule of thumb for the shape of numpy.ndarray: Frequency Dimension (Nf) -> Time Dimension (Nt) -> Different Wavelets (Ns, if exist) -> Component/Channel Dimension (3, or 3 x 3 / 6 x 3 for spectral maxtrix)
"""

def fractional_boxcar_kernel(w: float) -> np.ndarray:
    if w <= 0:
        raise ValueError("`w` must be positive.")
    

    # Fractional part to be split between the two edges
    residual = ((w - 1) % 2) / 2.0   # 0 ≤ residual < 1

    # Ensure an odd kernel length so the window is centred
    length = int(2 * ((w + 1) // 2) + 1)

    # Ones in the middle, fractional weights at the ends
    kernel = np.ones(length, dtype=float)
    kernel[0] = kernel[-1] = residual

    # Normalise to unit area
    kernel /= w
    return kernel

def smooth_spec_gaussian(spec, scales, a = 1, axis = 0):
    _spec = np.copy(spec)
    for i, s in enumerate(scales):
        _spec[i] = scipy.ndimage.gaussian_filter1d(spec[i], sigma = s * a, axis = axis, mode = 'constant', cval = 0.0)

    return _spec

def smooth_spec_box(spec: np.ndarray, scales: np.ndarray, w: float, axis = 0):    
    dj = np.abs(np.diff(np.log2(scales))[0])
    print(dj)
    w = w / dj
    _kernel = fractional_boxcar_kernel(w)
    kshape = [1] * spec.ndim
    kshape[axis] = _kernel.size
    _kernel = _kernel.reshape(kshape)
    _spec = scipy.signal.convolve(spec * scales[:, np.newaxis, np.newaxis, np.newaxis], _kernel, mode = 'same')
    return _spec

def wavelet_coef_psd(time: np.ndarray, signal: np.ndarray, scales: np.ndarray, bandwidth: float = 6.0, downsample: int = 1, downsample_signal: bool = True):
    """
    Compute complex Morlet wavelet coefficients and power spectral density (PSD).
    Complex Morlet wavelet can be written as:

    Parameters:
    ----------
    time : np.ndarray
        Time vector of shape (Nt,). Can be in seconds or `np.datetime64`.
    signal : np.ndarray
        Input signal of shape (Nt,).
    scales : np.ndarray
        Wavelet scales corresponding to frequencies of shape (Nf,).
    bandwidth : float, optional
        bandwidthandwidth parameter for the Morlet wavelet (default is 6.0).
    downsample : int, optional
        Downsampling factor for the output (default is 1, meaning no downsampling).
    downsample_signal : bool, optional
        If True, downsample the input signal before computing wavelet coefficients (default is False).
        If False, the signal is not downsampled, but the output coefficients and PSD are downsampled.

    Returns:
    -------
    time : np.ndarray
        Time vector after optional downsampling.
    frequency : np.ndarray
        Frequencies corresponding to the wavelet scales.
    coef : np.ndarray
        Complex wavelet coefficients of shape (Nf, Nt), where Nf is the number of frequencies.
    psd : np.ndarray
        Power spectral density of shape (Nf, Nt).
    signal : np.ndarray
        Downsampled (moving average) signal if downsampling is applied.

    Notes:
    -----
    - The function uses `scipy.signal.cwt` with the Morlet wavelet by default.
    - `scipy.signal.cwt` is deprecated in SciPy 1.12 and will be removed in SciPy 1.15. Alternatives like PyWavelets or ssqueezepy can be used.
    - The power spectral density (PSD) is computed as the squared magnitude of the wavelet coefficients, scaled by `2 * dt`.
    - Downsampling is applied to the time, coefficients, PSD, and signal if `downsample > 1`.
    """

    if downsample_signal:
        time = time[::downsample]
        signal = signal[::downsample]
        downsample = 1

    if isinstance(time[0], np.datetime64):
        elapsed_time = np.array(time).astype('datetime64[ns]').astype('float') / 1e9
    else:
        elapsed_time = np.array(time)

    dt = elapsed_time[1] - elapsed_time[0]

    # === Option 1: scipy.signal implementation ===
    # bandwidthut, scipy.signal.cwt is deprecated in SciPy 1.12 and will be removed in SciPy 1.15. 
    # They recommend using PyWavelets instead.
    # https://docs.scipy.org/doc/scipy-1.12.0/reference/generated/scipy.signal.cwt.html
    # However, as you can see in the bottom code block, pywt.cwt is actually problematic for the Morlet wavelet.

    # widths = bandwidth * scales / (2 * np.pi)
    # coef = scipy.signal.cwt(
    #     signal,
    #     scipy.signal.morlet2,
    #     widths = widths,
    #     w = bandwidth,
    #     dtype = np.complex128
    # )
    # frequency = 1 / dt / scales

    # === Option 2: pywt implementation ===
    # Unlike scipy.signal.cwt, pywt.cwt is not L2-normalized.
    # You need to multiply the coefficients by a factor to make them L2-normalized.
    # Another two problems with pywt.cwt: The precision of the wavelet may be unenough while you can still not adjust the precision yourself.
    # Check https://github.com/PyWavelets/pywt/issues/531
    # This defect has been proposed at 2019 but still not fixed yet.

    # central_frequency = 1.0
    # wavelet = 'cmor%.1f-%.1f' % (bandwidth, central_frequency)
    # coef, frequency = pywt.cwt(signal, scales, wavelet, dt, method = 'fft')
    # coef *= np.sqrt(np.sqrt(bandwidth) * np.sqrt(2 * np.pi))  # amplitude normalization for Morlet

    # === Option 3: ssqueezepy implementation (default used here) ===
    # This one is accurate and fast. They claim that this package is the fastest implementation of the wavelet transform in Python.
    # bandwidthut, this implementation is not elegantly designed and the input parameters are not well documented.
    # Also, it requires package numba, which may raise some compatibility issues.
    coef, scales = ssqueezepy.cwt(signal - np.mean(signal), ('morlet', {'mu': bandwidth}), scales = bandwidth / (2 * np.pi) * scales.astype(np.float32), fs = 1 / dt, l1_norm = False, padtype = 'zero')
    frequency = bandwidth / (2 * np.pi) / dt / scales

    # coef, _, frequency, _ = ssqueezepy.ssq_cwt(signal, ('morlet', {'mu': bandwidth}), scales = bandwidth / (2 * np.pi) * scales.astype(np.float32), fs = 1 / dt)
    # coef = (coef.T * bandwidth * np.sqrt(scales / np.pi)).T

    psd = (np.abs(coef) ** 2) * (2 * dt)
    
    return time[::downsample], frequency, coef[:, ::downsample], psd[:, ::downsample], signal[::downsample]

def wfft_coef_psd(time: np.ndarray, signal: np.ndarray, step: int = 1, window: int = 120):
    """
    Compute short-time Fourier transform (STFT) coefficients and power spectral density (PSD) using a sliding Hanning window.

    Parameters:
    ----------
    time : np.ndarray
        Time vector of shape (Nt,). Can be in seconds or `np.datetime64`.
    signal : np.ndarray
        Input signal of shape (Nt,).
    step : int, optional
        Step size for sliding the window (default is 1).
    window : int, optional
        Window length in samples (default is 120).

    Returns:
    -------
    wtime : np.ndarray
        Center time for each window after sliding.
    freq : np.ndarray
        Frequency vector corresponding to the FFT.
    coef : np.ndarray
        Complex FFT coefficients of shape (Nf, Nt), where Nf is the number of frequencies.
    psd : np.ndarray
        Power spectral density of shape (Nf, Nt).
    wsignal : np.ndarray
        Window-averaged signal of shape (Nt,).

    Notes:
    -----
    - The Hanning window is applied to each segment, and normalization is performed based on Parseval's theorem.
    - The PSD is computed as the squared magnitude of the FFT coefficients, scaled by `2 * dt / window`.
    """

    if isinstance(time[0], np.datetime64):
        elapsed_time = np.array(time).astype('datetime64[ns]').astype('float') / 1e9
    else:
        elapsed_time = np.array(time)

    dt = elapsed_time[1] - elapsed_time[0]
    # Apply sliding window view
    wtime = np.lib.stride_tricks.sliding_window_view(elapsed_time, window)[::step][:, 0] + dt * window / 2
    freq = np.fft.fftfreq(window, dt)[:window // 2]
    freq = np.fft.rfftfreq(window, dt)
    freq = np.abs(freq)
    wsignal = np.lib.stride_tricks.sliding_window_view(signal, window)[::step]
    wsignal = wsignal - np.nanmean(wsignal, axis=-1, keepdims=True)  # Remove mean to avoid DC component
    # Apply Hanning window and normalize based on Parseval's theorem
    wsignal = wsignal * np.sqrt(8 / 3) * np.hanning(window + 1)[:-1]

    coef = np.fft.fft(wsignal, axis=-1)[:, :window // 2].T
    coef = np.fft.rfft(wsignal, axis=-1).T


    psd = (np.abs(coef) ** 2) * dt / window
    
    # Double the PSD values except for the DC component and Nyquist frequency
    if window % 2 == 0:
        psd[1:-1] *= 2

    else:
        psd[1:] *= 2

    if isinstance(time[0], np.datetime64):
        wtime = (np.array(wtime) * 1e9).astype('datetime64[ns]')

    wsignal = np.mean(wsignal, axis=-1)

    return wtime, freq, coef, psd, wsignal


def svd_wave_analysis(spec: np.ndarray, freq_window: int = 5, time_window: int = 5):
    """
    Perform SVD-based wave polarization analysis to compute planarity, ellipticity, and coherence.

    Parameters:
    ----------
    coef : np.ndarray
        Complex coefficient tensor of shape (Nf, Nt, 3), where Nf is the number of frequencies, Nt is the number of time points, and 3 represents the 3 components.
    freq_window : int, optional
        Frequency-domain smoothing window size (default is 5).
    time_window : int, optional
        Time-domain smoothing window size (default is 5).

    Returns:
    -------
    planarity : np.ndarray
        Planarity of the wave of shape (Nf, Nt), defined as `1 - sqrt(s3 / s1)`, where s1 and s3 are the largest and smallest singular values.
    ellipticity_along_k: np.ndarray
        Ellipticity along the wave vector direction of shape (Nf, Nt), defined as the ratio of the second to the first singular value.
    coherence : np.ndarray
        Coherence between the first and second principal components (along vh1 and vh2) of shape (Nf, Nt), computed from the smoothed wavefield spectrum.
    degree_of_polarization : np.ndarray
        3D Degree of polarization of shape (Nf, Nt), defined as `sqrt[3 / 2 * tr(J^2) / tr^2(J) - 1 / 2]`, computed using the wavefield spectrum.
    vh : np.ndarray
        Right singular vectors of shape (Nf, Nt, 3, 3), representing the polarization basis.

    Notes:
    -----
    - The input coefficients are smoothed in both frequency and time domains before performing SVD.
    - The coherence is computed from the wavefield spectrum in the transformed basis.
    """

    spec_63 = np.concatenate([spec.real, spec.imag], axis=-2)
    u, s, vh = np.linalg.svd(spec_63, full_matrices=False)
    planarity = 1 - np.sqrt(s[:, :, 2] / s[:, :, 0])
    ellipticity_along_k = s[:, :, 1] / s[:, :, 0]

    # ellipticity_along_k = (s[:, :, 1] - s[:, :, 2]) / (s[:, :, 0] - s[:, :, 2])

    # Rotate the coefficients to the wave frame, in which the third component is the least significant
    spec_wf = np.einsum('ftia,ftab,ftjb->ftij', vh, spec, vh.conj())

    # There are two ways to compute the degree of polarization:
    # To see the difference in theory, please refer to the paper by Taubenschuss and Santonlik (2019).
    # Equation (28) in Taubenschuss and Santonlik 2019: 
    # degree_of_polarization = np.sqrt(3 / 2 * np.abs(np.trace(np.matmul(spec, spec), axis1 = 2, axis2 = 3) / (np.trace(spec, axis1 = 2, axis2 = 3) ** 2)) - 1 / 2)

    # Equation (74) in Taubenschuss and Santonlik 2019:
    # Be careful about np.linalg.eigh, which returns the eigenvalues in ascending order
    # While np.linalg.svd returns the singular values in descending order  
    w, v = np.linalg.eigh(spec)
    degree_of_polarization = (w[:, :, 2] - w[:, :, 1]) / np.sum(w, axis = -1)

    coherence = np.abs(spec_wf[:, :, 0, 1]) / np.sqrt(np.abs(spec_wf[:, :, 0, 0] * spec_wf[:, :, 1, 1]))

    # eigenvalues, _ = np.linalg.eig(spec_wf[:, :, :2, :2])
    # eigenvalues_r, _ = np.linalg.eig(spec_wf[:, :, :2, :2].real)
    # ellipticity_along_k = np.sqrt((np.min(eigenvalues_r[:, :, :].real, axis = -1) - np.min(eigenvalues[:, :, :].real, axis = -1)) \
    #                               / (np.max(eigenvalues_r[:, :, :].real, axis = -1) - np.min(eigenvalues[:, :, :].real, axis = -1)))

    eigenvalues_r, _ = np.linalg.eigh(spec_wf[:, :, :2, :2].real) # Ascending
    eigenvalues, _ = np.linalg.eigh(spec_wf[:, :, :2, :2]) # Ascending

    ellipticity_along_k = np.sqrt((eigenvalues_r[:, :, 0] - eigenvalues[:, :, 0]) \
                                  / (eigenvalues_r[:, :, 1] - eigenvalues[:, :, 0]))

    return planarity, ellipticity_along_k, coherence, degree_of_polarization, vh

def fac_wave_analysis(spec: np.ndarray, magf: np.ndarray):
    """
    Field-aligned-coordinate (FAC) wave analysis using the power-spectral tensor.

    Parameters
    ----------
    spec : np.ndarray
        Power (cross-power) tensor of shape (Nf, Nt, 3, 3); Hermitian and ≥ 0.
    magf : np.ndarray
        Background magnetic-field vector, shape (Nt, 3).

    Returns
    -------
    compressibility : np.ndarray
        Parallel power divided by the total power, shape (Nf, Nt).
    ellipticity_along_b : np.ndarray
        (RH amplitude – LH amplitude) / (RH + LH), shape (Nf, Nt).

    Notes
    -----
    • Power along any unit vector e is P = e† S e, where S is the spectral tensor.  
    • Left-hand (LH) and right-hand (RH) unit vectors in the ⟂-plane are  
        e_lh = (e1 – i e2)/√2, e_rh = (e1 + i e2)/√2.
    """

    # ----- 1. Build the orthonormal FAC basis -----
    dir_para = magf / np.linalg.norm(magf, axis=1, keepdims=True)           # (Nt, 3)

    # Choose a reference axis least parallel to B to guarantee a non-zero cross product
    dir_ref = np.eye(3)[np.argmin(np.abs(dir_para), axis=1)]               # (Nt, 3)

    dir_perp1 = np.cross(dir_para, dir_ref)
    dir_perp1 /= np.linalg.norm(dir_perp1, axis=1, keepdims=True)          # (Nt, 3)

    dir_perp2 = np.cross(dir_para, dir_perp1)
    dir_perp2 /= np.linalg.norm(dir_perp2, axis=1, keepdims=True)          # (Nt, 3)

    # ----- 2. Define LH / RH complex unit vectors -----
    dir_lh = (dir_perp1 - 1j * dir_perp2) / np.sqrt(2)                       # (Nt, 3)
    dir_rh = (dir_perp1 + 1j * dir_perp2) / np.sqrt(2)                       # (Nt, 3)

    # ----- 3. Compute power:  e† S e  -----
    # np.einsum broadcasts (Nt,3) → (Nf,Nt,3) automatically
    P_para = np.einsum('tj,ftkj,tk->ft', dir_para.conj(), spec, dir_para)
    P_lh   = np.einsum('tj,ftkj,tk->ft',   dir_lh.conj(), spec, dir_lh)
    P_rh   = np.einsum('tj,ftkj,tk->ft',   dir_rh.conj(), spec, dir_rh)

    # ----- 4. Compressibility & ellipticity -----
    # Amplitude is √power
    amp_para = np.sqrt(P_para.real, dtype=P_para.dtype)
    amp_lh   = np.sqrt(P_lh.real,   dtype=P_lh.dtype)
    amp_rh   = np.sqrt(P_rh.real,   dtype=P_rh.dtype)

    eps = np.finfo(float).eps  # avoid 0/0
    compressibility = P_para.real / (P_para.real + P_lh.real + P_rh.real + eps)
    ellipticity_along_b = ((amp_rh - amp_lh) / (amp_rh + amp_lh + eps)).real

    return compressibility, ellipticity_along_b
# Interpretation of the signal

## Hilbert Transform [`scipy.signal.hilbert`]

The Hilbert transform is a fundamental tool for analyzing the instantaneous amplitude and phase of a signal. By constructing the analytic signal, it enables us to extract the envelope and instantaneous frequency, which are essential in the study of modulated waves and transient phenomena. This section demonstrates how to implement the Hilbert transform in Python and interpret its results in both physical and engineering contexts.

<p align = 'center'>
<img src="Figure/figure_hilbert.png" width="100%"/>
</p>



```python
omega = 2 * np.pi * 8.0
time = np.linspace(0, 1, 2 ** 7, endpoint=False)
# Modulate the Sine Wave with a offseted Hanning Window
signal = np.sin(omega * time) * (0.1 + np.hanning(time.size))
signal_ht = scipy.signal.hilbert(signal)

signal_ht.real, sighal_ht.imag, np.abs(signal_ht)
```

## Digital Filter [`scipy.interpolate`]

Digital filters are fundamental tools for shaping, extracting, or suppressing specific features in time series data. In essence, a digital filter is a mathematical algorithm that modifies the amplitude and/or phase of certain frequency components of a discrete signal. Filters can be designed to remove noise, isolate trends, block out-of-band interference, or even simulate the response of a physical system. Anti-aliasing is also a common application, where filters are used to prevent high-frequency components from distorting the signal before downsampling.

Digital filters are divided into two main types:

- **Finite Impulse Response (FIR):** The output depends only on the current and <u>**a finite number of past input samples**</u>. FIR filters are always stable and can have exactly linear phase.  A general FIR digital filter is implemented as:
  $$
  y[n] = \sum_{k=0}^{M} h[k]\, x[n-k]
  $$

  - $x[n], y[n]$: Input, Output signal
  - $h[k]$: Filter coefficients (impulse response), length $M+1$
  - $M$: Filter order

- **Infinite Impulse Response (IIR):** The output depends on both current and **<u>past input samples *and* past outputs</u>**. IIR filters can achieve sharp cutoffs with fewer coefficients, but may be unstable and generally do not preserve linear phase. A general IIR filter has both input and output recursion:
  $$
  y[n] = \sum_{k=0}^{M} b[k]\, x[n-k] - \sum_{l=1}^{N} a[l]\, y[n-l]
  $$

  - $b[k]$: Feedforward (input) coefficients
  - $a[l]$: Feedback (output) coefficients, usually $a[0] = 1$
  - $M, N$: Orders for input and output

### Example: Moving Average

The **moving average filter** is actually a simple FIR filter. For a window length $L$, the coefficients are:
$$
h[k] = \frac{1}{L},\quad k = 0, 1, \ldots, L-1
$$
So the output is:
$$
y[n] = \frac{1}{L} \sum_{k=0}^{L-1} x[n-k]
$$
That is, the output is the **average of the most recent $L$ input samples**. **Therefore, the moving average filter is an FIR filter whose coefficients are all equal.**

### Example: Low-pass FIR Filtering

Suppose we want to smooth a time series by attenuating frequencies above a certain threshold (e.g., removing noise above 50 Hz). This can be accomplished with an FIR filter designed using `scipy.signal.firwin`:

```python
import numpy as np
from scipy.signal import firwin, lfilter

fs = 200.0  # Sampling frequency (Hz)
nyq = fs / 2.0
cutoff = 50.0  # Desired cutoff frequency (Hz)
numtaps = 101  # Filter length (number of coefficients)

# Design FIR low-pass filter
fir_coeff = firwin(numtaps, cutoff / nyq)
# Apply to data
filtered_signal = lfilter(fir_coeff, 1.0, raw_signal)
```

For IIR filters (such as Butterworth, Chebyshev), the `scipy.signal.butter` function is commonly used. **Note:** Filtering can introduce edge effects—always inspect the beginning and end of the filtered signal.

<p align = 'center'>
<img src="Figure/figure_filters.png" width="100%"/>
</p>




### Frequency Response and Interpretation

The effect of a digital filter can be fully characterized by its *frequency response*, i.e., how it amplifies or suppresses each frequency. Use `scipy.signal.freqz` to plot the amplitude and phase response of your filter, and check that it matches your physical requirements (e.g., minimal ripple in the passband, sufficient attenuation in the stopband).

| Filter Type    | Function (`scipy.signal`) | Main Features                | Use Case                  |
| -------------- | ------------------------- | ---------------------------- | ------------------------- |
| FIR (window)   | `firwin`, `firwin2`       | Stable, linear phase         | Smoothing, band selection |
| Butterworth    | `butter`                  | Smooth, monotonic response   | General purpose           |
| Chebyshev I/II | `cheby1`, `cheby2`        | Sharper cutoff, ripples      | Strong suppression        |
| Elliptic       | `ellip`                   | Fastest cutoff, both ripples | Selective, small band     |
| Bessel         | `bessel`                  | Linear phase, slow rolloff   | Transient preservation    |
| Median         | `medfilt`, `medfilt1d`    | Nonlinear, preserves edges   | Spike removal             |



<p align = 'center'>
<img src="Figure/figure_filters_response.png" width="100%"/>
</p>


It is not suggested to apply a filter with an over-narrow bandwidth unless you are already confident about the central frequency and waveform of the signal. Like, if you apply a 5-Butterworth 6-10 Hz bandpass filter to a pure white noise, you are going to get a filtered signal that looks like a sine wave with a frequency around 8 Hz. Thus, you are introducing artificial waves into the signal.

<p align = 'center'>
<img src="Figure/figure_filtered_noise.png" width="100%"/>
</p>


### Practical Tips

- **Zero-phase Filtering:** Use `scipy.signal.filtfilt` for zero-phase filtering to avoid phase distortion, especially for waveform analysis.
- **Edge Effects:** Discard a small number of samples at both ends after filtering, or pad the signal before filtering to reduce transient effects.
- **Causality:** Standard filters are causal (output depends only on current and past inputs). Non-causal (zero-phase) filtering requires processing both forward and backward in time, and is not physically realizable in real-time applications.

## Interpolation

**Interpolation** is the process of estimating unknown values between discrete data points. In scientific data analysis, especially in signal processing and time series studies, interpolation plays a vital role in resampling, aligning datasets, filling gaps, and reconstructing higher-resolution signals from coarse measurements.

### Why Do We Need Interpolation?

- **Resampling:** Convert irregularly sampled data to a regular time grid for spectral analysis.
- **Filling Gaps:** Restore missing or corrupted data in a time series.
- **Temporal Alignment:** Synchronize data from different sources with differing sampling rates.
- **Upsampling/Downsampling:** Increase or decrease data resolution, e.g., for visualization or model input.

## Common Interpolation Methods

##### 1. Nearest-Neighbor and Linear Interpolation [`interp1d`]

Selects the value of the nearest known data point. Simple and fast, but produces a “blocky” or step-like signal.

```python
scipy.interpolate.interp1d(
    x: array_like,           # 1D array of independent variable (e.g., time)
    y: array_like,           # N-D array of dependent data values
    kind: str = 'linear',    # 'linear', 'nearest', 'cubic', etc.
    axis: int = -1,          
    fill_value: Union[str, float, tuple] = 'extrapolate',
    bounds_error: bool = False
)
```

##### 2. Spline Interpolation [`CubicSpline`]

Fits smooth polynomial curves (usually cubic) through the data. Produces smooth and visually appealing results, but can introduce overshoot or ringing near sharp transitions.

```python
scipy.interpolate.CubicSpline(
    x: array_like,                      # 1D array of increasing x-values
    y: array_like,                      # 1D or 2D array of values to interpolate
    bc_type: Union[str, tuple] = 'not-a-knot',  # Boundary condition type
    extrapolate: Union[bool, str] = True
)
```

##### 3. Akima Interpolation [`Akima1DInterpolator`]

Akima interpolation is a piecewise method based on fitting local polynomials between data points using adaptive slopes that depend on the trends of neighboring intervals. Unlike cubic splines, it does not enforce global smoothness but instead focuses on avoiding oscillations and overshoots near sharp transitions. This makes Akima interpolation particularly effective for datasets with non-uniform behavior or outliers, where traditional spline methods may produce unwanted ringing. It maintains a good balance between smoothness and stability and is especially useful in applications requiring visually reliable curve fitting without excessive global influence.

Assumes that the data is periodic and uniformly sampled. The signal is extended using its discrete Fourier transform (DFT), and interpolation is performed in the frequency domain by zero-padding and inverse transforming. Fourier interpolation is ideal for band-limited signals and preserves the frequency content, but it may introduce artifacts if the periodicity assumption is violated.

```python
scipy.interpolate.Akima1DInterpolator(
    x: array_like,       # 1D array of x data
    y: array_like        # 1D array of y data
)
```

##### 4. Fourier Interpolation [`resample`]

Assumes data is periodic and uses the Fourier series for reconstruction. Ideal for band-limited signals with uniform sampling, preserves frequency content.

```python
scipy.signal.resample(
    x: array_like,           # Input 1D array (signal)
    num: int,                # Number of samples in output
    t: Optional[array_like] = None,   # Original time vector (optional)
    axis: int = 0,
    window: Union[str, tuple, array_like, callable] = None,
    domain: str = 'time'     # 'time' or 'freq'
)
```

##### 5. Sinc Interpolation

**Sinc interpolation** is the theoretical ideal method for reconstructing a uniformly sampled, band-limited signal from its discrete samples. According to the Shannon sampling theorem, a continuous signal with no frequency components above the Nyquist frequency can be perfectly reconstructed from its samples using a sinc function as the interpolation kernel:
$$
x(t) = \sum_{n=-\infty}^{+\infty} x[n]\, \mathrm{sinc}\left(\frac{t - nT_s}{T_s}\right)
$$
where $T_s$ is the sampling interval, and $\mathrm{sinc}(x) = {\sin(\pi x)}/{\pi x}$.

- **Perfect for band-limited, uniformly sampled signals** (theoretical limit).
- **Preserves all frequency content up to the Nyquist frequency.**
- The interpolation kernel is infinitely wide (non-local), so true sinc interpolation is not practically achievable (requires truncation or windowing).
- In practice, *windowed sinc* or a finite sum is used.

```python
def sinc_interp(t_interp, sig, t):
    dt = t[1] - t[0]
    weight = np.sinc((t_interp[:, None] - t[None, :]) / dt)
    return weight @ sig

sinc_interp(t_interp, sig, t)
```

<p align = 'center'>
<img src="Figure/figure_interpolation.png" width="100%"/>
</p>


### Interpolation and the Frequency Domain

**Interpolation in the time domain directly impacts the signal’s frequency content:**

- **Nearest-neighbor** acts as a zero-order hold, introducing high-frequency artifacts.
- **Linear** acts as a convolution with a triangle ($\mathrm{sinc^2}$ in frequency), attenuating high-frequency components.
- **Spline/Polynomial** offers smoother spectra but can still introduce artifacts at sharp features.
- **Sinc interpolation** and **Fourier interpolation** (i.e., zero-padding in the frequency domain) yield the most faithful reconstruction for band-limited signals.

> **Tip:**
>  For spectral analysis, prefer linear, sinc, or Fourier interpolation for uniformly sampled, band-limited signals. Spline interpolation is suitable for smooth, low-noise signals, but beware of overshoot and frequency artifacts.



<div STYLE="page-break-after: always;"></div>
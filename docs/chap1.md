# Initialize Your Data

## Signals and Time Series

Many familiar phenomena can naturally be described as time series, including:

- **Meteorology:** e.g. El-Nino ENSO Index

<p align = 'center'>
<img src="Figure/figure_meiv2.png" alt="Multivariate ENSO Index (MEI)" width="60%"/>
</p>


- **Geophysics:** e.g. Seismic Waves

<p align = 'center'>
<img src="Figure/figure_seismic_waves.png" width="60%"/>
</p>

- **Solar Physics:** e.g. Sunspot Number

<p align = 'center'><img src="Figure/figure_sunspot.png" width="60%"/></p>

Each of these examples is described in the **time domain**—meaning we specify a physical quantity (such as amplitude, voltage, or magnetic field strength) explicitly **as a function of time**.

While understanding the time-domain behavior of a system is fundamental, it can sometimes be challenging to discern **underlying patterns, periodicities, or oscillatory features** directly from time-domain data. This is precisely where **spectral analysis** proves invaluable, as it allows us to examine signals in the frequency domain, clearly identifying and characterizing these hidden structures.

## Sampling

> **Nyquist-Shannon Sampling Theorem:** A band-limited continuous-time signal $x(t)$ containing no frequency components higher than $f_{max}$,  can be perfectly reconstructed from its samples if it is sampled at a rate:
$$
f_{max} \le f_s/2
$$

The frequency upper limitation $f_s/2$ is also called ***Nyquist Frequency***.

When you measure a high frequency signal with a low cadence instrument, you will not only miss the high frequency component, **<u>but also measure an erroneous signal</u>**, so called ***Aliasing***.

<p align = 'center'>
<img src="Figure/figure_aliasing.png" width="60%"/>
</p>
Such a phenomenon is essentially unrelated to the Fourier transform as its frequency range ends up to $f_s/2$ and can be directly observed by naked eye. In real life, aliasing can be visualized by recording the running car wheel (or helicopter propeller) and television (or computer screen) with your smart phone. 

It should be noted that the Nyquist frequency is only the lowest frequency to reconstruct the signal and its frequency spectrum. It doesn't mean that you can fully capture the waveform by sampling at this frequency. For example, if you sample a sine wave at a frequency of 2.01 times its frequency, you are going to get a distorted signal that does not represent the original sine wave accurately. 

This effect always happens when you (down-)sampling the signal, a common way to avoid it is to apply a low pass filter so that the high frequency component doesn't contribute to the unreal signal. In the instrumental implementation, that filter typically consists of a set of resistor, inductor, and capacity and is putted before the analog-digital converter.

- **Generate the Timestamps**

  ```python
  import numpy as np
  
  # Length of signal, N
  n = 100
  fs = 10
  T = 1.0
  
  dt = 1 / fs
  
  # Way 1: Given sampling frequency, fs
  t =	 np.arange(0, n) / fs 
  # t = 0.0, 0.1, 0.2, ..., 9.9
  
  # Way 2: Given sampling period, dt
  t =	 np.arange(0, n) / fs 
  
  # Way 3: Given Signal Duration, T
  t = np.linspace(0, T, n, endpoint = False)
  
  ```

  - Tips: Set `endpoint = False` so that the last point is not included in the time array, which ensures that the sampling frequency is equal to $f_s$.

- **Generate the Signal**

  ```python
  # Assign wave (angular) frequency (omega = 2 * np.pi * frequency)
  f0 = 3
  f1 = 6
  amp1 = 1
  amp2 = 2
  omega0 = 2 * np.pi * f0
  omega1 = 2 * np.pi * f1
  
  # Way 1: Directly generate the signal
  sig = amp1 * np.sin(omega0 * t) + amp2 * np.sin(omega1 * t)
  
  # Way 2: Generate a complex function and then take the real (cosine) or imaginary (sine) part
  sig = (amp1 * np.exp(1j * omega0 * t) + amp2 * np.exp(1j * omega1 * t)).imag
  
  # Way 3: Use a function or anonymous function to generate the signal
  def sig_func(t):
      return amp1 * np.sin(omega0 * t) + amp2 * np.sin(omega1 * t)
  # or
  sig_func = lambda t: amp1 * np.sin(omega0 * t) + amp2 * np.sin(omega1 * t)
  
  sig = sig_func(t)
  ```



## Read Signals From Data

| Interface                              | Purpose                                  | Supported Formats                         | Dependencies         | Pros                                                         | Cons                                                         |
| -------------------------------------- | ---------------------------------------- | ----------------------------------------- | -------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| **scipy.io.loadmat**                   | Load MATLAB `.mat` files (read-only)     | `.mat` v4, v6, v7–7.2                     | SciPy                | Reads structs, nested arrays, sparse matrices :contentReference[oaicite:1]{index=1} | No support for v7.3 HDF5 format, needs external libs for those files |
| **scipy.io.readsav**                   | Read IDL `.sav` files                    | IDL `.sav`                                | SciPy                | Pure‑Python, preserves structure, case-insensitive dict :contentReference[oaicite:2]{index=2} | Read-only; limited support for advanced types                |
| **scipy.io.netcdf** *(deprecated)*     | Read NetCDF3 files                       | NetCDF3                                   | SciPy                | Easy legacy I/O setup                                        | Deprecated—users are advised to use `netCDF4` or `xarray` instead |
| **numpy.load**                         | Load NumPy’s `.npy` or `.npz` files      | `.npy`, `.npz`                            | NumPy                | Fast, supports memory mapping (`mmap_mode`) :contentReference[oaicite:3]{index=3} | Only supports NumPy’s native formats                         |
| **numpy.fromfile**                     | Load raw binary or simple text files     | Raw binary or delimited text files `.dat` | NumPy                | Very fast for known-format binary/text; structured dtypes; no overhead :contentReference[oaicite:4]{index=4} | No metadata; brittle if dtype misaligned; non-portable; poor error handling |
| **pandas.read_csv/excel/json/parquet** | Load structured or semi-structured data  | CSV, Excel, JSON, Parquet, etc.           | Pandas + format libs | Rich parsing; powerful DataFrame tools                       | Memory-heavy for large datasets; extra backend dependencies  |
| **spacepy.pycdf.CDF**                  | Read NASA Common Data Format (CDF) files | `.cdf` (NASA CDF)                         | SpacePy + NASA C lib | Dict-like, lazy loading, metadata support :contentReference[oaicite:5]{index=5} | Requires NASA CDF C library; installation complexity         |



## Generate Signals Artificially

### Common Waveforms

In addition to sine waves, there are some other commonly used waveforms built-in functions that are provided by `scipy.signal` module. These functions can be used to generate various types of signals for testing, simulation, and analysis purposes. Below is a brief overview of some typical waveforms:

<p align = 'center'>
<img src="Figure/figure_typical_signals.png" width="60%"/>
</p>

#### **Chirp Waveform (`chirp`)**

Generates a swept-frequency (chirp) signal, which is often used in radar, sonar, and frequency response analysis.

```python
scipy.signal.chirp(
    t: array_like,                 # Time array
    f0: float,                     # Initial frequency at time t=0
    t1: float,                     # Final time for frequency sweep
    f1: float,                     # Final frequency at time t1
    method: str = 'linear',       # Frequency sweep type: 'linear', 'quadratic', 'logarithmic', 'hyperbolic'
    phi: float = 0                # Initial phase in degrees
)
```

------

#### **Gaussian Pulse (`gausspulse`)**

Generates a Gaussian-modulated sinusoidal pulse, commonly used in ultrasound and narrow-band radar simulations.

```python
scipy.signal.gausspulse(
    t: array_like,             # Time array
    fc: float = 1000,          # Center frequency (Hz)
    bw: float = 0.5,           # Fractional bandwidth in frequency domain
    bwr: float = -6,           # Reference level (dB) relative to peak for bandwidth measurement
    tpr: float = -60,          # Truncation level (dB) for pulse duration
    retquad: bool = False,     # If True, return quadrature (analytic) signal
    retenv: bool = False       # If True, return envelope of the signal
)
```

------

#### **Square Wave (`square`)**

Generates a square waveform, useful in digital signal simulations, PWM applications, and modulation experiments.

```python
scipy.signal.square(
    t: array_like,             # Time array
    duty: float = 0.5          # Duty cycle (fraction of period signal is high)
)
```

------

#### **Sawtooth Wave (`sawtooth`)**

Generates a sawtooth waveform, widely used in signal synthesis and electronics simulations.

```python
scipy.signal.sawtooth(
    t: array_like,             # Time array
    width: float = 1           # Width of rising ramp, from 0 to 1 (1 = triangle wave)
)
```

------

#### **Unit Impulse (`unit_impulse`)**

Generates a discrete-time impulse (Dirac delta function), fundamental for impulse response analysis.

```python
scipy.signal.unit_impulse(
    shape: int or tuple[int],  # Output shape
    idx: int or tuple[int] = None,  # Index at which the impulse is 1 (default: center)
    dtype: type = float        # Data type of output array
)
```

<div STYLE="page-break-after: always;"></div>
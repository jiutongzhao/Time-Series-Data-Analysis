# Initialize Your Data

## Sampling

> **Nyquist-Shannon Sampling Theorem:** A band-limited continuous-time signal $x(t)$ containing no frequency components higher than $f_{max}$,  can be perfectly reconstructed from its samples if it is sampled at a rate:
$$
f_{max} \le f_s/2
$$

The frequency upper limitation $f_s/2$ is also called ***Nyquist Frequency***.

When you measure a high frequency signal with a low cadence instrument, you will not only miss the high frequency component, **<u>but also measure an erroneous signal</u>**, so called ***Aliasing***.

<p align = 'center'>
<img src="Figure/figure_aliasing.png" width="100%"/>
</p>
Such a phenomenon is essentially unrelated to the Fourier transform as its frequency range ends up to $f_s/2$ and can be directly observed by naked eye. In real life, aliasing can be visualized by recording the running car wheel (or helicopter propeller) and television (or computer screen) with your smart phone. 

<u>**A sampling frequency of two times of the wave frequency can not guarantee fully capturing the waveform.**</u> This fact is even true for pure sine waves. 



This effect always happens when you (down-)sampling the signal, a common way to avoid it is to apply a low pass filter so that the high frequency component doesn't contribute to the unreal signal. In the instrumental implementation, that filter typically consists of a set of resistor, inductor, and capacity and is putted before the analog-digital converter.

## Read Signals From Data

[Numpy](https://numpy.org/doc/2.2/reference/routines.io.html), [Scipy.io](https://docs.scipy.org/doc/scipy/tutorial/io.html), and [Pandas](https://pandas.pydata.org/docs/reference/io.html) provide several input/output interfaces for reading the commonly used data format.

- MATLAB `.mat` (v7.2 and below): **`scipy.io.loadmat('*.mat')`**—handles structs, nested arrays; HDF5 v7.3 requires HDF5 libraries.
- IDL `.sav`: **`scipy.io.readsav('*.sav')`** gives you a dict-like structure.
- NetCDF3 (`.nc`): **`scipy.io.netcdf.NetCDFFile('*.nc','r')`**, though it’s deprecated—prefer **`netCDF4`/`xarray`**.
- NASA CDF (`.cdf`): **`spacepy.pycdf.CDF('*.cdf');`**—dict-like API with lazy loading, but needs the NASA CDF C library
- Raw Binary, e.g., `.dat`: **`numpy.fromfile('*.dat', dtype = [np.int8, np.float64, ...])`**
- Text, e.g., `.txt, .csv, .TAB`: **`pandas.read_csv('*.txt', sep = '\s+')`**

## Generate Signals Artificially

- **Generate the ** ***Relative Timestamps*** **given two of signal length ($N$), total duration ($T$), and sampling frequency ($f_s$).**

    ```python
    N, T = 100, 1
    t = np.linspace(0, T, N, endpoint = False)
    ```

    Correspondingly, the ungiven parameter can be derived uniquely.

    ```python
    dt = T / N
    fs = 1 / dt
    ```

- **Generate the signal**

  - Real sine wave

  ```python
  omega0 = 2 * np.pi * 20
  sig = np.sin(omega0 * t)
  ```

  - Complex sine wave

  ```python
  sig = np.exp(j * omega * t)
  ```

  <p align = 'center'>
  <img src="Figure/figure_typical_signals.png" width="100%"/>
  </p>

  - **<u>Using these `scipy.signal` built-in functions</u>** helps to **<u>improve your code readability and reduce your chances of creating bugs: </u>**    

    1. **[Chirp Waveform](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.chirp.html) (`chirp`)**
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

    2. **[Gaussian Pulse](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.gausspulse.html#scipy.signal.gausspulse) (`gausspulse`)**

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

    3. **[Square Wave](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.square.html#scipy.signal.square) (`square`)**

    Generates a square waveform, useful in digital signal simulations, PWM applications, and modulation experiments.

    ```python
    scipy.signal.square(
        t: array_like,             # Time array
        duty: float = 0.5          # Duty cycle (fraction of period signal is high)
    )
    ```

    ------

    4. **[Sawtooth Wave](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.sawtooth.html#scipy.signal.sawtooth) (`sawtooth`)**

    Generates a sawtooth waveform, widely used in signal synthesis and electronics simulations.

    ```python
    scipy.signal.sawtooth(
        t: array_like,             # Time array
        width: float = 1           # Width of rising ramp, from 0 to 1 (1 = triangle wave)
    )
    ```

    ------

    5. **[Unit Impulse](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.unit_impulse.html#scipy.signal.unit_impulse) (`unit_impulse`)**

    Generates a discrete-time impulse (Dirac delta function), fundamental for impulse response analysis.

    ```python
    scipy.signal.unit_impulse(
        shape: int or tuple[int],  # Output shape
        idx: int or tuple[int] = None,  # Index at which the impulse is 1 (default: center)
        dtype: type = float        # Data type of output array
    )
    ```

- From ***Unix Timestamps*** to `np.datetime`

  - A Unix timestamp is the number of seconds that have elapsed since January 1, 1970 (midnight UTC/GMT), not counting leap seconds.

      ```python
      t = (t * 1e9).astype('datetime64[ns]')
      t = t.astype('datetime64[s]')
      ```
  
  - If you need to consider leap seconds in your investigation, use `astropy`.
  

<div STYLE="page-break-after: always;"></div>
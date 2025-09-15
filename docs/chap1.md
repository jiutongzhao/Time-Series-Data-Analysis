# Initialize Your Data

## Sampling

All the data that await analysis are yield from sampling, no matter it originates from the real-world observation or a simulation program. When you measure a high frequency signal with a low cadence instrument, you will not only miss the high frequency component, **<u>but also measure an erroneous signal</u>**, so called ***<u>Aliasing</u>***.

> **[Nyquist-Shannon Sampling Theorem](https://en.wikipedia.org/wiki/Nyquist%E2%80%93Shannon_sampling_theorem):** A band-limited continuous-time signal $x(t)$ containing no frequency components higher than $f_{max}$,  can be perfectly reconstructed from its samples if it is sampled at a rate:
$$
f_s > 2f_{max}
$$
> The frequency upper limitation $f_s/2$ is also called ***<u>Nyquist Frequency</u>***.



<p align = 'center'><img src="Figure/figure_aliasing.png" width="100%"/></p><p align = 'center'>
    <i>Aliasing effcet in a virtual signal sampling.</i>
</p>

Such a phenomenon is essentially unrelated to the Fourier transform as its frequency range ends up to $f_s/2$ and can be directly observed by naked eye. In real life, aliasing can be visualized by compressing the image with grid structure or recording the running helicopter propeller/car wheel.

<p align = 'center'><img src="Figure/figure_moire_pattern.png" width="45%"/> <img src="Figure/figure_helicopter.gif" width="38%"/></p>
<p align = 'center'><i>Aliasing effect in daily life. You can also zoom-in the left-most panel to see the difference before/after compression.</i><p>
<u>**Even a sampling frequency of two times of the wave frequency can not guarantee fully capturing the waveform.**</u> This fact is even true for pure sine waves. Ideally, you can capture the complete wave properties when you got a long enough samples when the sampling frequency is slightly higher than the Nyquist frequency. However, every realistic sample has a finite length. The higher sampling frequency you have, the shorter sample length is required.


<p align = 'center'>
<img src="Figure/figure_aliasing_nyquist.png" width="100%"/>
</p><p align = 'center'>
    <i>A sampling frequency (128 Hz) that slightly higher than the Nyquist frequency (2 × 63.4=126.8 Hz). The sampled signal is shown as a wave packet.</i>
</p>

This effect always happens when you (down-)sampling the signal, a common way to avoid it is to apply a low pass filter so that the high frequency component doesn't contribute to the unreal signal. After this low pass filter, the high frequency variations will be suppressed and not shown in the sampled signal. 

This technique and its implementation, `scipy.signal.decimate`, will be introduced in Chapter 5 of this document. 

In the instrumental implementation, that filter typically consists of a set of resistor, inductor, and capacity and is putted before the analog-digital converter.

<p align = 'center'>
<img src="Figure/figure_anti_aliasing_filter_design.jpg" width="100%"/>
</p><p align = 'center'>
    <i>An example circuit diagram of anti-aliasing filter.</i>
</p>


In brief, **Do Not Interpret Your Data with A Target Frequency near or even above Nyquist Frequency**

## Read Signals From Data

[Numpy](https://numpy.org/doc/stable/reference/routines.io.html), [Scipy.io](https://docs.scipy.org/doc/scipy/tutorial/io.html), and [Pandas](https://pandas.pydata.org/docs/reference/io.html) provide several input/output interfaces for reading the commonly used data format.

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
    # Ensure parameter endpoint is set to False, or use
    # t = np.linspace(0, T, N + 1, endpoint = False)
    ```

    Correspondingly, the ungiven parameter can be derived uniquely.

    ```python
    dt = T / N
    fs = 1 / dt
    ```

- **Generate the signal**

  - Real sine wave

  ```python
  omega0 = 2 * np.pi * 6.0
  
  N = 2 ** 10
  t = np.linspace(0, 1, N, endpoint=False)
  dt = t[1] - t[0]
  fs = 1 / dt  # Sampling frequency
  
  sig_sin = np.sin(omega0 * t)
  ```
  
  - Complex sine wave
  
  ```python
  sig = np.exp(j * omega * t)
  ```
  
  The usage of complex signal aids in shortening the code sometime, but also reduce the readability.
  
  <p align = 'center'>
  <img src="Figure/figure_typical_signals.png" width="100%"/>
  </p>

  - **<u>Using these `scipy.signal` built-in functions</u>** helps to **<u>improve your code readability and reduce your chances of creating bugs: </u>**    
  
    ```python
    sig_cos = np.cos(omega0 * t)
    
    # f0: Frequency at t = 0
    # f1: Frequency at t = t1
    sig_chrip = scipy.signal.chirp(t, f0 = omega0 / 2 / np.pi, t1 = 1, f1 = omega0 / 2 / np.pi * 3)
    
    # fc: central frequency
    # bw: bandwidth
    sig_gauss_pulse = scipy.signal.gausspulse((t - np.mean(t)), fc = omega0 / 1 / np.pi, bw = 0.5)
    
    sig_square = scipy.signal.square(omega0 * t)
    
    sig_sawtooth = scipy.signal.sawtooth(omega0 * t)
    
    # idx : None or int or tuple of int or 'mid', optional
    sig_unit_impulse = scipy.signal.unit_impulse(t.size, idx = 'mid')
    ```



## Management of Timestamps

There are several packages in Python for managing timestamps, and the choice depends on your project's specific requirements—such as high precision, time zone support, or compatibility with other libraries. Each package has its own advantages and disadvantages. In general, selecting a well-maintained package with a large user community can simplify finding help and resources when needed. **`datetime` as a built-in package, is always supported by different third-party package, therefore is suitable to be used during conversion**

| Library  |                          Advantages                          |
| :------: | :----------------------------------------------------------: |
|  Pandas  | Rich time series functionality, easy integration with DataFrames **[Use in Reading Table Files]** |
|  NumPy   | Efficient for large arrays and vectorized time operations **[Use in Data Analysis]** |
| datetime | Part of Python standard library for basic date/time operations **[Use in Timestamps Format Conversion]** |
| Astropy  | High precision handling of special time formats and leap seconds **[Use in Time System Conversion]** |

- **<u>Use `numpy.datetime64` to represent the timestamps</u>**. It is a fixed-size data type that represents dates and times in a variety of formats, including year-month-day, hour-minute-second, and nanoseconds since the Unix epoch. It is also compatible with NumPy's array operations, making it easy to perform calculations and manipulations on large datasets.

    ```python
    import numpy as np
    import pandas as pd
    import astropy.time
    import datetime
    
    # Using pandas
    t_pd = pd.Timestamp('2000-01-01T00:00:00')
    
    # Using astropy
    t_astropy = astropy.time.Time('2000-01-01T00:00:00')
    
    # Using datetime
    t_datetime = datetime.datetime(2000, 1, 1, 0, 0)
    
    # Convert to numpy.datetime64
    t_np = np.datetime64(t_pd)
    t_np = np.datetime64(t_astropy.to_datetime())
    t_np = np.datetime64(t_datetime)
    ```


- convert the timestamps to milliseconds or nanoseconds, respectively. The default unit is nanoseconds.

    ```python
    t_np = np.datetime64('2000-01-01T00:00:00')
    t_np_ms = t_np.astype('datetime64[ms]') # np.datetime64('2000-01-01T00:00:00.000')
    t_np_ns = t_np.astype('datetime64[ns]') # np.datetime64('2000-01-01T00:00:00.000000000')
    ```

- A Unix timestamp is the number of seconds that have elapsed since January 1, 1970 (midnight UTC/GMT), not counting leap seconds. Use `*.astype(float)` to convert `numpy.datetime64[s]` to Unix timestamp. When the unit is not seconds, the result is the number of time units since the epoch.

    ```python
    # np.datetime64 to Unix timestamp
    t_np.astype(float) # np.float64(946684800.0)
    t_np.astype('datetime64[ns]').astype(float) # np.float64(9.466848e+17)
    ```
    
- If you need to consider leap seconds in your investigation, use `astropy`.

    ```mermaid
    ---
    config:
      theme: 'base'
    ---
    
    flowchart LR
    A[*pandas.Timestamp*] --*.to_datetime64()--> B[*datetime.datetime*] --np.datetime64(*)--> C[*numpy.datetime64*] --astropy.time.Time(*)--> D[astropy.time.Time]
    B --pd.Timestamp(*)--> A
    C --.astype(object)--> B
    A --astropy.time.Time(*)--> D
    B --astropy.time.Time(*)--> D
    D --.to_datetime(*)--> B
    A --np.datetime64(*)--> C
    C --pd.Timestamp(*)--> A
    ```

- `np.timedelta64(timedelta, unit)` is a fixed-size data type that represents a duration, the difference between two dates or times. It can be used to perform arithmetic operations on dates and times, such as adding or subtracting a certain number of days or hours. The input `timedelta` can only be integer. But you can use `*.astype('timedelta64[s]')` to convert integer seconds to `np.timedelta64`.
  
  
    ```python
    # Generate numpy.datetime64 array using year_array, month_array, and day_array
    year = np.array([2023, 2023, 2023])
    month = np.array([1, 2, 3])
    day = np.array([1, 2, 3])
    t_array = np.array([np.datetime64(f"{y:04d}-{m:02d}-{d:02d}") for y, m, d in zip(year, month, day)])
    
    # Generate numpy.datetime64 array using year_array and doy_array (day of year)
    doy = np.array([1, 32, 60])  # Day of year for Jan 1, Feb 1, Mar 1 in a non-leap year
    t_array = np.array([np.datetime64(f"{y:04d}-01-01") + np.timedelta64(d - 1, 'D') for y, d in zip(year, doy)])
    
    # Generate numpy.datetime64 array using year and doy (day of year)
    year = 2023
    doy = np.array([1, 32, 60])
    t_array = np.datetime64(f"{year:04d}-01-01") + (doy - 1).astype('timedelta64[D]')
    ```

- `pd.Timestamp()` and `np.datetime64()` can not take **array** as input. So ***list comprehension*** is required in the conversion

    ```python
    # Convert datetime_array to pandas.Timestamp array
    datetime_array_to_pd = [pd.Timestamp(dt) for dt in datetime_array]
    print("datetime_array to pandas.Timestamp array:", datetime_array_to_pd)
    ```
    

## Sampling Method

After get your data, you should know that **what does each timestamps represent?** Is it a **Snapshot** Sample or an **Integrated/Average** Sample?

<p align = 'center'>
<img src="Figure/figure_sampling_methods.png" width="100%"/>
</p>

Due to the limitations of electronic circuitry and measurement principles, the acquisition of real-world signals always requires a finite sampling time. As a result, strictly perfect instantaneous sampling does not exist. However, if the sampling time is much shorter than the interval between data points or the timescale of interest, the measurement can be reasonably approximated as instantaneous. Average sampling itself also acts as a low-pass filter so has the advantage of anti-aliasing.

While artificial signal can be sampled in either ways, just be aware at what you are actually doing if you want to make the comparison between the observations and simulations.



<div STYLE="page-break-after: always;"></div>
# Noise



> **[Lindeberg–Lévy Central Limit Theorem (CLT)](https://en.wikipedia.org/wiki/Central_limit_theorem):**  Suppose $X_1, X_2, X_3,...,X_n$ is a sequence of independent and identically distributed random variables with $\mathbb{E}(X)=\mu$ and $\mathbb{D}(X)=\sigma^2<\infty$. Then, as $n$ approaches infinity, the random variables $\sqrt{n}(\bar{X}_n-\mu)$ converge in distribution to a **normal (or Gaussian)** distribution:
$$
\sqrt{n}(\bar{X}_n-\mu) \mathrel {\overset {d}{\longrightarrow }} \mathcal{N}(0,\sigma^2)
$$
>  

**<u>*Central Limit Theorem*</u>** states that the sum of a large number of **independent random variables**, regardless of their individual distributions, will tend to follow a ***<u>Normal (Gaussian)</u>*** distribution. This principle underpins much of statistical analysis and signal processing. The requirement of independent and random variables is intrinsically related to the definition of noise.

## Description of Noise

Every real‑world measurement—whether from a spacecraft magnetometer, a seismometer, or a studio microphone—contains the desired signal plus unwanted fluctuations we lump together as **<u>*noise*</u>**. These fluctuations can originate from the sensor (thermal agitation, quantization), the environment (vibrations, electromagnetic interference) or the intrinsic randomness of the source.

A compact way to describe noise is by the **power-law index** of the $PSD(f) \propto 1/f^\beta$, $\mathbf{\beta}$, traditionally labeled with colors:

| Noise Color / Info | White | Pink | Brown (Red) | Blue | Violet |
|:------------------:|:-----:|:----:|:-----------:|:----:|:------:|
| $\beta$            | $\beta=0$ | $\beta=1$ | $\beta=2$ | $\beta=-1$ | $\beta=-2$ |
| Examples           | Thermal/electronic background | Music, biological rhythms | Brownian motion, accumulated random walks | Halftoning dither, mask design | GPS acceleration errors |


These color names condense how energy spreads across frequencies and guide filter choice, window length, and averaging strategies.

<p align = 'center'>
<img src="Figure/figure_noise.png" alt="An example of DFT." width="100%"/>
</p>
<p align = 'center'><i> Three types of colored noise.</i></p>




## Signal-to-Noise Ratio (SNR) and Decibels

**SNR** measures how clearly the signal emerges from noise:
$$
\mathrm{SNR} = \frac{P_{\text{signal}}}{P_{\text{noise}}}
$$
where $P_{\text{signal}}$ and $P_{\text{noise}}$ are the average powers of the signal and noise, respectively.

Because the ratio can span orders of magnitude, it is usually expressed in **decibels (dB)**:
$$
\mathrm{SNR}_{\mathrm{dB}} = 10 \log_{10} \left( \frac{P_{\text{signal}}}{P_{\text{noise}}} \right)
$$
For amplitude‑based measurements (root‑mean‑square values):
$$
\mathrm{SNR}_{\mathrm{dB}} = 20 \log_{10} \left( \frac{A_{\text{signal}}}{A_{\text{noise}}} \right)
$$

***Decibel Quick Reference***


|     Decibel     | -20  | -10    | -6                   | -3                           | -1    |  0   |   1   |             3              |         6         |  10   |  20  |
| :-------------: | ---- | ------ | -------------------- | ---------------------------- | ----- | :--: | :---: | :------------------------: | :---------------: | :---: | :--: |
| Amplitude Ratio | 0.1  | 0.3162 | 0.501 $\approx$ 1/2  | 0.708 $\approx$ $1/\sqrt{2}$ | 0.891 |  1   | 1.122 | 1.413 $\approx$ $\sqrt{2}$ | 1.995 $\approx$ 2 | 3.162 |  10  |
|   Power Ratio   | 0.01 | 0.1    | 0.251 $\approx$  1/4 | 0.501                        | 0.794 |  1   | 1.259 |           1.995            | 3.981 $\approx$ 4 |  10   | 100  |

Due to the fact that $2^{10}\approx10^3$, 3 dB corresponds to a energy ratio of $10^{3/10}=\sqrt[10]{1000}\approx \sqrt[10]{1024}=2$.

The adoption of decibel instead of the conventional physical unit has three advantage:

- It allows the directive addition when compare the amplitude of the signal.
- When you are not confident about the magnitude of the uncalibrated data, you can just use dB to describe the ambiguous intensity.
- The [***Weber–Fechner law***](https://en.wikipedia.org/wiki/Weber-Fechner_law) states that human perception of stimulus intensity follows a logarithmic scale, which is why decibels—being logarithmic units—are used to align physical measurements with human sensory sensitivity, such as in sound and signal strength.

## [Quantization Noise](https://en.wikipedia.org/wiki/Quantization_(signal_processing))

Quantization noise arises when a continuous-valued signal is mapped to discrete levels in an analog-to-digital converter (or by fixed-point rounding in digital signal processing). The difference between the true value $x$ and the quantized value $\Delta\cdot\lfloor {x}/{\Delta} + {1}/{2}\rfloor$  is the ***<u>quantization error</u>***, where $\Delta$ is the quantization step size.

For a uniform mid-tread quantizer with step size $\Delta$, and when the input is sufficiently “busy” (uncorrelated with the quantizer and exercising many codes), the error can be modeled as **additive white noise**, uniformly distributed on $[-\Delta/2,\,+\Delta/2]$. Such a addictive noise is white with a total power of $\sigma_q^2=\Delta^2/12$.  Therefore, the signal-to-quantization-noise ratio (SQNR) increases by **6.02 dB ($\approx 10\mathrm{log}_{10}2$)** for each additional bit of resolution.

The uniformly distributed requirement can be fulfilled by the full-amplitude triangle waves or sawtooth waves, or a fast temporal variation. For small or slowly varying signals, the error becomes correlated with the signal, leading to distortion and spurious tones.

<p align = 'center'><img src="Figure/figure_quantization_noise.png" width="100%"/></p>
<p align = 'center'><i> Signals and power spectra for the sine waves with incresing background with different number of quantization bits.</i></p>

At lower amplitudes the quantization error becomes dependent on the input signal, resulting in distortion. This distortion is created after the anti-aliasing filter, and if these distortions are above 1/2 the sample rate they will alias back into the band of interest. In order to make the quantization error **independent** of the input signal, the signal is <u>***dithered***</u> by adding noise to the signal. This slightly reduces signal-to-noise ratio, but can completely eliminate the distortion.

## [Shot Noise](https://en.wikipedia.org/wiki/Shot_noise)

The uncertainty in the counting measurement create intrinsic randomness called **<u>*Shot noise (Schottky noise)*</u>**. Even with constant signal to be measured by counting the number of events, arrivals are [**<u>*Poisson-distributed*</u>**](https://en.wikipedia.org/wiki/Poisson_distribution), producing fluctuations.

$$
P(x;\lambda)=\frac{\lambda^x e^{-\lambda}}{x!}
$$
where $\lambda$ is the expected number of events in the interval, and $x$ is the actual number of events that occur. The mean and variance of the Poisson distribution are both equal to $\lambda$.

This noise is neither additive or multiplicative, but rather **signal-dependent**: its variance scales with the signal level. For a constant average rate $\bar{x}$, the shot noise power spectrum is white with total power $\sigma^2=\eta\bar{x}^2$ when using a detector with an one-count-level, $\eta$:
$$
\eta:=\frac{\mathrm{Signal\ Intensity}}{\mathbb{E}(\mathrm{Number\ of\ Events})}
$$


Thus, the $PSD$ noise level increases 10 dB for each 10-fold increase in signal intensity, corresponding to a 10-fold increase in the Poisson mean $\lambda$. 

<p align = 'center'><img src="Figure/figure_shot_noise.png" width="100%"/></p>
<p align = 'center'><i> Shot noise for a constant signal with different one-count-level.</i></p>

A more general signal model consists of a constant background $a$ plus a time-varying component $b\cdot\mathrm{sin}(\omega t)$:

$$
x(t)=a+b\cdot\mathrm{sin}(\omega t)
$$
In such case, the high frequency part of the power spectrum is still white, but with a total power $\sigma^2=\eta(a+b^2/2)$ that depends on the average signal level. While, the power at the frequency $\omega$ is boosted by the sinusoidal component, standing out from the white noise floor.



------

## Artificial Noise Generation `np.random.randn`

Colored noise can be generated by manipulating the Fourier coefficients of white noise. The following code snippet demonstrates how to create Brownian (red) noise and Violet noise from white noise using NumPy:

```python
white_noise = np.random.randn(time.size)

brownian_noise_fft = np.fft.rfft(white_noise)
brownian_noise_fft[1:] /= freq[1:] ** 1
brownian_noise_fft[0] = 0
brownian_noise = np.fft.irfft(brownian_noise_fft)

violet_noise_fft = np.fft.rfft(white_noise)
violet_noise_fft[1:] /= freq[1:] ** -1
violet_noise_fft[0] = 0
violet_noise = np.fft.irfft(violet_noise_fft)
```

Besides these two methods, one can also get colored noise by filtering white noise. A colored noise that accurately follows its expected power spectrum requires the order of the filter to be high enough. Even though this method is not as straightforward as the Fourier method, it is more efficient in terms of computation time and memory usage.

------

## [Autoregressive (AR) Model](https://en.wikipedia.org/wiki/Autoregressive_model)

Apart from the completely randomness, many real-world disturbances carry just a hint of “inertia”: the next value mostly echoes the present one, plus a fresh random jolt. Such behavior is well captured by a first-order **<u>*autoregressive model*</u>**:
$$
x[n]=\varphi x[n-1] + \mathcal{N}(0,1)
$$
whose single coefficient $\varphi$ sets the memory length. When $\varphi=0$ the series is white noise; as $\varphi\to1^{-}$ it approaches an integrated (Brownian) path with power piling up at low frequencies. 

<p align = 'center'><img src="Figure/figure_ar1.png" width="100%"/></p>
<p align = 'center'><i> AR1 Time series with different AR1 coefficient (alpha = 0.10, 0.90, and 0.99).</i></p>

The spectrum makes this clear:
$$
\mathbb{E}[PSD(f)]=\frac{\sigma^2(1-\varphi^2)}{1+\varphi^2-2\varphi \mathrm{cos}(2\pi f/f_s)}
$$

which is flat for $\varphi=0$ and climbs like $1/f^{2}$ near $f=0$ for $\varphi\approx1$. Thus AR1 offers the simplest realistic noise model—white at one extreme, red at the other—while remaining easy to simulate and fit.

The $$\varphi$$ parameter can be estimated by the lag-1 autocorrelation of the time series. 
$$
\varphi = \frac{\sum_{n=0}^{N-2}(x[n]-\bar{x})(x[n+1]-\bar{x})}{\sum_{n=0}^{N-1}(x[n]-\bar{x})^2}
$$
The AR model is naturally linked to the high-order dynamic system: For example, a second order ordinary differential equation can be discretized as:
$$
\frac{\mathrm{d}^2x}{\mathrm{d}t^2}=kx\to \frac{x[n]-2x[n-1]+x[n-2]}{\Delta t^2}=kx[n-1]
$$
The right-hand side can be re-arranged as
$$
x[n]=(2+k\Delta t^2)x[n-1]-x[n-2]
$$

$$
x[n]=\varphi_1 x[n-1]+\varphi_2 x[n-2]+\mathcal{N}(0,1)
$$


## Null Hypothesis In the Signal Interpretation

When you see a peak in the Fourier spectrum, you may want to know whether it is a real signal or just a random fluctuation of the noise. To answer this question, you need to set up a null hypothesis that the time series is purely noise, and then calculate the probability of observing such a peak (or higher) under this hypothesis. If this probability (the p-value) is very low, you can reject the null hypothesis and conclude that the peak is likely a real signal.

## "Noise"(Uncertainty) of Noise
From the power spectra of noises, one can see that the PSD of the generated noise may randomly deviates from the theoretical expectation, i.e., the exactly power-law PSD. 

The Fourier coefficient computed as 
$$
\begin{align}

X[k]:=\sum_0^{N-1}x[n]\mathrm{e}^{\mathit{i}2\pi  n k}

\end{align}
$$
can be deemed as a <u>**weighted summation**</u> of the signal $x[n]$. When $x[n]$ are independent identically distributed (*i.i.d*) random variables, their weighted summation approaches the Normal distribution when *N* is large enough, according to the ***Central Limit Theorem***. Thus, the *PSD*, defined as the square sum of the real and imaginary parts, naturally follows the *Kappa* Distribution with the freedom of 2. The above statement requires that he real and imaginary parts are independent to each other, which can be proved by calculating their covariance.



<p align = 'center'>
<img src="Figure/figure_noise_hist.png" width="100%"/>
</p>
### Significance Level

A Fourier spectrum always has some peaks no matter the signal is really periodic or totally random. It is of course not appreciated if you interpret a random Fourier peak as a sign of periodicity. 

The ***significance level*** comes out for the assessment of these Fourier power peak. It uses the hypothesis testing to testify whether the peak is significance. The null hypothesis is that 
$$
x[n]\sim N(\mu, \sigma)
$$
Thus, the power spectral density follows the exponential distribution $PSD/S(f)\sim \chi_2^2$. The $95\%$ and $99\%$ percentile of this distribution is $2.995$ and $4.605$, respectively. Therefore, 95% and 99% significance level of the PSD is $2.995\sigma$ and $4.605\sigma$.

### Welch Method [`scipy.signal.welch`]
Welch proposed that the averaging the power spectral density instead of the coefficient can largely reduce the flutuation levels of the spectrum. Therefore, we may just get a.

The averaging operation must be taken after the conversion from coefficient to power other wise the averaged coefficients are actually unchanged.

This method can be implemented by `scipy.signal.welch` function:

```python
freq, psd = scipy.signal.welch(sig, fs, window = 'hann', nperseg=2 ** 10)
```

Except for averaging, one can also  choose the median of the PSD across different segements and obtain a less disturbed PSD. This choice can be implemented by `scipy.signal.welch(signal, fs, average = 'median')`. The default parameter for `average` is `mean`, corresponding to the normal Welch method.

For each segement, you can also chose the window function to reduce the spectral leakage. The result of this method is shown below:

<p align = 'center'>
<img src="Figure/figure_noise_welch.png" alt="From Wikipedia [Gamma Distribution]." width="100%"/>
</p>

One can also verify that the distribution of the PSD convert to *Gamma* Distribution, which has a ***Probability Density Function (PDF)*** of:

$$
\begin{align}
PDF(x; \alpha, \lambda)=\frac{\lambda^\alpha}{\Gamma(\alpha)} x ^{\alpha - 1} e^{-\lambda x}
\end{align}
$$

The mean and variance of this distribution is $\alpha/\lambda$ and $\alpha / \lambda^2$. When the number of segments ($\alpha$) decrease/increase to 1/$+\infty$, the Gamma distribution degenerate to exponential/normal distribution.

<p align = 'center'>
<img src="Figure/figure_gamma_distribution.png" width="100%"/>
</p>


In ***Bartlett Method***, the ratio of ``N_STEP`` and ``N_PER_SEG`` is fixed at unity, which means every segement has no overlapping with each other. It can be regarded as a special case of the *Welch Method* while it is actually proposed earlier.

### Blackman-Tukey Method

***Blackman-Tukey method*** gives another approach to a high SNR estimation of *PSD* based on the *W.S.S* properties of the signal and *Wiener–Khinchin theorem*. This method consists of three steps:

1. Calculate the (***double-sided***) ACF of the signal
2. Apply a window function to the ACF
3. Do DFT to the windowed ACF.

```mermaid
graph LR;
    A[Signal] --> B[ACF];
    B --> C[Window];
    C --> D[DFT];
    D --> E[PSD];
```



It should be keep in mind that these methods are all build based on the assumption of wide-sense stationarity of the signal.[Explain WSS here]. A noise signal, no matter its color, is wide-sense stationary. However, a real time series of a physics quantity cannot gurantee its wide-sense stationarity. Since W.S.S is the only presumption of these method, they are also termed ***Nonparametric Estimator***.

Apart from splitting the signal into several segments, one can also downsample the signal and get multiple sub-signal with different startup time. However, the maximum frequency of the yield spectrum will also be reduced by a factor of ``N_DOWNSAMPLE``. At the same time, the frequency resolution remains to be $(N\Delta t)^{-1}$. 

<p align = 'center'>
<img src="Figure/figure_noise_blackman_tukey.png" width="100%"/>
</p>


## Signal Over Noise

A signal composed of a deterministic sinusoidal component $s(t)$ and additive noise $n(t)$ can be written as:
$$
x(t) = s(t) + n(t)
$$
Correspondingly, the Fourier coefficient at frequency $f$ is the sum of the signal and noise components::
$$
{X}(f) = {S}(f) + {N}(f)
$$
where ${S}(f)$ is the deterministic signal component (a fixed complex number), and ${N}(f)$ is the Fourier transform of the noise. If the noise $n(t)$ is zero-mean wide-sense stationary, then:
$$
\tilde{N}(f) \sim \mathcal{CN}(0, \sigma_n^2)
$$
That is, ${X}(f)$ is a complex Gaussian random variable:
$$
{X}(f) \sim \mathcal{CN}(\mu, \sigma_n^2), \quad \mu = {S}(f)
$$
The power spectrum estimate is:
$$
{S}_x(f) = |{X}(f)|^2
$$
Since $|{X}(f)|^2$ is the sum of squares of two independent Gaussian variables (real and imaginary parts), it strictly follows a non-central chi-squared distribution:
$$
{S}_x(f) \sim \sigma_n^2 \cdot \chi^2(2, \lambda), \quad \lambda = \frac{|\mu|^2}{\sigma_n^2}
$$
In other words, the deterministic signal provides a **complex offset** (mean $\mu$), and the noise determines the **variance** $\sigma_n^2$. The resulting power spectrum estimate is exactly a **<u>non-central chi-squared distribution</u>** with 2 degrees of freedom.

<p align = 'center'>
<img src="Figure/figure_signal_over_noise_hist.png" width="100%"/>
</p>




<div STYLE="page-break-after: always;"></div>

[https://en.wikipedia.org/wiki/Quantization_(signal_processing)]: https://en.wikipedia.org/wiki/Quantization_(signal_processing)
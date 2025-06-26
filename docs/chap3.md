# Power Spectral Density

## Parseval's Theorem and Energy Conservation
In the physical world, the square power of the amplitude often refers to some kind of ***energy*** or ***power***. For example, the square of the displacement ($x$) of a spring, $x^2$ is proportional to the elastic potential energy ($kx^2/2$, where $k$ describes the stiffness). In plasma physics, electromagnetic field contains the energy density ($u$) written as 

$$
u=u_E + u_B=\frac{1}{2}(\varepsilon_0 \mathit{E}^2 + \frac{1}{\mu_0}\mathit{B}^2)
$$



In this case, the ***energy*** of the signal naturally linked with the ***energy*** of the electromagnetic field. Nevertheless, the energy of a signal is an extensive property as it linearly increases with the length of the sample. In the ordinary investigation, the signal energy is always further converted as signal ***power***, which is an intensive property that describe the amplitude and is independent of signal length. The defition of power, *P*, can be written as:

$$
P= \frac{1}{T}\int_{-T/2}^{T/2}|x(t)|^2 \mathrm{d}t
$$



or 
$$
\begin{align}
P&=\frac{1}{N\Delta t}\sum_{n=0}^{N-1}|x(n\Delta t)|^2 \Delta t\\
&=\frac{1}{N^2\Delta f}\sum_{k=0}^{N-1}|X(k\Delta f)|^2 \Delta f \\
&=\sum_{k=0}^{N-1} \boxed{\frac{1}{Nf_s} |X(k\Delta f)|^2}\, \Delta f\\
&=\sum_{k=0}^{N-1}PSD[k]\Delta f
\end{align}
$$

$$
\begin{align}
\int_{-\infty}^\infty x^2(t)\, dt = \int_{-\infty}^\infty X^2(f)\, df
\end{align}
$$

$$
\sum_{n=0}^{N-1}|x[n]|^2 = \frac{1}{N}\sum_{k=0}^{N-1}|X[k]|^2
$$

for DFT. Considering that DFT yields both positive and negative frequency, we typically fold the DFT result. Naturally, the definition of *power spectral density (PSD)* is given as:

$$
\begin{align}
&\sum_{k=0}^{N-1} PSD[k\Delta f] \Delta f =\\
\mathrm{For\ Even \ }N:\ &\Delta f \left[PSD[0] + \sum_{k=1}^{{N}/{2}-1} 2\cdot PSD[k\Delta f] + PSD[f_{N/2}]\right]\\
\mathrm{For\ Odd \ }N:\ &\Delta f \left[PSD[0] + \sum_{k=1}^{{(N-1)}/{2}} 2\cdot PSD[k\Delta f]\right]
\end{align}
$$



$PSD[0]$ represents the DC component and is ignored in the spectral analysis for the most(but not all) time.

```python
N = coef.size
fs = 1 / dt
psd = (np.abs(coef) ** 2) / (N * fs)

if N % 2 == 0:
    psd[1:-1] *= 2
else:
    psd[1:] *= 2
```

According to the lineairty of $\mathcal{F}$, $X[k]$ should also be proportional to the signal amplitude. Easily catch that the coefficient at the exact wave frequency has the form of 

$$
\begin{align}
|X[k]| = \frac{1}{2}A_k \cdot f_s \cdot T 
\end{align}
$$



1/2 in this equation arises from the fact that $\int_0^{2\pi}\mathrm{sin^2}x \mathrm{d}x=1/2$.

## Wiener–Khinchin Theorem

For a wide-sense stationary (WSS) random process $x(t)$, the **autocorrelation function** depends only on the time difference $\tau$, not on absolute time:
$$
R_x(\tau) = \mathbb{E}[x(t)\,x(t + \tau)]
$$
The **power spectral density** is defined as the **Fourier transform** of the autocorrelation function:
$$
S_x(f) = \int_{-\infty}^{\infty} R_x(\tau)\,e^{-j 2\pi f \tau}\,d\tau
$$
This is known as the **Wiener–Khinchin theorem**, and it is valid *only* under the assumption of WSS. The PSD $S_x(f)$ then describes how the total power of the signal is distributed across different frequency components. The relationship between PSD and the Fourier coefficients has been introduced in the previous section.

This theorem tells the intrinsic relationship between the *PSD* and *ACF*. Its contra-position claims that if the PSD doesn't equal to the Fourier transform of the ACF, the signal is not a *w.s.s* signal. The difference between them signify the nature of the solar wind parameters —— They are different from the NOISE! But, for some specific frequency range, they agree with each other well. It should be noticed that the closeness between them doesn't gurantee the signal to be *w.s.s*.

## Wide-Sense Stationarity

Without WSS, the autocorrelation $R_x(t_1, t_2)$ becomes a function of two independent time variables rather than just the lag $\tau$. In such cases, the expectation of the instantaneous wave power $\mathbb{E}[{x^2(t)}]=R_x(t, t)\neq R_x(0=t-t)$ is not independent on $t$. Hence, the Fourier transform of the autocorrelation no longer represents a meaningful or consistent frequency-domain power measure.

> **Therefore, only stationary processes have a well-defined power spectral density, and only then can the spectrum be interpreted as the distribution of power over frequency.**



This condition separates **deterministic Fourier transforms** (which apply to individual signals) from **statistical spectral analysis** (which applies to ensembles of signals or realizations of random processes).

$$
x(t) = A_1 \mathrm{sin}(\omega_1 t) + A_2 \mathrm{sin} \left (\omega_2t + \frac{1}{2}\beta t^2\right )
$$

$$
\begin{align}

R_x(t_1, t_2) & = \mathbb{E}[{x(t_1)x(t_2)}]\\
& = \mathbb{E} \left [A_1^2 \cdot \mathrm{sin}(\omega_1 t_1) \cdot \mathrm{sin}(\omega_1 t_2) + A_2^2 \cdot \mathrm{sin}\left(\omega_2 t_1 + \frac{1}{2}\beta t_1^2 \right ) \cdot \mathrm{sin}\left(\omega_2 t_2 + \frac{1}{2}\beta t_2^2 \right ) \right] \\
 & + A_0 A_1 \left\{ \mathbb{E}\left[\mathrm{sin}(\omega_1 t_1)\mathrm{sin}\left(\omega_2 t_2 + \frac{1}{2}\beta t_2^2 \right ) + \mathrm{sin}(\omega_1 t_2)\mathrm{sin}\left(\omega_2 t_1 + \frac{1}{2}\beta t_1^2 \right ) \right] \right\}

\end{align}
$$

Assuming the window length $T\gt \tau=t_2-t_1\gg 1/f_0$, the square term can be converted to an univariate function of $\tau = t_2-t_1$ by product-to-sum identity. The cross terms can be treated in a similar way ***unless*** $\omega_1 \approx \omega_2+\beta t$, in which condition the 

$$
\mathrm{sin}(\omega_1 t_1)\mathrm{sin}\left(\omega_2 t_2 + \frac{1}{2}\beta t_2^2 \right )=\frac{1}{2} \left[ \mathrm{cos}\left( \omega_1 t_1 -\omega_2 t_2-\frac{1}{2}\beta t_2^2 \right) - \mathrm{cos}\left( \omega_1 t_1 + \omega_2 t_2+\frac{1}{2}\beta t_2^2 \right) \right]
$$

The second terms traverse the whole wave phase from $0$ to $2\pi$ therefore has a expectation of zero. As $\omega_1 t_1-\omega_2 t_2-\frac{1}{2}\beta t_2^2\approx \omega_1 (t_1-t_2) + \frac{1}{2}\beta t_2^2$, the first term can neither written as a function of $\tau$ nor converge to zero in the statistical sense. Thus, it is not wide-sense stationary. When $\omega_1$ is well separated with $\omega_2+\beta t$, the first term again traverse the whole wave phase and the signal return to wide-sense stationary. 

## What If the Signal Is Not Stationary?

For nonstationary signals, the PSD is ill-defined or misleading. In such cases, time-frequency analysis techniques such as:

- **Short-Time Fourier Transform (STFT)**: analyzes local frequency content assuming approximate stationarity within short windows;
- **Wavelet Transform**: offers multi-scale, adaptive analysis of transient and time-varying features;

can be used to track how the spectrum evolves over time, even though no stationary PSD exists.

## Constant and Linear Detrend

To prevent the periodic signals of interest being drowns in the boring background variation, It is always suggested to detrend the signal before the Fourier transform. Constant and Linear detrends are commonly used but polynomial detrend is also not rare.

```python
scipy.signal.welch(x, fs, detrend = 'constants' (default) | 'linear' | Fasle)
```

<p align = 'center'>
<img src="Figure/figure_detrend.png" width="100%"/>
</p>

A slowing change signal is intrinsically not *w.s.s.* as its mean values varies with time. It turns into *w.s.s* after an appropriate detrend.

## More Properties of Fourier Transform

A super powerful property of Fourier transform is that:
$$
\mathcal{F}\left[\frac{\mathrm{d}}{\mathrm{d}t}x(t)\right]=(i2\pi f)\cdot X(f)
$$



which can be easily proved by doing derivative to the both sides of the inverse Fourier transform:
$$
\begin{align}
\frac{\mathrm{d}}{\mathrm{d}t}[x(t)]&=\int_{-\infty}^{+\infty} X(f) (i2\pi f)e^{i 2 \pi f t} \mathrm{d}f\\
&=\int_{-\infty}^{+\infty} \left[(i2\pi f)\cdot X(f)\right] e^{i 2 \pi f t} \mathrm{d}f\\
&=\mathcal{F}^{-1}\left[(i2\pi f)\cdot X(f)\right]
\end{align}
$$



It can be denoted as 
$$
{{\mathrm{d}}/{\mathrm{d}t}}\leftrightarrow i 2\pi f
$$



One can also extend this property to
$$
({{\mathrm{d}/}{\mathrm{d}t}})^n\leftrightarrow (i 2\pi f)^n
$$



In plasma physics, the conventional way to express the electromagnetic field.

It should be noted that this derivation property change a little bit for discrete Fourier transform:

$$
\begin{align}
\frac{\Delta x(t)}{\Delta  t}&=\int_{-\infty}^{+\infty}X(f) \frac{e^{2\pi i f (t+\Delta t)}-e^{2\pi i f t}}{\Delta t} \mathrm{d}f\\
&=\mathcal{F}^{-1}[\frac{e^{2\pi if \Delta t} - 1}{\Delta t}\cdot X(f)]
\end{align}
$$





| Property                   | Continuous-Time Fourier Transform (FT)                 | Discrete Fourier Transform (DFT/FFT)                         |
| -------------------------- | ------------------------------------------------------ | :----------------------------------------------------------- |
| Linearity                  | $\mathcal{F}\{a\,x(t) + b\,y(t)\} = a\,X(f) + b\,Y(f)$ | $\mathrm{DFT}\{a\,x[n] + b\,y[n]\} = a\,X[k] + b\,Y[k]$      |
| Time Shift                 | $x(t - t_0) \rightarrow X(f)\,e^{-j 2\pi f t_0}$       | $x[n - n_0] \rightarrow X[k]\,e^{-j 2\pi k n_0 / N}$         |
| Frequency Shift            | $x(t)\,e^{j 2\pi f_0 t} \rightarrow X(f - f_0)$        | $x[n]\,e^{j 2\pi k_0 n / N} \rightarrow X[(k - k_0)\bmod N]$ |
| Time-Domain Convolution    | $x(t) * h(t) \rightarrow X(f)\,H(f)$                   | $x[n] \circledast h[n] \rightarrow X[k]\,H[k]$ (circular)    |
| Time-Domain Multiplication | $x(t)\,h(t) \rightarrow X(f) * H(f)$                   | $x[n]\,h[n] \rightarrow X[k] * H[k] / N$                     |
| Derivative (Time Domain)   | $\frac{d^n x(t)}{dt^n} \rightarrow (j 2 \pi f)^n X(f)$ | $\Delta^n x[n] \rightarrow X[k] \cdot (e^{j 2 \pi k / N} - 1)^n$ |
| Conjugate Symmetry         | $x(t) \in \mathbb{R} \Rightarrow X(-f) = X^*(f)$       | $x[n] \in \mathbb{R} \Rightarrow X[N - k] = X^*[k]$          |
| Parseval's Theorem         | $\int |x(t)|^2\,dt = \int |X(f)|^2\,df$                | $\sum |x[n]|^2 = \frac{1}{N} \sum |X[k]|^2$                  |
| Frequency Resolution       | Continuous (infinitesimal)                             | $\Delta f = \frac{f_s}{N}$                                   |
| Spectral Periodicity       | $X(f)$ not periodic                                    | $X[k]$ is periodic with period $N$                           |
| Periodic Input Duality     | Periodic $x(t) \Rightarrow$ discrete $X(f)$            | Periodic $x[n] \Rightarrow$ sparse $X[k]$                    |


<div STYLE="page-break-after: always;"></div>
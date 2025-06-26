# What Else Should You Know About the DFT and FFT?

## Gibbs Phenomenon

Although the Fourier transform can perfectly reconstruct **<u>discrete</u>** signal points, it struggles to represent analytic functions with sharp discontinuities. This is because a finite number of sinusoids cannot fully capture the infinite-frequency behavior near a jump, leading to **the Gibbs phenomenon**—persistent overshoot and ringing near discontinuities. 

<p align = 'center'>
<img src="Figure/figure_gibbs.png" width="100%"/>
</p>
Despite increasing the number of harmonics, the overshoot remains (~9%), becoming narrower but not disappearing. This is a fundamental limitation of the Fourier basis, not a numerical flaw, and in practice, techniques like windowing, filtering, or using alternative bases such as wavelets are employed to mitigate its effects.

## Uncertainty Principle

The Uncertainty Principle in signal processing states that a function cannot be simultaneously localized in both time and frequency: the more precisely you know a signal’s timing, the less precisely you can know its frequency content, and vice versa. 

<p align = 'center'><img src="Figure/figure_uncertainty_principle.png" width="100%"/></p>
<p align = 'center'><i>Top: Gaussian pulses with narrow (purple) to wide (blue) time span; Bottom: Their Fourier coefficient with wide to narrow frequency bandwidth. </i></p>


This trade-off is a fundamental limitation of the Fourier transform and is mathematically expressed through time–bandwidth products. It implies that short-duration signals must occupy a wide frequency range, while narrow-band signals cannot be sharply confined in time—a concept closely related to Heisenberg’s uncertainty principle in quantum mechanics.

## The performance of `numpy.fft.fft` and `scipy.signal.fft`

The invention of the ***(Cooley–Tukey) Fast Fourier Transform (FFT) algorithm*** reduced the time complexity of DFT from $\mathcal{O}(N^2)$ to $\mathcal{O}(N\mathrm{log}N)$ by efficiently decomposing the DFT into smaller computations, i.e., [divide-and-conquer](https://en.wikipedia.org/wiki/Divide-and-conquer_algorithm).  

Most tutorials introduce the ***radix-2*** FFT, which splits the signal into ***two*** sub-signals with exactly the same length and requires the length of the signal to be an integer power of ***2***. This requirement is hard to satisfy in common applications without zero-padding, which actually includes unwanted modification of the original signal. To overcome that, ***radix-3*** and ***radix-5*** FFTs are developed and implemented. 

Still, the divide-and-conquer strategy fails when the signal length *N* consists of at least one big prime number factor (e.g, 10007) as the signal is hard to split. In that situation, the ***Bluestein's algorithm***, which is essentially a ***Chirping-Z transform***, is used. This algorithm takes the $\mathcal{F}$ operation as a convolution and then uses the *convolution theorem* in the calculation of DFT coefficients. The convolution property allows us to extend the signal length to a proper, highly composite number with zero-padding (denoted as *M*), but the coefficients and frequency resolution remain unchanged. The final time complexity of *Bluestein's algorithm* goes to $\mathcal{O}(N+M\mathrm{log}M)$, where the first term originates from the iterate all the frequency component.

<p align = 'center'>
<img src="Figure/figure_numpy_fft_performance.png" width="100%"/>
</p>



From the performance test, we observe that signals with prime-number lengths (dark red dots) often incur higher computational costs. For example:
$$
\begin{align}
N&=181=182-1=2^1\times\boxed{7^1\times13^1}-1\\
N&=197=198-1=2^1\times3^2\times\boxed{11^1}-1\\
\end{align}
$$
In contrast, signals with highly composite number lengths (dark blue dots), such as those with lengths being integer powers of 2, usually have the lowest computation time.

However, some prime numbers like: 
$$
\begin{align}
N&=191=192-1=2^6\times3^1-1\\
N&=199=200-1 = 2^3\times5^2-1
\end{align}
$$
can also exhibit relatively efficient performance due to their proximity to highly factorable numbers.

Modern implementation of the FFT algorithm, such as `pocketfft`, combines the above two methods (*Cooley–Tukey* and *Bluestein*). This *C++* package is used in both `numpy` and `scipy(1.4.0+)`  for their FFT implementation. Besides, `fftw`, which stands for the somewhat whimsical title of *"Fastest Fourier Transform in the West"*, is also very popular and used in the `fft/ifft` functions of *MATLAB*. Its *Python* implementation  can be found in the `pyfftw` package.

The `scipy.signal.fft` additionally provides an input parameter `workers:` *`int, optional`* to assign the maximum number of workers to use for parallel computation. If negative, the value wraps around from `os.cpu_count()`. For parallel computation, you need to input a batch of signals with shape of $N\times K$.

***<u>Reference</u>**:*

1. Cooley, James W., and John W. Tukey, 1965, “An algorithm for the machine calculation of complex Fourier series,” Math. Comput. 19: 297-301.
2. Bluestein, L., 1970, “A linear filtering approach to the computation of discrete Fourier transform”. IEEE Transactions on Audio and Electroacoustics. 18 (4): 451-455.
3. https://dsp.stackexchange.com/questions/24375/fastest-implementation-of-fft-in-c
4. https://www.fftw.org/

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

## Polynomial Trend as Seen by DFT


<p align = 'center'>

<div STYLE="page-break-after: always;"></div>
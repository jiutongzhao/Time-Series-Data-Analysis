"""
Append the three to-do examples (gear-box cepstrum, CMB WK, DWT image
compression) to the right notebooks as new (markdown, code) cell pairs,
without disturbing the existing cells/outputs.

Run it from the repo root.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path


def md_cell(text: str) -> dict:
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": text.splitlines(keepends=True),
    }


def code_cell(text: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": text.splitlines(keepends=True),
    }


def append_cells(nb_path: Path, cells: list[dict]) -> None:
    nb = json.loads(nb_path.read_text(encoding="utf-8"))
    nb["cells"].extend(cells)
    nb_path.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")


# ---------------------------------------------------------------------------
# 1. Cepstrum: gear-box fault detection  (-> chap3)
# ---------------------------------------------------------------------------
GEARBOX_MD = r"""### Worked Example: Gear-Box Fault Detection

A textbook industrial application of cepstral analysis is **gear-mesh
diagnostics**.  A healthy gear-box vibration spectrum is dominated by the
*gear-mesh frequency* $f_m = N_t \, f_r$ (number of teeth $\times$ shaft
rotation frequency) and its harmonics.  A localized fault — a chipped tooth,
say — modulates the mesh signal at the *shaft-rotation frequency* $f_r$,
and an entire family of side-bands $f_m \pm k f_r$ appears around every
mesh harmonic.

In the spectrum these side-bands form a **comb** with comb spacing $f_r$.
The cepstrum collapses each such comb into a single peak (a *rahmonic*) at
the **quefrency** $\tau_r = 1/f_r$ — i.e., the period of one shaft
revolution.  This makes the fault signature far easier to spot than reading
side-bands from the log-spectrum directly.

The synthetic signal below imitates a gear-box with $N_t = 21$ teeth running
at $f_r = 25\,\mathrm{Hz}$.  We modulate the mesh and its first harmonic by
a once-per-revolution amplitude envelope, add a faint tooth-impact impulse
train at $N_t f_r$, and bury everything in white noise.  The cepstrum
clearly shows a peak at $\tau_r = 40\,\mathrm{ms}$ and its rahmonics, while
$f_r$ itself is invisible in the raw signal.
"""

GEARBOX_CODE = r'''# === Cepstrum example: gear-box fault detection ===
# Synthetic gear-mesh signal modulated at the shaft-rotation frequency.
fs   = 20_000.0                     # Hz
T    = 2.0                          # s
t    = np.arange(0, T, 1.0 / fs)
fr   = 25.0                         # shaft rotation frequency (Hz)
Nt   = 21                           # number of teeth
fm   = Nt * fr                      # gear-mesh frequency  -> 525 Hz

# Once-per-revolution amplitude modulation -> sidebands fm +/- k*fr
mod = 1.0 + 0.6 * np.cos(2 * np.pi * fr * t)
gear_signal  = mod * np.sin(2 * np.pi * fm * t)
gear_signal += 0.5 * mod * np.sin(2 * np.pi * 2 * fm * t)

# Faint tooth-impact impulse train at fm = Nt * fr
impulse_train = np.zeros_like(t)
impulse_train[::int(fs / fm)] = 1.0
gear_signal += 0.05 * impulse_train

# Add white noise (SNR ~ 5)
rng = np.random.default_rng(0)
gear_signal += 0.2 * rng.standard_normal(t.size)

# --- Real cepstrum ---
freq        = np.fft.rfftfreq(gear_signal.size, 1.0 / fs)
log_mag     = np.log(np.abs(np.fft.rfft(gear_signal)) + 1e-12)
log_mag    -= log_mag.mean()                       # drop DC -> cleaner cepstrum
cepstrum    = np.fft.irfft(log_mag).real
quefrency   = np.arange(cepstrum.size) / fs        # seconds

# --- Plot ---
fig, axes = plt.subplots(3, 1, figsize=(8, 7), constrained_layout=True)

axes[0].plot(t[:int(0.1 * fs)], gear_signal[:int(0.1 * fs)], lw=0.8)
axes[0].set(xlabel='Time [s]', ylabel='Amplitude',
            title='(a) Gear-box vibration (first 100 ms)')

axes[1].semilogy(freq, np.abs(np.fft.rfft(gear_signal)) / gear_signal.size)
axes[1].axvline(fm,     color='C3', ls='--', label=fr'$f_m={fm:.0f}$ Hz')
axes[1].axvline(2 * fm, color='C3', ls=':')
axes[1].set(xlim=(0, 1500), xlabel='Frequency [Hz]',
            ylabel='|X(f)|', title='(b) Spectrum: side-bands form combs around $f_m$')
axes[1].legend()

mask = (quefrency > 1e-3) & (quefrency < 0.2)
axes[2].plot(quefrency[mask] * 1e3, cepstrum[mask])
for k in range(1, 5):
    axes[2].axvline(k * 1000.0 / fr, color='C3', ls='--', alpha=0.6,
                    label=fr'$k/f_r$' if k == 1 else None)
axes[2].set(xlabel='Quefrency [ms]', ylabel='Cepstrum',
            title=fr'(c) Cepstrum: rahmonic at $1/f_r={1000/fr:.0f}$ ms reveals shaft-rate modulation')
axes[2].legend()

plt.savefig(save_path + 'figure_cepstrum_gearbox.png', dpi=150)
plt.show()
'''

# ---------------------------------------------------------------------------
# 2. Wiener-Khinchin: angular power spectrum of CMB-like signal  (-> chap3)
# ---------------------------------------------------------------------------
CMB_MD = r"""The Wiener-Khinchin theorem is the workhorse behind almost every
estimator of the **angular power spectrum of the Cosmic Microwave
Background**.  Modern pipelines (e.g. Planck, ACT) operate on the
2-sphere, but the underlying logic is one-dimensional:

1. Observe a (nearly) stationary stochastic field $T(\hat n)$.
2. Estimate its two-point correlation function $C(\theta) =
   \langle T(\hat n_1)\,T(\hat n_2)\rangle$ for $\hat n_1\!\cdot\!\hat n_2
   = \cos\theta$.
3. Take the (spherical) Fourier transform — equivalent to expanding
   $C(\theta)$ on Legendre polynomials — to obtain $C_\ell$.

Wiener-Khinchin guarantees that this $C_\ell$ is the same object as
$\langle |a_{\ell m}|^2\rangle$ obtained by a direct harmonic transform of
$T$ — provided the field is statistically isotropic (the angular analogue
of WSS).

To keep the code tutorial-friendly we work in **one dimension**: build a
synthetic temperature strip whose theoretical PSD is a CMB-style
power-law with an *acoustic-peak* bump, draw a realisation in the Fourier
domain, then verify the WK theorem by:

* estimating the PSD directly with Welch's method, and
* estimating it via the FT of the autocorrelation function.

The two estimators agree to within sample-variance, exactly as
Wiener-Khinchin predicts.
"""

CMB_CODE = r'''# === Wiener-Khinchin example: 1-D synthetic "CMB" strip ===
N      = 2 ** 14
dx     = 1.0                # arbitrary "pixel" units
fs     = 1.0 / dx
freq   = np.fft.rfftfreq(N, dx)
freq_pos = np.maximum(freq, freq[1])     # avoid divide-by-zero at DC

# Toy CMB-like target PSD: red-tilted continuum + acoustic-peak bump
ell0, sigma_ell, amp_peak = 0.05, 0.012, 4.0
C_target  = (1.0 / freq_pos) ** 1.2
C_target += amp_peak * np.exp(-0.5 * ((freq_pos - ell0) / sigma_ell) ** 2)
C_target[0] = 0.0                         # zero-mean field

# Draw a Gaussian realisation with this PSD
rng = np.random.default_rng(42)
phase = rng.uniform(0, 2 * np.pi, freq.size)
amp   = np.sqrt(C_target * fs * N / 2.0)
amp[0] = 0.0
T_k   = amp * np.exp(1j * phase)
T_k[-1] = T_k[-1].real if N % 2 == 0 else T_k[-1]
T_x   = np.fft.irfft(T_k, n=N)

# --- Estimate PSD two ways ---
# (1) Welch's method (direct spectral estimator)
f_welch, psd_welch = scipy.signal.welch(
    T_x, fs=fs, nperseg=1024, window='hann', detrend=False,
)

# (2) FT of the autocorrelation function (Wiener-Khinchin)
acf = np.fft.irfft(np.abs(np.fft.rfft(T_x)) ** 2, n=N) / N
acf = np.concatenate([acf[-N // 2:], acf[:N // 2]])  # symmetric about lag 0
lag = (np.arange(N) - N // 2) * dx

# Window the ACF (Blackman-Tukey) before FT-ing it
M = 1024
window = np.blackman(2 * M + 1)
acf_win = acf[N // 2 - M : N // 2 + M + 1] * window
psd_wk  = np.abs(np.fft.rfft(acf_win)) / fs
freq_wk = np.fft.rfftfreq(acf_win.size, dx)

# --- Plot ---
fig, axes = plt.subplots(2, 1, figsize=(8, 7), constrained_layout=True)

axes[0].plot(np.arange(N) * dx, T_x, lw=0.5)
axes[0].set(xlim=(0, 2000), xlabel='x', ylabel=r'$T(x)$',
            title='(a) Synthetic CMB-style temperature strip (one realisation)')

axes[1].loglog(freq_pos, C_target,         'k-',  lw=2.0, label='target $C(f)$')
axes[1].loglog(f_welch[1:], psd_welch[1:], 'C0-', alpha=0.8, label="Welch")
axes[1].loglog(freq_wk[1:], psd_wk[1:],    'C3--', alpha=0.8,
               label=r'$\mathcal{F}[R_x(\tau)]$  (Wiener-Khinchin)')
axes[1].axvline(ell0, color='gray', ls=':', alpha=0.7,
                label=f'acoustic peak @ {ell0:.3f}')
axes[1].set(xlabel='Frequency (analogue of multipole $\\ell$)',
            ylabel='PSD',
            title='(b) PSD: direct estimator vs. FT of ACF agree')
axes[1].legend()

plt.savefig(save_path + 'figure_wk_cmb.png', dpi=150)
plt.show()
'''

# ---------------------------------------------------------------------------
# 3. DWT image compression  (-> chap7)
# ---------------------------------------------------------------------------
DWT_IMG_MD = r"""**Worked example.**  We compress an 8-bit grey-scale photograph by
replacing every detail coefficient whose magnitude is below a
percentile-based threshold with zero, then reconstructing.  Sweeping the
threshold from 0 % (lossless) to 99 % (only the top 1 % of coefficients
kept) traces out a *rate-distortion* curve: how many coefficients we keep
versus the resulting peak signal-to-noise ratio (PSNR).

The Daubechies-4 (`db4`) basis is a good default — compactly supported,
four vanishing moments, and orthogonal so the energy at each level is
preserved.  With 95 % of the coefficients zeroed the picture is still
visually faithful, illustrating the strong sparsity that wavelet bases
expose in natural images.
"""

DWT_IMG_CODE = r'''# === DWT image compression: rate-distortion sweep ===
import pywt

# Use a sample image bundled with the repo if available, otherwise
# fall back to scipy's classic "ascent" 8-bit picture.
try:
    image = plt.imread(read_path + 'figure_indian_cuckoo.png')
    if image.ndim == 3:
        image = image[..., :3].mean(axis=-1)
    image = (image * 255).astype(float) if image.max() <= 1.0 else image.astype(float)
except FileNotFoundError:
    image = scipy.datasets.ascent().astype(float)

wavelet, level = 'db4', 4
coeffs    = pywt.wavedec2(image, wavelet=wavelet, level=level)
arr, slices = pywt.coeffs_to_array(coeffs)

ratios = [0.0, 0.5, 0.9, 0.95, 0.99]
fig, axes = plt.subplots(2, len(ratios), figsize=(3 * len(ratios), 6),
                         constrained_layout=True)
psnrs = []
for col, r in enumerate(ratios):
    if r == 0.0:
        arr_c = arr.copy()
    else:
        thresh = np.quantile(np.abs(arr), r)
        arr_c  = np.where(np.abs(arr) >= thresh, arr, 0.0)

    coeffs_c    = pywt.array_to_coeffs(arr_c, slices, output_format='wavedec2')
    image_recon = pywt.waverec2(coeffs_c, wavelet=wavelet)
    image_recon = image_recon[: image.shape[0], : image.shape[1]]

    mse  = np.mean((image - image_recon) ** 2)
    psnr = 10.0 * np.log10((image.max() ** 2) / max(mse, 1e-12))
    psnrs.append(psnr)
    kept = float((arr_c != 0.0).mean())

    axes[0, col].imshow(image_recon, cmap='gray')
    axes[0, col].set_title(f'{r * 100:.0f}% zeroed\nkeep={kept * 100:.1f}%, '
                            f'PSNR={psnr:.1f} dB')
    axes[0, col].set_axis_off()

    axes[1, col].imshow(np.abs(image - image_recon), cmap='magma')
    axes[1, col].set_title('|residual|')
    axes[1, col].set_axis_off()

plt.savefig(save_path + 'figure_dwt_image_compression.png', dpi=150)
plt.show()

print('keep ratio  PSNR (dB)')
for r, p in zip(ratios, psnrs):
    print(f'{(1 - r) * 100:7.1f}% : {p:6.2f}')
'''


def main(repo_root: Path) -> None:
    chap3 = repo_root / "chap3.ipynb"
    chap7 = repo_root / "chap7.ipynb"

    append_cells(chap3, [
        md_cell("## Worked example: gear-box fault detection (cepstrum)\n\n" + GEARBOX_MD),
        code_cell(GEARBOX_CODE),
        md_cell("## Worked example: Wiener-Khinchin and the CMB angular power spectrum\n\n" + CMB_MD),
        code_cell(CMB_CODE),
    ])
    append_cells(chap7, [
        md_cell("## Worked example: DWT image compression\n\n" + DWT_IMG_MD),
        code_cell(DWT_IMG_CODE),
    ])
    print("Appended example cells.")


if __name__ == "__main__":
    repo = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    main(repo)

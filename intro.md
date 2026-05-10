# Practical Spectral Analysis with Python

Welcome to the Jupyter Book build of *Practical Spectral Analysis with
Python*. Each chapter has two views:

- **Prose chapter** — the readable explanation, lifted from the original
  Docsify pages in `docs/chap*.md`.
- **Notebook** — the same content with every figure-generating code cell
  inline, runnable end-to-end in the local Anaconda environment.

The book covers, in order:

1. **Data initialization** — sampling, units, loading common time-series
   formats.
2. **Frequency domain** — the Fourier transform, frequency resolution,
   leakage, zero-padding.
3. **Power spectral density** — Parseval's theorem, the Wiener–Khinchin
   theorem (with a 1-D synthetic CMB worked example), the cepstrum (with
   a gear-box fault-detection example), the Lomb–Scargle periodogram.
4. **More about DFT/FFT** — the uncertainty principle, FFT-implementation
   performance, angular power spectra on a sphere.
5. **Noise** — colour, generation, autoregressive models, significance
   levels, Welch / Blackman–Tukey estimators.
6. **Signal processing** — interpolation, digital filters, the Hilbert
   transform.
7. **Wavelet & time–frequency** — STFT, continuous and discrete wavelet
   transforms, denoising, image compression.
8. **Multi-channel analysis** — PCA / minimum variance analysis,
   cross-spectral density, coherence, SVD wave analysis.

> **Why a new framework?** The old Docsify site rendered LaTeX through a
> client-side MathJax shim that was slow and inconsistent across
> browsers. Jupyter Book uses MyST + MathJax 3 with the AMS package
> preloaded, so `\begin{align}` and friends typeset exactly the way they
> do in a paper. It also renders `chap*.ipynb` natively, so the runnable
> notebooks become first-class pages in the site.

## Reproducing the book

```bash
# one-time install
pip install -r requirements-jupyterbook.txt

# build the HTML site into _build/html/
jupyter-book build .

# open the result
open _build/html/index.html       # macOS
xdg-open _build/html/index.html   # Linux
start _build\html\index.html      # Windows PowerShell
```

The full build instructions, including how to publish to GitHub Pages,
are in [`tools/README.md`](tools/README.md).

## Source

The book is open-source on
[GitHub](https://github.com/jiutongzhao/Time-Series-Data-Analysis) — pull
requests, typo fixes and worked-example contributions are very welcome.

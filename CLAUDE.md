# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Is

**Practical Spectral Analysis with Python** — an interactive educational book published as a Quarto website. It covers DFT/FFT, PSD, noise, filtering, time-frequency analysis, and multi-channel signals.

Live site: https://jiutongzhao.github.io/Time-Series-Data-Analysis/

## Build Commands

```bash
# First build all pages, then start the dev server (pages served instantly, no per-visit render)
quarto render && quarto preview

# Full static build only → _build/html/
quarto render

# Clear cached notebook outputs and rebuild from scratch
quarto clean && quarto render

# Single chapter (faster iteration)
quarto render chap3.qmd
```

Build outputs go to `_build/html/` (gitignored). GitHub Actions deploys the site to `gh-pages` on every push to `main`.

## Architecture

### Content Layer

Each chapter has two paired files:

| File | Purpose |
|------|---------|
| `chap*.qmd` | Authoritative source — prose, math (KaTeX), and embedded code |
| `chap*.ipynb` | Companion notebook — self-contained figure-generating code cells |

The `.qmd` files pull code snippets from the notebooks into collapsible tabs (`### Code <small>chap*.ipynb · cell N</small>`), so edits to displayed code must be made in both places.

### Assets

- `docs/Figure/` — all PNG/SVG figures referenced by `.qmd` files (canonical location, ~120 files)
- `docs/Data/` — audio/video/data files used in chapter content
- `Figure/` — older root-level copies; prefer `docs/Figure/` for new content
- `_includes/` — HTML snippets injected into every page (`mermaid-init.html`, `figure-modal-scripts.html`)
- `styles.css` — custom CSS (loaded via `_quarto.yml`)

All asset paths in `.qmd` files use `docs/Figure/` and `docs/Data/` (relative to project root), not `../Figure/`.

### Freeze Cache

The `.qmd` files contain no executable `{python}` cells — all code blocks are static Markdown. Therefore `_freeze/` is never created and `execute: freeze: auto` has no effect currently. `quarto preview` renders each page on first visit via Pandoc; run `quarto render` first to pre-build all pages so preview is instant.

### Adding a New Figure

1. Generate the figure in the relevant `chap*.ipynb`
2. Save the PNG to `docs/Figure/figure_<name>.png`
3. Reference it in the `.qmd` as `<img src="docs/Figure/figure_<name>.png" .../>`
4. Commit both the notebook and the PNG

### Adding New Audio/Video to the Site

Audio/video files must be explicitly listed under `resources:` in `_quarto.yml` (the `Data/` folder is too large to copy wholesale):

```yaml
resources:
  - docs/Data/new_file.wav
```

## Python Environment

Core dependencies (from CI and `practical_spectral_analysis_plot.ipynb`):

```
numpy  matplotlib  pandas  scipy  pywavelets  bottleneck  scienceplots  ssqueezepy
jupyter  ipykernel  nbformat
```

## Legacy Docsify Site

The `docs/` folder also contains `.md` files (`chap0_preface.md`, etc.) and `_sidebar.md` — these are the old Docsify source. They are no longer the source of truth; the `.qmd` files at the project root are. Do not edit the `.md` files in `docs/`.

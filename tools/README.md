# tools/

Helper scripts for the **Quarto** site of *Practical Spectral Analysis
with Python*.

## Build

The site is built with [Quarto](https://quarto.org/).  Configuration lives
in `_quarto.yml`; the source pages are `index.qmd` + `chap*.qmd` at the
repo root; assets are in `docs/Figure/` and a small whitelist of files
inside `docs/Data/`.

```bash
bash tools/build_book.sh              # static render -> _build/html/
bash tools/build_book.sh --preview    # dev server, live reload
bash tools/build_book.sh --clean      # wipe _build/ + freeze cache
```

## Scripts (kept)

| Script | What it does |
| --- | --- |
| `build_book.sh` | One-line wrapper around `quarto render` / `quarto preview`. |
| `add_examples.py` | One-shot: appends the three worked-example cells (gear-box cepstrum, WK CMB, DWT image compression) to `chap3.ipynb` / `chap7.ipynb`. |
| `reorganize_notebooks.py` | Idempotent: rewrites every `chap*.ipynb` into a tutorial with section headers and a "Further reading" footer. |
| `inject_figure_tabsets.py` | Idempotent: wraps every `figure_*.png` reference in a Quarto `::: {.panel-tabset}` Figure / Code block, pulling generating code from chap*.ipynb. |
| `cleanup.sh` | Interactive deletion of stale duplicates / unused media. Frees ~2.7 GB. Each phase asks for confirmation. |

## Updating the figure-code tabsets

Quarto's `panel-tabset` is rendered at build time, so the Figure/Code
toggle on every `figure_*.png` is purely static HTML — no JS runtime
needed.  Workflow when you change a notebook cell:

```bash
python3 tools/inject_figure_tabsets.py .   # idempotent re-wrap
bash    tools/build_book.sh --preview      # see the new code
```

# tools/

Helper scripts and the MyST-MD / Jupyter Book **v2** scaffold for
*Practical Spectral Analysis with Python*.

## MyST-MD / Jupyter Book v2

The site is built with **MyST-MD** (the Node-based "Jupyter Book v2"),
which renders LaTeX through KaTeX, ships a polished book theme, and
gives you a live-reload dev server while you write.

### One-time setup

```bash
# 1. Install Node 18+ if you don't already have it.
nvm install --lts && nvm use --lts        # https://github.com/nvm-sh/nvm

# 2. Install the MyST-MD CLI.
npm install -g mystmd
```

`mystmd` puts both `myst` and `jupyter book` on your PATH; they are
aliases for the same binary.

### Daily use

```bash
bash tools/build_book.sh --start    # dev server with live reload
bash tools/build_book.sh            # static build -> _build/html/
bash tools/build_book.sh --clean    # wipe _build/ first
```

The static site lives in `_build/html/` and is what the GitHub Pages
deploy uploads.

### Configuration files

| File | Purpose |
| --- | --- |
| `myst.yml` (repo root) | Project + site config: title, authors, license, math macros, two-level TOC, theme. |
| `intro.md` (repo root) | Landing page. |
| `.github/workflows/jupyter-book.yml` | CI workflow — installs `mystmd` via npm, runs `myst build --html`, publishes `_build/html` to the `gh-pages` branch on every push to `main`. |
| `tools/build_book.sh` | One-line build wrapper with helpful errors when the CLI is missing. |
| `tools/legacy_v1/` | Old v1 config (`_config.yml`, `_toc.yml`, Python requirements). Kept for reference. |

### Don't confuse it with v1

There are two completely different things called "Jupyter Book":

| | v1 | v2 (this repo) |
|---|---|---|
| Language | Python | Node.js |
| Install | `pip install jupyter-book` | `npm install -g mystmd` |
| Build | `jupyter-book build .` (hyphen) | `myst build --html` or `jupyter book build` (space) |
| Config | `_config.yml` + `_toc.yml` | `myst.yml` |

If you accidentally call `jupyter-book` (hyphen, the v1 Python tool)
in this repo it'll abort because there is no `_config.yml` at the root.

### Math, mermaid, raw HTML

`myst.yml` already wires up:

- LaTeX macros (`\FT`, `\iFT`, `\expect`, `\abs`, `\norm`) usable from
  any chapter.
- KaTeX rendering with the AMS package, so `\begin{align}…\end{align}`
  works out of the box.
- Mermaid diagrams (the block in `docs/chap3.md`) render natively in
  v2 — no extra plugin needed.
- `<img src=…>`, `<u>`, `<center>`, `<p>` and other inline HTML in the
  existing `docs/chap*.md` files are passed through.

### Notebook execution

`myst.yml` sets `jupyter.execute_notebooks: cache`, which means the
build re-uses the outputs already stored in `chap*.ipynb`. Re-execute
locally first if you want fresh outputs:

```bash
jupyter nbconvert --to notebook --inplace --execute chap*.ipynb
```

## Notebook helpers

### `add_examples.py`

One-shot, **non-idempotent** — appends three new (markdown + code) cell
pairs to the relevant notebooks:

* `chap3.ipynb` — gear-box fault detection (cepstrum)
* `chap3.ipynb` — Wiener–Khinchin & 1-D CMB-style power spectrum
* `chap7.ipynb` — DWT image compression rate-distortion sweep

Run only after restoring the notebooks from a clean state, otherwise
the cells will be appended again.

```bash
python3 tools/add_examples.py .
```

### `reorganize_notebooks.py`

Idempotent: every `chap*.ipynb` is rewritten so it begins with a
chapter title cell, contains a labelled `## Setup` block, and has a
short markdown header before each code cell. Re-running the script is
a no-op once a notebook has been reorganised.

```bash
python3 tools/reorganize_notebooks.py .
```

The label heuristic prefers, in order:

1. Banner-style comments such as `# === Foo ===`.
2. The first descriptive leading comment in the cell.
3. The first informative statement, ignoring imports, magics, and
   matplotlib boilerplate (`plt.close()`, `plt.show()`, …).

Per-chapter manual overrides live in the `CODE_LABEL_OVERRIDES` dict at
the top of the script. Edit it whenever a generated label is unhelpful.

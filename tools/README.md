# tools/

Helper scripts used to reorganise the per-chapter notebooks and inject the
worked-example sections (cepstrum, Wiener–Khinchin/CMB, DWT image
compression).

## `add_examples.py`

One-shot, **non-idempotent** — appends three new (markdown + code) cell
pairs to the relevant notebooks:

* `chap3.ipynb` — gear-box fault detection (cepstrum)
* `chap3.ipynb` — Wiener–Khinchin & 1-D CMB-style power spectrum
* `chap7.ipynb` — DWT image compression rate-distortion sweep

It does **not** edit the corresponding `docs/chap*.md` files (those were
updated by hand). Run only after restoring the notebooks from a clean
state, otherwise the cells will be appended again.

```bash
python3 tools/add_examples.py .
```

## `reorganize_notebooks.py`

Idempotent: every `chap*.ipynb` is rewritten so it begins with a chapter
title cell, contains a labelled `## Setup` block, and has a short markdown
header before each code cell. Re-running the script is a no-op once a
notebook has been reorganised.

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

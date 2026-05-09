"""Reorganise chap*.ipynb files into self-contained tutorials.

Idempotent. Run from repo root:
    python3 reorganize_notebooks.py /path/to/Time-Series-Data-Analysis
"""
from __future__ import annotations
import json
import re
import sys
from pathlib import Path


CHAPTERS = [
    ("chap1.ipynb", 1, "chap1.md",
     "Data initialization, sampling, units, and how to load common time-series formats."),
    ("chap2.ipynb", 2, "chap2.md",
     "From the time domain to the frequency domain: the Fourier transform, frequency resolution, leakage and zero-padding."),
    ("chap3.ipynb", 3, "chap3.md",
     "Power spectral density, Parseval's theorem, the Wiener-Khinchin theorem, the cepstrum, and Lomb-Scargle periodograms - including worked examples (CMB, gear-box)."),
    ("chap4.ipynb", 4, "chap4.md",
     "Practical considerations of the FFT: the uncertainty principle, performance of numpy.fft vs scipy.signal.fft, and angular-power-spectrum computations on a sphere."),
    ("chap5.ipynb", 5, "chap5.md",
     "Noise: colour, generation, autoregressive modelling, significance levels, and PSD estimators (Welch, Blackman-Tukey)."),
    ("chap6.ipynb", 6, "chap6.md",
     "Signal processing: interpolation/resampling, digital filters, the Hilbert transform."),
    ("chap7.ipynb", 7, "chap7.md",
     "Time-frequency analysis: STFT, continuous and discrete wavelet transforms, denoising, image compression."),
    ("chap8.ipynb", 8, "chap8.md",
     "Multi-channel signals: PCA / minimum variance analysis, cross-spectral density, coherence, SVD wave analysis."),
]

CODE_LABEL_OVERRIDES = {
    "chap3.ipynb": {
        20: "### Worked example: gear-box fault detection (cepstrum)",
        22: "### Worked example: Wiener-Khinchin and the CMB angular power spectrum",
    },
    "chap7.ipynb": {
        20: "### Worked example: DWT image compression",
    },
}

SKIP_PREFIXES = ("import ", "from ", "%", "#")
BOILERPLATE = ("plt.close", "plt.show", "plt.tight_layout",
               "fig.tight_layout", "fig.savefig", "plt.savefig")


def md_cell(text):
    return {"cell_type": "markdown", "metadata": {},
            "source": [ln + "\n" for ln in text.rstrip("\n").splitlines()] or [text]}


def is_already_reorganized(nb):
    if not nb["cells"]:
        return False
    c0 = nb["cells"][0]
    if c0["cell_type"] != "markdown":
        return False
    src = "".join(c0["source"]) if isinstance(c0["source"], list) else c0["source"]
    return src.lstrip().startswith("# Chapter ")


def label_for_code_cell(src):
    lines = [ln for ln in src.splitlines() if ln.strip()]
    if not lines:
        return "### (empty cell)"
    for ln in lines:
        m = re.match(r"\s*#\s*[=\-]{2,}\s*(.+?)\s*[=\-]{2,}\s*$", ln)
        if m:
            return "### " + m.group(1).strip()
    for ln in lines:
        s = ln.lstrip()
        if not s.startswith("#") or s.startswith("#!"):
            break
        text = s.lstrip("# ").strip()
        low = text.lower()
        if not text or low.startswith(("import", "from", "todo", "fixme")):
            continue
        if 4 < len(text) < 110:
            return "### " + text
    for ln in lines:
        s = ln.strip()
        if s.startswith(SKIP_PREFIXES):
            continue
        if any(s.startswith(b) for b in BOILERPLATE):
            continue
        if any(s.startswith(p) for p in ("plt.rcParams", "plt.style")):
            continue
        snippet = s[:80] + ("..." if len(s) > 80 else "")
        return "### `" + snippet + "`"
    s = lines[0].strip()
    return "### `" + s[:80] + ("...`" if len(s) > 80 else "`")


def chapter_title_block(num, title, subtitle, doc_name):
    return ("# Chapter {} - {}\n\n{}\n\n"
            "> The narrative chapter is in [`docs/{}`](docs/{}). "
            "This notebook reproduces every figure and table from that chapter, "
            "organised section-by-section so the code is easy to read alongside "
            "the prose.\n").format(num, title, subtitle, doc_name, doc_name)


def doc_h1(doc_path):
    if not doc_path.exists():
        return doc_path.stem
    for ln in doc_path.read_text(encoding="utf-8").splitlines():
        if ln.lstrip().startswith("# "):
            t = ln.lstrip("# ").strip()
            t = re.sub(r"</?[^>]+>", "", t).strip()
            return t or doc_path.stem
    return doc_path.stem


def reorganise(nb_path, num, doc_name, subtitle, docs_dir):
    nb = json.loads(nb_path.read_text(encoding="utf-8"))
    if is_already_reorganized(nb):
        print("  [skip] " + nb_path.name + " already reorganised")
        return False

    title = doc_h1(docs_dir / doc_name)
    new_cells = [md_cell(chapter_title_block(num, title, subtitle, doc_name))]

    overrides = CODE_LABEL_OVERRIDES.get(nb_path.name, {})
    cells = nb["cells"]
    seen_setup = False

    for idx, cell in enumerate(cells):
        if cell["cell_type"] == "markdown":
            new_cells.append(cell)
            continue

        src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]

        if not seen_setup and ("import " in src or "matplotlib" in src):
            new_cells.append(md_cell(
                "## Setup\n\nCommon imports, plot styling and data paths used "
                "throughout this notebook."))
            seen_setup = True
        else:
            label = overrides.get(idx) or label_for_code_cell(src)
            new_cells.append(md_cell(label))

        new_cells.append(cell)

    new_cells.append(md_cell(
        "## Further reading\n\n"
        "The narrative chapter, complete with derivations and references, is in "
        "[`docs/" + doc_name + "`](docs/" + doc_name + "). "
        "For the full table of contents see [`docs/_sidebar.md`](docs/_sidebar.md)."))

    nb["cells"] = new_cells
    nb_path.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")
    print("  [done] " + nb_path.name + " -> " + str(len(new_cells)) + " cells")
    return True


def main(repo_root):
    docs_dir = repo_root / "docs"
    for nb_name, num, doc_name, subtitle in CHAPTERS:
        nb_path = repo_root / nb_name
        if not nb_path.exists():
            print("  [missing] " + nb_name)
            continue
        reorganise(nb_path, num, doc_name, subtitle, docs_dir)


if __name__ == "__main__":
    repo = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    main(repo)

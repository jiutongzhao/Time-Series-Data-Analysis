# coding: utf-8
r"""Rewrite mermaid LaTeX in chap*.qmd / docs/chap*.md to Unicode at the
source level.  Mermaid doesn't render LaTeX inside node labels; storing
clean Unicode in the source means mermaid simply gets readable text.

Idempotent.  Run from repo root:

    python3 tools/delatex_mermaid.py .
"""
from __future__ import annotations
import re, sys
from pathlib import Path
from typing import List, Tuple

# `\b` treats `_` as a word character, so `\omega_0` would NOT trigger a
# rule ending in `\b`.  We use a negative look-ahead for any letter
# instead, which fires on `\omega_`, `\omega^`, `\omega(`, `\omega,` etc.
NB = r"(?![A-Za-z])"

PLACEHOLDER_OPEN  = ""
PLACEHOLDER_CLOSE = ""

GREEK = {
    "alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ", "epsilon": "ε",
    "zeta": "ζ", "eta": "η", "theta": "θ", "iota": "ι", "kappa": "κ",
    "lambda": "λ", "mu": "μ", "nu": "ν", "xi": "ξ", "pi": "π",
    "rho": "ρ", "sigma": "σ", "tau": "τ", "phi": "φ", "chi": "χ",
    "psi": "ψ", "omega": "ω",
    "varepsilon": "ε", "vartheta": "ϑ", "varphi": "ϕ",
    "Gamma": "Γ", "Delta": "Δ", "Theta": "Θ", "Lambda": "Λ", "Xi": "Ξ",
    "Pi": "Π", "Sigma": "Σ", "Phi": "Φ", "Psi": "Ψ", "Omega": "Ω",
}

REPLACEMENTS: List[Tuple[str, str]] = [
    # math delimiters
    (r"\$\$([^$]*?)\$\$", r"\1"),
    (r"\$([^$\n]+?)\$",   r"\1"),

    # escaped braces -> placeholders
    (r"\\\{", PLACEHOLDER_OPEN),
    (r"\\\}", PLACEHOLDER_CLOSE),

    # wrapped commands
    (r"\\hat\s*\{([^{}]+)\}", "\\1̂"),
    (r"\\bar\s*\{([^{}]+)\}", "\\1̄"),
    (r"\\vec\s*\{([^{}]+)\}", "\\1⃗"),
    (r"\\hat\s+(\S)",          "\\1̂"),
    (r"\\bar\s+(\S)",          "\\1̄"),
    (r"\\vec\s+(\S)",          "\\1⃗"),

    (r"\\(?:mathrm|mathbf|mathcal|mathbb|text|boldsymbol|operatorname)\s*\{([^{}]+)\}",
     r"\1"),

    # operators (use NB so trailing _ doesn't block the match)
    (r"\\langle\s*",     "⟨"),
    (r"\\rangle" + NB,   "⟩"),
    (r"\\Re"      + NB,  "Re"),
    (r"\\Im"      + NB,  "Im"),
    (r"\\times"   + NB,  "×"),
    (r"\\cdot"    + NB,  "·"),
    (r"\\parallel"+ NB,  "∥"),
    (r"\\perp"    + NB,  "⊥"),
    (r"\\to"      + NB,  "→"),
    (r"\\rightarrow" + NB, "→"),
    (r"\\leftarrow"  + NB, "←"),
    (r"\\Rightarrow" + NB, "⇒"),
    (r"\\Leftarrow"  + NB, "⇐"),
    (r"\\sum"        + NB, "Σ"),
    (r"\\prod"       + NB, "Π"),
    (r"\\int"        + NB, "∫"),
    (r"\\partial"    + NB, "∂"),
    (r"\\nabla"      + NB, "∇"),
    (r"\\infty"      + NB, "∞"),
    (r"\\pm"         + NB, "±"),
    (r"\\mp"         + NB, "∓"),
    (r"\\approx"     + NB, "≈"),
    (r"\\sim"        + NB, "∼"),
    (r"\\ne"         + NB, "≠"),
    (r"\\neq"        + NB, "≠"),
    (r"\\le"         + NB, "≤"),
    (r"\\leq"        + NB, "≤"),
    (r"\\ge"         + NB, "≥"),
    (r"\\geq"        + NB, "≥"),

    # \frac{A}{B} -> A/B
    (r"\\frac\s*\{([^{}]+)\}\s*\{([^{}]+)\}", r"\1/\2"),

    # spacing / line breaks
    (r"\\,",  " "),
    (r"\\;",  " "),
    (r"\\:",  " "),
    (r"\\ ",  " "),
    (r"\\\\", " "),
]

_GREEK_RE = re.compile(r"\\(" + "|".join(GREEK) + r")" + NB)


def delatex(s: str) -> str:
    for pat, rep in REPLACEMENTS:
        s = re.sub(pat, rep, s)
    s = _GREEK_RE.sub(lambda m: GREEK[m.group(1)], s)
    prev = None
    while prev != s:
        prev = s
        s = re.sub(r"\{([^{}]*)\}", r"\1", s)
    s = re.sub(r"\\([A-Za-z]+)", r"\1", s)
    s = s.replace(PLACEHOLDER_OPEN, "{").replace(PLACEHOLDER_CLOSE, "}")
    return s


MERMAID_BLOCK = re.compile(
    r"(^```mermaid\s*\n)(.*?)(^```\s*$)",
    re.MULTILINE | re.DOTALL,
)


def rewrite_qmd(text: str) -> tuple[str, int]:
    counter = [0]
    def sub(m: re.Match) -> str:
        counter[0] += 1
        return m.group(1) + delatex(m.group(2)) + m.group(3)
    return MERMAID_BLOCK.sub(sub, text), counter[0]


def main(repo_root: Path) -> None:
    total = 0
    for fp in sorted(list(repo_root.glob("chap*.qmd")) +
                     list((repo_root / "docs").glob("chap*.md"))):
        text = fp.read_text(encoding="utf-8")
        new, n = rewrite_qmd(text)
        if new != text:
            fp.write_text(new, encoding="utf-8")
            print(f"  {fp.relative_to(repo_root)!s:35s}  re-cleaned {n} mermaid block(s)")
            total += n
    print(f"Total mermaid blocks touched: {total}")


if __name__ == "__main__":
    main(Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve())

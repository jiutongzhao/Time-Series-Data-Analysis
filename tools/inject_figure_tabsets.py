"""Wrap every `figure_*.png` reference in a chap*.qmd inside a Quarto
panel-tabset (Figure / Code) using figure_code_map.js as the source of
generating code.

Design notes
------------
* Line-by-line state machine so we can:
    - skip code fences (```...```)
    - skip inside any pre-existing `::: {.panel-tabset}` block (avoids
      nesting tabsets inside the user's own examples on the index page)
* Look-ahead for the multi-line  <p ...>\n<img>\n</p>  wrapper so the
  whole 3-line block is consumed and replaced atomically.  Single-line
  <p ...><img></p> and bare <img> are also handled.
* Idempotent self-healing: any earlier `<!-- fig-tabset:NAME -->` block
  is detected and stripped before re-wrapping, even if the previous run
  left dangling <p>/</p> bits around it.

Run from the repo root:

    python3 tools/build_figure_code_map.py .
    python3 tools/inject_figure_tabsets.py .
"""
from __future__ import annotations
import json, re, sys
from pathlib import Path
from typing import Dict, List


# ---- 1. load figure -> code map -------------------------------------

def load_map(p: Path) -> Dict[str, dict]:
    if not p.exists():
        raise FileNotFoundError(p)
    text = p.read_text(encoding="utf-8")
    m = re.search(r"window\.__figureCodeMap\s*=\s*(\{.*\});", text, re.DOTALL)
    if not m:
        raise ValueError("figure_code_map.js missing payload")
    return json.loads(m.group(1))


# ---- 2. unwrap any previous injection (self-heal) -------------------
# Matches the entire injected block PLUS any orphan <p ...> / </p> tags
# left dangling around it by older buggy runs.

UNWRAP_RE = re.compile(
    r'(?:<p[^>]*>[ \t]*)?'                     # optional dangling <p ...> just before
    r'<!--\s*fig-tabset:(?P<name>[^\s>]+)\s*-->[ \t]*\n'
    r'.*?'
    r'\n[ \t]*:::[ \t]*\n'                     # closing :::
    r'(?:[ \t]*</p>[ \t]*\n)?',                # optional dangling </p> just after
    re.DOTALL,
)

def unwrap(text: str) -> tuple[str, int]:
    n = 0
    def sub(m: re.Match) -> str:
        nonlocal n
        n += 1
        # We don't try to perfectly restore the original `<p>...<img>...</p>`
        # markup; instead we leave a bare <img>.  The next wrap pass will
        # rebuild a fresh tabset around it.
        name = m.group("name")
        return f'<img src="docs/Figure/{name}" width="100%"/>\n'
    return UNWRAP_RE.sub(sub, text), n


# ---- 3. wrap fresh ---------------------------------------------------

IMG_RE = re.compile(
    r"<img[^>]+src=['\"](?:[^'\"]*?/)?"
    r"(?P<name>figure_[A-Za-z0-9_\-.]+\.(?:png|jpg|jpeg|gif|svg))"
    r"['\"][^>]*/?>",
    re.IGNORECASE,
)

P_OPEN_RE  = re.compile(r"^\s*<p\b[^>]*>\s*$", re.IGNORECASE)
P_CLOSE_RE = re.compile(r"^\s*</p>\s*$",       re.IGNORECASE)
INLINE_PIMG_RE = re.compile(
    r"^\s*<p\b[^>]*>\s*<img[^>]+>\s*</p>\s*$",
    re.IGNORECASE,
)
TABSET_OPEN_RE  = re.compile(r"^\s*:::+\s*\{\s*\.panel-tabset")
TABSET_CLOSE_RE = re.compile(r"^\s*:::+\s*$")
CODE_FENCE_RE   = re.compile(r"^\s*```")


def make_block(orig_text: str, name: str, entry: dict) -> str:
    code = entry["code"].rstrip()
    nb = entry.get("notebook", "?")
    cell = entry.get("cell", "?")
    return (
        f"<!-- fig-tabset:{name} -->\n"
        "::: {.panel-tabset}\n"
        "\n"
        "### Figure\n"
        "\n"
        f"{orig_text.strip()}\n"
        "\n"
        f"### Code  <small>{nb} · cell {cell}</small>\n"
        "\n"
        "```python\n"
        f"{code}\n"
        "```\n"
        "\n"
        ":::\n"
    )


def wrap_qmd(text: str, code_map: Dict[str, dict]) -> tuple[str, int]:
    lines = text.splitlines(keepends=True)
    out: List[str] = []
    in_tabset_depth = 0
    in_code_fence = False
    i = 0
    rewrites = 0

    while i < len(lines):
        raw = lines[i]
        line = raw.rstrip("\n")
        stripped = line.strip()

        # --- code fence tracking ---
        if CODE_FENCE_RE.match(line):
            in_code_fence = not in_code_fence
            out.append(raw); i += 1; continue
        if in_code_fence:
            out.append(raw); i += 1; continue

        # --- existing tabset tracking ---
        if TABSET_OPEN_RE.match(line):
            in_tabset_depth += 1
            out.append(raw); i += 1; continue
        if in_tabset_depth > 0 and TABSET_CLOSE_RE.match(line):
            in_tabset_depth -= 1
            out.append(raw); i += 1; continue
        if in_tabset_depth > 0:
            out.append(raw); i += 1; continue

        # --- pattern A: multi-line  <p ...>\n<img>\n</p> ---
        if (i + 2 < len(lines)
                and P_OPEN_RE.match(lines[i].rstrip("\n"))
                and IMG_RE.search(lines[i+1])
                and P_CLOSE_RE.match(lines[i+2].rstrip("\n"))):
            m = IMG_RE.search(lines[i+1])
            name = m.group("name")
            entry = code_map.get(name)
            if entry:
                full = lines[i] + lines[i+1] + lines[i+2]
                out.append(make_block(full, name, entry))
                i += 3; rewrites += 1; continue

        # --- pattern B: single-line  <p ...><img></p> ---
        if INLINE_PIMG_RE.match(line):
            m = IMG_RE.search(line)
            if m:
                name = m.group("name")
                entry = code_map.get(name)
                if entry:
                    out.append(make_block(line, name, entry))
                    i += 1; rewrites += 1; continue

        # --- pattern C: bare <img> on its own line ---
        if IMG_RE.search(line) and stripped.startswith("<img"):
            m = IMG_RE.search(line)
            name = m.group("name")
            entry = code_map.get(name)
            if entry:
                out.append(make_block(line, name, entry))
                i += 1; rewrites += 1; continue

        out.append(raw); i += 1

    return "".join(out), rewrites


# ---- driver ----------------------------------------------------------

def main(repo_root: Path) -> None:
    code_map = load_map(repo_root / "figure_code_map.js")
    print(f"Loaded {len(code_map)} figure -> code entries")

    files = sorted(repo_root.glob("chap*.qmd")) + [repo_root / "index.qmd"]
    grand = 0
    for q in files:
        if not q.exists():
            continue
        original = q.read_text(encoding="utf-8")
        cleaned, removed = unwrap(original)
        wrapped, added = wrap_qmd(cleaned, code_map)
        if removed or added:
            q.write_text(wrapped, encoding="utf-8")
            print(f"  {q.name:25s}  unwrapped {removed:3d}, wrapped {added:3d}")
            grand += added
        else:
            print(f"  {q.name:25s}  no changes")
    print(f"Total fresh tab-sets: {grand}")


if __name__ == "__main__":
    main(Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve())

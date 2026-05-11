"""Rewrite *.qmd asset references at the repo root.

* `Figure/...`  -> `docs/Figure/...`
* `Data/...`    -> `docs/Data/...`

Skips any path that is already prefixed with `docs/` and any URL.
Does NOT touch files in subdirectories (docs/, tools/...).
"""
from __future__ import annotations
import re
import sys
from pathlib import Path

# Match `Figure/...` or `Data/...` only when (a) at start of attribute,
# (b) just after src=" or ![]( or similar.
ATTR = re.compile(r'(?P<lead>\b(?:src|href|poster)\s*=\s*["\'])(?P<path>(?:Figure|Data)/[^"\']+)')
MD_LINK = re.compile(r'(?P<lead>!\[[^\]]*\]\()(?P<path>(?:Figure|Data)/[^)\s]+)')

def fix(text: str) -> tuple[str, int]:
    n = 0
    def repl(m):
        nonlocal n
        n += 1
        return m.group('lead') + 'docs/' + m.group('path')
    text2 = ATTR.sub(repl, text)
    text2 = MD_LINK.sub(repl, text2)
    return text2, n

def main(repo_root: Path) -> None:
    total = 0
    for qmd in sorted(repo_root.glob("*.qmd")):
        old = qmd.read_text(encoding="utf-8")
        new, n = fix(old)
        if n:
            qmd.write_text(new, encoding="utf-8")
            print(f"  {qmd.name:25s} -> {n} reference(s) rewritten")
            total += n
        else:
            print(f"  {qmd.name:25s}    (no changes)")
    print(f"Total references rewritten: {total}")

if __name__ == "__main__":
    main(Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve())

#!/usr/bin/env bash
# Build the MyST-MD / Jupyter Book v2 site.
#
#   bash tools/build_book.sh             # static site build  -> _build/html/
#   bash tools/build_book.sh --start     # dev server with live reload
#   bash tools/build_book.sh --clean     # wipe _build/ first
#
# Requires the Node-based MyST-MD CLI:
#     npm install -g mystmd
#
# (The package is `mystmd` on npm; it installs `myst` and `jupyter book`
# binaries on your PATH.)

set -euo pipefail
cd "$(dirname "$0")/.."

if [[ "${1-}" == "--clean" ]]; then
    rm -rf _build
fi

if ! command -v myst >/dev/null 2>&1; then
    cat <<'MSG'
`myst` (MyST-MD CLI) is not on PATH.

This repo is configured for MyST-MD / Jupyter Book v2. Install with:

    npm install -g mystmd

If you need Node first, on macOS/Linux:

    nvm install --lts && nvm use --lts          # via https://nvm.sh

then re-run this script.

Heads-up: there is also a Python "Jupyter Book v1" (the `pip install
jupyter-book` one). It uses _config.yml + _toc.yml, **not** myst.yml.
The legacy v1 config has been moved to tools/legacy_v1/ for reference.
MSG
    exit 1
fi

if [[ "${1-}" == "--start" ]]; then
    exec myst start
fi

myst build --html
echo
echo "Built site: $(pwd)/_build/html/index.html"

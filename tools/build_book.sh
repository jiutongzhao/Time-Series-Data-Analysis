#!/usr/bin/env bash
# Build the Quarto site (configured by ../_quarto.yml).
#
#   bash tools/build_book.sh                       # static render -> _build/html/
#   bash tools/build_book.sh --preview             # full-site dev server, live reload
#   bash tools/build_book.sh --preview chap3.qmd   # preview ONE file (fast)
#   bash tools/build_book.sh --clean               # wipe _build/ and freeze cache
#
# Requires the Quarto CLI:  https://quarto.org/docs/get-started/

set -euo pipefail
cd "$(dirname "$0")/.."

if ! command -v quarto >/dev/null 2>&1; then
    cat <<'MSG'
`quarto` is not on PATH.

Install the Quarto CLI from
    https://quarto.org/docs/get-started/

Legacy Jupyter Book v1 / v2 configs are archived in tools/legacy_v1/.
The live config is _quarto.yml.
MSG
    exit 1
fi

case "${1-}" in
    --clean)
        rm -rf _build .quarto/_freeze _freeze
        echo "Wiped _build/ and freeze cache."
        ;;
    --preview)
        # `--no-browser` keeps Quarto from invoking xdg-open, which is
        # noisy on WSL where there's no GUI browser. Open the URL
        # manually in your host browser.
        if [[ -n "${2-}" ]]; then
            # Single-file preview is dramatically faster: only re-renders
            # the file you're editing, ignores the other 11 chapters.
            exec quarto preview "$2" --no-browser
        else
            exec quarto preview --no-browser
        fi
        ;;
esac

quarto render
echo
echo "Built site: $(pwd)/_build/html/index.html"

#!/usr/bin/env bash
# Aggressive slim-down for the Time-Series-Data-Analysis repo.
#
# Frees roughly 2.7 GB of stale duplicates + unused media + obsolete
# scaffolding.  Each phase is gated behind a `read -p` prompt so you can
# bail out at any step.
#
# Run from the repo root:
#     bash tools/cleanup.sh           # interactive, asks before each phase
#     bash tools/cleanup.sh --yes     # non-interactive, do all phases
set -u
cd "$(dirname "$0")/.."

YES=0
[[ "${1-}" == "--yes" ]] && YES=1
ask() {
    [[ $YES -eq 1 ]] && return 0
    read -r -p "  ${1} (y/N) " a
    [[ "$a" =~ ^[Yy] ]]
}

du_safe() { du -sh "$@" 2>/dev/null | tail -1 ; }

echo "===== 1. caches / build artefacts ====="
du_safe _build .quarto _freeze 2>/dev/null
ask "delete _build/, .quarto/, _freeze/?" && rm -rf _build .quarto _freeze && echo "  done"

echo
echo "===== 2. obsolete JS modal (replaced by panel-tabset) ====="
echo "  (keeps _includes/mermaid-init.html which the site still uses)"
ls figure_modal.js figure_code_map.js _includes/figure-modal-scripts.html 2>/dev/null
ask "delete figure_modal.js, figure_code_map.js, _includes/figure-modal-scripts.html?" \
    && rm -f figure_modal.js figure_code_map.js _includes/figure-modal-scripts.html \
    && echo "  done"

echo
echo "===== 3. obsolete tools ====="
ls tools/legacy_v1 tools/delatex_mermaid.py tools/build_figure_code_map.py tools/fix_qmd_paths.py 2>/dev/null
ask "delete tools/{legacy_v1,delatex_mermaid.py,build_figure_code_map.py,fix_qmd_paths.py}?" \
    && rm -rf tools/legacy_v1 tools/delatex_mermaid.py \
              tools/build_figure_code_map.py tools/fix_qmd_paths.py && echo "  done"

echo
echo "===== 4. legacy / scratch notebooks ====="
ls readme.ipynb a_guide_to_wave_spectral_analysis_v2.ipynb \
   practical_spectral_analysis_plot.ipynb psatd.ipynb test_util.ipynb \
   Spectral_Analysis_Tutorial.ipynb 2>/dev/null
ask "delete the 6 legacy notebooks above?" && \
    rm -f readme.ipynb a_guide_to_wave_spectral_analysis_v2.ipynb \
          practical_spectral_analysis_plot.ipynb psatd.ipynb test_util.ipynb \
          Spectral_Analysis_Tutorial.ipynb && echo "  done"

echo
echo "===== 5. old single-file monolithic versions ====="
ls practical_spectral_analysis_in_plasma_physics.* 2>/dev/null
ask "delete practical_spectral_analysis_in_plasma_physics.*?" && \
    rm -f practical_spectral_analysis_in_plasma_physics.md \
          practical_spectral_analysis_in_plasma_physics.html \
          practical_spectral_analysis_in_plasma_physics.pdf \
          practical_spectral_analysis_in_plasma_physics_tmp.html && echo "  done"

echo
echo "===== 6. stale root-level Figure/ (duplicate of docs/Figure/) ====="
du_safe Figure
ask "delete root Figure/?" && rm -rf Figure && echo "  done"

echo
echo "===== 7. stale root-level Data/ (2.2 GB duplicate of docs/Data/) ====="
du_safe Data
ask "delete root Data/ (~2.2 GB)?" && rm -rf Data && echo "  done"

echo
echo "===== 8. heavy unused audio in docs/Data/ (315 + 21 + 37 MB) ====="
du_safe docs/Data/headphone_comparison.wav docs/Data/orpheus.wav \
        docs/Data/hd800.wav docs/Data/hd800s.wav \
        docs/Data/audio.zip docs/Data/audio 2>/dev/null
ask "delete the 4 heavy WAVs + audio.zip + audio/?" && \
    rm -f docs/Data/headphone_comparison.wav docs/Data/orpheus.wav \
          docs/Data/hd800.wav docs/Data/hd800s.wav docs/Data/audio.zip && \
    rm -rf docs/Data/audio && echo "  done"

echo
echo "===== 9. Docsify scaffolding (no longer used by Quarto) ====="
ls docs/_sidebar.md docs/_coverpage.md docs/index.html \
   docs/theme-simple.css docs/test.css docs/_test.css \
   docs/README.md docs/README\ -\ Copy.md docs/README\ -\ Copy\ \(2\).md docs/README\ -\ Copy\ \(3\).md 2>/dev/null
ask "delete Docsify scaffolding files in docs/?" && \
    rm -f docs/_sidebar.md docs/_coverpage.md docs/index.html \
          docs/theme-simple.css docs/test.css docs/_test.css \
          docs/README.md "docs/README - Copy.md" \
          "docs/README - Copy (2).md" "docs/README - Copy (3).md" \
          docs/chap0_preface\ -\ Copy.md 2>/dev/null && echo "  done"

echo
echo "===== 10. duplicate docs/chap*.md (canonical now = chap*.qmd at root) ====="
ls docs/chap*.md 2>/dev/null
ask "delete docs/chap*.md (8 files)?" && rm -f docs/chap*.md && echo "  done"

echo
echo "===== 11. miscellaneous status-report markdowns at root ====="
ls BUILD_GUIDE.md BUILD_NOW.md FINAL_STATUS.md PROJECT_STRUCTURE.md \
   READY_TO_BUILD.md README.pdf README_watermark.pdf 2>/dev/null
ask "delete those status/build markdowns and the pdf reports?" && \
    rm -f BUILD_GUIDE.md BUILD_NOW.md FINAL_STATUS.md PROJECT_STRUCTURE.md \
          READY_TO_BUILD.md README.pdf README_watermark.pdf && echo "  done"

echo
echo "===== final size report ====="
du -sh . 2>/dev/null
echo "Done.  Run 'git status' to review the deletions before committing."

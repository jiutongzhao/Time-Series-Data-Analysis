/* ============================================================
 *  Click-to-show-code on figure images.
 *
 *  Behaviour
 *  ---------
 *  Every  <img>  whose filename matches /^figure_.+\.(png|jpg|gif|svg)$/i
 *  becomes clickable.
 *
 *    - If the filename appears in window.__figureCodeMap, the modal
 *      shows the source cell that generated it (notebook + cell).
 *    - Otherwise we still pop up a modal explaining that no source
 *      is currently in the repo, so the click is never silent.
 *
 *  A short banner is logged to the console on attach so you can
 *  confirm in DevTools that the script ran.
 * ============================================================ */
(function () {
  "use strict";

  function basenameOf(src) {
    if (!src) return "";
    var clean = src.split("?")[0].split("#")[0];
    var parts = clean.split("/");
    return parts[parts.length - 1];
  }

  function escapeHtml(s) {
    return String(s).replace(/[&<>]/g, function (ch) {
      return ({"&": "&amp;", "<": "&lt;", ">": "&gt;"})[ch];
    });
  }

  function buildModal() {
    if (document.getElementById("fig-code-modal")) return document.getElementById("fig-code-modal");
    var modal = document.createElement("div");
    modal.id = "fig-code-modal";
    modal.className = "fig-code-modal hidden";
    modal.innerHTML = [
      '<div class="fig-code-backdrop"></div>',
      '<div class="fig-code-panel" role="dialog" aria-modal="true">',
      '  <header class="fig-code-header">',
      '    <div class="fig-code-meta">',
      '      <span class="fig-code-figname"></span>',
      '      <span class="fig-code-source"></span>',
      '    </div>',
      '    <div class="fig-code-actions">',
      '      <button type="button" class="fig-code-copy"  title="Copy code">Copy</button>',
      '      <button type="button" class="fig-code-close" title="Close (Esc)">×</button>',
      '    </div>',
      '  </header>',
      '  <pre class="fig-code-body"><code class="language-python"></code></pre>',
      '</div>',
    ].join("");
    document.body.appendChild(modal);

    function close() { modal.classList.add("hidden"); }
    modal.querySelector(".fig-code-backdrop").addEventListener("click", close);
    modal.querySelector(".fig-code-close").addEventListener("click", close);
    document.addEventListener("keydown", function (e) {
      if (e.key === "Escape" && !modal.classList.contains("hidden")) close();
    });
    modal.querySelector(".fig-code-copy").addEventListener("click", function () {
      var txt = modal.querySelector("code").textContent;
      if (!navigator.clipboard) return;
      navigator.clipboard.writeText(txt).then(function () {
        var btn = modal.querySelector(".fig-code-copy");
        var old = btn.textContent;
        btn.textContent = "Copied!";
        setTimeout(function () { btn.textContent = old; }, 1200);
      });
    });
    return modal;
  }

  function openModal(filename, entry) {
    var modal = buildModal();
    modal.querySelector(".fig-code-figname").textContent = filename;
    var source = modal.querySelector(".fig-code-source");
    var code   = modal.querySelector("code");
    var copyBtn = modal.querySelector(".fig-code-copy");

    if (entry) {
      source.textContent = entry.notebook + "  · cell " + entry.cell;
      code.innerHTML = escapeHtml(entry.code);
      copyBtn.style.display = "";
    } else {
      source.textContent = "no source recorded";
      code.innerHTML = '<span style="color:#888;font-style:italic;">' +
        '# This figure is in the repo, but its generating code is not\n' +
        '# currently part of any chap*.ipynb notebook.\n' +
        '# Add the corresponding cell to the right notebook and re-run\n' +
        '#   python3 tools/build_figure_code_map.py .' +
        '</span>';
      copyBtn.style.display = "none";
    }
    modal.classList.remove("hidden");
  }

  function attach() {
    var map = window.__figureCodeMap || {};
    var imgs = document.querySelectorAll("img");
    var clickable = 0;
    var unmapped  = 0;
    imgs.forEach(function (img) {
      var name = basenameOf(img.getAttribute("src") || "");
      if (!/^figure_.+\.(png|jpg|jpeg|gif|svg)$/i.test(name)) return;
      var entry = map[name];
      img.classList.add(entry ? "fig-clickable" : "fig-clickable fig-clickable-empty");
      img.title = entry
        ? "Click to view generating code (" + entry.notebook + ")"
        : "Click for details (no source recorded)";
      img.addEventListener("click", function (e) {
        // Don't follow links if the image is wrapped in an <a>
        e.preventDefault();
        e.stopPropagation();
        openModal(name, entry || null);
      });
      if (entry) clickable++; else unmapped++;
    });
    console.info(
      "[figure-modal] attached.  clickable:", clickable,
      "  unmapped (still clickable for explainer):", unmapped,
      "  map size:", Object.keys(map).length
    );
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", attach);
  } else {
    attach();
  }
})();

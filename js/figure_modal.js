/* ============================================================
 *  Click-to-show-code on figure images.
 *
 *  figure_code_map.js (198 KB) is loaded lazily — only when the
 *  user first clicks a figure, not on every page load.
 * ============================================================ */
(function () {
  "use strict";

  // Capture script dir at parse time — document.currentScript is null in callbacks.
  var scriptDir = (document.currentScript && document.currentScript.src)
    ? document.currentScript.src.replace(/[^/]*$/, "")
    : "";

  var mapLoaded = false;
  var mapLoading = false;
  var pendingCallbacks = [];

  function loadMap(cb) {
    if (mapLoaded) { cb(); return; }
    pendingCallbacks.push(cb);
    if (mapLoading) return;
    mapLoading = true;
    var s = document.createElement("script");
    s.src = scriptDir + "figure_code_map.js";
    s.onload = function () {
      mapLoaded = true;
      pendingCallbacks.forEach(function (fn) { fn(); });
      pendingCallbacks = [];
    };
    s.onerror = function () {
      mapLoaded = true; // treat as empty map
      pendingCallbacks.forEach(function (fn) { fn(); });
      pendingCallbacks = [];
    };
    document.head.appendChild(s);
  }

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
    var source  = modal.querySelector(".fig-code-source");
    var code    = modal.querySelector("code");
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
    var imgs = document.querySelectorAll("img");
    imgs.forEach(function (img) {
      var name = basenameOf(img.getAttribute("src") || "");
      if (!/^figure_.+\.(png|jpg|jpeg|gif|svg)$/i.test(name)) return;

      // All matching figures get a clickable style; we decide which class
      // after the map loads (first click). Pre-mark as potentially clickable.
      img.classList.add("fig-clickable");
      img.title = "Click to view generating code";
      img.style.cursor = "zoom-in";

      img.addEventListener("click", function (e) {
        e.preventDefault();
        e.stopPropagation();
        loadMap(function () {
          var map = window.__figureCodeMap || {};
          var entry = map[name];
          // Correct the CSS class now that we know
          if (!entry) {
            img.classList.add("fig-clickable-empty");
            img.classList.remove("fig-clickable");
          }
          openModal(name, entry || null);
        });
      });
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", attach);
  } else {
    attach();
  }
})();

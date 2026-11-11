/* PROFETdb: right-click info popups for features & analyses. */
(function () {
  "use strict";
  var POPUP_ID = "profet-info-popup";

  function D() { return window.PROFET_DESCRIPTIONS || { features: {}, analyses: {} }; }

  // Title-Case datatable header -> raw feature id
  var titleToRaw = null;
  function buildTitleMap() {
    titleToRaw = {};
    Object.keys(D().features || {}).forEach(function (id) {
      var title = id.replace(/_/g, " ").replace(/\w\S*/g, function (t) {
        return t.charAt(0).toUpperCase() + t.slice(1).toLowerCase();
      });
      titleToRaw[title] = id;
    });
  }
  function normalize(name) {
    if (!name) return null;
    name = String(name).trim();
    if (name.endsWith("_profile")) name = name.slice(0, -8);
    return name;
  }
  function lookupFeature(name) {
    var feats = D().features || {};
    var n = normalize(name);
    if (n && feats[n]) return { info: feats[n], label: n };
    if (!titleToRaw) buildTitleMap();
    if (name && titleToRaw[name.trim()]) {
      var raw = titleToRaw[name.trim()];
      return { info: feats[raw], label: raw };
    }
    return null;
  }

  function closePopup() {
    var el = document.getElementById(POPUP_ID);
    if (el) el.remove();
    document.removeEventListener("mousedown", onOutside, true);
  }
  function onOutside(e) {
    var el = document.getElementById(POPUP_ID);
    if (el && !el.contains(e.target)) closePopup();
  }
  function showPopup(x, y, label, info) {
    closePopup();
    var div = document.createElement("div");
    div.id = POPUP_ID;
    div.style.cssText =
      "position:fixed;z-index:99999;max-width:320px;background:#fff;border:1px solid #ccc;" +
      "border-radius:6px;box-shadow:0 4px 16px rgba(0,0,0,.2);padding:10px 12px;" +
      "font-size:13px;color:#222;line-height:1.4;";
    var brief = (info && info.brief) ? info.brief : "No description available yet.";
    var url   = (info && info.url) ? info.url : "#";
    div.innerHTML =
      '<div style="font-weight:600;margin-bottom:4px;">' + label + "</div>" +
      '<div style="margin-bottom:8px;">' + brief + "</div>" +
      '<a href="' + url + '" target="_blank" rel="noopener" ' +
      'style="color:#2c7be5;text-decoration:none;">Read more &#8599;</a>';
    document.body.appendChild(div);
    var r = div.getBoundingClientRect();
    div.style.left = Math.max(8, Math.min(x, window.innerWidth  - r.width  - 8)) + "px";
    div.style.top  = Math.max(8, Math.min(y, window.innerHeight - r.height - 8)) + "px";
    setTimeout(function () { document.addEventListener("mousedown", onOutside, true); }, 0);
  }

  function plotHoverFeature(gd) {
    var hd = gd && gd._hoverdata;              // Plotly stores last hovered point(s)
    if (!hd || !hd.length) return null;
    var p = hd[0];
    return p.text || p.x || p.y || p.theta || null;   // pgc/spider=text/theta, heatmap/corr=x|y
  }

  document.addEventListener("contextmenu", function (e) {
    var d = D(), t = e.target;

    // 1) analysis button
    var btn = t.closest && t.closest("button[id]");
    if (btn && d.analyses && d.analyses[btn.id]) {
      e.preventDefault();
      return showPopup(e.clientX, e.clientY, btn.textContent.trim(), d.analyses[btn.id]);
    }
    // 2) datatable header/cell
    var cell = t.closest && t.closest("[data-dash-column]");
    if (cell) {
      var m = lookupFeature(cell.getAttribute("data-dash-column"));
      if (m && m.info) { e.preventDefault(); return showPopup(e.clientX, e.clientY, m.label, m.info); }
    }
    // 3) SVG tick / point text label (raw feature text)
    var txt = t.closest && t.closest("text");
    if (txt && txt.textContent) {
      var m2 = lookupFeature(txt.textContent);
      if (m2 && m2.info) { e.preventDefault(); return showPopup(e.clientX, e.clientY, m2.label, m2.info); }
    }
    // 4) plot marker/cell without DOM text -> Plotly hover state
    var gd = t.closest && t.closest(".js-plotly-plot");
    if (gd) {
      var f = plotHoverFeature(gd), m3 = f && lookupFeature(f);
      if (m3 && m3.info) { e.preventDefault(); return showPopup(e.clientX, e.clientY, m3.label, m3.info); }
    }
    // else: leave the browser's native menu
  }, true);

  // hover hint: add an SVG <title> to feature tick labels the first time they're hovered
  document.addEventListener("mouseover", function (e) {
    var t = e.target.closest && e.target.closest("text");
    if (!t || t.getAttribute("data-info-hinted")) return;
    t.setAttribute("data-info-hinted", "1");
    var m = lookupFeature(t.textContent);
    if (m && m.info) {
      var ttl = document.createElementNS("http://www.w3.org/2000/svg", "title");
      ttl.textContent = "Right-click to show info";
      t.appendChild(ttl);
    }
  }, true);

  document.addEventListener("keydown", function (e) { if (e.key === "Escape") closePopup(); });
  window.addEventListener("scroll", closePopup, true);
})();
/* Scroll to a citation's reference entry inside the docs modal.
   Native "#ref-N" hash navigation updates the URL but does not scroll the
   modal's overflow container, so we drive the scroll explicitly. */
(function () {
  "use strict";
  document.addEventListener("click", function (e) {
    var a = (e.target && e.target.closest) ? e.target.closest('a[href]') : null;
    if (!a) return;
    var href = a.getAttribute("href") || "";
    var i = href.indexOf("#ref-");
    if (i === -1) return;
    var el = document.getElementById(href.slice(i + 1));   // e.g. "ref-5"
    if (!el) return;
    e.preventDefault();
    el.scrollIntoView({ behavior: "smooth", block: "start" });
  }, true);
})();

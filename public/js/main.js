// Vulgaris — minimal progressive enhancement.
// Everything here is decorative; the site is fully usable with this disabled.
(function () {
  "use strict";

  var reduceMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  if (reduceMotion || !("IntersectionObserver" in window)) return;

  var headings = document.querySelectorAll(".section-heading");
  if (!headings.length) return;

  var observer = new IntersectionObserver(
    function (entries) {
      entries.forEach(function (entry) {
        if (entry.isIntersecting) {
          entry.target.classList.add("in-view");
          observer.unobserve(entry.target);
        }
      });
    },
    { threshold: 0.4 }
  );

  headings.forEach(function (el) { observer.observe(el); });
})();

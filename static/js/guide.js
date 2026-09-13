/* ==========================================================================
   Vulgaris — hover preview para [[wikilinks]]
   Pasas el mouse sobre un enlace de nota y se abre la nota entera flotando,
   sin salir de la página. Es el hover-preview de Obsidian, en estático:
   se hace fetch de la página destino, se recorta su .note-body y se muestra.
   ========================================================================== */
(function () {
  "use strict";

  var links = document.querySelectorAll("a.wikilink");
  if (!links.length) return;

  // Sin popovers en táctil ni en pantallas chicas: ahí el tap debe navegar.
  if (window.matchMedia("(hover: none)").matches) return;
  if (window.matchMedia("(max-width: 820px)").matches) return;

  var cache = new Map();
  var pop = null;
  var openTimer = null;
  var closeTimer = null;
  var current = null;

  function build() {
    if (pop) return pop;
    pop = document.createElement("div");
    pop.className = "wikipop";
    pop.setAttribute("role", "tooltip");
    pop.hidden = true;
    pop.addEventListener("mouseenter", function () { clearTimeout(closeTimer); });
    pop.addEventListener("mouseleave", scheduleClose);
    document.body.appendChild(pop);
    return pop;
  }

  function fetchNote(href) {
    if (cache.has(href)) return Promise.resolve(cache.get(href));
    return fetch(href, { credentials: "same-origin" })
      .then(function (r) { return r.ok ? r.text() : null; })
      .then(function (html) {
        if (!html) return null;
        var doc = new DOMParser().parseFromString(html, "text/html");
        var body = doc.querySelector(".note-body");
        var title = doc.querySelector("h1");
        if (!body) return null;
        // Los enlaces de adentro no deben ser navegables desde el popover.
        body.querySelectorAll("a").forEach(function (a) { a.removeAttribute("href"); });
        var payload = {
          title: title ? title.textContent : "",
          html: body.innerHTML
        };
        cache.set(href, payload);
        return payload;
      })
      .catch(function () { return null; });
  }

  function place(link) {
    var r = link.getBoundingClientRect();
    var p = pop.getBoundingClientRect();
    var margin = 10;
    var left = r.left + window.scrollX;
    var maxLeft = window.scrollX + document.documentElement.clientWidth - p.width - margin;
    if (left > maxLeft) left = maxLeft;
    if (left < window.scrollX + margin) left = window.scrollX + margin;

    var below = r.bottom + margin + window.scrollY;
    var above = r.top - p.height - margin + window.scrollY;
    var fitsBelow = r.bottom + p.height + margin < document.documentElement.clientHeight;
    pop.style.left = left + "px";
    pop.style.top = (fitsBelow || above < window.scrollY ? below : above) + "px";
  }

  function open(link) {
    var href = link.getAttribute("href").split("#")[0];
    fetchNote(href).then(function (note) {
      if (!note || current !== link) return;
      build();
      pop.innerHTML =
        '<div class="wikipop-head"><span class="prompt">~/</span>' +
        escapeHtml(note.title) + "</div>" +
        '<div class="wikipop-body">' + note.html + "</div>" +
        '<div class="wikipop-foot">click para abrir la nota completa</div>';
      pop.hidden = false;
      pop.classList.add("is-open");
      place(link);
    });
  }

  function scheduleClose() {
    clearTimeout(openTimer);
    closeTimer = setTimeout(function () {
      if (!pop) return;
      pop.classList.remove("is-open");
      pop.hidden = true;
      current = null;
    }, 180);
  }

  function escapeHtml(s) {
    var d = document.createElement("div");
    d.textContent = s;
    return d.innerHTML;
  }

  links.forEach(function (link) {
    link.addEventListener("mouseenter", function () {
      clearTimeout(closeTimer);
      current = link;
      openTimer = setTimeout(function () { open(link); }, 220);
    });
    link.addEventListener("mouseleave", scheduleClose);
    link.addEventListener("focus", function () { current = link; open(link); });
    link.addEventListener("blur", scheduleClose);
  });

  document.addEventListener("keydown", function (e) {
    if (e.key === "Escape") scheduleClose();
  });
})();

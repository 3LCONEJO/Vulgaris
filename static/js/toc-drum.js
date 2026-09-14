/* ==========================================================================
   toc-drum.js — el índice del vault como rodillo cilíndrico.
   - Apaisado: panel lateral; la rueda/el arrastre giran el rodillo.
   - Retrato: el panel se esconde y queda el canto del libro en el borde
     derecho; el pulgar lo abre a pantalla completa.
   - El scroll de la nota gira el rodillo sobre los h2 de la página actual.
   Sin dependencias. Degrada a lista plana con prefers-reduced-motion.
   ========================================================================== */
(function () {
  "use strict";

  var root = document.querySelector("[data-tocd]");
  if (!root) return;

  var STEP_DESK = 24, STEP_MOB = 20;
  var H_DESK = 34, H_MOB = 58;
  var SNAP = 0.18;

  var panel   = root.querySelector(".tocd-panel");
  var windowEl= root.querySelector(".tocd-window");
  var stage   = root.querySelector(".tocd-stage");
  var drum    = root.querySelector(".tocd-drum");
  var gate    = root.querySelector(".tocd-gate");
  var card    = root.querySelector(".tocd-card");
  var cTitle  = root.querySelector(".tocd-card-title");
  var cBranch = root.querySelector(".tocd-card-branch");
  var cEstado = root.querySelector(".tocd-card-estado");
  var cBody   = root.querySelector(".tocd-card-body");
  var prog    = root.querySelector(".tocd-progress");
  var toggle  = root.querySelector(".tocd-toggle");
  var overlay = root.querySelector(".tocd-overlay");
  var rail    = root.querySelector(".tocd-rail");

  var flat = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  if (flat) root.classList.add("is-flat");

  /* --- filas: árbol del vault + los h2 de la nota que se está leyendo ------ */
  var current = drum.querySelector(".tocd-row.is-current");
  var body = document.querySelector(".note-body, .post-body");
  if (current && body) {
    var heads = [].slice.call(body.querySelectorAll("h2[id]"));
    var after = current.nextSibling;
    heads.forEach(function (h) {
      var a = document.createElement("a");
      a.className = "tocd-row tocd-row--head";
      a.href = "#" + h.id;
      a.dataset.depth = String(Math.min(2, (+current.dataset.depth || 0) + 1));
      a.dataset.branch = current.dataset.branch || "0";
      a.dataset.estado = "";
      a.dataset.body = "";
      a.style.setProperty("--c", current.style.getPropertyValue("--c") || "var(--accent)");
      a.innerHTML = '<span class="tocd-indent"></span><span class="tocd-label"></span><span class="tocd-dot"></span>';
      var full = h.textContent.trim();
      // "nf/modules/a_pycistopic/x.nf" -> "x.nf": la ruta entera no cabe en una fila
      var short = /^[\w./-]+\/[\w.-]+$/.test(full) ? full.split("/").pop() : full;
      a.querySelector(".tocd-label").textContent = short;
      if (short !== full) a.title = full;
      a._anchor = h;
      drum.insertBefore(a, after);
    });
  }

  var rows = [].slice.call(drum.querySelectorAll(".tocd-row"));
  if (!rows.length) return;
  var N = rows.length;
  var anchored = rows.filter(function (r) { return r._anchor; });

  /* --- canto: una marca por fila ------------------------------------------ */
  var ticks = rows.map(function (r) {
    var t = document.createElement("span");
    t.className = "tocd-tick";
    t.style.setProperty("--c", r.style.getPropertyValue("--c") || "var(--accent)");
    rail.appendChild(t);
    return t;
  });

  /* --- estado -------------------------------------------------------------- */
  var pos = rows.indexOf(current) >= 0 ? rows.indexOf(current) : 0;
  var target = pos, shown = -1, dragging = null, startY = 0, startPos = 0, moved = false;
  var portrait = false, open = false, railShown = -1;
  if (overlay && "inert" in overlay) overlay.inert = true;
  if (rail) {
    rail.setAttribute("aria-valuemin", "1");
    rail.setAttribute("aria-valuemax", String(N));
    rail.setAttribute("aria-orientation", "vertical");
  }

  var mq = window.matchMedia("(max-aspect-ratio: 1/1), (max-width: 820px)");
  function applyMode() {
    portrait = mq.matches;
    root.classList.toggle("is-portrait", portrait);
    if (!portrait) setOpen(false);
    moveStage();
    shown = -1;
  }
  var cardHome = card && card.parentNode;
  function moveStage() {
    var host = portrait && open ? overlay : windowEl;
    if (stage.parentNode !== host) host.appendChild(stage);
    // La tarjeta vive dentro de .tocd-panel, que en retrato es display:none —
    // ahi no puede verse nunca, por mucho que .is-open le ponga display:flex.
    // Se saca a la raiz (no al overlay: su mask difumina justo donde cae).
    if (card) {
      var chost = portrait && open ? root : cardHome;
      if (card.parentNode !== chost) chost.appendChild(card);
    }
  }
  function setOpen(v) {
    if (open === v) return;
    open = v;
    root.classList.toggle("is-open", v);
    moveStage();
    // el stage lleva enlaces enfocables: no puede quedar dentro de algo
    // marcado aria-hidden mientras esta abierto
    if (overlay) {
      overlay.setAttribute("aria-hidden", v ? "false" : "true");
      if ("inert" in overlay) overlay.inert = !v;
    }
  }
  mq.addEventListener ? mq.addEventListener("change", applyMode) : mq.addListener(applyMode);
  applyMode();

  /* --- posición de lectura -> rodillo -------------------------------------- */
  function indexFromScroll() {
    if (!anchored.length) return rows.indexOf(current) >= 0 ? rows.indexOf(current) : pos;
    var probe = window.scrollY + window.innerHeight * 0.3;
    for (var i = 0; i < anchored.length; i++) {
      var top = anchored[i]._anchor.getBoundingClientRect().top + window.scrollY;
      var nextEl = anchored[i + 1];
      var next = nextEl ? nextEl._anchor.getBoundingClientRect().top + window.scrollY : document.body.scrollHeight;
      if (probe < next) {
        var idx = rows.indexOf(anchored[i]);
        var frac = Math.max(0, Math.min(1, (probe - top) / Math.max(1, next - top)));
        return Math.max(0, Math.min(N - 1, idx + frac * (nextEl ? rows.indexOf(nextEl) - idx : 1)));
      }
    }
    return rows.indexOf(anchored[anchored.length - 1]);
  }

  function goTo(i) {
    var row = rows[i];
    if (!row) return;
    if (row._anchor) {
      window.scrollTo({ top: row._anchor.getBoundingClientRect().top + window.scrollY - window.innerHeight * 0.18, behavior: "smooth" });
    }
  }

  function clamp(v) { return Math.max(0, Math.min(N - 1, v)); }

  /* --- pintado ------------------------------------------------------------- */
  function paint() {
    var H = portrait && open ? H_MOB : H_DESK;
    var STEP = portrait && open ? STEP_MOB : STEP_DESK;
    var idx = Math.round(pos);

    if (!flat) {
      var rad = (H / 2) / Math.tan((STEP * Math.PI) / 360);
      drum.style.transform = "translateZ(" + -rad + "px)";
      rows.forEach(function (el, i) {
        var d = i - pos, ad = Math.abs(d);
        if (ad > 3.4) { el.style.visibility = "hidden"; return; }
        el.style.visibility = "visible";
        el.style.transform = "rotateX(" + -d * STEP + "deg) translateZ(" + rad + "px)";
        el.style.opacity = String(Math.max(0, 1 - ad * 0.3));
        el.style.filter = "blur(" + Math.min(2.2, Math.max(0, ad - 0.35) * 0.8) + "px)";
        el.classList.toggle("is-active", ad < 0.5);
        el.setAttribute("tabindex", ad < 0.5 ? "0" : "-1");
      });
    } else {
      rows.forEach(function (el, i) { el.classList.toggle("is-active", i === idx); });
    }

    ticks.forEach(function (t, i) {
      t.classList.toggle("is-active", Math.abs(i - pos) < 0.5);
      t.style.opacity = Math.abs(i - pos) < 0.5 ? "1" : String(Math.max(0.22, 0.62 - Math.abs(i - pos) * 0.08));
    });

    if (idx !== shown) {
      shown = idx;
      var r = rows[idx];
      var color = r.style.getPropertyValue("--c") || "var(--accent)";
      root.style.setProperty("--tocd-active", color);
      if (cTitle) { cTitle.textContent = r.querySelector(".tocd-label").textContent; cTitle.style.color = r._anchor ? color : ""; }
      if (cBranch) cBranch.textContent = r._anchor ? "en esta nota" : "nota del vault";
      if (cEstado) cEstado.textContent = r.dataset.estado || "";
      if (cBody) cBody.textContent = r.dataset.body || "";
      if (gate) gate.style.background = "linear-gradient(90deg, color-mix(in srgb, " + color + " 12%, transparent), transparent 78%)";
    }
    if (prog) prog.textContent = String(idx + 1).padStart(2, "0") + "/" + String(N).padStart(2, "0");
    if (rail && idx !== railShown) {
      railShown = idx;
      rail.setAttribute("aria-valuenow", String(idx + 1));
      rail.setAttribute("aria-valuetext", rows[idx].querySelector(".tocd-label").textContent);
    }
  }

  function loop() {
    pos += (target - pos) * SNAP;
    if (Math.abs(target - pos) < 0.0015) pos = target;
    paint();
    requestAnimationFrame(loop);
  }
  requestAnimationFrame(loop);

  window.addEventListener("scroll", function () { if (!dragging) target = indexFromScroll(); }, { passive: true });

  /* --- gestos -------------------------------------------------------------- */
  if (!flat) {
    stage.addEventListener("wheel", function (e) {
      if (portrait) return;
      e.preventDefault();
      target = clamp(Math.round(pos) + (e.deltaY > 0 ? 1 : -1));
      goTo(target);
    }, { passive: false });

    stage.addEventListener("pointerdown", function (e) {
      if (portrait) return;
      dragging = "desk"; startY = e.clientY; startPos = pos; moved = false;
      stage.classList.add("is-grabbing");
    });

  }

  /* el canto y su arrastre van fuera de `if (!flat)`: en retrato son la unica
     forma de abrir el indice, con o sin reduced-motion */
  {
    rail.addEventListener("pointerdown", function (e) {
      dragging = "mob"; startY = e.clientY; startPos = pos; moved = false;
      setOpen(true);
      if (navigator.vibrate) navigator.vibrate(4);
      e.preventDefault();
    });

    window.addEventListener("pointermove", function (e) {
      if (!dragging) return;
      var dy = e.clientY - startY;
      if (Math.abs(dy) > 3) moved = true;
      if (dragging === "desk") { pos = target = clamp(startPos - dy / H_DESK); return; }
      e.preventDefault();
      var next = clamp(startPos + dy / (window.innerHeight / (N + 2)));
      if (Math.round(next) !== Math.round(pos) && navigator.vibrate) navigator.vibrate(3);
      pos = target = next;
    }, { passive: false });

    window.addEventListener("pointerup", function () {
      if (!dragging) return;
      var mode = dragging;
      dragging = null;
      stage.classList.remove("is-grabbing");
      var i = clamp(Math.round(pos));
      target = i;
      goTo(i);
      if (mode === "mob") setTimeout(function () { setOpen(false); }, moved ? 260 : 900);
    });
  }

  /* clic en una fila: si es un h2 de esta nota hace scroll, si no navega */
  rows.forEach(function (r, i) {
    r.addEventListener("click", function (e) {
      target = i;
      if (r._anchor) { e.preventDefault(); goTo(i); }
      if (portrait) setOpen(false);
    });
  });

  /* teclado sobre el canto */
  rail.addEventListener("keydown", function (e) {
    if (e.key === "ArrowDown" || e.key === "ArrowUp") {
      e.preventDefault();
      target = clamp(Math.round(pos) + (e.key === "ArrowDown" ? 1 : -1));
      setOpen(true);
      goTo(target);
      clearTimeout(rail._t);
      rail._t = setTimeout(function () { setOpen(false); }, 1200);
    }
    if (e.key === "Enter") rows[Math.round(pos)].click();
  });

  /* botón panel-left */
  if (toggle) {
    toggle.addEventListener("click", function () {
      var collapsed = root.classList.toggle("is-collapsed");
      var label = collapsed ? "Mostrar índice" : "Ocultar índice";
      toggle.setAttribute("aria-label", label);
      toggle.setAttribute("title", label);
      try { localStorage.setItem("tocd-collapsed", collapsed ? "1" : "0"); } catch (err) {}
    });
    try { if (localStorage.getItem("tocd-collapsed") === "1") toggle.click(); } catch (err) {}
  }
})();

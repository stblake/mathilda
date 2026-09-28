/* Interactive builtin-network graph for the Mathilda home page.
 *
 * Renders every public built-in function as a node (coloured by category) with
 * cross-reference edges mined from the reference pages (see site/graph_data.py,
 * which emits assets/builtins_graph.json). Uses the vendored Cytoscape.js +
 * fcose layout. Interactions: hover -> tooltip, click -> open the builtin's doc
 * page, legend -> toggle categories, search -> highlight/centre matches.
 *
 * MkDocs Material has `navigation.instant` on, so initialise through the
 * `document$` observable (fires on first load and every instant navigation)
 * rather than a one-shot DOMContentLoaded, and guard against double-init.
 */
(function () {
  "use strict";

  function initBuiltinGraph() {
    var el = document.getElementById("builtin-graph");
    if (!el || el.dataset.cyInited === "1") return;
    if (typeof cytoscape === "undefined") return; // vendored libs not loaded
    el.dataset.cyInited = "1";

    // Register the fcose layout once (loading twice throws — ignore that).
    try {
      if (window.cytoscapeFcose) cytoscape.use(window.cytoscapeFcose);
    } catch (e) { /* already registered */ }

    var loading = document.getElementById("bg-loading");
    var legendEl = document.getElementById("bg-legend");
    var searchEl = document.getElementById("bg-search");
    var resetEl = document.getElementById("bg-reset");
    var countEl = document.getElementById("bg-count");
    var tip = document.getElementById("bg-tooltip");

    function theme() {
      var s = document.body.getAttribute("data-md-color-scheme");
      return s === "slate" ? "dark" : "light";
    }
    function palette() {
      return theme() === "dark"
        ? { label: "#e6e6e6", edge: "rgba(255,255,255,0.16)", edgeHi: "#ffd479", faded: 0.06 }
        : { label: "#1b1b1f", edge: "rgba(0,0,0,0.12)", edgeHi: "#c77b00", faded: 0.07 };
    }

    var dataUrl = new URL("assets/builtins_graph.json", document.baseURI).href;

    fetch(dataUrl)
      .then(function (r) {
        if (!r.ok) throw new Error("HTTP " + r.status);
        return r.json();
      })
      .then(function (data) { build(data); })
      .catch(function (err) {
        if (loading) loading.textContent = "Could not load the function graph.";
        console.error("[builtin-graph]", err);
      });

    function build(data) {
      var pal = palette();
      var maxDeg = 1;
      data.nodes.forEach(function (n) { if (n.degree > maxDeg) maxDeg = n.degree; });

      var elements = [];
      data.nodes.forEach(function (n) {
        elements.push({ data: {
          id: n.id, label: n.label, category: n.category, color: n.color,
          url: n.url, status: n.status, summary: n.summary, degree: n.degree
        }});
      });
      data.edges.forEach(function (e, i) {
        elements.push({ data: {
          id: "e" + i, source: e.source, target: e.target, weight: e.weight
        }});
      });

      var cy = cytoscape({
        container: el,
        elements: elements,
        wheelSensitivity: 0.25,
        minZoom: 0.05,
        maxZoom: 4,
        pixelRatio: 1,
        textureOnViewport: true,
        hideEdgesOnViewport: true,
        style: [
          { selector: "node", style: {
            "background-color": "data(color)",
            "width": "mapData(degree, 0, " + maxDeg + ", 9, 46)",
            "height": "mapData(degree, 0, " + maxDeg + ", 9, 46)",
            "label": "data(label)",
            "font-size": "mapData(degree, 0, " + maxDeg + ", 10, 26)",
            "color": pal.label,
            // Only the larger (higher-degree) labels show when zoomed out, so the
            // initial view reads as a coloured map, not a wall of text; the rest
            // reveal as you zoom in, and every node still has a hover tooltip.
            "min-zoomed-font-size": 13,
            "text-outline-width": 2,
            "text-outline-color": theme() === "dark" ? "#0b0b0d" : "#ffffff",
            "text-outline-opacity": 0.8,
            "border-width": 0,
            "z-index": 10
          }},
          { selector: "edge", style: {
            "width": "mapData(weight, 1, 12, 0.5, 3)",
            "line-color": pal.edge,
            "curve-style": "haystack",
            "haystack-radius": 0,
            "opacity": 1,
            "z-index": 1
          }},
          { selector: "node.hl", style: {
            "border-width": 3,
            "border-color": pal.edgeHi,
            "z-index": 30
          }},
          { selector: "edge.hl", style: {
            "line-color": pal.edgeHi, "width": 2.2, "z-index": 25
          }},
          { selector: ".faded", style: { "opacity": pal.faded } },
          { selector: "node.match", style: {
            "border-width": 4, "border-color": pal.edgeHi, "z-index": 40
          }},
          { selector: ".hidden", style: { "display": "none" } }
        ],
        layout: { name: "preset" } // real layout runs below with an overlay
      });

      // Run fcose after the overlay has had a chance to paint.
      setTimeout(function () {
        var layout = cy.layout({
          name: "fcose",
          quality: "default",
          randomize: true,
          animate: false,
          // ~1k nodes: a spectral start + a bounded iteration count keeps the
          // one-off layout compute to a second or two rather than freezing.
          sampleSize: 25,
          nodeRepulsion: 13000,
          idealEdgeLength: 78,
          gravity: 0.22,
          numIter: 1500,
          nodeSeparation: 130,
          packComponents: true
        });
        layout.one("layoutstop", function () {
          if (loading) loading.style.display = "none";
          cy.fit(cy.elements(":visible"), 30);
        });
        layout.run();
      }, 30);

      // ---- tooltip ---------------------------------------------------------
      function showTip(node, evt) {
        if (!tip) return;
        var d = node.data();
        var status = d.status ? '<span class="bg-tip-status">' + d.status + "</span>" : "";
        tip.innerHTML =
          '<div class="bg-tip-name">' + d.label + "</div>" +
          '<div class="bg-tip-cat"><span class="bg-swatch" style="background:' +
            d.color + '"></span>' + d.category + status + "</div>" +
          (d.summary && d.summary !== d.label
            ? '<div class="bg-tip-sum">' + escapeHtml(d.summary) + "</div>" : "");
        tip.style.display = "block";
        moveTip(evt);
      }
      function moveTip(evt) {
        if (!tip || tip.style.display !== "block") return;
        var rect = el.getBoundingClientRect();
        var oe = evt.originalEvent || {};
        var x = (oe.clientX || 0) - rect.left + 14;
        var y = (oe.clientY || 0) - rect.top + 14;
        var maxX = rect.width - tip.offsetWidth - 8;
        var maxY = rect.height - tip.offsetHeight - 8;
        tip.style.left = Math.max(6, Math.min(x, maxX)) + "px";
        tip.style.top = Math.max(6, Math.min(y, maxY)) + "px";
      }
      function hideTip() { if (tip) tip.style.display = "none"; }

      cy.on("mouseover", "node", function (evt) {
        var n = evt.target;
        n.addClass("hl");
        n.connectedEdges().addClass("hl");
        n.connectedEdges().connectedNodes().addClass("hl");
        showTip(n, evt);
        el.style.cursor = "pointer";
      });
      cy.on("mousemove", "node", moveTip);
      cy.on("mouseout", "node", function (evt) {
        cy.elements().removeClass("hl");
        hideTip();
        el.style.cursor = "";
      });
      cy.on("tap", "node", function (evt) {
        var url = evt.target.data("url");
        if (url) window.location.href = new URL(url, document.baseURI).href;
      });
      cy.on("pan zoom", hideTip);

      // ---- legend (category toggles) --------------------------------------
      var hiddenCats = {};
      if (legendEl) {
        legendEl.innerHTML = "";
        data.categories.forEach(function (c) {
          var chip = document.createElement("button");
          chip.className = "bg-chip";
          chip.type = "button";
          chip.setAttribute("data-cat", c.name);
          chip.innerHTML =
            '<span class="bg-swatch" style="background:' + c.color + '"></span>' +
            escapeHtml(c.name) + ' <span class="bg-chip-n">' + c.count + "</span>";
          chip.addEventListener("click", function () {
            var on = !hiddenCats[c.name];
            hiddenCats[c.name] = on;
            chip.classList.toggle("off", on);
            var sel = cy.nodes().filter(function (n) {
              return n.data("category") === c.name;
            });
            if (on) sel.addClass("hidden"); else sel.removeClass("hidden");
          });
          legendEl.appendChild(chip);
        });
      }

      // ---- search ----------------------------------------------------------
      function runSearch() {
        var q = (searchEl ? searchEl.value : "").trim().toLowerCase();
        cy.nodes().removeClass("match");
        if (!q) { cy.elements().removeClass("faded"); return; }
        var matches = cy.nodes().filter(function (n) {
          return n.data("label").toLowerCase().indexOf(q) !== -1;
        });
        if (matches.length === 0) { cy.elements().removeClass("faded"); return; }
        cy.elements().addClass("faded");
        matches.removeClass("faded").addClass("match");
        var nbrs = matches.connectedEdges();
        nbrs.removeClass("faded");
        nbrs.connectedNodes().removeClass("faded");
      }
      function gotoTop() {
        var q = (searchEl ? searchEl.value : "").trim().toLowerCase();
        if (!q) return;
        var exact = cy.nodes().filter(function (n) {
          return n.data("label").toLowerCase() === q;
        });
        var pref = exact.length ? exact : cy.nodes().filter(function (n) {
          return n.data("label").toLowerCase().indexOf(q) === 0;
        });
        var target = (pref.length ? pref : cy.nodes().filter(function (n) {
          return n.data("label").toLowerCase().indexOf(q) !== -1;
        }))[0];
        if (target) cy.animate({ center: { eles: target }, zoom: 1.4 }, { duration: 350 });
      }
      if (searchEl) {
        searchEl.addEventListener("input", runSearch);
        searchEl.addEventListener("keydown", function (e) {
          if (e.key === "Enter") { e.preventDefault(); gotoTop(); }
        });
      }
      if (resetEl) {
        resetEl.addEventListener("click", function () {
          if (searchEl) searchEl.value = "";
          hiddenCats = {};
          cy.elements().removeClass("faded match hidden");
          if (legendEl) Array.prototype.forEach.call(
            legendEl.querySelectorAll(".bg-chip"),
            function (c) { c.classList.remove("off"); });
          cy.fit(cy.elements(), 30);
        });
      }
      if (countEl) {
        countEl.textContent =
          data.nodes.length + " functions · " +
          data.edges.length + " links · " +
          data.categories.length + " categories";
      }

      // ---- restyle on light/dark toggle -----------------------------------
      var obs = new MutationObserver(function () {
        var p = palette();
        cy.style().selector("node").style({
          "color": p.label,
          "text-outline-color": theme() === "dark" ? "#0b0b0d" : "#ffffff"
        }).selector("edge").style({ "line-color": p.edge })
          .selector("node.hl").style({ "border-color": p.edgeHi })
          .selector("edge.hl").style({ "line-color": p.edgeHi })
          .selector("node.match").style({ "border-color": p.edgeHi })
          .update();
      });
      obs.observe(document.body, { attributes: true, attributeFilter: ["data-md-color-scheme"] });
    }

    function escapeHtml(s) {
      return String(s).replace(/[&<>"']/g, function (c) {
        return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c];
      });
    }
  }

  if (typeof window.document$ !== "undefined" && window.document$.subscribe) {
    window.document$.subscribe(function () { initBuiltinGraph(); });
  } else {
    document.addEventListener("DOMContentLoaded", initBuiltinGraph);
  }
})();

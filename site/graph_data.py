#!/usr/bin/env python3
"""Build ``site/docs/assets/builtins_graph.json`` for the landing-page network graph.

The home page renders an interactive map of every public built-in function
(Cytoscape.js + fcose): nodes are builtins coloured by category, and edges are
cross-references mined from the generated reference pages. This keeps the graph in
perfect sync with the reference because it consumes ``generate.py``'s own output
(``assets/builtins.json`` + the per-builtin Markdown), so ``make docs`` runs it
straight after ``generate.py``.

Edges (undirected, weighted) come from two signals, strongest first:

  * **Cross-reference links** — every relative link on a page that points at
    another builtin's page (the curated "See also" line and any inline links),
    e.g. ``[Plus](../../arithmetic/Plus/)``. These are hand/spec-curated
    relatedness and get the higher weight.
  * **Code-span mentions** — a backtick code span whose leading identifier is
    another builtin name (``Sin[x]`` -> ``Sin``). These enrich connectivity but
    are noisier, so they are weighted lower and are dropped when the target is a
    "mega-hub" mentioned on a large fraction of all pages (``List``, ``Plus`` …),
    which would otherwise wire that node to everything.

To keep the graph legible each node keeps only its strongest ``PER_NODE_CAP``
edges. Isolated nodes are still emitted so the map shows *every* builtin.

Stdlib only — CI never runs this (the generated JSON is committed); it runs
locally alongside ``generate.py``.
"""
from __future__ import annotations

import colorsys
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
DOCS = HERE / "docs"
ASSETS = DOCS / "assets"
DOC_OUT = DOCS / "documentation"

BUILTINS_JSON = ASSETS / "builtins.json"
OUT_JSON = ASSETS / "builtins_graph.json"

# --- tuning ---------------------------------------------------------------
LINK_WEIGHT = 3          # a curated cross-reference link to another builtin page
MENTION_WEIGHT = 1       # a backtick code-span mentioning another builtin
PER_NODE_CAP = 8         # keep at most this many (strongest) edges per node
ANCHOR_WEIGHT = 1        # faint edge tying an isolated node to its category hub
HUB_MENTION_FRACTION = 0.20  # mention-edges to names cited on > this fraction of
                             # all pages are dropped (they are structural hubs)

# A relative link that targets another builtin's page: `](../../<cat>/<Name>/)`.
# Category index links (`.../index.md)`) and absolute GitHub links don't match
# (they lack the trailing `Name/)` shape), so only real builtin cross-refs pass.
_LINK_RE = re.compile(r"\]\((?:\.\./)+[a-z0-9][a-z0-9-]*/([^/)]+)/\)")
# The leading builtin identifier of a backtick code span, e.g. `Sin[x]` -> Sin,
# `$Assumptions` -> $Assumptions, `PolyGamma` -> PolyGamma.
_CODE_SPAN_RE = re.compile(r"`([^`\n]+)`")
_IDENT_RE = re.compile(r"^\$?[A-Za-z][A-Za-z0-9]*")
_H1_RE = re.compile(r"^#\s+(.+?)\s*$", re.M)


def load_builtins():
    if not BUILTINS_JSON.exists():
        raise SystemExit(
            f"{BUILTINS_JSON} not found — run `make docs` (generate.py) first.")
    return json.loads(BUILTINS_JSON.read_text())


def category_colors(categories):
    """Deterministic, evenly spread colour per category (stable across runs)."""
    cats = sorted(categories)
    n = max(len(cats), 1)
    colors = {}
    for i, cat in enumerate(cats):
        hue = i / n
        # Alternate lightness/saturation slightly so neighbouring hues, of which
        # there are ~38, stay distinguishable.
        light = 0.52 if i % 2 == 0 else 0.62
        sat = 0.68 if i % 3 else 0.60
        r, g, b = colorsys.hls_to_rgb(hue, light, sat)
        colors[cat] = "#%02x%02x%02x" % (round(r * 255), round(g * 255), round(b * 255))
    return colors


def find_page(entry, by_title):
    """Locate the Markdown page for a builtins.json entry.

    The reliable key is the page's H1, which generate.py sets to the builtin
    name; the on-disk filename stem can differ (e.g. FLINT` entries), so we index
    every page by its H1 title and look the name up there, falling back to the
    conventional ``<slug>/<name>.md`` path."""
    hit = by_title.get(entry["name"])
    if hit is not None:
        return hit
    guess = DOC_OUT / entry["slug"] / f'{entry["name"]}.md'
    return guess if guess.exists() else None


def main():
    data = load_builtins()
    names = {e["name"] for e in data}
    name_to_cat = {e["name"]: e["category"] for e in data}

    # Index every generated page by its H1 title (== builtin name).
    by_title = {}
    page_text = {}
    for md in DOC_OUT.rglob("*.md"):
        if md.name == "index.md":
            continue
        text = md.read_text(errors="replace")
        m = _H1_RE.search(text)
        if not m:
            continue
        title = m.group(1).strip()
        by_title[title] = md
        page_text[title] = text

    # Accumulate link and mention weights separately (so mention-edges to hubs
    # can be dropped without touching curated links), plus per-target citation
    # counts used to detect the hubs.
    link_w = {}     # (src, dst) -> number of cross-reference links
    mention_w = {}  # (src, dst) -> 1 if src's page mentions dst in a code span
    cited_by = {}   # dst -> number of distinct pages mentioning it

    for entry in data:
        name = entry["name"]
        page = find_page(entry, by_title)
        if page is None:
            continue
        text = page_text.get(name) or page.read_text(errors="replace")

        # 1) curated cross-reference links
        for target in _LINK_RE.findall(text):
            if target != name and target in names:
                link_w[(name, target)] = link_w.get((name, target), 0) + 1

        # 2) code-span mentions (leading identifier of each span), once per page
        mentioned = set()
        for span in _CODE_SPAN_RE.findall(text):
            im = _IDENT_RE.match(span.strip())
            if not im:
                continue
            ident = im.group(0)
            if ident in names and ident != name:
                mentioned.add(ident)
        for ident in mentioned:
            mention_w[(name, ident)] = 1
            cited_by[ident] = cited_by.get(ident, 0) + 1

    # Structural hubs: names cited on a large fraction of pages. Their *mention*
    # edges are dropped (they would wire that node to everything); curated links
    # to them survive.
    n_pages = max(len(page_text), 1)
    hub_cut = HUB_MENTION_FRACTION * n_pages
    hubs = {d for d, c in cited_by.items() if c > hub_cut}

    # Combine into directed weights, then collapse to undirected.
    directed = {}
    for (src, dst), c in link_w.items():
        directed[(src, dst)] = directed.get((src, dst), 0) + c * LINK_WEIGHT
    for (src, dst), _ in mention_w.items():
        if dst in hubs:
            continue
        directed[(src, dst)] = directed.get((src, dst), 0) + MENTION_WEIGHT

    undirected = {}
    for (src, dst), w in directed.items():
        if w <= 0:
            continue
        key = (src, dst) if src < dst else (dst, src)
        undirected[key] = undirected.get(key, 0) + w

    # Per-node cap: keep each node's strongest edges, then union.
    incident = {}
    for (a, b), w in undirected.items():
        incident.setdefault(a, []).append((w, b))
        incident.setdefault(b, []).append((w, a))
    keep = set()
    for node, lst in incident.items():
        lst.sort(reverse=True)
        for w, other in lst[:PER_NODE_CAP]:
            key = (node, other) if node < other else (other, node)
            keep.add(key)
    edge_weight = {k: undirected[k] for k in keep}

    # Degree so far (used to pick each category's hub).
    deg = {}
    for (a, b) in edge_weight:
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1

    # Attach otherwise-isolated functions to their category's hub (its
    # highest-degree member) so every node joins its colour group — fcose would
    # otherwise scatter disconnected singletons away from their category.
    by_cat_nodes = {}
    for e in data:
        by_cat_nodes.setdefault(e["category"], []).append(e["name"])
    cat_hub = {c: sorted(m, key=lambda nm: (-deg.get(nm, 0), nm))[0]
               for c, m in by_cat_nodes.items()}
    for e in data:
        nm = e["name"]
        if deg.get(nm, 0) != 0:
            continue
        hub = cat_hub[e["category"]]
        if hub == nm:
            continue
        key = (nm, hub) if nm < hub else (hub, nm)
        edge_weight.setdefault(key, ANCHOR_WEIGHT)

    edges = [{"source": a, "target": b, "weight": w}
             for (a, b), w in sorted(edge_weight.items())]

    # Final degrees for node sizing.
    deg = {}
    for e in edges:
        deg[e["source"]] = deg.get(e["source"], 0) + 1
        deg[e["target"]] = deg.get(e["target"], 0) + 1

    colors = category_colors(set(name_to_cat.values()))
    cat_counts = {}
    for c in name_to_cat.values():
        cat_counts[c] = cat_counts.get(c, 0) + 1

    nodes = [{
        "id": e["name"],
        "label": e["name"],
        "category": e["category"],
        "color": colors[e["category"]],
        "url": e["url"],
        "status": e.get("status", ""),
        "summary": e.get("summary", ""),
        "degree": deg.get(e["name"], 0),
    } for e in sorted(data, key=lambda e: e["name"])]

    categories = [{
        "name": c,
        "color": colors[c],
        "count": cat_counts[c],
    } for c in sorted(cat_counts)]

    out = {
        "nodes": nodes,
        "edges": edges,
        "categories": categories,
    }
    OUT_JSON.write_text(json.dumps(out))
    print(f"[graph_data] {len(nodes)} nodes, {len(edges)} edges, "
          f"{len(categories)} categories -> {OUT_JSON.relative_to(HERE)}")


if __name__ == "__main__":
    main()

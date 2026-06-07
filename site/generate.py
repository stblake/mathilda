#!/usr/bin/env python3
"""
Mathilda documentation-site generator.

Produces one MkDocs page per public built-in function under
``site/docs/documentation/<category>/<Name>.md``, plus the documentation-centre
landing page, per-category index pages, awesome-pages ``.pages`` nav files, and
``assets/builtins.json`` for client-side use.

Pipeline (see site/README or the project plan for rationale):

  1. Discover public builtins by parsing ``symtab_set_docstring("Name", ...)``
     C string literals across ``src/**/*.c`` (deterministic; no REPL quoting).
  2. Pull attributes by driving the compiled ``./Mathilda`` binary once.
  3. Map function -> category from the headings of ``docs/spec/builtins/*.md``.
  4. Mine ``In[]:=/Out[]=`` example blocks from each function's spec section and
     RE-VERIFY every example by feeding the inputs back through ``./Mathilda``
     and capturing the real output (the verification gate).
  5. Derive a Stable / Partial / Experimental status heuristically.
  6. Build baseline references (source module + spec page).
  7. Merge optional hand-curated overlays from ``site/overlays/<Name>.md``.
  8. Emit pages, indexes, nav, and the JSON index.

Run from anywhere:  ``python3 site/generate.py``  (needs the built ``./Mathilda``).
The generated Markdown is committed so CI only needs MkDocs, not the C toolchain.
"""

import json
import re
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent          # repo root
SITE = ROOT / "site"
SRC = ROOT / "src"
SPEC_DIR = ROOT / "docs" / "spec" / "builtins"
LIMITATIONS = ROOT / "docs" / "spec" / "limitations.md"
TESTS_DIR = ROOT / "tests"
DOC_OUT = SITE / "docs" / "documentation"
ASSETS = SITE / "docs" / "assets"
OVERLAYS = SITE / "overlays"
IMPL = SITE / "impl"
MATHILDA = ROOT / "Mathilda"

GITHUB_BLOB = "https://github.com/stblake/mathilda/blob/main"

# A built-in name: an uppercase-initial identifier, or a $-system symbol.
NAME_RE = re.compile(r"^\$?[A-Z][A-Za-z0-9]*$")

# Status taxonomy ----------------------------------------------------------
STATUS_BADGE = {
    "Stable": ("success", "Stable",
               "documented, exercised by the test suite and/or worked examples, "
               "with no known limitations recorded."),
    "Partial": ("warning", "Partial",
                "implemented with documented limitations or caveats; some argument "
                "forms fall through to symbolic/unevaluated output."),
    "Experimental": ("note", "Experimental",
                     "present and registered, but lightly documented and not yet "
                     "covered by dedicated tests."),
}
# Strong signals that a function is only partially implemented. We deliberately
# avoid generic phrases like "leave X unevaluated", which describe *correct*
# symbolic passthrough rather than a limitation.
PARTIAL_KEYWORDS = (
    "not implemented", "not yet supported", "not yet implemented",
    "experimental", "approximation only", "currently limited", "todo", "stub",
    "partial support", "is a stub",
)


# ===========================================================================
# 1. Discover builtins + parse docstrings from C source
# ===========================================================================
def _read_c_string_literal(text, i):
    """Parse one C string literal starting at text[i] == '"'. Return (value, end)."""
    assert text[i] == '"'
    i += 1
    out = []
    escapes = {"n": "\n", "t": "\t", "r": "\r", '"': '"', "\\": "\\", "0": "\0"}
    while i < len(text):
        c = text[i]
        if c == "\\":
            nxt = text[i + 1] if i + 1 < len(text) else ""
            out.append(escapes.get(nxt, nxt))
            i += 2
        elif c == '"':
            return "".join(out), i + 1
        else:
            out.append(c)
            i += 1
    return "".join(out), i


def _parse_concatenated_literals(text, i):
    """From index i (just after the comma), read adjacent C string literals
    separated only by whitespace/comments, until the closing ')'."""
    pieces = []
    while i < len(text):
        c = text[i]
        if c.isspace():
            i += 1
        elif c == '"':
            val, i = _read_c_string_literal(text, i)
            pieces.append(val)
        elif text[i:i + 2] == "/*":
            end = text.find("*/", i)
            i = len(text) if end == -1 else end + 2
        elif text[i:i + 2] == "//":
            end = text.find("\n", i)
            i = len(text) if end == -1 else end + 1
        elif c == ")":
            break
        else:
            break
    return "".join(pieces)


DOC_CALL_RE = re.compile(r'symtab_set_docstring\s*\(\s*"((?:[^"\\]|\\.)*)"\s*,')


def discover_builtins():
    """Return {name: {"doc": str, "module": "src/xxx.c"}} for public builtins."""
    found = {}
    for cfile in sorted(SRC.rglob("*.c")):
        if "external" in cfile.parts:
            continue
        text = cfile.read_text(errors="replace")
        for m in DOC_CALL_RE.finditer(text):
            raw_name = m.group(1)
            # Internal context-qualified helpers contain a backtick.
            if "`" in raw_name:
                continue
            name = raw_name.encode().decode("unicode_escape")
            doc = _parse_concatenated_literals(text, m.end())
            rel = cfile.relative_to(ROOT).as_posix()
            # First registration wins (info.c is the canonical hub and is scanned
            # in sorted order; prefer a longer docstring if a later one is richer).
            if name not in found or len(doc) > len(found[name]["doc"]):
                found[name] = {"doc": doc.strip(), "module": rel}
    return found


# ===========================================================================
# 2. Drive the binary
# ===========================================================================
OUT_RE = re.compile(r"^Out\[\d+\]=\s?(.*)$")


def run_session(lines, timeout=60):
    """Feed `lines` to ./Mathilda, return the list of Out[] values in order.
    Each evaluated input yields exactly one Out[] line (verified)."""
    inp = "\n".join(lines) + "\n"
    try:
        proc = subprocess.run([str(MATHILDA)], input=inp, capture_output=True,
                              text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None
    outs = []
    sl = proc.stdout.splitlines()
    i = 0
    while i < len(sl):
        m = OUT_RE.match(sl[i])
        if m:
            buf = [m.group(1)]
            i += 1
            while (i < len(sl) and sl[i].strip() != ""
                   and not sl[i].startswith("Out[") and not sl[i].startswith("In[")):
                buf.append(sl[i])
                i += 1
            outs.append("\n".join(buf).rstrip())
        else:
            i += 1
    return outs


def get_attributes(names):
    """Return {name: [attr, ...]} by querying Attributes[name] for each name."""
    outs = run_session([f"Attributes[{n}]" for n in names], timeout=120)
    attrs = {}
    for name, out in zip(names, outs or []):
        out = (out or "").strip()
        if out.startswith("{") and out.endswith("}"):
            inner = out[1:-1].strip()
            attrs[name] = [a.strip() for a in inner.split(",")] if inner else []
        else:
            attrs[name] = []
    return attrs


# ===========================================================================
# 3. Category mapping + 4. example mining (from docs/spec/builtins/*.md)
# ===========================================================================
def slugify(title):
    return re.sub(r"[^a-z0-9]+", "-", title.lower()).strip("-")


def parse_spec_files():
    """Return (categories, sections).
    categories: ordered list of (slug, title, spec_rel_path).
    sections:   {name: {"category": slug, "body": str}}  (function -> spec section)."""
    categories = []
    sections = {}
    for md in sorted(SPEC_DIR.glob("*.md")):
        text = md.read_text(errors="replace")
        h1 = re.search(r"^#\s+(.+)$", text, re.M)
        title = h1.group(1).strip() if h1 else md.stem.replace("-", " ").title()
        slug = md.stem
        spec_rel = md.relative_to(ROOT).as_posix()
        categories.append((slug, title, spec_rel))

        # Split on H2 headings; each heading may list several names.
        parts = re.split(r"^##\s+(.+)$", text, flags=re.M)
        # parts = [preamble, heading1, body1, heading2, body2, ...]
        for h, body in zip(parts[1::2], parts[2::2]):
            head = re.sub(r"\(.*?\)", "", h).strip()        # drop "(+)", "(parser-level)"
            for token in re.split(r"[,/]", head):
                token = token.strip().strip("`")
                if NAME_RE.match(token) and token not in sections:
                    sections[token] = {"category": slug, "body": body}
    return categories, sections


IN_RE = re.compile(r"^In\[\d+\]:=\s?(.*)$")
OUTLINE_RE = re.compile(r"^Out\[\d+\]=")
FENCE_RE = re.compile(r"^```")


def mine_example_inputs(body, limit=8):
    """Extract ordered lists of input expressions, grouped per code block.
    Returns a list of blocks; each block is a list of input strings."""
    blocks = []
    in_fence = False
    cur = None
    for line in body.splitlines():
        if FENCE_RE.match(line):
            if in_fence:
                if cur:
                    blocks.append(cur)
                cur = None
                in_fence = False
            else:
                in_fence = True
                cur = []
            continue
        if in_fence:
            m = IN_RE.match(line)
            if m and m.group(1).strip():
                cur.append(m.group(1).strip())
    total = sum(len(b) for b in blocks)
    # Trim to `limit` inputs total, keeping whole blocks where possible.
    trimmed, count = [], 0
    for b in blocks:
        if count >= limit:
            break
        room = limit - count
        trimmed.append(b[:room])
        count += min(len(b), room)
    return trimmed


def verify_block(inputs):
    """Run a block's inputs through the binary; return [(in, out), ...] for the
    ones that produced a non-trivial, evaluated result."""
    outs = run_session(inputs)
    if outs is None or len(outs) != len(inputs):
        return []
    pairs = []
    for expr, out in zip(inputs, outs):
        out = out.strip()
        if not out or out == "Null":
            continue
        pairs.append((expr, out))
    return pairs


# ===========================================================================
# 5. Status + 6. references
# ===========================================================================
def load_tests_text():
    text = []
    if TESTS_DIR.exists():
        for t in TESTS_DIR.rglob("*.c"):
            text.append(t.read_text(errors="replace"))
    return "\n".join(text)


def derive_status(name, doc, has_examples, tests_text, lim_text):
    low = doc.lower()
    if re.search(rf"\b{re.escape(name)}\b", lim_text) or any(k in low for k in PARTIAL_KEYWORDS):
        return "Partial"
    tested = (f"{name}[" in tests_text) or (f'"{name}"' in tests_text)
    if tested or has_examples:
        return "Stable"
    return "Experimental"


def build_references(name, module, category_slug, spec_rel):
    refs = [f"Source: [`{module}`]({GITHUB_BLOB}/{module})"]
    if spec_rel:
        refs.append(f"Specification: [`{spec_rel}`]({GITHUB_BLOB}/{spec_rel})")
    else:
        refs.append(f"Specification index: "
                    f"[`Mathilda_spec.md`]({GITHUB_BLOB}/Mathilda_spec.md)")
    return refs


# ===========================================================================
# 7. Overlays + implementation notes
# ===========================================================================
def _parse_front_matter(text):
    """Split a curated Markdown file into ({key: value-or-list}, body).

    Accepts an optional leading ``---\\n … \\n---`` YAML-ish block where each
    line is ``key: value`` and indented ``- item`` lines extend the previous
    key into a list. Values may be single- or double-quoted."""
    meta, body = {}, text
    fm = re.match(r"^---\n(.*?)\n---\n?(.*)$", text, re.S)
    if fm:
        body = fm.group(2)

        def _unquote(s):
            s = s.strip()
            if len(s) >= 2 and s[0] == s[-1] and s[0] in "\"'":
                q = s[0]
                s = s[1:-1].replace("\\" + q, q).replace("\\\\", "\\")
            return s

        cur_key = None
        for line in fm.group(1).splitlines():
            if re.match(r"^\s*-\s+", line) and cur_key:
                meta.setdefault(cur_key, []).append(_unquote(line.split("-", 1)[1]))
            elif ":" in line:
                k, v = line.split(":", 1)
                cur_key = k.strip()
                v = v.strip()
                meta[cur_key] = _unquote(v) if v else []
    return meta, body.strip()


def load_overlay(name):
    """Parse site/overlays/<Name>.md: optional front matter (status:,
    references: list) followed by a Markdown body shown under
    "Notes & additional examples"."""
    f = OVERLAYS / f"{name}.md"
    if not f.exists():
        return None
    meta, body = _parse_front_matter(f.read_text())
    return {"meta": meta, "body": body}


def load_impl(name):
    """Parse site/impl/<Name>.md: the source-grounded implementation summary
    rendered into the "Implementation notes" section. Optional front matter:
    ``references:`` (literature list, merged into References) and ``source:``
    (the real implementation file, overriding the info.c docstring-hub link)."""
    f = IMPL / f"{name}.md"
    if not f.exists():
        return None
    meta, body = _parse_front_matter(f.read_text())
    return {"meta": meta, "body": body}


# ===========================================================================
# 8. Page rendering
# ===========================================================================
def render_docstring(doc):
    """Render the raw docstring verbatim in a fenced block — exactly the text a
    user sees from `?Name` in the REPL. A fenced block is faithful and immune to
    Markdown injection from arbitrary `[`, `]`, `*`, `_` characters in docstrings."""
    body = doc.replace("\t", "    ").rstrip()
    # Guard against the (vanishingly unlikely) docstring containing a fence.
    fence = "```"
    while fence in body:
        fence += "`"
    return f"{fence}text\n{body}\n{fence}"


def render_examples(blocks):
    out = []
    for pairs in blocks:
        if not pairs:
            continue
        lines = ["```mathematica"]
        for i, (expr, res) in enumerate(pairs, 1):
            lines.append(f"In[{i}]:= {expr}")
            lines.append(f"Out[{i}]= {res}")
            lines.append("")
        if lines[-1] == "":
            lines.pop()
        lines.append("```")
        out.append("\n".join(lines))
    return "\n\n".join(out)


def render_impl_notes(name, attrs, section_body, impl_body=None):
    """Assemble the "Implementation notes" section.

    When a hand-authored, source-grounded summary exists in ``site/impl/<Name>.md``
    it leads the section. The ``**Features**`` bullets mined from the spec (which
    describe user-facing capabilities, complementary to the algorithm prose) follow,
    and the attribute list closes it. Without an impl file the behaviour is the
    original Features + Attributes rendering."""
    notes = []
    if impl_body:
        notes.append(impl_body.strip())
    feat = re.search(r"\*\*Features\*\*:?\s*\n((?:[ \t]*-.*\n?)+)", section_body or "")
    if feat:
        notes.append(feat.group(1).rstrip())
    if attrs:
        notes.append(f"**Attributes:** {', '.join(f'`{a}`' for a in attrs)}.")
    else:
        notes.append("**Attributes:** none registered.")
    return "\n\n".join(notes) if notes else "_No additional implementation notes._"


def render_page(info):
    name = info["name"]
    color, label, rationale = STATUS_BADGE[info["status"]]
    lines = [f"# {name}", ""]
    lines.append(f'!!! {color} "Status: {label}"')
    lines.append(f"    {rationale}")
    lines.append("")
    lines.append("## Description")
    lines.append("")
    lines.append(render_docstring(info["doc"]) if info["doc"]
                 else "_No description available._")
    lines.append("")

    lines.append("## Examples")
    lines.append("")
    ex = render_examples(info["examples"])
    if ex:
        lines.append("All examples below are verified against the current "
                     "Mathilda build.")
        lines.append("")
        lines.append(ex)
    else:
        lines.append("_No verified examples yet for this function._")
    lines.append("")

    lines.append("## Implementation notes")
    lines.append("")
    lines.append(info["notes"])
    lines.append("")

    lines.append("## Implementation status")
    lines.append("")
    lines.append(f"**{label}** — {rationale}")
    lines.append("")

    lines.append("## References")
    lines.append("")
    for r in info["references"]:
        lines.append(f"- {r}")
    lines.append("")

    if info.get("overlay_body"):
        lines.append("## Notes & additional examples")
        lines.append("")
        lines.append(info["overlay_body"])
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


# ===========================================================================
# Main
# ===========================================================================
def main():
    if not MATHILDA.exists():
        sys.exit(f"error: {MATHILDA} not found — run `make` first.")

    print("Discovering builtins from source ...")
    builtins = discover_builtins()
    names = sorted(builtins)
    print(f"  {len(names)} public builtins discovered.")

    print("Parsing category specs ...")
    categories, sections = parse_spec_files()
    cat_titles = {slug: title for slug, title, _ in categories}
    cat_spec = {slug: rel for slug, _, rel in categories}

    print("Querying attributes from ./Mathilda ...")
    attrs = get_attributes(names)

    tests_text = load_tests_text()
    lim_text = LIMITATIONS.read_text(errors="replace") if LIMITATIONS.exists() else ""

    OTHER_SLUG = "other-advanced"
    pages = {}
    verified_examples = 0

    print("Building + verifying pages ...")
    for name in names:
        sec = sections.get(name)
        cat = sec["category"] if sec else OTHER_SLUG
        body = sec["body"] if sec else ""
        spec_rel = cat_spec.get(cat, "")

        blocks = []
        for inputs in mine_example_inputs(body):
            pairs = verify_block(inputs)
            if pairs:
                blocks.append(pairs)
                verified_examples += len(pairs)

        status = derive_status(name, builtins[name]["doc"], bool(blocks),
                               tests_text, lim_text)

        # Source-grounded implementation notes (site/impl/<Name>.md). An optional
        # `source:` overrides the docstring-hub module for the Source: reference.
        impl = load_impl(name)
        impl_body = None
        impl_refs = []
        source_module = builtins[name]["module"]
        if impl:
            im = impl["meta"]
            if isinstance(im.get("source"), str) and im["source"].strip():
                source_module = im["source"].strip()
            if isinstance(im.get("references"), list) and im["references"]:
                impl_refs = im["references"]
            impl_body = impl["body"] or None

        refs = build_references(name, source_module, cat, spec_rel)
        if impl_refs:
            refs = impl_refs + refs
        notes = render_impl_notes(name, attrs.get(name, []), body, impl_body)

        overlay = load_overlay(name)
        overlay_body = None
        if overlay:
            m = overlay["meta"]
            if m.get("status") in STATUS_BADGE:
                status = m["status"]
            if isinstance(m.get("references"), list) and m["references"]:
                refs = m["references"] + refs
            overlay_body = overlay["body"] or None

        # Dedupe references, preserving first-seen order (impl, then overlay, then base).
        seen_refs, deduped = set(), []
        for r in refs:
            if r not in seen_refs:
                seen_refs.add(r)
                deduped.append(r)
        refs = deduped

        pages[name] = {
            "name": name, "category": cat, "doc": builtins[name]["doc"],
            "attrs": attrs.get(name, []), "examples": blocks, "status": status,
            "references": refs, "notes": notes, "module": builtins[name]["module"],
            "overlay_body": overlay_body,
        }

    # ---- emit per-builtin pages, grouped by category ----------------------
    # Build ordered category list (spec order, then Other).
    cat_order = [slug for slug, _, _ in categories]
    cat_titles[OTHER_SLUG] = "Other & Advanced"
    if any(p["category"] == OTHER_SLUG for p in pages.values()):
        cat_order.append(OTHER_SLUG)

    by_cat = {slug: [] for slug in cat_order}
    for name in sorted(pages):
        by_cat[pages[name]["category"]].append(name)

    # Wipe previously-generated category dirs (idempotent regeneration).
    if DOC_OUT.exists():
        for child in DOC_OUT.iterdir():
            if child.is_dir():
                for f in child.glob("*.md"):
                    f.unlink()
                pages_file = child / ".pages"
                if pages_file.exists():
                    pages_file.unlink()

    npages = 0
    for slug in cat_order:
        members = by_cat.get(slug, [])
        if not members:
            continue
        cdir = DOC_OUT / slug
        cdir.mkdir(parents=True, exist_ok=True)
        title = cat_titles.get(slug, slug)
        # category index page
        idx = [f"# {title}", "",
               f"{len(members)} built-in function(s) in this category.", ""]
        for name in members:
            st = pages[name]["status"]
            first = pages[name]["doc"].split("\n", 1)[0] if pages[name]["doc"] else ""
            idx.append(f"- [`{name}`]({name}.md) — {first}  _({st})_" if first
                       else f"- [`{name}`]({name}.md)  _({st})_")
        (cdir / "index.md").write_text("\n".join(idx) + "\n")
        # awesome-pages: title + order (index first, then alphabetical)
        (cdir / ".pages").write_text(
            f"title: {title}\nnav:\n  - index.md\n  - ...\n")
        for name in members:
            (cdir / f"{name}.md").write_text(render_page(pages[name]))
            npages += 1

    # ---- documentation centre landing page --------------------------------
    DOC_OUT.mkdir(parents=True, exist_ok=True)
    landing = ["# Documentation Center", "",
               "Every public built-in function in Mathilda, grouped by category. "
               "Each page follows the same shape: **Description** (the function's "
               "docstring), **Examples** (verified against the current build), "
               "**Implementation notes**, **Implementation status**, and "
               "**References**.", "",
               f"_{len(pages)} functions across {len([s for s in cat_order if by_cat.get(s)])} "
               "categories. Use the search box (press `/`) to jump to any function._",
               "", "## Categories", ""]
    for slug in cat_order:
        members = by_cat.get(slug, [])
        if not members:
            continue
        title = cat_titles.get(slug, slug)
        landing.append(f"### [{title}]({slug}/index.md)")
        landing.append("")
        landing.append("  ".join(f"[`{n}`]({slug}/{n}.md)" for n in members))
        landing.append("")
    landing.append("## Alphabetical index")
    landing.append("")
    for name in sorted(pages):
        slug = pages[name]["category"]
        landing.append(f"- [`{name}`]({slug}/{name}.md)")
    (DOC_OUT / "index.md").write_text("\n".join(landing) + "\n")

    # nav ordering for the documentation/ folder
    (DOC_OUT / ".pages").write_text("title: Documentation\nnav:\n  - index.md\n  - ...\n")

    # ---- JSON index -------------------------------------------------------
    ASSETS.mkdir(parents=True, exist_ok=True)
    index = [{
        "name": n, "category": cat_titles.get(pages[n]["category"]),
        "slug": pages[n]["category"], "status": pages[n]["status"],
        "summary": (pages[n]["doc"].split("\n", 1)[0] if pages[n]["doc"] else ""),
        "url": f"documentation/{pages[n]['category']}/{n}/",
    } for n in sorted(pages)]
    (ASSETS / "builtins.json").write_text(json.dumps(index, indent=2))

    # ---- summary ----------------------------------------------------------
    statuses = {}
    for p in pages.values():
        statuses[p["status"]] = statuses.get(p["status"], 0) + 1
    print(f"\nDone. {npages} pages written.")
    print(f"  Verified examples: {verified_examples}")
    print(f"  Status breakdown:  {statuses}")
    overlays = list(OVERLAYS.glob("*.md")) if OVERLAYS.exists() else []
    print(f"  Overlays merged:   {len(overlays)}")
    impls = list(IMPL.glob("*.md")) if IMPL.exists() else []
    print(f"  Impl notes merged: {len(impls)}")
    uncategorized = len(by_cat.get(OTHER_SLUG, []))
    print(f"  Other & Advanced:  {uncategorized}")


if __name__ == "__main__":
    main()

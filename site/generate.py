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

# The ``FLINT` `` context. Its routines are documented in a dedicated section
# rather than the per-``System``` category machinery: they are backtick-qualified
# (so excluded from ``discover_builtins`` and the spec-heading heuristic) and are
# grouped together as the direct-access wrappers over the FLINT-backed kernels.
FLINT_SLUG = "flint"
FLINT_TITLE = "FLINT context"
FLINT_SPEC = SPEC_DIR / "flint.md"

# A built-in name: an uppercase-initial identifier, or a $-system symbol.
# Context-qualified internal helpers (e.g. ``Solve`SolveLinearSystem``) contain a
# backtick and therefore fail this pattern; they are intentionally excluded from
# the public documentation centre.
NAME_RE = re.compile(r"^\$?[A-Z][A-Za-z0-9]*$")

# ---------------------------------------------------------------------------
# Curated category map (authoritative grouping)
# ---------------------------------------------------------------------------
# The spec-file heuristic (parse_spec_files) only assigns a category when a
# function's name literally appears in an ``## H2`` heading of a
# ``docs/spec/builtins/*.md`` file. Functions documented under *thematic*
# headings ("## Trig Functions"), or not yet documented at all, fall through to
# "Other & Advanced". This table is the single, reviewable place that pins such
# functions — and a few mis-filed ones — to the right category. The value is a
# spec-file slug (the ``*.md`` stem). Entries here override the heuristic;
# functions already placed correctly by their spec file need no entry.
CATEGORY_OVERRIDES = {
    # Elementary functions documented under thematic headings -----------------
    "Sin": "elementary-functions", "Cos": "elementary-functions",
    "Tan": "elementary-functions", "Cot": "elementary-functions",
    "Sec": "elementary-functions", "Csc": "elementary-functions",
    "ArcSin": "elementary-functions", "ArcCos": "elementary-functions",
    "ArcTan": "elementary-functions", "Sinh": "elementary-functions",
    "Cosh": "elementary-functions", "Tanh": "elementary-functions",
    "Exp": "elementary-functions", "Log": "elementary-functions",
    # Arithmetic (rounding, basic operators, complex parts, precision) --------
    "Floor": "arithmetic", "Ceiling": "arithmetic", "Round": "arithmetic",
    "Rationalize": "arithmetic", "Complex": "arithmetic", "Divide": "arithmetic",
    "Subtract": "arithmetic", "Factorial2": "arithmetic",
    "Accuracy": "arithmetic", "Precision": "arithmetic",
    "SetAccuracy": "arithmetic", "SetPrecision": "arithmetic",
    # Pattern matching --------------------------------------------------------
    "Blank": "pattern-matching", "BlankSequence": "pattern-matching",
    "BlankNullSequence": "pattern-matching", "Optional": "pattern-matching",
    "Repeated": "pattern-matching", "RepeatedNull": "pattern-matching",
    "Longest": "pattern-matching", "Shortest": "pattern-matching",
    "MatchQ": "pattern-matching", "Default": "pattern-matching",
    "HoldPattern": "pattern-matching",
    # Assignment and rules ----------------------------------------------------
    "Rule": "assignment-and-rules", "RuleDelayed": "assignment-and-rules",
    "Clear": "assignment-and-rules", "AddTo": "assignment-and-rules",
    "SubtractFrom": "assignment-and-rules", "Increment": "assignment-and-rules",
    "Decrement": "assignment-and-rules", "PreIncrement": "assignment-and-rules",
    "PreDecrement": "assignment-and-rules", "DownValues": "assignment-and-rules",
    "OwnValues": "assignment-and-rules", "ReplaceList": "assignment-and-rules",
    "ReplaceRepeated": "assignment-and-rules",
    # Expression information (predicates, attributes, hold/forms) -------------
    "NumericQ": "expression-information", "SetAttributes": "expression-information",
    "TeXForm": "expression-information", "Hold": "expression-information",
    "Flat": "expression-information", "Orderless": "expression-information",
    "OneIdentity": "expression-information",
    # Functional programming --------------------------------------------------
    "Function": "functional-programming", "Slot": "functional-programming",
    "SlotSequence": "functional-programming", "FixedPoint": "functional-programming",
    # Scoping constructs / contexts ------------------------------------------
    "Begin": "scoping-constructs", "BeginPackage": "scoping-constructs",
    "End": "scoping-constructs", "EndPackage": "scoping-constructs",
    "Context": "scoping-constructs", "$Context": "scoping-constructs",
    "$ContextPath": "scoping-constructs",
    # Control flow / session-evaluation hooks --------------------------------
    "$RecursionLimit": "control-flow", "$Pre": "control-flow",
    "$Post": "control-flow", "$PrePrint": "control-flow",
    "$PreRead": "control-flow", "$Epilog": "control-flow",
    # Calculus ----------------------------------------------------------------
    "Limit": "calculus",
    # Statistics --------------------------------------------------------------
    "Quartiles": "statistics", "Commonest": "statistics",
    # Linear algebra ----------------------------------------------------------
    "ConjugateTranspose": "linear-algebra",
    # Structural manipulation -------------------------------------------------
    "UpTo": "structural-manipulation",
    # Solutions of Equations (Other-rescues; the core solver lives in its own
    # spec file solutions-of-equations.md) ------------------------------------
    "RootSum": "solutions-of-equations",
    "GeneratedParameters": "solutions-of-equations",
    "InverseFunctions": "solutions-of-equations",
    "VerifySolutions": "solutions-of-equations",
    "Eliminate": "solutions-of-equations",
    # Algebra (Other-rescues + moved out of Structural Manipulation) ----------
    "Decompose": "algebra",
    "Coefficient": "algebra", "CoefficientList": "algebra", "Collect": "algebra",
    "Discriminant": "algebra", "Expand": "algebra",
    "ExpandDenominator": "algebra", "ExpandNumerator": "algebra",
    "Factor": "algebra", "FactorSquareFree": "algebra", "FactorTerms": "algebra",
    "FactorTermsList": "algebra", "GroebnerBasis": "algebra",
    "HornerForm": "algebra", "PolynomialExtendedGCD": "algebra",
    "PolynomialGCD": "algebra", "PolynomialLCM": "algebra",
    "PolynomialMod": "algebra", "PolynomialQ": "algebra",
    "PolynomialQuotient": "algebra", "PolynomialRemainder": "algebra",
    "PowerExpand": "algebra", "Resultant": "algebra", "Variables": "algebra",
}

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


DOC_CALL_ANY_RE = re.compile(r'symtab_set_docstring\s*\(')
FLINT_NAME_RE = re.compile(r"^FLINT`[A-Za-z0-9]+")


def discover_flint_builtins():
    """Return {full_name: {"doc": str, "module": "src/xxx.c"}} for the
    ``FLINT` `` context routines.

    Unlike the public ``System``` builtins, the FLINT wrappers register their
    docstring against a *symbol constant* first argument
    (``symtab_set_docstring(SYM_FLINT_Zeta, "FLINT`Zeta ...")``) rather than a
    string literal, so ``DOC_CALL_RE`` misses them. Here the routine name is
    recovered from the leading ``FLINT`Name`` token of the docstring text."""
    found = {}
    for cfile in sorted(SRC.rglob("*.c")):
        if "external" in cfile.parts:
            continue
        text = cfile.read_text(errors="replace")
        for m in DOC_CALL_ANY_RE.finditer(text):
            # Skip the first argument (a symbol or string) up to its comma, then
            # parse the concatenated docstring literals that follow.
            comma = text.find(",", m.end())
            if comma == -1:
                continue
            doc = _parse_concatenated_literals(text, comma + 1).strip()
            nm = FLINT_NAME_RE.match(doc)
            if not nm:
                continue
            name = nm.group(0)
            rel = cfile.relative_to(ROOT).as_posix()
            if name not in found or len(doc) > len(found[name]["doc"]):
                found[name] = {"doc": doc, "module": rel}
    return found


def parse_flint_spec():
    """Parse ``docs/spec/builtins/flint.md`` into an ordered list of
    ``(full_name, body)`` pairs, one per ``## FLINT`Name`` H2 section, preserving
    spec order for the nav. Section bodies feed the shared example miner."""
    if not FLINT_SPEC.exists():
        return []
    text = FLINT_SPEC.read_text(errors="replace")
    out = []
    parts = re.split(r"^##\s+(.+)$", text, flags=re.M)
    for h, body in zip(parts[1::2], parts[2::2]):
        head = h.strip().strip("`")
        if head.startswith("FLINT`"):
            out.append((head, body))
    return out


# ===========================================================================
# 2. Drive the binary
# ===========================================================================
OUT_RE = re.compile(r"^Out\[\d+\]=\s?(.*)$")


def run_session(lines, timeout=60):
    """Feed `lines` to ./Mathilda and return the printed output value for each
    input line, in order (empty string when an input produced no value).
    Returns None on timeout.

    When stdin is not a terminal the REPL speaks a line-based NDJSON protocol
    (see src/repl.c): each request ``{"id": N, "expr": "..."}`` yields one or
    more response objects, the value carried by ``{"id": N, "type": "expr",
    "payload": "..."}``. Aligning replies by id is exact even for inputs that
    print nothing, so this is more robust than scraping ``Out[]=`` banners
    (which the pipe protocol no longer emits)."""
    reqs = [json.dumps({"id": i, "expr": line}) for i, line in enumerate(lines, 1)]
    reqs.append(json.dumps({"type": "quit"}))
    inp = "\n".join(reqs) + "\n"
    try:
        proc = subprocess.run([str(MATHILDA)], input=inp, capture_output=True,
                              text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None
    payloads = {}
    for raw in proc.stdout.splitlines():
        raw = raw.strip()
        if not raw or raw[0] != "{":
            continue
        try:
            msg = json.loads(raw)
        except ValueError:
            continue
        if isinstance(msg, dict) and msg.get("type") == "expr" and "payload" in msg:
            payloads[msg.get("id")] = msg["payload"]
    return [payloads.get(i, "") for i in range(1, len(lines) + 1)]


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
        # The FLINT spec is consumed by the dedicated FLINT pass, not the
        # per-System-category machinery (its backtick names have no System page).
        if md.name == FLINT_SPEC.name:
            continue
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
# FLINT context section
# ===========================================================================
def emit_flint_section(tests_text, lim_text):
    """Build and write the dedicated ``FLINT` `` documentation section.

    Discovers the backtick-qualified FLINT routines, queries their attributes,
    verifies the examples mined from ``docs/spec/builtins/flint.md``, and writes
    ``documentation/flint/<Suffix>.md`` (one page per routine) plus the category
    index and ``.pages`` nav. Returns an ordered list of entry dicts
    ``{name, slug, status, summary, url}`` for the landing page and JSON index."""
    flint = discover_flint_builtins()
    if not flint:
        return []

    spec_order = parse_flint_spec()
    bodies = {name: body for name, body in spec_order}
    # Spec order first, then any registered routine missing from the spec.
    names = [n for n, _ in spec_order if n in flint]
    names += [n for n in sorted(flint) if n not in bodies]

    attrs = get_attributes(names)
    spec_rel = FLINT_SPEC.relative_to(ROOT).as_posix()

    cdir = DOC_OUT / FLINT_SLUG
    cdir.mkdir(parents=True, exist_ok=True)

    entries = []
    verified = 0
    for name in names:
        slug = name.split("`", 1)[1]            # filesystem/URL-safe suffix
        body = bodies.get(name, "")
        blocks = []
        for inputs in mine_example_inputs(body):
            pairs = verify_block(inputs)
            if pairs:
                blocks.append(pairs)
                verified += len(pairs)
        status = derive_status(name, flint[name]["doc"], bool(blocks),
                               tests_text, lim_text)
        refs = build_references(name, flint[name]["module"], FLINT_SLUG, spec_rel)
        notes = render_impl_notes(name, attrs.get(name, []), body)
        info = {
            "name": name, "doc": flint[name]["doc"], "attrs": attrs.get(name, []),
            "examples": blocks, "status": status, "references": refs,
            "notes": notes, "overlay_body": None,
        }
        (cdir / f"{slug}.md").write_text(render_page(info))
        summary = flint[name]["doc"].split("\n", 1)[0] if flint[name]["doc"] else ""
        entries.append({"name": name, "slug": slug, "status": status,
                        "summary": summary,
                        "url": f"documentation/{FLINT_SLUG}/{slug}/"})

    # category index + awesome-pages nav
    idx = [f"# {FLINT_TITLE}", "",
           f"{len(entries)} routine(s) in the ``FLINT` `` context — direct "
           "access to the FLINT-backed kernels.", "",
           "## Building with FLINT", "",
           "These routines require Mathilda to be compiled against "
           "[FLINT](https://flintlib.org/) **>= 3.0** (the release that merged "
           "ANTIC for number-field arithmetic). FLINT is **optional** and "
           "**auto-detected**: the makefile enables it (`USE_FLINT=1`, the "
           "default) whenever `pkg-config` reports `flint >= 3.0`, and otherwise "
           "prints a warning and falls back to `USE_FLINT=0` — the classical, "
           "still-rigorous algebraic-extension and numeric paths.", "",
           "```bash",
           "# Install FLINT >= 3.0",
           "brew install flint                    # macOS (Homebrew)",
           "sudo apt install libflint-dev         # Ubuntu 24.04+/Debian Bookworm+",
           "",
           "# Build — FLINT is picked up automatically",
           "make -j",
           "",
           "# Force it off (classical fallback):",
           "make -j USE_FLINT=0",
           "```", "",
           "Confirm the installed version with `pkg-config --modversion flint`. "
           "When FLINT is unavailable these `` `FLINT` `` routines are not "
           "registered, and the public builtins that delegate to them (`Factor`, "
           "`PolynomialGCD`, `Cancel`, `Together`, `Zeta`, …) transparently use "
           "the classical implementations.", ""]
    for e in entries:
        first = e["summary"]
        idx.append(f"- [`{e['name']}`]({e['slug']}.md) — {first}  _({e['status']})_"
                   if first else f"- [`{e['name']}`]({e['slug']}.md)  _({e['status']})_")
    (cdir / "index.md").write_text("\n".join(idx) + "\n")
    (cdir / ".pages").write_text(
        f"title: {FLINT_TITLE}\nnav:\n  - index.md\n  - ...\n")

    print(f"  FLINT section: {len(entries)} routines, {verified} verified examples.")
    return entries


# ===========================================================================
# Main
# ===========================================================================
def main():
    if not MATHILDA.exists():
        sys.exit(f"error: {MATHILDA} not found — run `make` first.")

    print("Discovering builtins from source ...")
    builtins = discover_builtins()
    # Drop context-qualified internal helpers (``Head`Helper``) — they carry a
    # docstring for ?-lookup but are not part of the public documentation centre.
    names = sorted(n for n in builtins if NAME_RE.match(n))
    dropped = len(builtins) - len(names)
    print(f"  {len(names)} public builtins discovered "
          f"({dropped} internal helper(s) excluded).")

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
        doc_slug = sec["category"] if sec else None          # where it's documented
        cat = CATEGORY_OVERRIDES.get(name) or doc_slug or OTHER_SLUG
        body = sec["body"] if sec else ""
        # The "Specification" reference points at the file that actually
        # documents the function, which may differ from its (overridden)
        # display category.
        spec_rel = cat_spec.get(doc_slug or cat, "")

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

    # ---- FLINT` context section (its own generation path) -----------------
    flint_entries = emit_flint_section(tests_text, lim_text)
    npages += len(flint_entries)
    n_cats = len([s for s in cat_order if by_cat.get(s)]) + (1 if flint_entries else 0)

    # ---- documentation centre landing page --------------------------------
    DOC_OUT.mkdir(parents=True, exist_ok=True)
    landing = ["# Documentation Center", "",
               "Every public built-in function in Mathilda, grouped by category. "
               "Each page follows the same shape: **Description** (the function's "
               "docstring), **Examples** (verified against the current build), "
               "**Implementation notes**, **Implementation status**, and "
               "**References**.", "",
               f"_{len(pages) + len(flint_entries)} functions across {n_cats} "
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
    if flint_entries:
        landing.append(f"### [{FLINT_TITLE}]({FLINT_SLUG}/index.md)")
        landing.append("")
        landing.append("  ".join(
            f"[`{e['name']}`]({FLINT_SLUG}/{e['slug']}.md)" for e in flint_entries))
        landing.append("")
    landing.append("## Alphabetical index")
    landing.append("")
    alpha = [(name, pages[name]["category"], name) for name in pages]
    alpha += [(e["name"], FLINT_SLUG, e["slug"]) for e in flint_entries]
    for name, slug, stem in sorted(alpha):
        landing.append(f"- [`{name}`]({slug}/{stem}.md)")
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
    index += [{
        "name": e["name"], "category": FLINT_TITLE, "slug": FLINT_SLUG,
        "status": e["status"], "summary": e["summary"], "url": e["url"],
    } for e in flint_entries]
    index.sort(key=lambda d: d["name"])
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

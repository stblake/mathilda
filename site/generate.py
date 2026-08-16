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
import collections
import os
import base64
import io
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
PLOTS_OUT = SITE / "docs" / "assets" / "plots"
DOC_OUT = SITE / "docs" / "documentation"
ASSETS = SITE / "docs" / "assets"
OVERLAYS = SITE / "overlays"
IMPL = SITE / "impl"
MATHILDA = ROOT / "Mathilda"

# Every documented example is re-verified by running it, and some of them call Plot or Manipulate --
# which opened a Raylib window each, dozens of them, over whatever the user was doing. The binary
# honours this by evaluating and printing as usual while withholding the window.
os.environ.setdefault("MATHILDA_NO_WINDOW", "1")

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
    payloads, plots = {}, {}
    for raw in proc.stdout.splitlines():
        raw = raw.strip()
        if not raw or raw[0] != "{":
            continue
        try:
            msg = json.loads(raw)
        except ValueError:
            continue
        if not isinstance(msg, dict):
            continue
        if msg.get("type") == "expr" and "payload" in msg:
            payloads[msg.get("id")] = msg["payload"]
        elif msg.get("type") == "image":
            # An image result arrives as its own message type, not as an expr. Collecting only
            # "expr" recorded nothing for it, so every image-valued example was discarded as
            # unevaluated -- which is why the Image page showed one example about dimensions and
            # none showing an actual picture. "-Image-" stands in for the printed form, and the
            # RGBA payload is kept so the page can SHOW the result.
            payloads.setdefault(msg.get("id"), "-Image-")
            pay = msg.get("payload")
            if isinstance(pay, dict):
                plots.setdefault(msg.get("id"), {"__image__": pay})
        elif msg.get("type") == "plot":
            # A graphics result arrives as a plot payload, not an expr, so
            # collecting only "expr" recorded nothing for it and every graphics
            # example was then discarded as unevaluated. Plot's page was left
            # with its one symbolic form, Plot[Sin[x], {x, a, b}], which is a
            # signature rather than a worked example. "-Graphics-" is what the
            # REPL prints for the same result.
            payloads.setdefault(msg.get("id"), "-Graphics-")
            # Keep the figure itself, not just the fact that there was one, so a
            # page can show the plot without the reader running anything.
            plots.setdefault(msg.get("id"), msg.get("payload"))
    return ([payloads.get(i, "") for i in range(1, len(lines) + 1)],
            {i: plots[i] for i in plots})


def get_attributes(names):
    """Return {name: [attr, ...]} by querying Attributes[name] for each name."""
    outs, _ = run_session([f"Attributes[{n}]" for n in names], timeout=120)
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


# Relative links mined from docs/spec/builtins/<slug>.md assume those files are
# siblings (e.g. `[packed list](packed-arrays.md)`), which holds in the spec tree
# but not on the site, where each category is a folder of per-builtin pages. Left
# unrewritten they dangle and abort `mkdocs build --strict`, blocking the whole
# deploy. Match only bare, same-directory `<slug>.md[#anchor]` targets — absolute
# URLs and already-relative (`../`) links are excluded by the lookahead.
_SPEC_LINK_RE = re.compile(
    r"\]\((?!https?://|/|\.{1,2}/)([a-z0-9][a-z0-9-]*)\.md(#[A-Za-z0-9._-]+)?\)")


def rewrite_spec_links(text, cat_slugs, anchor_to_cat):
    """Rewrite spec-relative category links into site-relative links.

    Every generated page lives at ``documentation/<category>/<Name>.md``, so a
    category index is ``../<slug>/index.md`` and a specific builtin's page is
    ``../<its-category>/<Name>.md``. An anchor (``#mapat``) that names a builtin
    is resolved to that builtin's page — the spec keeps a whole category in one
    file with per-function headings, whereas the site splits them out — otherwise
    it degrades to the category index. Targets that are not a known category file
    (``cat_slugs``) are left untouched."""
    def repl(m):
        slug, anchor = m.group(1), m.group(2)
        if slug not in cat_slugs:
            return m.group(0)
        if anchor:
            hit = anchor_to_cat.get(anchor[1:].lower())
            if hit:
                cat, disp = hit
                return f"](../{cat}/{disp}.md)"
        return f"](../{slug}/index.md)"
    return _SPEC_LINK_RE.sub(repl, text)


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
            # A heading naming several functions is a curated grouping; record
            # the whole set so each member can point at the others.
            group = [t.strip().strip("`") for t in re.split(r"[,/]", head)]
            group = [t for t in group if NAME_RE.match(t)]
            for token in group:
                if token not in sections:
                    sections[token] = {"category": slug, "body": body,
                                       "siblings": group}
    return categories, sections


IN_RE = re.compile(r"^In\[\d+\]:=\s?(.*)$")
OUTLINE_RE = re.compile(r"^Out\[\d+\]=")
FENCE_RE = re.compile(r"^```")


EXAMPLE_LIMIT = 60


GROUP_RE = re.compile(r"^#{3,5}\s+(.+?)\s*$")

# The subsection names a page may declare, in the order a reader wants them: the simple cases,
# then the range of inputs, then what the options do, then a real use, then the algebra, then the
# one that is just nice to look at. `Options: <Name>` declares a sub-block of Options, which is
# how one option with eight examples stays legible.
GROUP_ORDER = ["Basic Examples", "Scope", "Options", "Applications",
               "Properties & Relations", "Neat Examples"]


def mine_example_inputs(body, limit=EXAMPLE_LIMIT):
    """Extract ordered lists of input expressions, grouped per code block.
    Returns a list of ``(group, inputs)`` -- the group is the nearest preceding
    ``###``-``#####`` heading when it names one of ``GROUP_ORDER`` (or an
    ``Options: <Name>`` sub-block), and ``None`` otherwise.

    DECLARED, NOT INFERRED. The grouping used to be guessed from position and
    syntax -- first block is basic, anything containing ``Name ->`` is an
    option -- which meant a page could not say what its own examples were for.
    An identity like ``Dilation[img, 0] === img`` is a *relation*, not a scope
    case, and no amount of pattern-matching on the input text can know that.
    Pages that declare nothing keep the old inference, so this is additive.

    The limit was 8, which silently capped every page: a section could carry forty examples and the
    page would show the first eight, so writing more had no effect and the omission was invisible --
    there is no message, the extra examples simply are not there. 60 is above what any section
    currently holds, so the cap no longer decides what a reader sees.""" 
    blocks = []
    in_fence = False
    cur = None
    group = None
    for line in body.splitlines():
        if not in_fence:
            mh = GROUP_RE.match(line)
            if mh:
                title = mh.group(1).strip()
                base = title.split(":", 1)[0].strip()
                group = title if base in GROUP_ORDER else None
        if FENCE_RE.match(line):
            if in_fence:
                if cur:
                    blocks.append((group, cur))
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
            elif (cur and line.strip() and not OUTLINE_RE.match(line)
                  and not line.startswith("(*")):
                # A continuation of the previous input. Without this an example
                # wrapped across two lines was truncated at the break and then
                # fed to the binary as a syntax error, so it silently vanished
                # from the page.
                cur[-1] = cur[-1] + " " + line.strip()
    # Trim to `limit` inputs total, keeping whole blocks where possible.
    trimmed, count = [], 0
    for grp, b in blocks:
        if count >= limit:
            break
        room = limit - count
        trimmed.append((grp, b[:room]))
        count += min(len(b), room)
    return trimmed


# A saved figure is a convenience, not the record -- the example is runnable
# either way -- so it is not worth shipping megabytes for one page. Surface and
# field plots (Plot3D, StreamPlot) carry tens of thousands of samples and blow
# past this; line plots land around 10 KB each.
FIGURE_MAX_BYTES = 150_000


def compact_figure(payload):
    """Round coordinates and drop the figure if it is still too large.

    An image payload passes through untouched: it is base64 pixels, not coordinates, and rounding
    would corrupt it.

    Six significant digits is far beyond what a few hundred pixels can show,
    and trimming the tail of every float is most of the saving available
    without resampling the curve, which would misrepresent it."""
    def shrink(v):
        if isinstance(v, float):
            return round(v, 6)
        if isinstance(v, list):
            return [shrink(x) for x in v]
        if isinstance(v, dict):
            return {k: shrink(x) for k, x in v.items()}
        return v

    small = shrink(payload)
    return small if len(json.dumps(small)) <= FIGURE_MAX_BYTES else None


# A worked example written as prose: `Expr[...]` -> `result`, usually in a
# semicolon-separated run under a "Worked examples" lead-in. 108 of these exist
# across 19 spec files.
PROSE_EX_RE = re.compile(
    r"`([A-Z][A-Za-z0-9]*\[[^`]{3,200}?)`\s*(?:\u2192|->)\s*`[^`\n]{1,120}`")


def mine_prose_examples(body, limit=10):
    """Input expressions from prose worked examples.

    Only the INPUT is taken. The written result is prose -- `\u221a\u03c0/2`,
    `\u0393[s]`, `\u03c0 a/2` -- not Mathilda syntax, so it cannot be trusted as
    an expected value and is not shown. Running the inputs turns a claim in prose
    into an example the reader can execute, with whatever this build actually
    produces."""
    seen, out = set(), []
    for expr in PROSE_EX_RE.findall(body or ""):
        e = " ".join(expr.split())
        if e not in seen:
            seen.add(e)
            out.append(e)
        if len(out) >= limit:
            break
    return out


def strip_prose_examples(text):
    """Remove the prose runs whose inputs were promoted to real examples, so the
    page does not state them twice."""
    if not text:
        return text
    out = PROSE_EX_RE.sub("", text)
    out = re.sub(r"(?im)^\s*Worked examples[^\n]*:\s*$", "", out)
    out = re.sub(r"[ \t]*;[ \t]*(?=\n)", "", out)
    out = re.sub(r"\n{3,}", "\n\n", out)
    return out.strip() or None


def verify_block(inputs):
    """Run a block's inputs; return ([(in, out), ...], {input_expr: figure}).

    The figure map carries the Plotly payload for any input that produced
    graphics, keyed by the input text so it survives the later trimming and
    renumbering of examples."""
    result = run_session(inputs)
    if result is None:
        return [], {}
    outs, plots = result
    if len(outs) != len(inputs):
        return [], {}
    pairs, figures = [], {}
    for idx, (expr, out) in enumerate(zip(inputs, outs), 1):
        out = out.strip()
        if not out:
            continue
        if out == "Null":
            # A SETUP LINE, not a dead example. `v = Image3D[...];` yields Null, and dropping it left
            # every later example in the block referring to a variable the page never defines -- the
            # Image3D page showed `ImageDimensions[v]` and `Part[ImageData[v], 2, 3, 4]` with `v`
            # undefined, so both printed unevaluated. It is kept with no Out line, which is how a
            # notebook shows a cell that assigns and prints nothing.
            pairs.append((expr, ""))
            continue
        pairs.append((expr, out))
        if idx in plots and plots[idx] is not None:
            fig = compact_figure(plots[idx])
            if fig is not None:
                figures[expr] = fig
    return pairs, figures


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
SIG_RE = re.compile(r"^[A-Z$][A-Za-z0-9$]*(\[|\s*$)")


def render_docstring(doc):
    """Render the docstring as a signature list, with trailing notes collapsed.

    A fenced block was faithful to `?Name` but wrong for a page: a fence never
    wraps, so a docstring written as one long sentence ran off the right edge.

    Docstrings here follow "signature, indented explanation", but many then add
    general prose -- "Sin is Listable. Numeric inputs are evaluated via libm..."
    Treating that as another signature would set a paragraph in code style, so
    a line only counts as a signature if it actually looks like a call. The rest
    becomes notes, in a collapsed block under the signatures, where it is
    available without pushing the examples off the screen.

    Signatures stay in code style so `[`, `*` and `_` inside them cannot be read
    as Markdown; prose is escaped for the same reason."""
    body = doc.replace("\t", "    ").rstrip()
    if not body:
        return "_No description available._"

    def esc(text):
        return re.sub(r"([*_`\[\]<>])", r"\\\1", text)

    forms, notes = [], []
    sig, expl = None, []

    def flush():
        if sig is not None:
            forms.append((sig, " ".join(expl)))

    for raw in body.split("\n"):
        line = raw.rstrip()
        if not line.strip():
            continue
        if line.startswith("    "):
            (expl if sig is not None else notes).append(line.strip())
        elif SIG_RE.match(line.strip()):
            flush()
            sig, expl = line.strip(), []
        else:
            # Prose, not a call: everything from here is a note.
            flush()
            sig, expl = None, []
            notes.append(line.strip())
    flush()

    out = []
    for form, text in forms:
        out.append(f"**`{form}`**")
        out.append("")
        if text:
            out.extend([esc(text), ""])
    if not forms:
        out.append(esc(" ".join(notes)))
        notes = []
    if notes:
        out.extend(["<details>", "<summary>Notes</summary>", "",
                    esc(" ".join(notes)), "", "</details>", ""])
    return "\n".join(out).strip()


OPTION_RE = re.compile(r"\b[A-Z][A-Za-z0-9]*\s*->")


# A one-line explanation attached to an example, written as a Mathilda comment
# in the spec or overlay: `Expand[(x+1)^3]  (* binomial coefficients *)`.
EX_NOTE_RE = re.compile(r"\s*\(\*(.+?)\*\)\s*")


def _example_fence(pairs, start=1, figures=None):
    """Render examples, each preceded by its explanation when it has one.

    An example carrying a `(* ... *)` comment gets that line as prose above its
    own fenced block, so the notebook shows a sentence and then the input it
    describes. Examples without one are grouped into a single fence as before --
    the alternative would be inventing an explanation, and a documentation page
    is the last place to put made-up prose.

    The comment is stripped from the input it annotates: it has already been
    said above the cell, and repeating it inside the expression the reader is
    about to run is noise."""
    out, n = [], start
    run = []

    def flush_run():
        if not run:
            return
        block = ["```mathematica"]
        for expr, res, _ in run:
            block.append(expr)
            if res:
                block.append(res)
            block.append("")
        if block[-1] == "":
            block.pop()
        block.append("```")
        out.append("\n".join(block))
        run.clear()

    for expr, res in pairs:
        note = EX_NOTE_RE.search(expr)
        clean = EX_NOTE_RE.sub(" ", expr).strip() if note else expr
        inp = f"In[{n}]:= {clean}"
        outp = f"Out[{n}]= {res}" if res else None
        n += 1
        pic = image_data_uri((figures or {}).get(expr))
        if note or pic:
            # An image result gets its own fence so the picture sits directly under the input that
            # produced it. Grouped into a shared fence it would follow some later, unrelated example.
            flush_run()
            if note:
                text = note.group(1).strip()
                out.append(text[0].upper() + text[1:] if text else text)
            out.append("```mathematica\n" + inp
                       + (("\n" + outp) if outp else "") + "\n```")
            if pic:
                out.append(pic)
        else:
            run.append((inp, outp, None))
    flush_run()
    return "\n\n".join(out), n



def image_data_uri(fig):
    """An image figure as an inline `![](data:image/png;base64,...)`, or None.

    Inline rather than a file on disk: an example image is a few kilobytes, and a data URI needs no
    change to the page-sync script, no path resolution inside the notebook app, and cannot go stale
    against the page that references it.

    Small results are scaled up with NEAREST neighbour so a 2x2 example reads as four squares rather
    than a four-pixel speck -- the same reason the notebook renders its canvases pixelated.
    """
    if not isinstance(fig, dict) or "__image__" not in fig:
        return None
    pay = fig["__image__"]
    try:
        from PIL import Image as PILImage
        w, h = int(pay["w"]), int(pay["h"])
        raw = base64.b64decode(pay["data"])
        if w <= 0 or h <= 0 or len(raw) < w * h * 4:
            return None
        im = PILImage.frombytes("RGBA", (w, h), raw[: w * h * 4])
        scale = max(1, min(24, 96 // max(w, h))) if max(w, h) < 96 else 1
        if scale > 1:
            im = im.resize((w * scale, h * scale), PILImage.NEAREST)
        buf = io.BytesIO()
        im.save(buf, format="PNG", optimize=True)
        b64 = base64.b64encode(buf.getvalue()).decode("ascii")
    except Exception:
        return None
    note = ""
    if pay.get("depth"):
        note = " (slice %s of %s)" % (pay.get("slice"), pay["depth"])
    return "![%sx%s result%s](data:image/png;base64,%s)" % (pay["w"], pay["h"], note, b64)

def render_examples(blocks, applications=None, worked=None, figures=None):
    """Group the verified examples into subsections, each with its count.

    Two modes. When the spec DECLARES its groups (a `#### Scope` heading above a
    fenced block, say), those declarations are honoured verbatim and in the
    canonical reader order, and an `Options: <Name>` label becomes a nested
    block so one option carrying eight examples stays legible. When it declares
    nothing -- 700-odd pages -- the old inference still applies: an example that
    passes an option is an options example, the first block is the basic set,
    later blocks are scope.

    Inference cannot tell a relation from a scope case: `Dilation[img, 0] ===
    img` is an algebraic identity and reads as neither. That is the whole reason
    a page is allowed to say so itself.

    A flat list tells the reader nothing about which examples are the simple
    cases and which demonstrate an option; the counts make the shape of the page
    visible before scrolling it."""
    declared = [b for b in blocks if b[0]]
    groups = []                      # [(title, [pairs...], depth)]

    if declared:
        buckets, order = {}, []
        for label, pairs in blocks:
            key = label or "Scope"   # an untagged block among tagged ones is scope
            if key not in buckets:
                buckets[key] = []
                order.append(key)
            buckets[key].extend(pairs)

        def rank(key):
            base = key.split(":", 1)[0].strip()
            return (GROUP_ORDER.index(base) if base in GROUP_ORDER else len(GROUP_ORDER),
                    order.index(key))

        open_parent = None       # the `Options` heading currently accepting children
        for key in sorted(order, key=rank):
            base, _, sub = key.partition(":")
            base, sub = base.strip(), sub.strip()
            if sub:
                # Options with sub-blocks: the parent heading is emitted once and
                # carries the total of its children, so `Options (21)` above
                # `ColorSpace (8)` and `ImageSize (3)` adds up.
                if open_parent != base:
                    groups.append((base, [], 0))
                    open_parent = base
                groups.append((sub, buckets[key], 1))
            else:
                groups.append((base, buckets[key], 0))
                open_parent = None
    else:
        basic, scope, options = [], [], []
        for idx, (_, pairs) in enumerate(blocks):
            for pair in pairs:
                if OPTION_RE.search(pair[0]):
                    options.append(pair)
                elif idx == 0:
                    basic.append(pair)
                else:
                    scope.append(pair)
        groups = [("Basic examples", basic, 0), ("Scope", scope, 0),
                  ("Options", options, 0)]

    if worked:
        groups.append(("Worked examples", worked, 0))
    total = sum(len(g[1]) for g in groups)
    app_pairs = applications[1] if applications and applications[0] == "pairs" else None
    if applications:
        total += len(app_pairs) if app_pairs is not None else applications[1]

    if not total:
        return None, 0

    # A parent heading counts the children beneath it.
    counts = []
    for i, (title, pairs, depth) in enumerate(groups):
        n = len(pairs)
        if depth == 0 and not pairs:
            for t2, p2, d2 in groups[i + 1:]:
                if d2 == 0:
                    break
                n += len(p2)
        counts.append(n)

    out, n = [], 1
    for (title, pairs, depth), cnt in zip(groups, counts):
        if not cnt:
            continue
        hashes = "###" if depth == 0 else "####"
        out.append(f"{hashes} {title} ({cnt})")
        out.append("")
        if pairs:
            fence, n = _example_fence(pairs, n, figures)
            out.append(fence)
            out.append("")
    if applications:
        if app_pairs is not None:
            fence, n = _example_fence(app_pairs, n, figures)
            out.append(f"### Applications ({len(app_pairs)})")
            out.append("")
            out.append(fence)
        else:
            body, count = applications
            out.append(f"### Applications ({count})")
            out.append("")
            out.append(body)
        out.append("")
    return "\n".join(out).strip(), total


def render_impl_notes(name, attrs, section_body, impl_body=None,
                      cat_slugs=None, anchor_to_cat=None):
    """Assemble the "Implementation notes" section.

    When a hand-authored, source-grounded summary exists in ``site/impl/<Name>.md``
    it leads the section. The ``**Features**`` bullets mined from the spec (which
    describe user-facing capabilities, complementary to the algorithm prose) follow,
    and the attribute list closes it. Without an impl file the behaviour is the
    original Features + Attributes rendering.

    ``cat_slugs``/``anchor_to_cat`` (when supplied) rewrite the spec-relative
    category links in the mined bullets to site-relative ones — see
    ``rewrite_spec_links``. The mined prose is the only spec-authored content that
    reaches a page as live Markdown (docstrings are inert fenced text), so it is
    the only place those links need fixing."""
    notes = []
    if impl_body:
        notes.append(impl_body.strip())
    bullets = mine_features(section_body)
    if bullets:
        if cat_slugs is not None:
            bullets = rewrite_spec_links(bullets, cat_slugs, anchor_to_cat or {})
        notes.append(bullets)
    if attrs:
        notes.append(f"**Attributes:** {', '.join(f'`{a}`' for a in attrs)}.")
    else:
        notes.append("**Attributes:** none registered.")
    return "\n\n".join(notes) if notes else "_No additional implementation notes._"


def mine_features(section_body):
    """The ``**Features**`` bullet block, INCLUDING wrapped continuation lines.

    The spec files are hard-wrapped near 80 columns, so a bullet routinely spills
    onto an indented continuation line that does not itself start with ``-``. A
    regex of the form ``(?:[ \t]*-.*\n?)+`` stops dead at the first such line,
    which silently truncated 168 of the 305 sections that have a Features block
    -- 3351 bullet lines in total, and often mid-sentence. ``TrigFactor`` kept 2
    lines of 42.

    Scanning line by line instead: a bullet, an indented continuation, or a blank
    line continues the block; anything flush-left that is not a bullet (the next
    ``**Bold**`` block, prose, a table) ends it."""
    if not section_body:
        return None
    fi = section_body.find("**Features**")
    if fi < 0:
        return None
    nl = section_body.find("\n", fi)
    if nl < 0:
        return None
    out = []
    for line in section_body[nl + 1:].split("\n"):
        if re.match(r"^[ \t]*-", line) or re.match(r"^[ \t]+\S", line) or line.strip() == "":
            out.append(line)
        else:
            break
    block = "\n".join(out).strip("\n").rstrip()
    return block if block.lstrip().startswith("-") else None


_FENCE_RE = re.compile(r"^```.*?^```[ \t]*$", re.M | re.S)


def mine_spec_detail(section_body):
    """Spec prose that follows the Features block -- option tables, algorithm
    notes, diagnostics -- which never reached a page at all.

    For FindClusters this is the entire Method compatibility matrix, the
    suboptions, and the DistanceFunction / CriterionFunction / PerformanceGoal
    documentation: 53 of the 622 sections carry such a block.

    Fenced blocks are dropped because they are the *source* the example miner
    reads, and those examples are already rendered under "Examples" after being
    re-verified against the binary. Keeping them here would print each twice, and
    the un-verified copy is the worse of the two.

    Only content after an existing Features block is taken. Without that anchor
    the remainder would start at the top of the section and duplicate the
    signature bullets that "Description" already shows."""
    feat = mine_features(section_body)
    if not feat:
        return None
    idx = section_body.find(feat)
    rest = section_body[idx + len(feat):]
    rest = _FENCE_RE.sub("", rest)
    rest = _promote_labels(rest)
    # A spec section nests four deep; the reference page has only two levels
    # below the page title, so H4 becomes H3 rather than rendering as literal
    # '#### text' inside a paragraph.
    rest = re.sub(r"^#{4,}\s+", "### ", rest, flags=re.M)
    rest = _drop_empty_headings(rest)
    rest = re.sub(r"\n{3,}", "\n\n", rest).strip()
    return rest or None


def _drop_empty_headings(text):
    """Remove headings left with nothing under them.

    The spec marks its example blocks with an '### Examples' heading, and those
    blocks are stripped here because the examples are mined, re-verified and
    shown in the Examples section instead. That left the heading behind: a page
    could carry three or four 'Examples' subsections in a row, each empty and
    each expanding to nothing."""
    lines = text.split("\n")
    keep = [True] * len(lines)
    for i, line in enumerate(lines):
        m = re.match(r"^(#{2,6})\s", line)
        if not m:
            continue
        level = len(m.group(1))
        has_content = False
        for j in range(i + 1, len(lines)):
            nxt = re.match(r"^(#{2,6})\s", lines[j])
            if nxt and len(nxt.group(1)) <= level:
                break
            if lines[j].strip():
                has_content = True
                break
        keep[i] = has_content
    return "\n".join(l for l, k in zip(lines, keep) if k)


def _promote_labels(text):
    """Turn the spec's ``**Method** — ...`` labelled blocks into ``### Method``.

    The spec marks a labelled block with a bold lead-in followed by an em dash or
    a colon. Rendered as-is these are bold paragraphs, so a long page reads as one
    undifferentiated run and the option documentation is invisible in the page
    structure.

    The em-dash/colon requirement is what makes this safe: bold used for emphasis
    mid-prose (``**infinite**``, ``**cliff**``) can also start a line, and
    promoting that would invent a heading out of a sentence fragment."""
    def repl(m):
        label, sep, tail = m.group(1), m.group(2), m.group(3)
        return f"### {label}\n\n{tail.lstrip()}" if tail.strip() else f"### {label}\n"
    return re.sub(r"^\*\*([A-Z][\w /`\"'-]{0,48}?)\*\*\s*(—|:)[ \t]*(.*)$",
                  repl, text, flags=re.M)


# ---------------------------------------------------------------------------
# Enrichment: algorithm, related symbols, tests, measured performance
# ---------------------------------------------------------------------------

PERF_CACHE = SITE / "perf.json"
BENCH_DIR  = ROOT / "benchmarks"
EXP_DIRS   = [BENCH_DIR, ROOT / "docs" / "experiments"]
TESTS_DIR  = ROOT / "tests"


def comment_to_markdown(body):
    """Turn a C block comment into Markdown, preserving what is load-bearing.

    These comments are hand-formatted and their layout carries meaning: aligned
    two-column tables, worked examples with `->`, ASCII diagrams. Reflowing
    everything into paragraphs destroys exactly the parts worth reading -- an
    early version merged two separate FindClusters examples onto one line and
    flattened a table of semantics into a sentence.

    So a line is treated as preformatted when it is indented, when it contains
    run-together spacing used for column alignment, or when it shows a result
    with `->`. Short ALL-CAPS lines are the section headings these comments use
    and become real headings. Everything else is ordinary prose and is
    reflowed."""
    lines = [re.sub(r"^\s*\* ?", "", ln) for ln in body.split("\n")]
    while lines and not lines[0].strip():
        lines.pop(0)
    while lines and not lines[-1].strip():
        lines.pop()
    if not lines:
        return None

    def preformatted(ln):
        if re.match(r"^(\s{2,}|\t)", ln):
            return True
        stripped = ln.strip()
        if "  " in stripped:            # column alignment
            return True
        if "->" in stripped and stripped.startswith(("FindClusters", "In[", "{")):
            return True
        return bool(re.match(r"^[|+`]", stripped))

    def is_heading(ln):
        t = ln.strip()
        return (2 <= len(t) <= 48 and t == t.upper()
                and re.match(r"^[A-Z][A-Z0-9 ,&/()'-]*$", t) is not None
                and any(c.isalpha() for c in t))

    out, para, code = [], [], []

    def flush_para():
        if para:
            out.extend([" ".join(x.strip() for x in para), ""])
            para.clear()

    def flush_code():
        if code:
            while code and not code[-1].strip():
                code.pop()
            if code:
                out.extend(["```text", *code, "```", ""])
            code.clear()

    for ln in lines:
        if not ln.strip():
            if code:
                code.append("")
            else:
                flush_para()
        elif is_heading(ln):
            flush_code()
            flush_para()
            out.extend([f"### {ln.strip().title()}", ""])
        elif preformatted(ln):
            flush_para()
            code.append(ln)
        else:
            flush_code()
            para.append(ln)
    flush_code()
    flush_para()
    return "\n".join(out).strip() or None


def module_builtin_counts(impl_files):
    """How many public builtins each source file implements."""
    c = {}
    for path in impl_files.values():
        c[path] = c.get(path, 0) + 1
    return c


def source_algorithm(name, module, counts):
    """The implementation's own account of how it works.

    The file header comment, used only when the file is effectively dedicated
    to this function (at most three builtins) and the comment is substantial.
    Those two conditions are what keep this honest: a shared file's header
    describes the module, not the function, and a one-line header is a title
    rather than an explanation. The modern one-symbol-per-file modules --
    list/, numbertheory/, stats/ -- are where the real algorithm prose lives."""
    if not module or counts.get(module, 0) > 3:
        return None
    path = ROOT / module
    if not path.exists():
        return None
    m = re.match(r"\s*/\*(.*?)\*/", path.read_text(errors="replace"), re.S)
    if not m:
        return None
    md = comment_to_markdown(m.group(1))
    if not md or len(md.split()) < 40:
        return None
    return md


def implementation_files():
    """symbol -> the file where its C function is DEFINED.

    `module` in the builtins map is where the function is *registered*, which
    for most of the tree is a shared init hub -- src/list/list_init.c alone
    registers dozens. That is the right answer for "where do I look this up"
    and the wrong one for "what does this file implement": the algorithm prose
    lives next to the code, in find_clusters.c, not in the registration line.

    Resolved in two hops: the registration names the C function, and a scan of
    the tree says where that function is defined. Exact, and free of guesswork
    about how a symbol name maps to an identifier."""
    reg = re.compile(r'symtab_add_builtin\s*\(\s*"((?:[^"\\]|\\.)*)"\s*,\s*(\w+)')
    define = re.compile(r"^\s*(?:static\s+)?Expr\s*\*\s*(\w+)\s*\(", re.M)
    fn_to_file, sym_to_fn = {}, {}
    for cfile in sorted(SRC.rglob("*.c")):
        if "external" in cfile.parts:
            continue
        text = cfile.read_text(errors="replace")
        rel = cfile.relative_to(ROOT).as_posix()
        for m in define.finditer(text):
            fn_to_file.setdefault(m.group(1), rel)
        for m in reg.finditer(text):
            sym = m.group(1).encode().decode("unicode_escape")
            sym_to_fn.setdefault(sym, m.group(2))
    return {sym: fn_to_file[fn] for sym, fn in sym_to_fn.items() if fn in fn_to_file}


def build_related(sections, valid):
    """symbol -> related symbols, strongest signal first.

    Two sources, both semantic rather than alphabetical. A spec heading that
    names several functions ("## AtomQ, NumberQ, IntegerQ") is a curated
    grouping, so its members are each other's closest relatives; after that,
    any function this one's own spec section names in backticks."""
    related = {}
    for name, info in sections.items():
        sibs = [s for s in info.get("siblings", []) if s != name and s in valid]
        body = info.get("body", "") or ""
        mentioned = []
        for m in re.findall(r"`([A-Z$][A-Za-z0-9$]*)`", body):
            if m != name and m in valid and m not in sibs and m not in mentioned:
                mentioned.append(m)
        picked = sibs + mentioned
        if picked:
            related[name] = picked[:8]
    return related


def tests_for(name):
    """Test files that exercise this function, as evidence rather than a claim.

    "Exercised by the test suite" appears in every Stable badge; naming the
    file lets a reader check."""
    if not TESTS_DIR.is_dir():
        return []
    hits = []
    pat = re.compile(r"\b%s\s*\[" % re.escape(name))
    for f in sorted(TESTS_DIR.glob("test_*.c")):
        try:
            if pat.search(f.read_text(errors="replace")):
                hits.append(f.relative_to(ROOT).as_posix())
        except OSError:
            continue
    return hits[:4]


def load_perf():
    """Measured timings from tools/docs_perf.py, if it has been run.

    Kept in a cache rather than measured inline: timing hundreds of functions
    would dominate generation, and a stale-but-labelled number beats making
    every doc build slow."""
    if not PERF_CACHE.exists():
        return {}
    try:
        return json.loads(PERF_CACHE.read_text())
    except (OSError, ValueError):
        return {}


def load_bench_rows():
    """The most recent complete benchmark run: real Mathilda / Wolfram / Python
    times for each case, already cross-checked for agreement."""
    res = sorted((BENCH_DIR / "results").glob("*.json")) if (BENCH_DIR / "results").is_dir() else []
    res = [p for p in res if "partial" not in p.name]
    if not res:
        return []
    try:
        return json.loads(res[-1].read_text()).get("rows", [])
    except (OSError, ValueError):
        return []


def experiment_attribution(valid):
    """symbol -> {experiment names}.

    An experiment counts for a function when its header names the function or
    the source file implementing it, or when the body calls it repeatedly. The
    header is the strong signal (it states what the experiment measures); the
    repeat-call rule catches experiments whose header is written in prose."""
    attr = {}
    for base in EXP_DIRS:
        if not base.is_dir():
            continue
        for d in sorted(base.iterdir()):
            if not d.is_dir() or not re.match(r"^\d\d-", d.name):
                continue
            for mfile in d.glob("*.m"):
                try:
                    txt = mfile.read_text(errors="replace")
                except OSError:
                    continue
                head = re.match(r"\s*\(\*(.*?)\*\)", txt, re.S)
                hdr = head.group(1) if head else ""
                names = set(re.findall(r"\b([A-Z][A-Za-z0-9]{2,})\b", hdr)) & valid
                for sym, cnt in collections.Counter(
                        re.findall(r"\b([A-Z][A-Za-z0-9]{2,})\[", txt)).items():
                    if sym in valid and cnt >= 3:
                        names.add(sym)
                for sym in names:
                    attr.setdefault(sym, set()).add(d.name)
    return attr


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

    ex, total = render_examples(info["examples"], info.get("applications"),
                                info.get("worked"), figures=info.get("figures"))
    lines.append(f"## Examples ({total})" if total else "## Examples")
    lines.append("")
    if ex:
        lines.append("Every input below was run against the current Mathilda "
                     "build and its output recorded.")
        lines.append("")
        lines.append(ex)
    else:
        lines.append("_No verified examples yet for this function._")
    lines.append("")

    if info.get("spec_detail"):
        # Options, algorithm behaviour and diagnostics, promoted out of the spec.
        # Placed before "Implementation notes" because it documents what the
        # function DOES; the notes below are about how it is built.
        lines.append("## Options & behaviour")
        lines.append("")
        lines.append(info["spec_detail"])
        lines.append("")

    if info.get("algorithm"):
        # How it actually works, taken from the implementation's own header
        # comment. This is the section a reference page usually cannot have:
        # it is written by whoever wrote the code, next to the code, so it goes
        # stale only when the file does.
        lines.append("## Algorithm")
        lines.append("")
        lines.append(info["algorithm"])
        lines.append("")

    perf = info.get("perf")
    bench = info.get("bench") or []
    if perf or bench:
        lines.append("## Performance")
        lines.append("")
        if perf:
            lines.append(f"Measured on {perf.get('host', 'the build machine')}"
                         f" at commit `{perf.get('commit', 'unknown')[:9]}`.")
            lines.append("")
            lines.append("| case | n | time |")
            lines.append("|---|---:|---:|")
            for row in perf.get("rows", []):
                lines.append(f"| {row['case']} | {row['n']:,} | {row['time']} |")
            lines.append("")
        if bench:
            lines.append("Against other systems, from the benchmark suite "
                         "(same input, results cross-checked for agreement):")
            lines.append("")
            lines.append("| case | Mathilda | Wolfram | Python |")
            lines.append("|---|---:|---:|---:|")
            for r in bench:
                t = r.get("times", {})
                def fmt(v):
                    return "--" if v is None else (f"{v:.3g} s" if isinstance(v, (int, float)) else str(v))
                lines.append(f"| {r.get('label', '')} | {fmt(t.get('mathilda'))} "
                             f"| {fmt(t.get('wolfram'))} | {fmt(t.get('python'))} |")
            lines.append("")

    lines.append("## Implementation notes")
    lines.append("")
    lines.append(info["notes"])
    lines.append("")

    lines.append("## References")
    lines.append("")
    if info.get("related"):
        # Related functions live here rather than in a section of their own:
        # both answer "where do I go from this page", and one heading for that
        # is enough. From the spec's own groupings, not an alphabetical
        # neighbourhood -- a heading naming several functions is a curated set,
        # and a function named in this one's prose is one the author thought
        # worth mentioning.
        urls = info.get("related_urls", {})
        lines.append("**See also:** " + ", ".join(
            f"[{r}](../../{urls[r]}/)" if r in urls else f"`{r}`"
            for r in info["related"]))
        lines.append("")
    for r in info["references"]:
        lines.append(f"- {r}")
    for t in info.get("tests", []):
        # Evidence for the "exercised by the test suite" claim in the badge.
        lines.append(f"- Tests: [`{t}`](https://github.com/stblake/mathilda/blob/main/{t})")
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
        for group, inputs in mine_example_inputs(body):
            pairs, _ = verify_block(inputs)   # FLINT routines emit no graphics
            if pairs:
                blocks.append((group, pairs))
                verified += len(pairs)
        status = derive_status(name, flint[name]["doc"], bool(blocks),
                               tests_text, lim_text)
        refs = build_references(name, flint[name]["module"], FLINT_SLUG, spec_rel)
        notes = render_impl_notes(name, attrs.get(name, []), body)
        detail = mine_spec_detail(body)
        info = {
            "name": name, "doc": flint[name]["doc"], "attrs": attrs.get(name, []),
            "examples": blocks, "status": status, "references": refs,
            "notes": notes, "overlay_body": None, "spec_detail": detail,
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

    # Link-rewrite context for the spec-mined bullets (see rewrite_spec_links).
    # A builtin's display category is fixed by CATEGORY_OVERRIDES / spec section /
    # OTHER_SLUG — computable up front, so `#anchor` links naming a builtin resolve
    # to that builtin's page even though `pages` is not fully built yet.
    cat_slugs = {slug for slug, _, _ in categories}

    def _cat_of(nm):
        sec = sections.get(nm)
        return CATEGORY_OVERRIDES.get(nm) or (sec["category"] if sec else None) or OTHER_SLUG

    anchor_to_cat = {nm.lower(): (_cat_of(nm), nm) for nm in names}

    # ---- enrichment indices, computed once for the whole run ----
    valid = set(names)
    impl_files = implementation_files()
    mod_counts = module_builtin_counts(impl_files)
    related = build_related(sections, valid)
    related_urls = {nm: f"{_cat_of(nm)}/{nm}" for nm in names}
    perf_data = load_perf()

    attribution = experiment_attribution(valid)
    bench_rows = load_bench_rows()
    bench_by_symbol = {}
    for sym, exps in attribution.items():
        rows = [r for r in bench_rows if r.get("experiment") in exps]
        # Keep the page readable: the slowest few cases say the most about
        # where a function actually costs something.
        rows.sort(key=lambda r: (r.get("times") or {}).get("mathilda") or 0, reverse=True)
        if rows:
            bench_by_symbol[sym] = rows[:6]
    print(f"  enrichment: {len(related)} with related symbols, "
          f"{len(bench_by_symbol)} with benchmark rows, "
          f"{len(perf_data)} with measured timings")

    plots_written = 0
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
        figures = {}
        for group, inputs in mine_example_inputs(body):
            pairs, figs = verify_block(inputs)
            if pairs:
                blocks.append((group, pairs))
                figures.update(figs)
                verified_examples += len(pairs)

        # Worked examples written as prose. Their inputs are real expressions;
        # their stated results are prose, so the binary supplies the output.
        worked_pairs = []
        prose_inputs = mine_prose_examples(body)
        if prose_inputs:
            worked_pairs, wfigs = verify_block(prose_inputs)
            figures.update(wfigs)
            verified_examples += len(worked_pairs)

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
        notes = render_impl_notes(name, attrs.get(name, []), body, impl_body,
                                  cat_slugs, anchor_to_cat)
        detail = mine_spec_detail(body)
        if detail:
            detail = strip_prose_examples(detail)
        if detail:
            detail = rewrite_spec_links(detail, cat_slugs, anchor_to_cat)

        overlay = load_overlay(name)
        overlay_body = None
        applications = None
        if overlay:
            m = overlay["meta"]
            if m.get("status") in STATUS_BADGE:
                status = m["status"]
            if isinstance(m.get("references"), list) and m["references"]:
                refs = m["references"] + refs
            overlay_body = overlay["body"] or None
            # The overlay's own "Worked examples" belong with the examples, not
            # in a footnote section below the references. Lift them out; the
            # remaining notes stay where they are.
            if overlay_body:
                m = re.search(r"^###\s*Worked examples\s*$(.*?)(?=^###\s|\Z)",
                              overlay_body, re.S | re.M)
                if m and m.group(1).strip():
                    app_body = m.group(1).strip()
                    # Parse into (input, output) pairs so these render through the
                    # same path as every other example group -- which is what
                    # lifts a `(* ... *)` note out as a sentence above the cell.
                    app_pairs = []
                    for blk in re.findall(r"```mathematica\n(.*?)```", app_body, re.S):
                        lines = blk.rstrip().split("\n")
                        for i, ln in enumerate(lines):
                            mi = re.match(r"^In\[\d+\]:=\s?(.*)$", ln)
                            if not mi:
                                continue
                            mo = (re.match(r"^Out\[\d+\]=\s?(.*)$", lines[i + 1])
                                  if i + 1 < len(lines) else None)
                            if mo:
                                app_pairs.append((mi.group(1).strip(), mo.group(1).strip()))
                    if app_pairs:
                        applications = ("pairs", app_pairs)
                    else:
                        applications = (app_body, len(re.findall(r"^In\[", app_body, re.M)) or 1)
                    overlay_body = (overlay_body[:m.start()] +
                                    overlay_body[m.end():]).strip() or None

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
            "overlay_body": overlay_body, "spec_detail": detail,
            "applications": applications, "worked": worked_pairs,
            "algorithm": source_algorithm(name, impl_files.get(name), mod_counts),
            "figures": figures,
            "related": related.get(name, []),
            "related_urls": related_urls,
            "tests": tests_for(name),
            "perf": perf_data.get(name),
            "bench": bench_by_symbol.get(name, []),
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
            # Saved figures, so a graphics example shows its plot before the
            # reader runs anything. Keyed by the input text rather than by
            # position, which survives example trimming and renumbering.
            figs = pages[name].get("figures") or {}
            if figs:
                pdir = PLOTS_OUT / slug
                pdir.mkdir(parents=True, exist_ok=True)
                (pdir / f"{name}.json").write_text(json.dumps(figs))
                plots_written += 1

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
    print(f"  Saved figures:     {plots_written} function(s)")
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

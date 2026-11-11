#!/usr/bin/env python3
"""import_abbasi_pde.py -- convert Nasser Abbasi's "PDE in CAS" corpus into
Mathilda DSolve PDE corpus files.

Inputs (both fetched once and checked in under DSolve_test_status/pde/raw/):

  Mathematica_PDE_IN_CAS_problems.txt   the 2012 problems in Wolfram syntax, one
                                         "(*PDE number N*)" block each
  lookup.tsv                             per-problem section / description /
                                         Mathematica 15 and Maple 2026.1 result,
                                         scraped from Ch1.S4.htm (--scrape)

Each block is a short script:

    k=1/100; L=5; f=x;
    pde = D[u[x, t], t] == k*D[u[x, t], {x, 2}];
    bc  = {Derivative[1, 0][u][0, t] == 0, Derivative[1, 0][u][L, t] == 0};
    ic  = u[x, 0] == f;
    sol = AbsoluteTiming[TimeConstrained[DSolve[{pde, bc, ic}, u[x, t], {x, t}], 600]];

The converter splits the block into top-level statements, keeps every
`symbol = value` assignment, finds the statement containing DSolve[...] (or
DSolveValue[...]), peels the AbsoluteTiming/TimeConstrained/Simplify wrappers,
splits the DSolve arguments at top level and inlines the assignments into them
(word-boundary safe, to a fixed point).  The output record is

    {"pde-N", "6.3.25.10", eqns, fn, vars, {opts...}, "mma", "maple"}

where eqns is the equation or a flat List of PDE + conditions, fn is the unknown
as given (u or u[x, t]), vars the independent-variable List, opts the remaining
DSolve options (Assumptions -> ...), and mma/maple are "Y" / "N" / "?".

The raw file is untrusted data: this script only reads it, and it is run with
`python3 -I`.  A block that cannot be converted is never silently dropped -- it
is reported and written to skipped.tsv with a reason.

Usage:
  python3 -I tools/pde_corpus/import_abbasi_pde.py --scrape Ch1.S4.htm   # writes raw/lookup.tsv
  python3 -I tools/pde_corpus/import_abbasi_pde.py                       # writes DE_PDE_*.m
"""
import argparse
import html
import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PDE_DIR = os.path.join(ROOT, "DSolve_test_status", "pde")
RAW = os.path.join(PDE_DIR, "raw", "Mathematica_PDE_IN_CAS_problems.txt")
LOOKUP = os.path.join(PDE_DIR, "raw", "lookup.tsv")

WRAPPERS = ("AbsoluteTiming", "TimeConstrained", "Simplify", "FullSimplify",
            "Quiet", "Timing")


# --------------------------------------------------------------------------
# lookup-table scraping
# --------------------------------------------------------------------------

def _cell_text(td):
    td = re.sub(r"<math.*?alttext=\"(.*?)\".*?</math>", lambda m: m.group(1), td, flags=re.S)
    td = re.sub(r"<[^>]+>", " ", td)
    return " ".join(html.unescape(td).split())


def scrape(path):
    """Ch1.S4.htm -> rows (N, section, category, description, mma, mma_t, maple, maple_t)."""
    h = open(path, encoding="utf-8").read()
    out = {}
    for row in re.findall(r"<tr.*?</tr>", h, re.S):
        m = re.search(r'<a href="[^"]*#X(\d+)" title="([0-9.]+)', row)
        if not m:
            continue
        n, section = int(m.group(1)), m.group(2).rstrip(".")
        cells = [_cell_text(c) for c in re.findall(r"<td.*?</td>", row, re.S)]
        if len(cells) < 7:
            continue

        def mark(s):
            return "Y" if "✓" in s else ("N" if ("✗" in s or "×" in s or "X" in s) else "?")
        out[n] = (n, section, cells[1], cells[2].replace("\t", " "),
                  mark(cells[3]), cells[4], mark(cells[5]), cells[6])
    return [out[k] for k in sorted(out)]


# --------------------------------------------------------------------------
# Wolfram-source splitting (strings, comments and brackets respected)
# --------------------------------------------------------------------------

def strip_comments(s):
    out, i, depth = [], 0, 0
    while i < len(s):
        if s.startswith("(*", i):
            depth += 1; i += 2; continue
        if depth and s.startswith("*)", i):
            depth -= 1; i += 2; continue
        if not depth:
            if s[i] == '"':                      # copy a string literal verbatim
                j = i + 1
                while j < len(s) and s[j] != '"':
                    j += 2 if s[j] == "\\" else 1
                out.append(s[i:j + 1]); i = j + 1; continue
            out.append(s[i])
        i += 1
    return "".join(out)


def split_top(s, sep):
    """Split s at top-level occurrences of the one-character separator sep."""
    parts, depth, cur, i = [], 0, [], 0
    while i < len(s):
        c = s[i]
        if c == '"':
            j = i + 1
            while j < len(s) and s[j] != '"':
                j += 2 if s[j] == "\\" else 1
            cur.append(s[i:j + 1]); i = j + 1; continue
        if c in "[{(":
            depth += 1
        elif c in "]})":
            depth -= 1
        if c == sep and depth == 0:
            parts.append("".join(cur)); cur = []
        else:
            cur.append(c)
        i += 1
    parts.append("".join(cur))
    return parts


def call_args(s, start):
    """s[start] == '['; return (args list, index after the matching ']')."""
    depth, i = 0, start
    while i < len(s):
        c = s[i]
        if c == '"':
            j = i + 1
            while j < len(s) and s[j] != '"':
                j += 2 if s[j] == "\\" else 1
            i = j + 1; continue
        if c in "[{(":
            depth += 1
        elif c in "]})":
            depth -= 1
            if depth == 0:
                return split_top(s[start + 1:i], ","), i + 1
        i += 1
    raise ValueError("unbalanced brackets")


# Names the source uses for whole equations / condition lists.  Filled from the
# corpus itself (every target some block binds to an equation) in main().
STRUCTURAL_NAMES = set()

ASSIGN = re.compile(r"^\s*([A-Za-z$][A-Za-z0-9$]*)\s*=(?![=!])\s*(.*)$", re.S)


def inline(text, env):
    """Substitute assigned symbols (word-boundary safe) to a fixed point."""
    for _ in range(12):
        new = text
        for name, val in env.items():
            new = re.sub(r"(?<![A-Za-z0-9$`])" + re.escape(name) + r"(?![A-Za-z0-9$`])",
                         lambda m, v=val: "(" + v + ")", new)
        if new == text:
            return text
        text = new
    raise ValueError("assignment inlining did not converge")


def convert_block(n, body):
    """One block -> (eqns, fn, vars, opts) as Wolfram source strings."""
    stmts = [st.strip() for st in split_top(strip_comments(body), ";")]
    env, call = {}, None
    for st in stmts:
        if not st:
            continue
        m = ASSIGN.match(st)
        if m and m.group(1) == "sol" or (not m and "DSolve" in st):
            rhs = m.group(2) if m else st
            k = re.search(r"\bDSolve(Value)?\[", rhs)
            if k and call is None:
                args, _ = call_args(rhs, k.end() - 1)
                call = args
            continue
        if m and "DSolve" not in m.group(2):
            env[m.group(1)] = m.group(2).strip()
    if call is None:
        raise ValueError("no DSolve call")
    args = [a.strip() for a in call]
    env.pop("sol", None)
    if len(args) < 3:
        raise ValueError("DSolve with fewer than 3 arguments")
    eqn = inline(args[0], env)
    fn = inline(args[1], env)
    rest = [inline(a, env) for a in args[2:]]
    # The source was run as ONE Wolfram session, so a structural name the block
    # never assigns (`bc` in pde-392) silently picked up the previous block's
    # value.  That is a defect of the source, not a problem to reproduce: refuse
    # it.  Parameter names (alpha, T0) are never bound to equations, so are allowed.
    for name in STRUCTURAL_NAMES:
        if re.search(r"(?<![A-Za-z0-9$`])" + name + r"(?![A-Za-z0-9$`\[])",
                     " ".join([eqn, fn] + rest)):
            raise ValueError("references `%s`, which this block never assigns" % name)
    opts = [a for a in rest if "->" in a or ":>" in a]
    pos = [a for a in rest if a not in opts]
    if len(pos) == 1:
        vars_ = pos[0]
    else:                                   # legacy DSolve[eq, u[x,t], x, t]
        vars_ = "{" + ", ".join(pos) + "}"
    # {pde, bc, ic} with bc itself a List is flattened by the prelude, which
    # evaluates the slot; the record keeps the source spelling.
    return eqn, fn, vars_, opts


def read_blocks(path):
    text = open(path, encoding="utf-8").read()
    pieces = re.split(r"\(\*\s*PDE number\s+(\d+)\s*\*\)", text)
    blocks = []
    for i in range(1, len(pieces), 2):
        body = pieces[i + 1].split("(*----------------------------*)")[0]
        blocks.append((int(pieces[i]), body))
    return blocks


def load_lookup():
    rows = {}
    if os.path.exists(LOOKUP):
        for line in open(LOOKUP, encoding="utf-8"):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 8:
                rows[int(f[0])] = f
    return rows


def file_key(section):
    """Group records into corpus files: chapter 6 by handbook chapter (6.k),
    every other chapter as one file."""
    parts = section.split(".")
    if parts[0] == "6" and len(parts) > 1:
        return "6_" + parts[1]
    return parts[0]


def wl_string(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scrape", metavar="Ch1.S4.htm")
    a = ap.parse_args()
    if a.scrape:
        rows = scrape(a.scrape)
        with open(LOOKUP, "w", encoding="utf-8") as f:
            f.write("# N\tsection\tcategory\tdescription\tmma\tmma_time\tmaple\tmaple_time\n")
            f.write("# scraped from https://www.12000.org/my_notes/pde_in_CAS/"
                    "maple_2026_1_and_mma_15/Ch1.S4.htm\n")
            for r in rows:
                f.write("\t".join(str(x) for x in r) + "\n")
        print(f"lookup: {len(rows)} rows -> {LOOKUP}")
        return 0

    lookup = load_lookup()
    groups, skipped = {}, []
    blocks = read_blocks(RAW)
    for _, body in blocks:
        for st in split_top(strip_comments(body), ";"):
            m = ASSIGN.match(st)
            # a name is structural when some block binds it to an equation
            # (or a List of them); parameters (alpha, beta, T0) never are
            if m and m.group(1) != "sol" and "==" in m.group(2):
                STRUCTURAL_NAMES.add(m.group(1))
    for n, body in blocks:
        try:
            eqn, fn, vars_, opts = convert_block(n, body)
        except ValueError as e:
            skipped.append((n, str(e)))
            continue
        lk = lookup.get(n)
        section = lk[1] if lk else "?"
        mma, maple = (lk[4], lk[6]) if lk else ("?", "?")
        rec = "{%s, %s, %s, %s, %s, {%s}, %s, %s}" % (
            wl_string("pde-%d" % n), wl_string(section), eqn, fn, vars_,
            ", ".join(opts), wl_string(mma), wl_string(maple))
        groups.setdefault(file_key(section), []).append(rec)

    for key, recs in sorted(groups.items()):
        out = os.path.join(PDE_DIR, "DE_PDE_%s.m" % key)
        with open(out, "w", encoding="utf-8") as f:
            f.write("(* GENERATED by tools/pde_corpus/import_abbasi_pde.py from\n"
                    " * raw/Mathematica_PDE_IN_CAS_problems.txt (N. Abbasi, PDE in CAS).\n"
                    " * Do not hand-edit: fix the converter and regenerate.\n"
                    " * Record: {label, section, eqns, fn, vars, opts, mma, maple}. *)\n{\n")
            f.write(",\n".join(recs))
            f.write("\n}\n")
    with open(os.path.join(PDE_DIR, "skipped.tsv"), "w", encoding="utf-8") as f:
        f.write("# N\treason\n")
        for n, why in skipped:
            f.write("%d\t%s\n" % (n, why))
    total = sum(len(v) for v in groups.values())
    print(f"blocks: {len(blocks)}  converted: {total}  skipped: {len(skipped)}")
    for key, recs in sorted(groups.items()):
        print(f"  DE_PDE_{key}.m: {len(recs)}")
    for n, why in skipped:
        print(f"  SKIP pde-{n}: {why}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

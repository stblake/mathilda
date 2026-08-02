#!/usr/bin/env python3
"""
nd_surface_audit.py -- the same operation on all THREE representations, timed.

THE FAILURE MODE THIS EXISTS FOR. Mathilda has two array representations that a
user can hand to a builtin, and they are not the same object:

    a PACKED List     invisible; `Head` is List, it prints as {1., 2., 3.}, and
                      src/eval.c's transparency gate MATERIALISES it to a nested
                      List for any head that has not opted in via pack.c's AWARE
                      list. Failing safe means failing silently. This is what
                      automatic packing produces, and also what `ToNDArray[...]`
                      returns -- ToNDArray does NOT build the visible form.
    a visible NDArray `NDArray[...]`, the head the user writes. The gate does
                      NOT touch it (`is_packed_list` is false), so it reaches
                      the builtin as a buffer whatever the AWARE list says --
                      and equally, reaches heads that cannot read one at all.

Those two facts have opposite signs, and the consequence is that **a head can be
fast on one surface and slow on the other, on identical data**. It has happened:

    Fourier    complete NDArray fast path, absent from AWARE -- so
               Fourier[packedList] ran 15.7x slower than Fourier[NDArray[...]]
               and every existing audit passed.

And the converse, which nothing in the tree looks for at all: a head that only
ever grew a *packed* path answers a visible NDArray by falling through to the
generic expression walker, one Expr per element.

The two static audits cannot see either case.
  * tools/check_packed_aware.py reads dispatch sites out of the source. A head
    with NO fast path on EITHER surface has no dispatch site, so it is invisible
    to it by construction -- which is exactly how `DeleteDuplicates` sat at 72x
    NumPy through four sweeps.
  * tools/numeric_coverage.py joins the registries. Registration is not speed.

Only running the same expression on all three representations separates them.
So that is what this does, reusing tools/numeric_sweep.py's probe catalogue --
the same 283 expressions, the same data, the same checksums -- over:

    plain     MATHILDA_NO_PACK=1, the one-Expr-per-element floor
    packed    default
    ndarray   every array in the preamble rebuilt as NDArray[x, DataType -> dt]

FINDINGS. Six kinds, and each is a defect of a different shape:

    DISAGREE   the surfaces do not compute the same thing. A representation may
               never change a value (pack.h), so this is always a bug. It found
               two: every real kernel truncating a visible int64 buffer, and
               MemberQ/Count/Position/Cases/FreeQ answering False/0/{} because a
               visible array exposes no args[] for the matcher to walk.
    ND-UNSUPPORTED  answers on plain and packed, UNEVALUATED on a visible array.
               Fails loudly, so it ranks below DISAGREE. 26 heads.
    NO-PATH    packed is no faster than plain -- the head has no buffer path
               being reached at all, whatever the registries claim. The
               DeleteDuplicates class.
    SKEW       packed is materially slower than the identical visible NDArray.
               A dispatch gap: the fast path exists and the gate ate the buffer
               before it could fire. The Fourier class.
    ND-SLOW    the visible NDArray is materially slower than the packed List.
               The converse: a head that only grew a packed path and treats a
               visible array as an opaque atom. `setop_i64` testing
               is_packed_list where is_ndarray was meant cost Union 145x.
    NO-ANSWER  answered on NO surface: the function does not exist, or declines
               these arguments everywhere. A COVERAGE finding, emphatically not
               a disagreement -- the first version of this classifier reported
               all 29 as DISAGREE, because numeric_sweep.agree() returns False
               for "UNEVAL" against itself by its own contract.

NumPy stays the ruler for absolute speed, and that column comes straight from
numeric_sweep.

Usage:
    python3 tools/nd_surface_audit.py                    # full audit
    python3 tools/nd_surface_audit.py --group struct     # one group
    python3 tools/nd_surface_audit.py --only tally,union
    python3 tools/nd_surface_audit.py --json out.json
    python3 tools/nd_surface_audit.py --findings         # only the flagged rows
"""
import argparse
import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numeric_sweep as ns  # noqa: E402

ROOT = ns.ROOT

# ---------------------------------------------------------------------------
# The array variables numeric_sweep's preamble defines. Rebuilding each as a
# visible NDArray gives the identical values on that surface -- identical, so
# a timing difference between the two runs is a dispatch difference and nothing
# else.
#
# tf/tg are deliberately absent: they are Boolean lists, there is no boolean
# dtype (performance.md §13, gap C.1), and ToNDArray declines them. Converting
# them would silently leave them as plain Lists and make the ndarray column a
# copy of the packed one for the mask group.
# ---------------------------------------------------------------------------
ARRAY_VARS = [
    "v", "w", "u", "p", "iv", "jv", "kv", "bv", "cv", "sv",
    "A1000", "A500", "A300", "A200", "A6",
    "b1000", "b500", "b6", "im", "k5", "k5v", "idx", "src",
]

# The conversion has to build a VISIBLE NDArray, and `ToNDArray` does not: it
# returns a PACKED List (`Head` is List, `present_as` is NDA_HEAD_LIST), which
# is a third thing again. Using it here would have made this column a
# force-repack of the packed column and left the visible surface — the one with
# no transparency gate in front of it — untested. `NDArray[...]` is the visible
# constructor, and it needs the dtype named explicitly: `NDArray[{1, 2, 3}]`
# defaults to float64 and would silently change every integer probe's dtype.
ND_PRE = ns.MPRE + """
ndof[x_] := Module[{t = ToNDArray[x]},
  If[NDArrayQ[t], NDArray[x, DataType -> DataType[t]], x]];
""" + "\n".join(f"{x} = ndof[{x}];" for x in ARRAY_VARS) + """
ndbad = Select[{%s}, Head[#] =!= NDArray &];
If[Length[ndbad] > 0, Print["!ND-CONVERT-FAILED\t", Length[ndbad]]];
Clear[ndbad];

(* The checksums have to survive a VISIBLE NDArray result, and numeric_sweep's
   do not: ck is Total[Flatten[{x}]], Flatten only descends into List heads, and
   a visible array's head is NDArray -- so Flatten[{NDArray[...]}] stays nested,
   Total does not reduce, and every correct row read as a value mismatch.
   FromNDArray is the identity on anything that is not an array, so redefining
   them here changes nothing about WHAT is summed, only whether it reduces. *)
ck[x_]  := N[Total[Flatten[{FromNDArray[x]}]]];
ckt[x_] := N[Count[Flatten[{FromNDArray[x]}], True]];
""" % ", ".join(ARRAY_VARS)


# Below this, wall-clock noise swamps the comparison and every ratio is
# meaningless. 50 us is ~20x the process-to-process jitter measured on this
# host and still far under any row a workload cares about.
FLOOR = 50e-6

# A surface is "materially" slower at this ratio. Chosen from the measured
# spread of repeated identical runs (< 1.2x), doubled.
SKEW = 1.5

# packed/plain below this means the buffer bought nothing.
NOGAIN = 1.3


def run_surface(probes, name, tmpdir):
    saved = ns.MPRE
    env = None
    if name == "plain":
        env = {"MATHILDA_NO_PACK": "1"}
    elif name == "ndarray":
        ns.MPRE = ND_PRE
    try:
        sys.stderr.write(f"surface {name}: {len(probes)} probes\n")
        return ns.run_mathilda(probes, env_extra=env, tmpdir=tmpdir)
    finally:
        ns.MPRE = saved


# ---------------------------------------------------------------------------
# Survival: does the RESULT come back packed?
#
# A separate question from "is this head fast", and the one that explains most
# of the slow rows downstream of it. Packing is a chain: `Mod[Range[10^6], 1000]`
# answering with a plain List of 10^6 Integers does not make Mod slow -- Mod is
# fine -- it makes Tally, Union, Sort, BinCounts and every other consumer of
# that value slow, because the buffer they were written for never arrives.
#
# The test is exact rather than heuristic: a result is a MISSED packing iff it
# is not packed AND ToNDArray would have packed it. That distinguishes "this
# head returns a String / a ragged list / an Association, and rightly does not
# pack" from "this head built the very array the representation exists for and
# handed back boxes".
#
# With one correction, which the first run needed: ToNDArray ignores the SIZE
# THRESHOLD (PACK_MIN_ELEMENTS = 250, src/pack.h) on purpose, so it packs a
# 3-element list that automatic packing correctly leaves alone. Reporting those
# put `Extract` (3 elements), `TakeLargest` (10), `MinMax` (2) and `Cross` (3)
# on the finding list, where nothing is wrong at all. A result shorter than the
# threshold is reported as "short", not as a miss.
# ---------------------------------------------------------------------------
PACK_MIN_ELEMENTS = 250

SURVIVAL_TMPL = """
Clear[r];
r = %s;
Print["@", "%s", "\\t", If[NDArrayQ[r], "packed", "plain"], "\\t",
      If[NDArrayQ[ToNDArray[r]], "packable", "unpackable"], "\\t",
      If[ListQ[r] || NDArrayQ[r], Length[Flatten[{FromNDArray[r]}]], 0]];
Clear[r];
"""


def survival_source(probes):
    out = [ns.MPRE]
    for pr in probes:
        out.append(SURVIVAL_TMPL % (pr["m"], pr["id"]))
    return "\n".join(out)


def run_survival(probes, tmpdir):
    saved = ns.m_source
    ns.m_source = survival_source
    try:
        raw = ns.run_mathilda(probes, tmpdir=tmpdir)
    finally:
        ns.m_source = saved
    return raw


def parse_survival(text):
    out = {}
    for line in text.splitlines():
        if not line.startswith("@"):
            continue
        parts = line[1:].split("\t")
        if len(parts) < 4:
            continue
        out[parts[0].strip()] = {
            "packed": parts[1].strip() == "packed",
            "packable": parts[2].strip() == "packable",
            "length": parts[3].strip(),
        }
    return out


def survival_main(probes):
    saved_parse = ns.parse_m
    ns.parse_m = parse_survival
    try:
        with tempfile.TemporaryDirectory() as td:
            got = run_survival(probes, td)
    finally:
        ns.parse_m = saved_parse

    hdr = (f"{'probe':<26} {'group':<13} {'result':>9} {'packable':>10} "
           f"{'elems':>9}  expr")
    print(hdr)
    print("-" * 100)
    missed, short = [], []
    for p in probes:
        g = got.get(p["id"])
        if not g:
            print(f"{p['id']:<26} {p['group']:<13} {'—':>9} {'—':>10} "
                  f"{'—':>9}  {p['m']}")
            continue
        try:
            n = int(float(g["length"]))
        except ValueError:
            n = 0
        flag = ""
        if not g["packed"] and g["packable"]:
            if n >= PACK_MIN_ELEMENTS:
                flag = "  <-- MISSED"
                missed.append((p["id"], p["group"], n))
            else:
                flag = "  (short)"
                short.append(p["id"])
        print(f"{p['id']:<26} {p['group']:<13} "
              f"{'packed' if g['packed'] else 'plain':>9} "
              f"{'yes' if g['packable'] else 'no':>10} {n:>9}  {p['m']}{flag}")
    print(f"\n{len(missed)} results are packable, at or over the {PACK_MIN_ELEMENTS}"
          f"-element threshold, and come back as plain Lists:")
    for pid, grp, n in missed:
        print(f"  {pid:<26} {grp:<13} {n:>9} elements")
    print(f"\n{len(short)} more are packable but under the threshold, where not "
          f"packing is correct:\n  " + ", ".join(short))
    return 0


def classify(r):
    """Every finding this row supports, worst first."""
    out = []
    pl, pk, nd = r["plain"], r["packed"], r["ndarray"]

    chks = {k: r[k + "_chk"] for k in ("plain", "packed", "ndarray")
            if r.get(k + "_chk") is not None}
    if chks:
        vals = list(chks.values())
        uneval = [v == "UNEVAL" for v in vals]
        if all(uneval):
            # The head answered on NO surface. That is a COVERAGE finding -- the
            # function does not exist, or declines these arguments everywhere --
            # and emphatically not a disagreement between representations. The
            # first version of this reported all 40-odd of them as DISAGREE,
            # because numeric_sweep.agree() returns False for "UNEVAL" against
            # itself, by design ("the head did not answer; never agreement").
            # Reading that wall of red cost real time; `numeric_coverage.py
            # --missing` names the 96 heads that simply are not defined.
            out.append("NO-ANSWER")
        elif any(uneval):
            # Answered on some surfaces and not others. When it is the visible
            # array that declines, this is the support gap the round was written
            # to find, so name it specifically.
            if r.get("ndarray_chk") == "UNEVAL":
                out.append("ND-UNSUPPORTED")
            else:
                out.append("DISAGREE")
        elif not all(ns.agree(vals[0], v) for v in vals[1:]):
            out.append("DISAGREE")

    if pk and nd and max(pk, nd) > FLOOR:
        if pk / nd > SKEW:
            out.append("SKEW")
        elif nd / pk > SKEW:
            out.append("ND-SLOW")

    if pl and pk and pl > FLOOR and pl / pk < NOGAIN:
        out.append("NO-PATH")

    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="comma-separated probe ids")
    ap.add_argument("--group", help="comma-separated groups: " +
                    ",".join(ns.GROUPS))
    ap.add_argument("--json", metavar="PATH")
    ap.add_argument("--findings", action="store_true",
                    help="print only the rows with a finding")
    ap.add_argument("--numpy", action="store_true",
                    help="also run the NumPy column (slower)")
    ap.add_argument("--survival", action="store_true",
                    help="do not time anything: report which probes RETURN a "
                         "packed array, and which hand back a plain List that "
                         "ToNDArray would have packed")
    args = ap.parse_args()

    probes = ns.PROBES
    if args.group:
        want = set(args.group.split(","))
        probes = [p for p in probes if p["group"] in want]
    if args.only:
        want = set(args.only.split(","))
        probes = [p for p in probes if p["id"] in want]

    if args.survival:
        return survival_main(probes)

    res, pyres = {}, {}
    with tempfile.TemporaryDirectory() as td:
        for surf in ("packed", "ndarray", "plain"):
            res[surf] = run_surface(probes, surf, td)
        if args.numpy:
            sys.stderr.write("numpy: running\n")
            pyres = ns.run_python(probes, td)

    rows = []
    for p in probes:
        row = {"id": p["id"], "group": p["group"], "expr": p["m"]}
        for surf in ("plain", "packed", "ndarray"):
            got = res[surf].get(p["id"], {})
            row[surf] = got.get("min")
            row[surf + "_chk"] = got.get("chk")
        q = pyres.get(p["id"], {})
        row["numpy"] = q.get("min")
        row["gain"] = (row["plain"] / row["packed"]
                       if row["plain"] and row["packed"] else None)
        row["skew"] = (row["packed"] / row["ndarray"]
                       if row["packed"] and row["ndarray"] else None)
        row["vs_numpy"] = (row["packed"] / row["numpy"]
                           if row["packed"] and row["numpy"] else None)
        row["findings"] = classify(row)
        rows.append(row)

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(rows, fh, indent=1)
        sys.stderr.write(f"wrote {args.json}\n")

    shown = [r for r in rows if r["findings"]] if args.findings else rows
    hdr = (f"{'probe':<26} {'group':<13} {'plain':>11} {'packed':>11} "
           f"{'ndarray':>11} {'gain':>7} {'skew':>7}  findings")
    print(hdr)
    print("-" * len(hdr))
    for r in shown:
        g = f"{r['gain']:.1f}x" if r["gain"] else "—"
        s = f"{r['skew']:.2f}x" if r["skew"] else "—"
        print(f"{r['id']:<26} {r['group']:<13} {ns.fmt_t(r['plain']):>11} "
              f"{ns.fmt_t(r['packed']):>11} {ns.fmt_t(r['ndarray']):>11} "
              f"{g:>7} {s:>7}  {','.join(r['findings'])}")

    print()
    for kind in ("DISAGREE", "ND-UNSUPPORTED", "SKEW", "ND-SLOW", "NO-PATH",
                 "NO-ANSWER"):
        hits = [r for r in rows if kind in r["findings"]]
        print(f"{kind:<10} {len(hits):>3}" +
              ("  " + ", ".join(h["id"] for h in hits) if hits else ""))
    nores = [r["id"] for r in rows if r["packed"] is None]
    if nores:
        print(f"{'NO-RESULT':<10} {len(nores):>3}  " + ", ".join(nores))
    print(f"\n{len(rows)} probes on 3 surfaces")
    return 0


if __name__ == "__main__":
    sys.exit(main())

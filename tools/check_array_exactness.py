#!/usr/bin/env python3
"""
check_array_exactness.py -- does any fast path return a two-headed array?

THE INVARIANT. A routine handed a packed array / NDArray must answer with a
scalar or with an array of ONE element head. It must never invent a mixture:

    DiagonalMatrix[{1., 2., 3.}]  ->  {{1., 0, 0}, {0, 2., 0}, {0, 0, 3.}}

is wrong twice over. It disagrees with Mathematica, which gives all Reals -- and
because a matrix of two heads fits no uniform buffer, it also silently drops
every consumer downstream off the fast path. That one cost 320x against NumPy
and was defended in four places as a correctness requirement before anyone
checked it.

WHY THIS NEEDS A TOOL. The defect is invisible to the three checks that already
exist. `tools/check_packed_aware.py` asks whether a head OPTED IN; these heads
had. The differential sweeps in tests/test_packed_list.c compare the packed and
unpacked paths against EACH OTHER, and a producer that invents an exact zero
does so on both, so they agree -- on the wrong answer. And no timing catches it,
because the answer is right enough to print. What is left is to look at the
element heads of the result, which is all this file does.

WHAT IT DOES NOT CHECK. This reads element HEADS, not values and not whether
the head answered at all: an UNEVALUATED result flattens to its own arguments,
so a routine that stopped working entirely still reports single-headed. That is
not hypothetical -- it is how the NullSpace regression of 2026-08-01 passed this
tool while returning `NullSpace[...]` unevaluated for every matrix past the
packing threshold. Pair it with value assertions; `test_no_two_headed_results`
in tests/test_packed_list.c is the companion.

WHAT IT REPORTS. Three states per probe:

    PACKED   the result is a buffer -- uniform by construction, nothing to check
    OK       an ordinary List, but its numeric leaves are all exact or all
             inexact
    MIXED    exact AND inexact numeric leaves in one result

MIXED IS A WORK LIST, NOT A FAILURE. Some heads are legitimately mixed and
Mathematica agrees -- `PadRight[v, n]` with the default exact `0` fill on Real
data is the documented example, and `Clip` with an exact bound is another. Those
live in the EXEMPT table below WITH THE MATHEMATICA OUTPUT THAT JUSTIFIES THEM,
because the whole reason this file exists is that "Mathematica does it too" was
once asserted without being checked. Anything not exempt is a bug until a
`wolframscript` run says otherwise.

Usage:
    python3 tools/check_array_exactness.py            # audit, table + summary
    python3 tools/check_array_exactness.py --all      # include OK/PACKED rows
    python3 tools/check_array_exactness.py --json out.json
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
MATHILDA = os.environ.get("MATHILDA_BIN", os.path.join(ROOT, "Mathilda"))
TIMEOUT = int(os.environ.get("EXACTNESS_TIMEOUT", "300"))
BATCH = 12

import numeric_sweep  # noqa: E402  -- reuses its probe corpus and preamble


# ---------------------------------------------------------------------------
# Extra corpus: the array-returning surface the timing sweep under-covers.
#
# numeric_sweep's probes are chosen to be TIMED, so they favour one big shape
# per head and skip the small/edge forms where an exactness bug actually shows.
# These are chosen to be READ: every one returns an array derived from machine
# numbers, at a size above the packing threshold so the fast path is live.
# ---------------------------------------------------------------------------
EXTRA = [
    # --- constructors: the zeros they invent are the whole question ---
    ("x_diagonalmatrix_real",  "DiagonalMatrix[Table[N[i], {i, 30}]]"),
    ("x_diagonalmatrix_int",   "DiagonalMatrix[Table[i, {i, 30}]]"),
    ("x_diagonalmatrix_mixed", "DiagonalMatrix[Table[If[i == 1, 1., i], {i, 30}]]"),
    ("x_diagonalmatrix_band",  "DiagonalMatrix[Table[N[i], {i, 30}], 2]"),
    ("x_diagonalmatrix_pad",   "DiagonalMatrix[Table[N[i], {i, 20}], 0, 40]"),
    ("x_identitymatrix",       "IdentityMatrix[30]"),
    ("x_unitvector",           "UnitVector[400, 7]"),
    ("x_unitvector_mach",      "UnitVector[400, 7, WorkingPrecision -> MachinePrecision]"),
    ("x_subdivide_real",       "Subdivide[0., 1., 400]"),
    ("x_subdivide_mixed",      "Subdivide[0, 1., 400]"),
    ("x_subdivide_mixed2",     "Subdivide[3., 1, 400]"),
    ("x_subdivide_rat",        "Subdivide[1/2, 1., 400]"),
    ("x_constantarray",        "ConstantArray[1., 400]"),
    ("x_range_real",           "Range[1., 400.]"),
    ("x_table_real",           "Table[N[i], {i, 400}]"),

    # --- linear algebra: every routine that rebuilds a matrix ---
    ("x_inverse",       "Inverse[Table[If[i == j, 30., 1./(1. + Abs[i - j])], {i, 30}, {j, 30}]]"),
    ("x_inverse_diag",  "Inverse[DiagonalMatrix[Table[N[i], {i, 30}]]]"),
    ("x_rowreduce",     "RowReduce[Table[If[i == j, 30., 1./(1. + Abs[i - j])], {i, 30}, {j, 30}]]"),
    ("x_rowreduce_diag","RowReduce[DiagonalMatrix[Table[N[i], {i, 30}]]]"),
    ("x_linearsolve",   "LinearSolve[Table[If[i == j, 30., 1./(1. + Abs[i - j])], {i, 30}, {j, 30}],"
                        " Range[1., 30.]]"),
    ("x_pseudoinverse", "PseudoInverse[Table[If[i == j, 30., 1./(1. + Abs[i - j])], {i, 30}, {j, 30}]]"),
    ("x_leastsquares",  "LeastSquares[Table[If[i == j, 30., 1./(1. + Abs[i - j])], {i, 30}, {j, 30}],"
                        " Range[1., 30.]]"),
    ("x_matrixpower",   "MatrixPower[Table[N[1./(1. + Abs[i - j])], {i, 30}, {j, 30}], 3]"),
    ("x_nullspace",     "NullSpace[Table[N[i j], {i, 30}, {j, 30}]]"),
    ("x_orthogonalize", "Orthogonalize[Table[If[i == j, 30., 1./(1. + Abs[i - j])], {i, 30}, {j, 30}]]"),
    ("x_transpose",     "Transpose[Table[N[i j], {i, 30}, {j, 30}]]"),
    ("x_conjtranspose", "ConjugateTranspose[Table[N[i j], {i, 30}, {j, 30}]]"),
    ("x_cholesky",      "CholeskyDecomposition[Table[If[i == j, 30., 1./(1. + Abs[i - j])], {i, 30}, {j, 30}]]"),
    ("x_hilbert",       "HilbertMatrix[30]"),
    ("x_toeplitz",      "ToeplitzMatrix[Table[N[i], {i, 30}]]"),
    ("x_hankel",        "HankelMatrix[Table[N[i], {i, 30}]]"),
    ("x_vandermonde",   "VandermondeMatrix[Table[N[i], {i, 30}]]"),
    ("x_cross",         "Cross[Range[1., 3.], Reverse[Range[1., 3.]]]"),
    ("x_normalize",     "Normalize[Range[1., 400.]]"),

    # --- structural: a fill or a pad is an invented element too ---
    ("x_padright_real", "PadRight[Range[1., 400.], 500, 0.]"),
    ("x_padright_dflt", "PadRight[Range[1., 400.], 500]"),
    ("x_padleft_dflt",  "PadLeft[Range[1., 400.], 500]"),
    ("x_arraypad",      "ArrayPad[Range[1., 400.], 10]"),
    ("x_riffle",        "Riffle[Range[1., 400.], 0.]"),
    ("x_riffle_exact",  "Riffle[Range[1., 400.], 0]"),
    ("x_join_mixed",    "Join[Range[1., 400.], Range[400]]"),
    ("x_replacepart",   "ReplacePart[Range[1., 400.], 1 -> 0]"),
    ("x_insert",        "Insert[Range[1., 400.], 0, 1]"),
    ("x_append",        "Append[Range[1., 400.], 0]"),
    ("x_clip_exact",    "Clip[Range[-2., 2., 0.01], {-1, 1}]"),
    ("x_clip_real",     "Clip[Range[-2., 2., 0.01], {-1., 1.}]"),

    # --- elementwise / reductions whose result type is the question ---
    ("x_floor",         "Floor[Range[0.5, 400.5]]"),
    ("x_round",         "Round[Range[0.5, 400.5]]"),
    ("x_sign",          "Sign[Range[-200., 200.]]"),
    ("x_unitstep",      "UnitStep[Range[-200., 200.]]"),
    ("x_unitbox",       "UnitBox[Range[-200., 200.]]"),
    ("x_quotient",      "Quotient[Range[1., 400.], 2]"),
    ("x_accumulate",    "Accumulate[Range[1., 400.]]"),
    ("x_differences",   "Differences[Range[1., 400.]]"),
    ("x_boole",         "Boole[Map[# > 0.5 &, Range[1., 400.]/400.]]"),
    ("x_chop",          "Chop[Range[1., 400.]/400.]"),
    ("x_n_of_int",      "N[Range[400]]"),
    ("x_re",            "Re[Range[1., 400.]]"),
    ("x_im",            "Im[Range[1., 400.]]"),
]


# ---------------------------------------------------------------------------
# Deliberate mixtures, each with the Mathematica output that justifies it.
#
# An entry here is a claim about another system, so it carries the evidence.
# Nothing goes in on the strength of "surely Mathematica does the same".
# ---------------------------------------------------------------------------
EXEMPT = {
    "padright": "PadRight[v, n] with the DEFAULT fill pads with the exact "
                "Integer 0, and Mathematica does the same: "
                "PadRight[{1., 2.}, 4] is {1., 2., 0, 0}. Giving an explicit "
                "0. fill produces a uniform Real result and packs.",
    "x_padright_dflt": "see padright",
    "x_padleft_dflt":  "see padright -- PadLeft[{1., 2.}, 4] is {0, 0, 1., 2.}",
    "arraypad": "ArrayPad with the DEFAULT fill pads with the exact Integer 0, "
                "and Mathematica does the same: ArrayPad[{1., 2., 3.}, 2] is "
                "{0, 0, 1., 2., 3., 0, 0}. An explicit real fill (ArrayPad[v, 2, "
                "0.]) gives a uniform Real result and packs. Same rule as padright.",
    "x_arraypad": "see arraypad -- ArrayPad[Range[1., 400.], 10] borders the "
                  "Real vector with exact Integer 0s, matching Mathematica.",
    "x_clip_exact": "Clip[x, {lo, hi}] returns the BOUND itself where it "
                    "clips, so an exact bound gives an exact element: "
                    "Mathematica's Clip[{-2., 0., 2.}, {-1, 1}] is "
                    "{-1, 0., 1}. Documented in packed-arrays.md.",
    "x_riffle_exact": "Riffle[v, 0] inserts the separator VERBATIM; "
                      "Mathematica's Riffle[{1., 2.}, 0] is {1., 0, 2.}.",
    "x_join_mixed": "Join concatenates without coercion, and Mathematica's "
                    "Join[{1.}, {1}] is {1., 1}. The mixture is the caller's.",
    "x_insert": "Insert places the element verbatim: Mathematica's "
                "Insert[{1., 2.}, 0, 1] is {0, 1., 2.}.",
    "x_append": "Append places the element verbatim, as Insert does.",
    "x_replacepart": "ReplacePart writes the given value verbatim: "
                     "Mathematica's ReplacePart[{1., 2.}, 1 -> 0] is {0, 2.}.",
    "tally": "Tally returns {value, count} pairs -- an inexact value beside an "
             "exact count is the shape of the answer, not a coercion bug.",
    "chop": "Chop replaces a small element with the EXACT 0, and Mathematica "
            "does the same: Chop[{1., 1.*^-12}] is {1., 0}. The mixture is the "
            "point of the operation -- 'this entry is exactly zero' -- not a "
            "lost head.",
    "lu": "LUDecomposition returns the TUPLE {lu, permutation, cond}. The "
          "exact elements are the permutation vector, which is a list of row "
          "indices: Mathematica's LUDecomposition[{{2., 1.}, {1., 3.}}] is "
          "{{{2., 1.}, {0.5, 2.5}}, {1, 2}, 3.2}. Three different things in a "
          "List, not one matrix of two heads.",
    "x_hilbert": "HilbertMatrix is exact Rationals by default; the 1 at [[1,1]] "
                 "is an Integer because 1/1 reduces. Uniformly EXACT, so it is "
                 "flagged only by the exact/exact split, not a real mixture.",
}


MPRE_SMALL = r"""
(* Head audit, not a timing run: the sizes only have to clear the 250-element
   packing threshold so the fast paths are live. *)
exq[x_] := MemberQ[{Integer, Rational}, Head[x]];
inq[x_] := SameQ[Head[x], Real];
audit[id_, r_] := Module[{f, ni, nr, tot},
  If[NDArrayQ[r],
     Print["@", id, "\t", "PACKED", "\t", DataType[r], "\t", 0, "\t", 0],
     f = Flatten[{r}];
     tot = Length[f];
     ni = Count[f, _?exq];
     nr = Count[f, _?inq];
     Print["@", id, "\t",
           Which[ni > 0 && nr > 0, "MIXED", tot == 0, "EMPTY", True, "OK"],
           "\t", tot, "\t", ni, "\t", nr]]];
"""


def m_source(probes):
    out = [numeric_sweep.MPRE, MPRE_SMALL]
    for pid, expr in probes:
        out.append(f"""
Clear[rr];
rr = {expr};
audit["{pid}", rr];
Clear[rr];
""")
    return "\n".join(out)


def parse(text):
    got = {}
    for line in text.splitlines():
        if not line.startswith("@"):
            continue
        parts = line[1:].split("\t")
        if len(parts) < 5:
            continue
        pid = parts[0].strip()
        got[pid] = {
            "state": parts[1].strip(),
            "detail": parts[2].strip(),
            "n_exact": parts[3].strip(),
            "n_inexact": parts[4].strip(),
        }
    return got


def run(probes, tmpdir):
    """Batched with a per-probe retry, so one hang costs one row not twelve."""
    results = {}

    def run_group(group):
        path = os.path.join(tmpdir, "audit.m")
        with open(path, "w") as fh:
            fh.write(m_source(group))
        try:
            cp = subprocess.run([MATHILDA, "-file", path], capture_output=True,
                                text=True, timeout=TIMEOUT)
            return parse(cp.stdout)
        except subprocess.TimeoutExpired:
            return {}

    for i in range(0, len(probes), BATCH):
        chunk = probes[i:i + BATCH]
        got = run_group(chunk)
        missing = [p for p in chunk if p[0] not in got]
        results.update(got)
        if missing and len(chunk) > 1:
            for p in missing:
                results.update(run_group([p]))
        sys.stderr.write(f"  {min(i + BATCH, len(probes))}/{len(probes)}\n")
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true",
                    help="include OK and PACKED rows, not just MIXED")
    ap.add_argument("--json", metavar="PATH")
    ap.add_argument("--only", help="comma-separated probe ids")
    args = ap.parse_args()

    probes = [(p["id"], p["m"]) for p in numeric_sweep.PROBES] + EXTRA
    if args.only:
        want = set(args.only.split(","))
        probes = [p for p in probes if p[0] in want]

    with tempfile.TemporaryDirectory() as td:
        got = run(probes, td)

    rows = []
    for pid, expr in probes:
        r = got.get(pid)
        rows.append({
            "id": pid, "expr": expr,
            "state": r["state"] if r else "NORESULT",
            "n": r["detail"] if r else "",
            "n_exact": r["n_exact"] if r else "",
            "n_inexact": r["n_inexact"] if r else "",
            "exempt": EXEMPT.get(pid),
        })

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(rows, fh, indent=1)

    mixed = [r for r in rows if r["state"] == "MIXED" and not r["exempt"]]
    excused = [r for r in rows if r["state"] == "MIXED" and r["exempt"]]
    show = rows if args.all else mixed

    if show:
        print(f"{'probe':<26} {'state':<9} {'exact':>7} {'inexact':>8}  expression")
        print("-" * 100)
        for r in sorted(show, key=lambda x: (x["state"] != "MIXED", x["id"])):
            print(f"{r['id']:<26} {r['state']:<9} {r['n_exact']:>7} "
                  f"{r['n_inexact']:>8}  {r['expr'][:52]}")
        print()

    packed = sum(1 for r in rows if r["state"] == "PACKED")
    ok = sum(1 for r in rows if r["state"] in ("OK", "EMPTY"))
    none = sum(1 for r in rows if r["state"] == "NORESULT")
    print(f"{len(rows)} probes: {packed} packed, {ok} single-headed, "
          f"{len(mixed)} MIXED, {len(excused)} mixed-but-exempt, {none} no result")
    if mixed:
        print("\nA two-headed result from a machine input is a bug until a "
              "wolframscript run says otherwise.\nCheck each, then fix it or "
              "add an EXEMPT entry carrying the Mathematica output.")
    return 1 if mixed else 0


if __name__ == "__main__":
    sys.exit(main())

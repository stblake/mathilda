#!/usr/bin/env python3
"""
compile_coverage.py -- does every head with a numeric fast path also compile?

THE INVARIANT THIS CHECKS. A head earns a numeric fast path by being something
a numeric workload runs over machine numbers. Compile[] and auto-compilation
exist to run exactly those workloads. So a head that has a buffer path in the
interpreter and no lowering in the compiler is a contradiction: the same
expression is fast when typed at the REPL and slow inside the Compile[] that
was supposed to speed it up.

It is worse than "no speedup", because the compilable subset is a CLIFF, not a
gradient (tasks/lessons.md, project_compile_subset_is_a_cliff). One unlowerable
head anywhere in the body makes the WHOLE body fall back to the interpreter,
so `Compile[{{v, _Real, 1}}, Total[v] / Mean[v]]` loses the compiled Total as
well as the uncompiled Mean, and says nothing about it. The failure is silent
by design: same answer, 10-40x slower, no symptom.

WHAT IT DOES. Joins three registries against the compiler:

    src/ndkernels.c    heads with an elementwise NDArray kernel
    src/pack.c AWARE   heads the transparency gate does not materialise for
    src/compile/       CompileDiagnostics[argspec, body], asked of the binary

For each head it tries the typed argument specs that head could plausibly take
and keeps the best answer. `CompileDiagnostics` is the system's own report --
it returns Compiled -> False with the innermost subexpression that could not be
lowered -- so this asks the compiler rather than guessing from the source.

Reading the output: `scalar` and `array` are separate columns because they are
separate lowerings, reached by different code. Sin has a scalar opcode AND an
array kernel; Mean had neither; Sort has an array delegation and no scalar form
(correctly -- Sort of a scalar is not a thing). A head that is `-` in a column
where its interpreter fast path is `y` is the finding. `array` is one family
across every rank and spelling on purpose: the question is whether the head has
a machine lowering at all, so "no rank-1 Inverse" would be noise.

Usage:
    python3 tools/compile_coverage.py              # gaps only
    python3 tools/compile_coverage.py --all        # every head
    python3 tools/compile_coverage.py --json out.json
    python3 tools/compile_coverage.py --only Mean,Median
"""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
MATHILDA = os.environ.get("MATHILDA_BIN", os.path.join(ROOT, "Mathilda"))

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numeric_coverage as nc          # noqa: E402  (registry readers)


# ---------------------------------------------------------------------------
# The typed shapes a compiled body can be asked for.
#
# Each entry is (id, argspec, body-template). The template's %s is the head.
# A head "compiles at" a shape if CompileDiagnostics says Compiled -> True for
# it; the report keeps the shape family (scalar / array) rather than the
# individual spec, because that is the granularity a gap is fixed at.
# ---------------------------------------------------------------------------

SHAPES = [
    ("s_r",   "scalar", "{{x, _Real}}",                    "%s[x]"),
    ("s_i",   "scalar", "{{x, _Integer}}",                 "%s[x]"),
    ("s_rr",  "scalar", "{{x, _Real}, {y, _Real}}",        "%s[x, y]"),
    ("s_ii",  "scalar", "{{x, _Integer}, {y, _Integer}}",  "%s[x, y]"),
    ("v_r",   "array",  "{{v, _Real, 1}}",                 "%s[v]"),
    ("v_i",   "array",  "{{v, _Integer, 1}}",              "%s[v]"),
    ("v_rk",  "array",  "{{v, _Real, 1}, {k, _Integer}}",  "%s[v, k]"),
    # An INTEGER array with an integer scalar.  The integer-only kernels (GCD,
    # LCM, DivisorSigma, IntegerLength) are defined on Z and correctly decline
    # a real buffer, so probing them only at {_Real, 1} reports a gap that is
    # the probe's, not the head's.
    ("v_ik",  "array",  "{{v, _Integer, 1}, {k, _Integer}}", "%s[v, k]"),
    ("v_rx",  "array",  "{{v, _Real, 1}, {y, _Real}}",     "%s[v, y]"),
    ("v_xr",  "array",  "{{v, _Real, 1}, {y, _Real}}",     "%s[y, v]"),
    ("v_rr",  "array",  "{{v, _Real, 1}, {w, _Real, 1}}",  "%s[v, w]"),
    ("m_r",   "array",  "{{m, _Real, 2}}",                 "%s[m]"),
    ("m_rr",  "array",  "{{m, _Real, 2}, {p, _Real, 2}}",  "%s[m, p]"),
    ("m_rv",  "array",  "{{m, _Real, 2}, {v, _Real, 1}}",  "%s[m, v]"),
    # A matrix with an integer scalar.  MatrixPower[m, k] has no valid 1-arg
    # spelling (MatrixPower[m] is unevaluated), so without this shape the tool
    # probes MatrixPower[m] and reports the probe's gap, not the head's — the
    # same failure the functional spellings below fix for Map / Fold.
    ("m_rk",  "array",  "{{m, _Real, 2}, {k, _Integer}}",  "%s[m, k]"),
    # Functional spellings.  Leaving these out is what made the first run of
    # this tool report Map, Select, Scan, Fold and Nest as gaps: every one of
    # them compiles, just not as H[v].  A probe list that cannot express a
    # head's real call shape reports the PROBE's gap as the head's.
    ("f_v",   "array",  "{{v, _Real, 1}}",                 "%s[2. # &, v]"),
    ("v_f",   "array",  "{{v, _Real, 1}}",                 "%s[v, # > 0.5 &]"),
    ("p_v",   "array",  "{{v, _Real, 1}}",                 "%s[# > 0.5 &, v]"),
    ("f_x_n", "array",  "{{x, _Real}, {k, _Integer}}",     "%s[2. # &, x, k]"),
    ("f_z_v", "array",  "{{v, _Real, 1}}",                 "%s[#1 + #2 &, 0., v]"),
]

FAMILIES = ["scalar", "array"]

PROBE = """
d = CompileDiagnostics[%(spec)s, %(body)s];
Print["@", "%(id)s", "\\t", Lookup[d, "Compiled", False], "\\t",
      Lookup[d, "Reason", ""], "\\t", Lookup[d, "Subexpression", ""]];
"""


def probe_source(items):
    return "\n".join(PROBE % {"id": i, "spec": s, "body": b}
                     for i, s, b in items)


def parse(text):
    out = {}
    for line in text.splitlines():
        if not line.startswith("@"):
            continue
        parts = line[1:].split("\t")
        if len(parts) < 2:
            continue
        out[parts[0].strip()] = {
            "compiled": parts[1].strip() == "True",
            "reason": parts[2].strip() if len(parts) > 2 else "",
            "sub": parts[3].strip() if len(parts) > 3 else "",
        }
    return out


def run(items, batch=40, timeout=90):
    got = {}
    with tempfile.TemporaryDirectory() as td:
        for i in range(0, len(items), batch):
            chunk = items[i:i + batch]
            fd, path = tempfile.mkstemp(suffix=".m", dir=td)
            with os.fdopen(fd, "w") as fh:
                fh.write(probe_source(chunk))
            try:
                p = subprocess.run([MATHILDA, "-file", path],
                                   capture_output=True, text=True,
                                   timeout=timeout)
                got.update(parse(p.stdout))
            except subprocess.TimeoutExpired:
                pass
            finally:
                os.unlink(path)
            sys.stderr.write(f"\r  probing {min(i + batch, len(items))}/"
                             f"{len(items)}   ")
            sys.stderr.flush()
    sys.stderr.write("\n")
    return got


# ---------------------------------------------------------------------------
# Which interpreter fast path does this head have?  Both halves matter: an
# elementwise kernel is what a SCALAR body needs lowered, and an AWARE buffer
# path is what an ARRAY body needs delegated.  A head can have one and not the
# other, and the two gaps are fixed in different places.
# ---------------------------------------------------------------------------

def fastpath_inventory():
    kernels = nc.kernel_heads()
    packs = nc.pack_lists()
    builtins = nc.all_builtins()
    inv = {}
    for h in sorted(set(kernels) | packs["aware"] | packs["int64_ok"]):
        if h not in builtins:
            continue                    # a kernel name that is not a head
        inv[h] = {
            "kernel": sorted(kernels.get(h, [])),
            "aware": h in packs["aware"],
            "int64": h in packs["int64_ok"],
            "file": builtins[h],
        }
    return inv


# Heads for which a compiled lowering is meaningless rather than missing, with
# the reason.  This is the file's only hand-maintained judgement, and it is
# deliberately short: everything else is derived, so a head that stops being
# exempt shows up as a gap by itself rather than staying quietly excused.
EXEMPT = {
    "Set": "assignment; the compiler has its own OP_MOVE lowering for locals",
    "SetDelayed": "assignment",
    "CompoundExpression": "lowered as a sequence, not a callable head",
    "If": "lowered as control flow (K_JZ)",
    "Which": "control flow",
    "Print": "side effect; no value to type",
    "Echo": "side effect",
    "Return": "control flow",
    "Module": "scoping; lowered as register allocation",
    "Block": "scoping",
    "With": "scoping; lowered by substitution",
    "N": "numericalisation is a no-op on a machine value",
    "TeXForm": "formatting",
    "Normal": "representation conversion",
    "ToNDArray": "representation conversion",
    "FromNDArray": "representation conversion",
    "ToPackedArray": "representation conversion",
    "FromPackedArray": "representation conversion",
    "NDArrayQ": "asks about the REPRESENTATION, which the compiler erases",
    "DataType": "asks about the representation",
    "ByteCount": "asks about the representation",
    "Head": "asks about the representation",
    "AtomQ": "asks about the representation",
    "Precision": "asks about the representation",
    "Accuracy": "asks about the representation",
    "MatrixQ": "shape test; the compiler already knows the rank statically",
    "VectorQ": "shape test",
    "Depth": "shape test",
    "Dimensions": "shape test; A_SIZE covers the scalar Length case",
    "RandomSample": "draws from the global generator; not a pure kernel",
    "RandomChoice": "draws from the global generator",
    "Counts": "returns an Association, which has no compiled type",
    "Tally": "returns ragged {value, count} pairs; no uniform compiled type",
    "Commonest": "count-keyed selection; result length is data-dependent",
    "Union": "result length is data-dependent",
    "Intersection": "result length is data-dependent",
    "Complement": "result length is data-dependent",
    "DeleteDuplicates": "result length is data-dependent",
    "Position": "result length is data-dependent",
    "LatticeReduce": "exact integer lattice; no machine type",
    "FindIntegerNullVector": "exact integer",
    "Prime": "exact integer table lookup, unbounded result",
    "IntegerDigits": "result length is data-dependent",
    "QuotientRemainder": "returns a pair, not a scalar",
    "MinMax": "returns a pair",
    "Quartiles": "returns a triple",
    "LUDecomposition": "returns a heterogeneous triple",
    "QRDecomposition": "returns a pair of matrices",
    "JordanDecomposition": "returns a pair of matrices",
    "SchurDecomposition": "returns a pair/tuple of matrices",
    "SingularValueDecomposition": "returns a triple of matrices",
    "Eigenvectors": "result shape is data-dependent (defective matrices)",
    "Eigenvalues": "may be complex for a real matrix; result type not static",
    "CharacteristicPolynomial": "returns a symbolic Plus in x, not a machine array",
    "PositiveDefiniteMatrixQ": "answers a Boolean; the matrix is the argument",
    "NegativeDefiniteMatrixQ": "answers a Boolean; the matrix is the argument",
    "HypergeometricPFQ": "parameter LISTS, not machine arrays",
    "Hypergeometric0F1": "rewrites to HypergeometricPFQ",
    "Hypergeometric1F1": "rewrites to HypergeometricPFQ",
    "Hypergeometric2F1": "rewrites to HypergeometricPFQ",
    "Interpolation": "returns an InterpolatingFunction object (no compiled machine "
                     "type); AWARE only to read a packed data table straight through",
    "ListInterpolation": "returns an InterpolatingFunction object (no compiled "
                         "machine type); AWARE only to read a packed value tensor",
    "InterpolatingPolynomial": "returns a symbolic polynomial in x, not a machine "
                               "array (cf. CharacteristicPolynomial)",
}


# ---------------------------------------------------------------------------
# The RATCHET.
#
# 52 heads are open as of 2026-08-02. A gate that fails on all of them from the
# day it lands is not a gate — nobody can tell a new regression from the
# standing backlog. So the exit status keys on the DIFFERENCE: a head that is a
# gap and is not listed here fails the build, and a head listed here that has
# started compiling is reported so the line can be deleted.
#
# The list is grouped by what closing it needs, because the groups close
# together rather than one head at a time. COMPILE_MISSING.md at the repo root
# is the same list written out, with the interpreter entry point to delegate to
# where one exists -- read that before picking one up.
# ---------------------------------------------------------------------------

BASELINE = {
    # Two ARRAY operands whose extra operand is a FUNCTION head, not a plain
    # array — Inner[f, a, b, g] / Outer[f, a, b].  The A_NDFN2 opcode (added with
    # COMPILE_MISSING.md §2/§3) carries two array registers; these also want the
    # callback lowering Map / Fold have, so they wait on that rather than on the
    # opcode.  (Dot, Cross, LinearSolve, LeastSquares, ListConvolve/Correlate and
    # Join now lower via A_NDFN2 / V_NDFN2; Inverse, Normalize, MatrixPower,
    # ReverseSort, ConjugateTranspose, PseudoInverse via ND_FNS / A_NDFN.  Tr,
    # Det, MatrixRank, Norm were the §1 rank-2 -> scalar reductions.)
    "Inner", "Outer",
    # No buffer fast path in the interpreter EITHER (COMPILE_MISSING.md §9): both
    # materialise on the first line and run the exact fraction-free path, so a
    # compiled form would delegate to a fast path that does not exist.  A LAPACK
    # path comes first, then a lowering.
    "RowReduce", "NullSpace",
    # ND_FNS entries whose extra argument is not an integer (a list, a real, or a
    # scalar separator), so the table's shape does not reach them yet.  Riffle is
    # array + SEPARATOR (a scalar or list), not two arrays, so A_NDFN2 does not
    # fit it — it belongs with this group, not §3.
    "Riffle", "Partition", "PadLeft", "PadRight", "Append", "Prepend",
    "Catenate", "Extract", "Rescale",
    "ExponentialMovingAverage",
    # ArrayPad amounts are a scalar/list spec and the fill a value/scheme; the
    # ArrayReshape dims argument is a List spec — neither is a CT_INT operand the
    # A_NDFN table can carry.
    "ArrayPad", "ArrayReshape",
    # A callback lowering, as Map / Select / Fold already have.
    "Scan", "MapAll", "MapAt", "MapIndexed", "MapThread", "SelectFirst",
    "NestWhile", "NestWhileList",
    # SCALAR forms of the integer-only kernels.  Their ARRAY forms landed on
    # 2026-08-02; the scalar side wants an int -> int kernel opcode that can
    # ABORT (EulerPhi of a negative has no answer), which the K_UN opcodes
    # cannot currently do — they write NaN.
    "MoebiusMu", "EulerPhi", "DivisorSigma", "PowerMod",
    # Ternary/n-ary special function with a machine kernel and no lowering.
    "LerchPhi",
    # The elementwise sign predicates.  Now that NDT_BOOL exists (§13 gap C.1
    # closed) they produce a packed bool array, so a compiled ARRAY form is
    # finally possible — an elementwise compare-to-zero writing CT_BOOL elements
    # (A_STORE_B).  No such lowering yet: they moved here from EXEMPT because the
    # dtype that blocked them is no longer missing, only the lowering is.
    "Positive", "Negative", "NonNegative", "NonPositive",
    # Boole: the CONVERSE bridge, a bool array -> an int64 one.  Joined pack.c's
    # AWARE with the NDT_BOOL work; a compiled array form (bool -> A_STORE_I) is
    # possible but likewise undone.  Its real call shape is Boole[boolArray],
    # which the probe list cannot express (no _Boolean array shape), so the tool
    # reports it at the real/int shapes it can — either way, no lowering yet.
    "Boole",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--only", default="")
    ap.add_argument("--json", default="")
    args = ap.parse_args()

    if not os.path.exists(MATHILDA):
        sys.exit(f"no binary at {MATHILDA}; run make first")

    inv = fastpath_inventory()
    if args.only:
        want = {h.strip() for h in args.only.split(",") if h.strip()}
        inv = {h: v for h, v in inv.items() if h in want}

    items = [(f"{h}/{sid}", spec, tmpl % h)
             for h in inv for sid, fam, spec, tmpl in SHAPES]
    sys.stderr.write(f"asking CompileDiagnostics {len(items)} questions "
                     f"over {len(inv)} heads\n")
    got = run(items)

    rows = []
    for h, v in sorted(inv.items()):
        fam_ok = {f: False for f in FAMILIES}
        reasons = {}
        for sid, fam, spec, tmpl in SHAPES:
            r = got.get(f"{h}/{sid}")
            if not r:
                continue
            if r["compiled"]:
                fam_ok[fam] = True
            elif fam not in reasons and r["sub"]:
                reasons[fam] = r["reason"]
        # What SHOULD compile.  An elementwise kernel is a SCALAR lowering and
        # an AWARE buffer path is an ARRAY one, and a head can want either or
        # both.  "array" is deliberately one family across every rank and
        # spelling: the question is whether the head has a machine lowering at
        # all, and reporting "no rank-1 Inverse" would be noise, not a gap.
        want_scalar = bool(v["kernel"])
        want_array = v["aware"] or bool(v["kernel"])
        gap = []
        if want_scalar and not fam_ok["scalar"]:
            gap.append("scalar")
        if want_array and not fam_ok["array"]:
            gap.append("array")
        rows.append(dict(head=h, file=v["file"], kernel=",".join(v["kernel"]),
                         aware=v["aware"], int64=v["int64"],
                         scalar=fam_ok["scalar"], array=fam_ok["array"],
                         gap=gap,
                         reason=reasons.get("array") or reasons.get("scalar", ""),
                         exempt=EXEMPT.get(h, "")))

    gaps = [r for r in rows if r["gap"] and not r["exempt"]]

    show = rows if args.all else gaps
    print(f"\n{'head':<26} {'kernel':<14} {'aware':>5} {'scalar':>7} "
          f"{'array':>6}  gap / reason")
    print("-" * 100)
    for r in sorted(show, key=lambda r: (r["head"])):
        print(f"{r['head']:<26} {r['kernel']:<14} "
              f"{'y' if r['aware'] else '-':>5} "
              f"{'y' if r['scalar'] else '-':>7} "
              f"{'y' if r['array'] else '-':>6}  {','.join(r['gap'])}"
              f"{('  [' + r['exempt'] + ']') if r['exempt'] else ''}")

    nex = sum(1 for r in rows if r["exempt"] and r["gap"])
    print(f"\n{len(inv)} heads with a numeric fast path; "
          f"{len(gaps)} do not compile at a shape they should; "
          f"{nex} exempt with a reason")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(rows, fh, indent=1)
        print(f"wrote {args.json}")

    # The ratchet.  Only a head that is newly unlowered fails; the standing
    # backlog is BASELINE and is reported, not enforced.  A --only run cannot
    # see the whole set, so it never fails on absences.
    names = {r["head"] for r in gaps}
    new = sorted(names - BASELINE)
    if new:
        print("\nNEW: these heads have a numeric fast path and no compiled "
              "lowering, and are not in BASELINE:")
        for h in new:
            print(f"  {h}")
        print("Add the lowering, or add the head to BASELINE / EXEMPT in "
              f"{os.path.relpath(__file__, ROOT)} with the reason.")
        return 1
    if not args.only:
        fixed = sorted(BASELINE - names - set(EXEMPT))
        if fixed:
            print("\nFIXED since the baseline was recorded — delete these "
                  "from BASELINE:")
            for h in fixed:
                print(f"  {h}")
        print(f"\nOK: no new gaps ({len(names & BASELINE)} known open, "
              f"tracked in BASELINE).")
    return 0


if __name__ == "__main__":
    sys.exit(main())

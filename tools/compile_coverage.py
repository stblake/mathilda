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
    # Three and four arguments: the n-ary kernels (LerchPhi, PolyLog[n, p, z],
    # WhittakerM/W, HypergeometricU, 1F1, 2F1, EllipticPi, ...). The array shapes
    # put the array FIRST, in the MIDDLE and LAST -- the compiled n-ary lowering
    # takes it in any one position (VKERNN's position byte), and a head that
    # threads over only some of them still shows as lowered.
    ("s_rrr",  "scalar", "{{x, _Real}, {y, _Real}, {w, _Real}}", "%s[x, y, w]"),
    ("s_rrrr", "scalar", "{{x, _Real}, {y, _Real}, {w, _Real}, {u, _Real}}",
     "%s[x, y, w, u]"),
    ("v_xxv",  "array",  "{{v, _Real, 1}, {x, _Real}, {y, _Real}}", "%s[x, y, v]"),
    ("v_vxx",  "array",  "{{v, _Real, 1}, {x, _Real}, {y, _Real}}", "%s[v, x, y]"),
    ("v_xvx",  "array",  "{{v, _Real, 1}, {x, _Real}, {y, _Real}}", "%s[x, v, y]"),
    ("v_xxxv", "array",  "{{v, _Real, 1}, {x, _Real}, {y, _Real}, {w, _Real}}",
     "%s[x, y, w, v]"),
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

# Heads whose real call shape no generic template can spell -- literal
# parameter LISTS around the numeric argument -- get their own probes, so the
# tool enforces their lowering instead of carrying them in BASELINE.
HEAD_SHAPES = {
    # PowerMod[x, e, m]: ternary, so no generic template spells it. The scalar
    # form lowers to OP_POWMOD_I, the array form delegates to ndint_powermod
    # through A_NDFN (two trailing integers).
    "PowerMod": [
        ("pm_s", "scalar", "{{x, _Integer}}",    "PowerMod[x, 5, 1009]"),
        ("pm_v", "array",  "{{v, _Integer, 1}}", "PowerMod[v, 5, 1009]"),
    ],
    # DirichletL[k, j, s]: integer labels broadcast, s scalar or a tile; a real
    # character (5, 3) types Real, a complex one (5, 2) Complex.
    "DirichletL": [
        ("dl_s",  "scalar", "{{x, _Real}}",    "DirichletL[5, 3, x]"),
        ("dl_sc", "scalar", "{{x, _Real}}",    "DirichletL[5, 2, x]"),
        ("dl_v",  "array",  "{{v, _Real, 1}}", "DirichletL[5, 3, v]"),
        ("dl_vc", "array",  "{{v, _Real, 1}}", "DirichletL[5, 2, v]"),
    ],
    # Boole: the bool -> int64 bridge.  Its real call shape is Boole[boolArray],
    # which no _Real/_Integer template spells, so it gets a _Boolean array
    # probe (and the composed Boole[Positive[v]] / scalar Boole[x > 0.5]).
    "Boole": [
        ("bo_s",  "scalar", "{{x, _Real}}",       "Boole[x > 0.5]"),
        ("bo_v",  "array",  "{{b, _Boolean, 1}}", "Boole[b]"),
        ("bo_pv", "array",  "{{v, _Real, 1}}",    "Boole[Positive[v]]"),
    ],
    "MeijerG": [
        ("mg_s",  "scalar", "{{x, _Real}}",
         "MeijerG[{{1/3}, {1/4}}, {{1/5, 1/6}, {1/7}}, x]"),
        ("mg_sr", "scalar", "{{x, _Real}}",
         "MeijerG[{{1/3}, {}}, {{1/2}, {1}}, x, 1/2]"),
        ("mg_v",  "array",  "{{v, _Real, 1}}",
         "MeijerG[{{1/3}, {1/4}}, {{1/5, 1/6}, {1/7}}, v]"),
    ],
}

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
    "FromCharacterCode": "reads an int64 code buffer but returns a String, "
                         "which no compiled register holds",
    "FromDigits": "reads an int64 digit buffer, but the answer is an exact "
                  "integer that outgrows a machine register past ~18 digits; "
                  "the buffer path exists for long digit lists",
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
    "DirichletCharacter": "exact roots of unity (0, +-1, unit*E^(I Pi r)); not a machine value",
    "FunctionExpand": "symbolic rewriter of special-function nodes; never computes a number",
    "IntegerDigits": "result length is data-dependent",
    "QuotientRemainder": "returns a pair, not a scalar",
    "MinMax": "returns a pair",
    "Quartiles": "returns a triple",
    "LUDecomposition": "returns a heterogeneous triple",
    "LUBackSubstitution": "consumes LUDecomposition's heterogeneous triple {lu, p, c}, which no compiled register carries",
    "QRDecomposition": "returns a pair of matrices",
    "JordanDecomposition": "returns a pair of matrices",
    "SchurDecomposition": "returns a pair/tuple of matrices",
    "SingularValueDecomposition": "returns a triple of matrices",
    "Eigenvectors": "result shape is data-dependent (defective matrices)",
    "Eigensystem": "returns the ragged pair {values, vectors} (a vector and a "
                   "matrix, complex or real by data); no single machine array",
    "Eigenvalues": "may be complex for a real matrix; result type not static",
    "MatrixFunction": "its first argument is a function value (f), which no compiled register or ND table operand carries; MatrixExp / MatrixLog / MatrixPower lower instead",
    "CharacteristicPolynomial": "returns a symbolic Plus in x, not a machine array",
    "PositiveDefiniteMatrixQ": "answers a Boolean; the matrix is the argument",
    "NegativeDefiniteMatrixQ": "answers a Boolean; the matrix is the argument",
    "HypergeometricPFQ": "parameter LISTS, not machine arrays (the scalar form with "
                         "literal Lists does lower, via OP_KERNN; the array form does not)",
    "Interpolation": "returns an InterpolatingFunction object (no compiled machine "
                     "type); AWARE only to read a packed data table straight through",
    "ListInterpolation": "returns an InterpolatingFunction object (no compiled "
                         "machine type); AWARE only to read a packed value tensor",
    "InterpolatingPolynomial": "returns a symbolic polynomial in x, not a machine "
                               "array (cf. CharacteristicPolynomial)",
    "FunctionInterpolation": "HoldAll; samples a held expression into an "
                             "InterpolatingFunction object (no compiled machine type)",
    "GaussianQuadratureWeights": "builds an n x 2 rule from an integer order; no "
                                 "element-wise numeric argument to lower",
    "GaussianQuadratureError": "returns a symbolic error term, not a number",
    # Image / graphics heads.  An Image (or Graphics) object has no compiled
    # machine type, so a head that constructs one, queries one, or takes one as
    # its argument cannot lower -- the AWARE entry only lets the transparency
    # gate read the image's packed pixel buffer straight through.
    "Image": "constructs an Image object; no compiled image type",
    "Image3D": "constructs an Image3D object; no compiled image type",
    "ImageQ": "predicate on an Image; the Image is the argument",
    "Image3DQ": "predicate on an Image3D; the Image is the argument",
    "ImageChannels": "reads Image metadata (channel count); the Image is the argument",
    "ImageDimensions": "reads Image metadata (dimensions); the Image is the argument",
    "ImageType": "reads an Image's pixel type (a string); the Image is the argument",
    "ImageData": "extracts an Image's pixel array; the Image is the argument",
    "ImageResize": "returns an Image; no compiled image type",
    "ImageConvolve": "returns an Image; no compiled image type",
    "Binarize": "returns a binarized Image; no compiled image type",
    "ColorConvert": "operates on an Image/color; no compiled image type",
    "GaussianFilter": "image-pipeline filter; declines a plain numeric list, "
                      "needs an Image",
    "ArrayPlot": "returns a Graphics object; no compiled type",
    "FindThreshold": "computes a threshold FROM an Image; the Image is the argument",
    # Representation predicate (cf. NDArrayQ).
    "PackedArrayQ": "asks about the REPRESENTATION, which the compiler erases",
    # Fitting / structural results that are not a machine array.
    "Fit": "returns a symbolic Plus in the fit variable, not a machine array "
           "(cf. CharacteristicPolynomial)",
    "DesignMatrix": "basis is a symbolic function list and a variable, not "
                    "compilable argument types",
    "HermiteDecomposition": "returns a pair of (exact-integer) matrices "
                            "(cf. QRDecomposition)",
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
    # PrimePowerQ: the SCALAR form compiles (OP_PPQ_I, int -> CT_BOOL, the
    # interpreter's own mth_prime_power_base_u64), so loops over PrimePowerQ[k]
    # and Table[PrimePowerQ[v[[i]]], ...] stay compiled.  The direct ARRAY form
    # PrimePowerQ[v] wants an int -> bool elementwise array lowering (A_STORE_B).
    # That lowering now exists -- emit_elementwise_over_array, which the sign
    # predicates and Boole use -- and wiring a unary predicate to it is the
    # remaining step (as for EvenQ and PrimeQ below).
    "PrimePowerQ",
    # CoprimeQ and Divisible (v0.471): the SCALAR forms compile (Divisible as
    # MOD_I == 0, CoprimeQ[a, b] as GCD_I == 1), and their packed buffers answer
    # straight off int64 with a packed bool List (ndint_binary_predicate).  The
    # direct ARRAY form is the int -> bool elementwise lowering PrimePowerQ
    # wants, but BINARY (array op scalar), which emit_elementwise_over_array
    # does not yet spell.
    "CoprimeQ", "Divisible",
    # PrimeQ (v0.472): same story -- the scalar form compiles (OP_PRIMEQ_I,
    # the interpreter's mth_is_prime_u64), the packed buffer answers off int64
    # (ndint_primeq), the direct array form wants the bool elementwise lowering.
    "PrimeQ",
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
    items += [(f"{h}/{sid}", spec, body)
              for h in inv for sid, fam, spec, body in HEAD_SHAPES.get(h, [])]
    sys.stderr.write(f"asking CompileDiagnostics {len(items)} questions "
                     f"over {len(inv)} heads\n")
    got = run(items)

    rows = []
    for h, v in sorted(inv.items()):
        fam_ok = {f: False for f in FAMILIES}
        reasons = {}
        for sid, fam, spec, tmpl in SHAPES + HEAD_SHAPES.get(h, []):
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

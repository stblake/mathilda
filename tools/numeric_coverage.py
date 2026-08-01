#!/usr/bin/env python3
"""
numeric_coverage.py -- which builtins have a numeric fast path, and which do not?

THE QUESTION THIS ANSWERS. The nineteen experiments in docs/experiments/ were
each driven by a *workload*: pick a kernel, measure it, find why it is slow. That
finds what the workloads happen to touch and is silent about everything else --
and the recurring finding across all nineteen was an operation that "had a
working buffer path and a working list path, and quietly took the second".
A workload-driven sweep can only find those one accident at a time.

This is the coverage-driven complement. It enumerates EVERY builtin the system
registers and joins it against the three registries that decide whether a numeric
argument reaches machine code:

    src/ndkernels.c     the elementwise kernel registry (REG_U / REG_B / REG_N)
    src/pack.c AWARE    heads the transparency gate does NOT materialise for
    src/pack.c INT64_OK heads whose buffer path is exact on an int64 buffer

A head that is in none of them and that a numeric workload can reach is, by
construction, running at one Expr node per element.

The static half cannot tell "no fast path" from "no fast path NEEDED" -- Names[]
and OpenWrite have no business on a float64 buffer. That judgement is the
CATEGORY table below, which is hand-assigned and is the part of this file worth
reviewing. Everything else is derived.

Usage:
    python3 tools/numeric_coverage.py                 # summary table
    python3 tools/numeric_coverage.py --json out.json # machine-readable
    python3 tools/numeric_coverage.py --gaps          # only the unbacked heads
"""
import argparse
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")


# ---------------------------------------------------------------------------
# Source extraction
# ---------------------------------------------------------------------------

def _walk_c():
    for dirpath, dirnames, filenames in os.walk(SRC):
        # src/external is third-party and off limits (CLAUDE.md).
        dirnames[:] = [d for d in dirnames if d != "external"]
        for fn in filenames:
            if fn.endswith((".c", ".h")):
                yield os.path.join(dirpath, fn)


def all_builtins():
    """Every name passed to symtab_add_builtin, with the file that registers it."""
    pat = re.compile(r'symtab_add_builtin\(\s*"([^"]+)"')
    out = {}
    for path in _walk_c():
        rel = os.path.relpath(path, ROOT)
        with open(path, encoding="utf-8", errors="replace") as fh:
            for m in pat.finditer(fh.read()):
                out.setdefault(m.group(1), rel)
    return out


def _string_list(text, name):
    """Pull the string literals out of `static const char* const NAME[] = {...}`."""
    m = re.search(r'\b' + re.escape(name) + r'\s*\[\s*\]\s*=\s*\{', text)
    if not m:
        return set()
    i = m.end()
    depth = 1
    start = i
    while i < len(text) and depth:
        if text[i] == '{':
            depth += 1
        elif text[i] == '}':
            depth -= 1
        i += 1
    body = text[start:i - 1]
    # Strip comments first: the lists carry long prose blocks that quote head
    # names, and those must not be read as entries.
    body = re.sub(r'/\*.*?\*/', '', body, flags=re.S)
    body = re.sub(r'//[^\n]*', '', body)
    return set(re.findall(r'"([^"]+)"', body))


def pack_lists():
    text = open(os.path.join(SRC, "pack.c"), encoding="utf-8").read()
    return {
        "aware": _string_list(text, "AWARE"),
        "not_aware": _string_list(text, "NOT_AWARE"),
        "int64_ok": _string_list(text, "INT64_OK"),
        "broadcast_ok": _string_list(text, "BROADCAST_OK"),
    }


def kernel_heads():
    """Heads with a registered elementwise NDArray kernel, by arity."""
    out = {}
    for path in _walk_c():
        with open(path, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
        for arity, macro in (("unary", "REG_U"), ("binary", "REG_B"), ("nary", "REG_N")):
            for m in re.finditer(r'\b' + macro + r'\(\s*([A-Za-z0-9_]+)', text):
                out.setdefault(m.group(1), set()).add(arity)
        for m in re.finditer(
                r'symtab_set_ndarray_(unary|binary|nary)_kernel\(\s*"([^"]+)"', text):
            out.setdefault(m.group(2), set()).add(m.group(1))
    return {k: sorted(v) for k, v in out.items()}


# ---------------------------------------------------------------------------
# The judgement call: does this head have any business on a numeric buffer?
#
# "numeric" -- a real workload can hand it a machine-number array, and the
#              answer should come back at machine speed. These are the heads
#              this sweep is about.
# "scalar"  -- numeric, but inherently one number at a time (RandomInteger's
#              scalar form, $MachineEpsilon). Cheap either way; listed so the
#              inventory is complete rather than because there is work here.
# "symbolic"-- operates on expression structure, exact arithmetic or the
#              symbolic algebra. A buffer path is meaningless or wrong.
# "io"      -- files, streams, sessions, graphics, formatting.
#
# Anything not named here defaults to "symbolic", which is the conservative
# direction: it keeps the head OUT of the gap list, so a wrong default under-
# reports rather than inventing work. The --unclassified flag lists them.
# ---------------------------------------------------------------------------

NUMERIC = {
    # --- elementwise transcendental / arithmetic -------------------------
    "Sin", "Cos", "Tan", "Cot", "Sec", "Csc",
    "Sinh", "Cosh", "Tanh", "Coth", "Sech", "Csch",
    "ArcSin", "ArcCos", "ArcTan", "ArcCot", "ArcSec", "ArcCsc",
    "ArcSinh", "ArcCosh", "ArcTanh", "ArcCoth", "ArcSech", "ArcCsch",
    "Exp", "Log", "Log2", "Log10", "Sqrt", "CubeRoot", "Power", "Surd",
    "Plus", "Times", "Subtract", "Divide", "Minus",
    "Abs", "Sign", "UnitStep", "Ramp", "Clip", "Rescale", "Chop", "Threshold",
    "Re", "Im", "Arg", "Conjugate",
    "Floor", "Ceiling", "Round", "IntegerPart", "FractionalPart",
    "Mod", "Quotient", "QuotientRemainder", 
    "Boole", "Unitize", "KroneckerDelta", "DiscreteDelta",
    "Max", "Min", "Gudermannian", "InverseGudermannian",
    "Haversine", "InverseHaversine", "Sinc", "LogisticSigmoid",

    # --- special functions ------------------------------------------------
    "Gamma", "LogGamma", "PolyGamma", "Beta", "Pochhammer", "Factorial",
    "Factorial2", "Binomial", "Multinomial", "HarmonicNumber",
    "Erf", "Erfc", "Erfi", "InverseErf", "InverseErfc",
    "ExpIntegralE", "ExpIntegralEi", "LogIntegral",
    "SinIntegral", "CosIntegral", "SinhIntegral", "CoshIntegral",
    "FresnelC", "FresnelS", "ProductLog",
    "BesselJ", "BesselY", "BesselI", "BesselK",
    "AiryAi", "AiryBi", "AiryAiPrime", "AiryBiPrime",
    "Zeta", "HurwitzZeta", "LerchPhi", "PolyLog", "StieltjesGamma",
    "Hypergeometric0F1", "Hypergeometric1F1", "Hypergeometric2F1",
    "HypergeometricPFQ", "QPochhammer",
    "LegendreP", "LegendreQ", "ChebyshevT", "ChebyshevU",
    "HermiteH", "LaguerreL", "GegenbauerC", "JacobiP",
    "SphericalHarmonicY", "SphericalBesselJ", "SphericalBesselY",
    "EllipticK", "EllipticE", "EllipticF", "EllipticPi", "EllipticTheta",
    "JacobiSN", "JacobiCN", "JacobiDN", "JacobiAmplitude",
    "WeierstrassP", "WeierstrassPPrime",
    "Fibonacci", "LucasL", "BernoulliB", "EulerE",
    "GammaRegularized", "BetaRegularized", "InverseGammaRegularized",
    "InverseBetaRegularized",

    # --- reductions and statistics ---------------------------------------
    "Total", "Mean", "Median", "Variance", "StandardDeviation",
    "RootMeanSquare", "Quartiles", "Quantile", "Skewness", "Kurtosis",
    "Correlation", "Covariance", "GeometricMean", "HarmonicMean",
    "MeanDeviation", "MedianDeviation", "InterquartileRange",
    "CentralMoment", "Moment", "Standardize", "Norm", "Normalize",
    "Accumulate", "Differences", "Ratios", "Tally", "Counts", "BinCounts",
    "MovingAverage", "MovingMedian", "ExponentialMovingAverage",
    "TotalVariation", "MaximalBy", "MinimalBy", "Ordering", "Commonest",

    # --- structural over arrays -------------------------------------------
    "Sort", "ReverseSort", "SortBy", "Reverse", "Flatten", "Transpose",
    "Take", "Drop", "Part", "Extract", "First", "Last", "Most", "Rest",
    "RotateLeft", "RotateRight", "Partition", "Riffle", "Join",
    "PadLeft", "PadRight", "Length", "Dimensions", "Depth", "ArrayReshape",
    "Union", "Intersection", "Complement", "DeleteDuplicates", "Split",
    "SplitBy", "Gather", "GatherBy", "Position", "Pick", "Select",
    "Insert", "Delete", "ReplacePart", "Append", "Prepend", "Range",
    "ConstantArray", "Array", "Table", "Subdivide", "Flatten", "Level",
    "Nearest", "TakeWhile", "LengthWhile", "SelectFirst", "Count",
    "FirstPosition", "ArrayFlatten", "ArrayPad", "Diagonal", "UpperTriangularize",
    "LowerTriangularize", "SparseArray", "Normal",

    # --- functional over arrays -------------------------------------------
    "Map", "MapIndexed", "MapThread", "MapAt", "MapAll", "Apply", "Scan",
    "Fold", "FoldList", "Nest", "NestList", "NestWhile", "NestWhileList",
    "FixedPoint", "FixedPointList", "Outer", "Inner", "Thread",
    "AllTrue", "AnyTrue", "NoneTrue", "Do", "Sum", "Product",

    # --- comparison / boolean over arrays ---------------------------------
    "Equal", "Unequal", "Less", "Greater", "LessEqual", "GreaterEqual",
    "Positive", "Negative", "NonNegative", "NonPositive",
    "SameQ", "UnsameQ", "Not", "And", "Or", "Xor",

    # --- linear algebra ----------------------------------------------------
    "Dot", "Cross", "Det", "Inverse", "PseudoInverse", "LinearSolve",
    "LeastSquares", "RowReduce", "NullSpace", "MatrixRank", "MatrixPower",
    "MatrixExp", "LUDecomposition", "QRDecomposition", "CholeskyDecomposition",
    "SchurDecomposition", "SingularValueDecomposition", "SingularValueList",
    "Eigenvalues", "Eigenvectors", "Eigensystem", "Tr", "IdentityMatrix",
    "DiagonalMatrix", "HankelMatrix", "ToeplitzMatrix", "VandermondeMatrix",
    "KroneckerProduct", "MatrixLog", "Adjugate", "Minors", 
    "PositiveDefiniteMatrixQ", "NegativeDefiniteMatrixQ", "HermitianMatrixQ",
    "SymmetricMatrixQ", "OrthogonalMatrixQ", "UnitaryMatrixQ",
    "ConjugateTranspose", "Projection", "Orthogonalize", 

    # --- transforms and signal --------------------------------------------
    "Fourier", "InverseFourier", "FourierDCT", "FourierDST",
    "FourierDCTMatrix", "FourierDSTMatrix", "FourierMatrix",
    "ListConvolve", "ListCorrelate", "PeriodogramArray", "Spectrogram",

    # --- numerics ----------------------------------------------------------
    "N", "NIntegrate", "NSum", "NProduct", "NDSolve", "NSolve", "NRoots",
    "FindRoot", "FindMinimum", "FindMaximum", "FindFit", "Fit",
    "Interpolation", "ListInterpolation", "InterpolatingPolynomial",
    "NMinimize", "NMaximize", "NLimit", "ND", "NResidue", "NSeries",
    "RandomReal", "RandomInteger", "RandomComplex", "RandomVariate",
    "RandomChoice", "RandomSample", "RandomPrime",

    # --- compiled surface ---------------------------------------------------
    "Compile", "NDArrayQ", "DataType", "ToNDArray",
    "FromNDArray", "ToPackedArray", "FromPackedArray", "ByteCount",

    # --- exact integer arithmetic that a numeric workload still reaches ----
    "GCD", "LCM", "PowerMod", "Mod", "PrimeQ", "NextPrime", "Prime",
    "PrimePi", "IntegerDigits", "IntegerLength", "BitAnd", "BitOr",
    "BitXor", "BitNot", "BitShiftLeft", "BitShiftRight", "IntegerExponent",
    "EvenQ", "OddQ", "DivisorSigma", "EulerPhi", "MoebiusMu", "JacobiSymbol",
    "Divisible", "CoprimeQ", "PrimeNu", "PrimeOmega", "LiouvilleLambda",
    "FromDigits", "DigitCount", "DigitSum", "RealDigits", "MantissaExponent",

    # --- second pass: heads that fell through to the symbolic default on the
    # first run of this script and are, on inspection, array operations a
    # numeric workload reaches. Kept separate so the two passes stay legible.
    "MinMax", "TakeLargest", "TakeSmallest", "TakeLargestBy", "TakeSmallestBy",
    "Catenate", "MemberQ", "Cases", "DeleteCases", "FirstCase",
    "ReIm", "Rationalize", "Numerator", "Denominator", "Complex",
    "UnitVector", "HilbertMatrix", "DiagonalMatrixQ", "SquareMatrixQ",
    "UpperTriangularMatrixQ", "MatrixQ", "VectorQ", "ListQ", "OrderedQ",
    "PositionIndex", "CountsBy", "GroupBy", "Lookup", "Keys", "Values",
    "ReverseSortBy", "DeleteDuplicatesBy", "Distribute", "DesignMatrix",
    "Histogram", "AdjacencyMatrix", "Subsets", "Tuples", "Permutations",
}

SCALAR = {
    "$MachineEpsilon", "$MachinePrecision", "$MaxMachineNumber",
    "$MinMachineNumber", "$MaxNumber", "$MinNumber",
    "Precision", "Accuracy", "SetPrecision", "SetAccuracy",
    "MachineNumberQ", "NumericQ", "NumberQ", "IntegerQ", "RealQ",
    "AbsoluteTiming", "Timing", "AbsoluteTime", "DateList", "SeedRandom",
    "Pi", "E", "EulerGamma", "GoldenRatio", "Catalan", "Degree", "I",
    "Infinity", "Indeterminate", "ComplexInfinity", "Khinchin", "Glaisher",
}


def category(name, kinds):
    if name in NUMERIC:
        return "numeric"
    if name in SCALAR:
        return "scalar"
    if any(k in name for k in ("Form", "Print", "Write", "Open", "Close",
                               "Stream", "File", "Export", "Import", "Plot",
                               "Graphics", "Show", "Style", "Color")):
        return "io"
    return "symbolic"


# ---------------------------------------------------------------------------

def live_names():
    """Names[] from the running system.

    The source scan can only classify heads that EXIST. A head that NumPy has
    and Mathilda does not never appears in it at all, so "no fast path" and "no
    function" look identical from the source -- and they are not the same
    finding. The first run of tools/numeric_sweep.py hit this: nine statistics
    heads (Quantile, Skewness, Kurtosis, Correlation, Covariance, Standardize,
    GeometricMean, HarmonicMean, TotalVariation) read as slow rows when in fact
    they return unevaluated. Asking the binary is the only way to tell.
    """
    import subprocess
    binpath = os.environ.get("MATHILDA_BIN", os.path.join(ROOT, "Mathilda"))
    if not os.path.exists(binpath):
        return None
    try:
        import tempfile
        with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as fh:
            fh.write('Scan[Print["#", #] &, Names["*"]];\n')
            path = fh.name
        cp = subprocess.run([binpath, "-file", path], capture_output=True,
                            text=True, timeout=120)
        os.unlink(path)
        return {ln[1:].strip() for ln in cp.stdout.splitlines()
                if ln.startswith("#")}
    except Exception:
        return None


def build():
    builtins = all_builtins()
    packs = pack_lists()
    kernels = kernel_heads()

    rows = []
    for name in sorted(builtins):
        k = kernels.get(name)
        row = {
            "name": name,
            "file": builtins[name],
            "category": category(name, k),
            "kernel": k or [],
            "aware": name in packs["aware"] or bool(k),
            "aware_explicit": name in packs["aware"],
            "not_aware": name in packs["not_aware"],
            "int64_ok": name in packs["int64_ok"],
            "broadcast_ok": name in packs["broadcast_ok"],
        }
        # A head cleared by NOT_AWARE is not aware however it got the bit.
        if row["not_aware"]:
            row["aware"] = False
        row["backed"] = bool(row["kernel"]) or row["aware"]
        rows.append(row)
    return rows, packs, kernels


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", metavar="PATH", help="write the full inventory as JSON")
    ap.add_argument("--gaps", action="store_true",
                    help="list only numeric heads with no fast path")
    ap.add_argument("--unclassified", action="store_true",
                    help="list heads that fell through to the symbolic default")
    ap.add_argument("--missing", action="store_true",
                    help="numeric heads this file names that the binary does "
                         "not define at all")
    args = ap.parse_args()

    if args.missing:
        live = live_names()
        if live is None:
            print("could not query the binary; build ./Mathilda first",
                  file=sys.stderr)
            return 2
        absent = sorted(n for n in NUMERIC if n not in live)
        for n in absent:
            print(n)
        print(f"\n{len(absent)} of {len(NUMERIC)} numeric heads are not "
              f"defined by the running system", file=sys.stderr)
        return 0

    rows, packs, kernels = build()

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(rows, fh, indent=1, sort_keys=True)
        print(f"wrote {args.json}: {len(rows)} builtins")

    by_cat = {}
    for r in rows:
        by_cat.setdefault(r["category"], []).append(r)

    if args.unclassified:
        for r in sorted(by_cat.get("symbolic", []), key=lambda r: r["name"]):
            print(f"{r['name']:<32} {r['file']}")
        return 0

    numeric = by_cat.get("numeric", [])
    gaps = [r for r in numeric if not r["backed"]]

    if args.gaps:
        print(f"{'head':<30} {'file':<34} note")
        print("-" * 82)
        for r in sorted(gaps, key=lambda r: (r["file"], r["name"])):
            note = "NOT_AWARE (deliberate)" if r["not_aware"] else ""
            print(f"{r['name']:<30} {r['file']:<34} {note}")
        print(f"\n{len(gaps)} of {len(numeric)} numeric heads have no registered path")
        return 0

    print(f"{'category':<12} {'total':>6} {'backed':>7} {'kernel':>7} "
          f"{'aware':>6} {'int64':>6}")
    print("-" * 50)
    for cat in ("numeric", "scalar", "symbolic", "io"):
        rs = by_cat.get(cat, [])
        print(f"{cat:<12} {len(rs):>6} "
              f"{sum(1 for r in rs if r['backed']):>7} "
              f"{sum(1 for r in rs if r['kernel']):>7} "
              f"{sum(1 for r in rs if r['aware']):>6} "
              f"{sum(1 for r in rs if r['int64_ok']):>6}")
    print("-" * 50)
    print(f"{'TOTAL':<12} {len(rows):>6}")
    print()
    print(f"registered kernels : {len(kernels)}")
    print(f"AWARE list         : {len(packs['aware'])}")
    print(f"INT64_OK list      : {len(packs['int64_ok'])}")
    print(f"numeric gaps       : {len(gaps)}  (--gaps to list)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
numeric_sweep.py -- every numeric builtin, timed against NumPy/SciPy on the
same data.

WHAT THIS IS FOR. tools/numeric_coverage.py reads the source and says which
heads are *registered* for a machine path. Registration is not speed: a head can
be on the AWARE list and still walk its buffer one element at a time, and a head
absent from every list can still be fast because it never touches elements. Only
a measurement separates the two, and only a second system says whether the
measurement is any good.

So this runs each numeric head over machine-number arrays and reports
nanoseconds per element beside the same operation in NumPy/SciPy. The rule the
nineteen experiments in docs/experiments/ arrived at is that NumPy is the ruler:
"slower than Mathematica" is a competitive gap, "slower than NumPy" is a gap to
the machine, and on this host both link the same Accelerate BLAS.

METHOD, following docs/experiments/README.md:
  * elapsed wall clock -- AbsoluteTiming, never Timing[], which sums CPU time
    over threads and overstates any threaded path by roughly the core count;
  * minimum of REPS repetitions, with the maximum recorded as a caching
    tripwire (a builtin that memoises answers repeat calls in ~0 s, and the
    minimum alone would make the row vacuous);
  * the same deterministic data in both systems, and a scalar checksum of every
    result compared between them -- a timing row means nothing until the two
    systems agree they computed the same thing;
  * --diag re-runs the flagged heads under MATHILDA_NO_PACK=1, which turns the
    packed representation off process-wide. A head whose packed and unpacked
    times are equal has no buffer path being reached, whatever the registry says.

Usage:
    python3 tools/numeric_sweep.py                     # full sweep, both systems
    python3 tools/numeric_sweep.py --only sin,cos,dot  # a subset, by id
    python3 tools/numeric_sweep.py --group elementwise # a whole group
    python3 tools/numeric_sweep.py --json out.json
    python3 tools/numeric_sweep.py --diag              # add the no-pack column
    python3 tools/numeric_sweep.py --system mathilda   # skip Python
"""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MATHILDA = os.environ.get("MATHILDA_BIN", os.path.join(ROOT, "Mathilda"))
PYTHON = os.environ.get("HPC_PYTHON", "/usr/local/bin/python3.11")

REPS = 4
BATCH = 12

# Per-BATCH wall-clock ceiling. A batch that blows it is retried one probe per
# process, so the ceiling is also the per-probe ceiling on the retry. Kept
# modest on purpose: a probe that cannot answer in this long is a finding
# ("does not complete"), and waiting longer does not make it a better one. The
# linear-algebra group is what set this -- RowReduce and NullSpace on a 300x300
# run the exact fraction-free path and do not finish.
TIMEOUT = int(os.environ.get("SWEEP_TIMEOUT", "180"))

# Where run_mathilda drops results as they arrive. Without it a run that is
# killed, or that stalls in one slow group, loses every measurement taken
# before it -- which is most of them.
CHECKPOINT = os.environ.get("SWEEP_CHECKPOINT", "")


# ---------------------------------------------------------------------------
# Shared data. The same values in both systems, built the same way, so the
# checksum comparison is meaningful rather than coincidental.
#
# The matrices are diagonally dominant on purpose: a matrix built out of Range
# is rank 2, so Det is 0, Inverse is singular and LinearSolve measures the error
# path rather than the solve. That is the kind of benchmark bug §11 of
# performance.md is about.
# ---------------------------------------------------------------------------

N = 1000000

MPRE = r"""
n = 1000000;
v  = Range[1., n]/n;
w  = Reverse[v];
u  = 2. v - 1.;
p  = 1. + v;
iv = Range[n];
jv = Mod[Range[n]*7919, 1000];
kv = Mod[Range[n]*7919, 999983];
bv = Take[v, 200000];
cv = Take[iv, 200000];
sv = Take[v, 100000];
A1000 = Table[If[i == j, 1000., 1./(1. + Abs[i - j])], {i, 1000}, {j, 1000}];
A500  = Table[If[i == j, 500.,  1./(1. + Abs[i - j])], {i, 500},  {j, 500}];
A300  = Table[If[i == j, 300.,  1./(1. + Abs[i - j])], {i, 300},  {j, 300}];
A200  = Table[If[i == j, 200.,  1./(1. + Abs[i - j])], {i, 200},  {j, 200}];
A6    = Table[If[i == j, 6.,    1./(1. + Abs[i - j])], {i, 6},    {j, 6}];
b1000 = Range[1., 1000.]/1000.;
b500  = Range[1., 500.]/500.;
b6    = Range[1., 6.]/6.;
im    = Table[N[Mod[i*j, 251]]/251., {i, 1024}, {j, 1024}];
k5    = Table[1./25., {5}, {5}];
k5v   = Table[1./5., {5}];
idx   = Mod[Range[n]*7919, 100000] + 1;
src   = Range[1., 100000.];
tf    = Map[# > 0.5 &, sv];
tg    = Map[# < 0.7 &, sv];
ck[x_] := N[Total[Flatten[{x}]]];
ckt[x_] := N[Count[Flatten[{x}], True]];
"""

PYPRE = r"""
import numpy as np, time, json, sys, math
try:
    import scipy.special as sp
except Exception:
    sp = None
try:
    import scipy.linalg as sla
except Exception:
    sla = None
try:
    import scipy.signal as ssig
except Exception:
    ssig = None
try:
    import scipy.stats as sst
except Exception:
    sst = None

n = 1000000
v  = np.arange(1.0, n + 1.0) / n
w  = v[::-1].copy()
u  = 2.0 * v - 1.0
p  = 1.0 + v
iv = np.arange(1, n + 1, dtype=np.int64)
jv = (iv * 7919) % 1000
kv = (iv * 7919) % 999983
bv = v[:200000].copy()
cv = iv[:200000].copy()
sv = v[:100000].copy()

def _dd(m, d):
    i = np.arange(1, m + 1)[:, None]
    j = np.arange(1, m + 1)[None, :]
    return np.where(i == j, float(d), 1.0 / (1.0 + np.abs(i - j)))

A1000 = _dd(1000, 1000)
A500  = _dd(500, 500)
A300  = _dd(300, 300)
A200  = _dd(200, 200)
A6    = _dd(6, 6)
b1000 = np.arange(1.0, 1001.0) / 1000.0
b500  = np.arange(1.0, 501.0) / 500.0
b6    = np.arange(1.0, 7.0) / 6.0
_i = np.arange(1, 1025)[:, None]; _j = np.arange(1, 1025)[None, :]
im = ((_i * _j) % 251).astype(np.float64) / 251.0
k5  = np.full((5, 5), 1.0 / 25.0)
k5v = np.full(5, 1.0 / 5.0)
idx = ((np.arange(1, n + 1) * 7919) % 100000)          # 0-based already
src = np.arange(1.0, 100001.0)
tf = sv > 0.5
tg = sv < 0.7

def ck(x):
    a = np.asarray(x)
    if a.dtype == object:
        return float(len(a))
    return float(np.sum(a))
"""


# ---------------------------------------------------------------------------
# The probe catalogue.
#
#   P(id, group, mathilda_expr, py=..., chk=..., pychk=..., note=...)
#
# `chk` defaults to ck[r] on both sides: the sum of the flattened result. It is
# computed OUTSIDE the timed region, so an expensive check costs wall clock but
# never contaminates a measurement. Heads whose result is boolean use ckt[r]
# (count of True) against NumPy's np.sum of a bool array, which is the same
# number.
#
# py=None means "no NumPy equivalent" and prints as an em dash. That is a real
# answer for arbitrary precision and symbolic dispatch, and it is NOT a licence
# to skip a head that merely needed thought.
# ---------------------------------------------------------------------------

PROBES = []


def P(pid, group, m, py=None, chk="ck[r]", pychk="ck(r)", note=""):
    PROBES.append(dict(id=pid, group=group, m=m, py=py, chk=chk,
                       pychk=pychk, note=note))


# --- elementwise unary, float64, n = 10^6 ---------------------------------
for name, pyf in [
    ("Sin", "np.sin(v)"), ("Cos", "np.cos(v)"), ("Tan", "np.tan(v)"),
    ("Sinh", "np.sinh(v)"), ("Cosh", "np.cosh(v)"), ("Tanh", "np.tanh(v)"),
    ("ArcSin", "np.arcsin(v)"), ("ArcCos", "np.arccos(v)"),
    ("ArcTan", "np.arctan(v)"), ("ArcSinh", "np.arcsinh(v)"),
    ("ArcTanh", "np.arctanh(v*0.99)"), ("Exp", "np.exp(v)"),
    ("Log", "np.log(v)"), ("Abs", "np.abs(u)"), ("Sqrt", "np.sqrt(v)"),
]:
    arg = "v"
    if name == "Abs":
        arg = "u"
    if name == "ArcTanh":
        arg = "v*0.99"
    P(name.lower(), "elementwise", f"{name}[{arg}]", pyf)

P("cot", "elementwise", "Cot[v]", "1.0/np.tan(v)")
P("sec", "elementwise", "Sec[v]", "1.0/np.cos(v)")
P("csc", "elementwise", "Csc[v]", "1.0/np.sin(v)")
P("coth", "elementwise", "Coth[v]", "1.0/np.tanh(v)")
P("sech", "elementwise", "Sech[v]", "1.0/np.cosh(v)")
P("csch", "elementwise", "Csch[v]", "1.0/np.sinh(v)")
P("arccot", "elementwise", "ArcCot[p]", "np.arctan(1.0/p)")
P("arcsec", "elementwise", "ArcSec[p]", "np.arccos(1.0/p)")
P("arccsc", "elementwise", "ArcCsc[p]", "np.arcsin(1.0/p)")
P("arccosh", "elementwise", "ArcCosh[p]", "np.arccosh(p)")
P("arccoth", "elementwise", "ArcCoth[1. + p]", "np.arctanh(1.0/(1.0+p))")
P("arcsech", "elementwise", "ArcSech[v]", "np.arccosh(1.0/v)")
P("arccsch", "elementwise", "ArcCsch[v]", "np.arcsinh(1.0/v)")
P("log10", "elementwise", "Log[10, v]", "np.log10(v)")
P("log2", "elementwise", "Log[2, v]", "np.log2(v)")
P("cuberoot", "elementwise", "CubeRoot[u]", "np.cbrt(u)")
P("sign", "elementwise", "Sign[u]", "np.sign(u)")
P("unitstep", "elementwise", "UnitStep[u]", "(u >= 0).astype(np.int64)")
P("ramp", "elementwise", "Ramp[u]", "np.maximum(u, 0.0)")
P("floor", "elementwise", "Floor[10. u]", "np.floor(10.0*u)")
P("ceiling", "elementwise", "Ceiling[10. u]", "np.ceil(10.0*u)")
P("round", "elementwise", "Round[10. u]", "np.round(10.0*u)")
P("integerpart", "elementwise", "IntegerPart[10. u]", "np.trunc(10.0*u)")
P("fractionalpart", "elementwise", "FractionalPart[10. u]",
  "10.0*u - np.trunc(10.0*u)")
P("re", "elementwise", "Re[v]", "np.real(v).copy()")
P("im", "elementwise", "Im[v]", "np.imag(v).copy()")
P("conjugate", "elementwise", "Conjugate[v]", "np.conj(v)")
P("clip", "elementwise", "Clip[u, {-0.5, 0.5}]", "np.clip(u, -0.5, 0.5)")
P("rescale", "elementwise", "Rescale[v]", "(v - v.min())/(v.max()-v.min())")
P("chop", "elementwise", "Chop[u]", "np.where(np.abs(u) < 1e-10, 0.0, u)")
P("boole", "elementwise", "Boole[Map[# > 0.5 &, sv]]",
  "(sv > 0.5).astype(np.int64)")
P("unitize", "elementwise", "Unitize[u]", "(u != 0).astype(np.int64)")
P("logisticsigmoid", "elementwise", "LogisticSigmoid[u]",
  "1.0/(1.0+np.exp(-u))")
P("gudermannian", "elementwise", "Gudermannian[u]", "np.arctan(np.sinh(u))")
P("haversine", "elementwise", "Haversine[u]", "np.sin(u/2.0)**2")
P("sinc", "elementwise", "Sinc[p]", "np.sin(p)/p")

# --- elementwise binary ----------------------------------------------------
P("plus", "elementwise", "v + w", "v + w")
P("times", "elementwise", "v w", "v * w")
P("subtract", "elementwise", "v - w", "v - w")
P("divide", "elementwise", "v/p", "v / p")
P("power_sq", "elementwise", "v^2", "v**2")
P("power_arr", "elementwise", "p^w", "p ** w")
P("triad", "elementwise", "v + 3. w", "v + 3.0*w")
P("mod_real", "elementwise", "Mod[10. v, 3.]", "np.mod(10.0*v, 3.0)")
P("max2", "elementwise", "MapThread[Max, {v, w}]", "np.maximum(v, w)")
P("min2", "elementwise", "MapThread[Min, {v, w}]", "np.minimum(v, w)")
P("arctan2", "elementwise", "ArcTan[v, w]", "np.arctan2(w, v)")

# --- comparison and boolean masks -----------------------------------------
#
# Greater and EvenQ are NOT Listable in the Wolfram Language, so `v > 0.5` and
# `EvenQ[iv]` stay unevaluated on a list in Mathematica exactly as they do here.
# The first version of this group probed them anyway and read five ~3000x rows
# that were measuring the evaluator declining to answer. The idiomatic mask is
# Map over a pure function, or arithmetic on UnitStep, and those are what a
# Wolfram programmer writes -- so those are what is timed.
P("greater", "mask", "Map[# > 0.5 &, v]", "v > 0.5", chk="ckt[r]")
P("less", "mask", "Map[# < 0.5 &, v]", "v < 0.5", chk="ckt[r]")
P("greaterequal", "mask", "Map[# >= 0.5 &, v]", "v >= 0.5", chk="ckt[r]")
P("greater_unitstep", "mask", "UnitStep[v - 0.5]",
  "(v >= 0.5).astype(np.int64)")
P("equal_arr", "mask", "MapThread[Equal, {v, w}]", "v == w", chk="ckt[r]")
P("unequal_arr", "mask", "MapThread[Unequal, {v, w}]", "v != w", chk="ckt[r]")
P("positive", "mask", "Positive[u]", "u > 0", chk="ckt[r]")
P("negative", "mask", "Negative[u]", "u < 0", chk="ckt[r]")
P("nonnegative", "mask", "NonNegative[u]", "u >= 0", chk="ckt[r]")
P("evenq", "mask", "Map[EvenQ, iv]", "iv % 2 == 0", chk="ckt[r]")
P("oddq", "mask", "Map[OddQ, iv]", "iv % 2 == 1", chk="ckt[r]")
P("and_mask", "mask", "MapThread[And, {tf, tg}]", "tf & tg", chk="ckt[r]")
P("or_mask", "mask", "MapThread[Or, {tf, tg}]", "tf | tg", chk="ckt[r]")
P("not_mask", "mask", "Not /@ tf", "~tf", chk="ckt[r]")
P("mask_select", "mask", "Pick[sv, tf]", "sv[tf]")
P("mask_count", "mask", "Count[Map[# > 0.5 &, sv], True]",
  "int(np.count_nonzero(sv > 0.5))")

# --- reductions -------------------------------------------------------------
P("total", "reduce", "Total[v]", "np.sum(v)")
P("mean", "reduce", "Mean[v]", "np.mean(v)")
P("median", "reduce", "Median[v]", "np.median(v)")
P("max", "reduce", "Max[v]", "np.max(v)")
P("min", "reduce", "Min[v]", "np.min(v)")
P("minmax", "reduce", "MinMax[v]", "np.array([v.min(), v.max()])")
P("variance", "reduce", "Variance[v]", "np.var(v, ddof=1)")
P("stddev", "reduce", "StandardDeviation[v]", "np.std(v, ddof=1)")
P("rms", "reduce", "RootMeanSquare[v]", "np.sqrt(np.mean(v*v))")
P("norm2", "reduce", "Norm[v]", "np.linalg.norm(v)")
P("normalize", "reduce", "Normalize[v]", "v/np.linalg.norm(v)")
P("standardize", "reduce", "Standardize[v]", "(v-v.mean())/v.std(ddof=1)")
P("quantile", "reduce", "Quantile[v, 0.9]", "np.quantile(v, 0.9)")
P("skewness", "reduce", "Skewness[v]", "sst.skew(v)")
P("kurtosis", "reduce", "Kurtosis[v]", "sst.kurtosis(v, fisher=False)")
P("geometricmean", "reduce", "GeometricMean[v]", "np.exp(np.mean(np.log(v)))")
P("harmonicmean", "reduce", "HarmonicMean[v]", "len(v)/np.sum(1.0/v)")
P("correlation", "reduce", "Correlation[v, w]", "np.corrcoef(v, w)[0,1]")
P("covariance", "reduce", "Covariance[v, w]", "np.cov(v, w, ddof=1)[0,1]")
P("total_int", "reduce", "Total[iv]", "int(np.sum(iv))")
P("max_int", "reduce", "Max[iv]", "int(np.max(iv))")
P("totalvariation", "reduce", "TotalVariation[v]",
  "np.sum(np.abs(np.diff(v)))")
P("alltrue", "reduce", "AllTrue[v, # > 0 &]", "bool(np.all(v > 0))",
  chk="ckt[r]", pychk="float(r)")
P("anytrue", "reduce", "AnyTrue[v, # > 0.9 &]", "bool(np.any(v > 0.9))",
  chk="ckt[r]", pychk="float(r)")

# --- scans -------------------------------------------------------------------
P("accumulate", "scan", "Accumulate[v]", "np.cumsum(v)")
P("accumulate_int", "scan", "Accumulate[cv]", "np.cumsum(cv)")
P("differences", "scan", "Differences[v]", "np.diff(v)")
P("ratios", "scan", "Ratios[v]", "v[1:]/v[:-1]")
P("foldlist_max", "scan", "FoldList[Max, -1., v]",
  "np.maximum.accumulate(np.concatenate(([-1.0], v)))")
P("foldlist_plus", "scan", "FoldList[Plus, 0., v]",
  "np.cumsum(np.concatenate(([0.0], v)))")
P("foldlist_ema", "scan", "FoldList[0.9 #1 + 0.1 #2 &, 0., bv]",
  "ssig.lfilter([0.1], [1.0, -0.9], bv, zi=np.array([0.9*0.0]))[0]")
P("movingaverage", "scan", "MovingAverage[bv, 5]",
  "np.convolve(bv, np.ones(5)/5.0, mode='valid')")
P("exponentialmovingaverage", "scan", "ExponentialMovingAverage[bv, 0.1]",
  None)

# --- structural --------------------------------------------------------------
P("sort", "struct", "Sort[v]", "np.sort(v)")
P("sort_int", "struct", "Sort[kv]", "np.sort(kv)")
P("reversesort", "struct", "ReverseSort[v]", "np.sort(v)[::-1].copy()")
P("sortby", "struct", "SortBy[bv, Minus]", "bv[np.argsort(-bv)]")
P("ordering", "struct", "Ordering[v]", "np.argsort(v)+1")
P("reverse", "struct", "Reverse[v]", "v[::-1].copy()")
P("rotateleft", "struct", "RotateLeft[v, 1]", "np.roll(v, -1)")
P("rotateright", "struct", "RotateRight[v, 1]", "np.roll(v, 1)")
P("take", "struct", "Take[v, 500000]", "v[:500000].copy()")
P("drop", "struct", "Drop[v, 500000]", "v[500000:].copy()")
P("part_span", "struct", "v[[200000 ;; 800000]]", "v[199999:800000].copy()")
P("part_gather", "struct", "src[[idx]]", "src[idx]")
P("first", "struct", "First[v]", "v[0]")
P("last", "struct", "Last[v]", "v[-1]")
P("most", "struct", "Most[v]", "v[:-1].copy()")
P("rest", "struct", "Rest[v]", "v[1:].copy()")
P("join", "struct", "Join[v, w]", "np.concatenate([v, w])")
P("riffle", "struct", "Riffle[v, 0.]", None)
P("partition_tile", "struct", "Partition[v, 2]", "v.reshape(-1, 2).copy()")
P("partition_slide", "struct", "Partition[bv, 5, 1]",
  "np.lib.stride_tricks.sliding_window_view(bv, 5)")
P("padright", "struct", "PadRight[v, 1100000, 0.]",
  "np.pad(v, (0, 100000))")
P("padleft", "struct", "PadLeft[v, 1100000, 0.]", "np.pad(v, (100000, 0))")
P("flatten", "struct", "Flatten[A1000]", "A1000.ravel().copy()")
P("transpose", "struct", "Transpose[A1000]", "A1000.T.copy()")
P("arrayreshape", "struct", "ArrayReshape[v, {1000, 1000}]",
  "v.reshape(1000, 1000).copy()")
P("union", "struct", "Union[kv]", "np.unique(kv)")
P("intersection", "struct", "Intersection[jv, kv]", "np.intersect1d(jv, kv)")
P("complement", "struct", "Complement[kv, jv]", "np.setdiff1d(kv, jv)")
P("deleteduplicates", "struct", "DeleteDuplicates[jv]", None)
P("split", "struct", "Split[Sort[jv]]", None, chk="N[Length[r]]",
  pychk="float(len(r))")
P("extract", "struct", "Extract[v, {{1}, {2}, {3}}]", "v[[0,1,2]]")
P("select", "struct", "Select[bv, # > 0.5 &]", "bv[bv > 0.5]")
P("takewhile", "struct", "TakeWhile[v, # < 0.5 &]", "v[:np.argmax(v >= 0.5)]")
P("takelargest", "struct", "TakeLargest[bv, 10]",
  "np.sort(bv)[-10:][::-1].copy()")
P("takesmallest", "struct", "TakeSmallest[bv, 10]", "np.sort(bv)[:10]")
P("catenate", "struct", "Catenate[{v, w}]", "np.concatenate([v, w])")
P("memberq", "struct", "MemberQ[v, 0.5]", "bool(np.any(v == 0.5))",
  chk="ckt[r]", pychk="float(r)")
P("position", "struct", "Position[jv, 0]", None, chk="N[Length[r]]",
  pychk="float(len(r))")
P("count", "struct", "Count[jv, 0]", "int(np.count_nonzero(jv == 0))")
P("tally", "struct", "Tally[jv]", None, chk="N[Length[r]]",
  pychk="float(len(r))")
P("counts", "struct", "Counts[jv]", None, chk="N[Length[r]]",
  pychk="float(len(r))")
P("bincounts", "struct", "BinCounts[jv, {0, 1000, 1}]",
  "np.bincount(jv, minlength=1000)")
P("gatherby", "struct", "GatherBy[cv, EvenQ]", None, chk="N[Length[r]]",
  pychk="float(len(r))")
P("arraypad", "struct", "ArrayPad[v, 10]", "np.pad(v, (10, 10))")
P("replacepart", "struct", "ReplacePart[v, 1 -> 0.]", None)
P("insert", "struct", "Insert[v, 0., 1]", "np.insert(v, 0, 0.0)")
P("delete", "struct", "Delete[v, 1]", "np.delete(v, 0)")
P("append", "struct", "Append[v, 0.]", "np.append(v, 0.0)")
P("prepend", "struct", "Prepend[v, 0.]", "np.insert(v, 0, 0.0)")

# --- construction ------------------------------------------------------------
P("range_real", "build", "Range[1., 1000000.]", "np.arange(1.0, 1000001.0)")
P("range_int", "build", "Range[1000000]", "np.arange(1, 1000001)")
P("constantarray", "build", "ConstantArray[1., 1000000]", "np.ones(1000000)")
P("table_num", "build", "Table[N[i], {i, 1000000}]",
  "np.arange(1.0, 1000001.0)")
P("array_fn", "build", "Array[N, 1000000]", "np.arange(1.0, 1000001.0)")
P("subdivide", "build", "Subdivide[0., 1., 999999]", "np.linspace(0.0,1.0,1000000)")
P("identitymatrix", "build", "IdentityMatrix[1000]", "np.eye(1000)")
P("diagonalmatrix", "build", "DiagonalMatrix[b1000]", "np.diag(b1000)")
P("unitvector", "build", "UnitVector[1000000, 1]",
  "np.eye(1, 1000000, 0).ravel()")
P("randomreal", "build", "RandomReal[{0, 1}, 1000000]", "np.random.random(1000000)",
  chk="N[Length[r]]", pychk="float(len(r))")
P("randominteger", "build", "RandomInteger[{0, 100}, 1000000]",
  "np.random.randint(0, 101, 1000000)", chk="N[Length[r]]",
  pychk="float(len(r))")
P("randomvariate", "build", "RandomVariate[NormalDistribution[0, 1], 1000000]",
  "np.random.standard_normal(1000000)", chk="N[Length[r]]",
  pychk="float(len(r))")
P("randomsample", "build", "RandomSample[cv]", "np.random.permutation(cv)",
  chk="N[Length[r]]", pychk="float(len(r))")
P("randomchoice", "build", "RandomChoice[cv, 200000]",
  "np.random.choice(cv, 200000)", chk="N[Length[r]]", pychk="float(len(r))")

# --- functional --------------------------------------------------------------
P("map_pure", "functional", "Map[Sin, v]", "np.sin(v)")
P("map_lambda", "functional", "Map[# 2. + 1. &, v]", "v*2.0 + 1.0")
P("mapindexed", "functional", "MapIndexed[#1 &, bv]", "bv.copy()")
P("mapthread_plus", "functional", "MapThread[Plus, {v, w}]", "v + w")
P("apply_plus", "functional", "Apply[Plus, v]", "np.sum(v)")
P("fold_plus", "functional", "Fold[Plus, 0., v]", "np.sum(v)")
P("nest_step", "functional", "Nest[# + 1. &, v, 10]", "v + 10.0")
P("nestlist_step", "functional", "NestList[# + 1. &, sv, 10]",
  "np.array([sv + k for k in range(11)])")
P("outer_sub", "functional", "Outer[Subtract, b1000, b1000]",
  "b1000[:,None] - b1000[None,:]")
P("inner_dot", "functional", "Inner[Times, b1000, b1000, Plus]",
  "np.dot(b1000, b1000)")
P("thread_plus", "functional", "Thread[Plus[v, w]]", "v + w")
P("do_loop", "functional", "Do[q = i^2, {i, 1000000}]", None,
  chk="0.", pychk="0.0")
P("sum_loop", "functional", "Sum[N[i], {i, 1000000}]",
  "float(np.sum(np.arange(1.0, 1000001.0)))")

# --- linear algebra ----------------------------------------------------------
P("dot_mm", "linalg", "A1000 . A1000", "A1000 @ A1000")
P("dot_mv", "linalg", "A1000 . b1000", "A1000 @ b1000")
P("dot_vv", "linalg", "v . w", "np.dot(v, w)")
P("det", "linalg", "Det[A500]", "np.linalg.det(A500)")
P("inverse", "linalg", "Inverse[A500]", "np.linalg.inv(A500)")
P("linearsolve", "linalg", "LinearSolve[A1000, b1000]",
  "np.linalg.solve(A1000, b1000)")
P("leastsquares", "linalg", "LeastSquares[A500, b500]",
  "np.linalg.lstsq(A500, b500, rcond=None)[0]")
P("matrixrank", "linalg", "MatrixRank[A300]", "np.linalg.matrix_rank(A300)")
P("rowreduce", "linalg", "RowReduce[A300]", None)
P("nullspace", "linalg", "NullSpace[A300]", None, chk="N[Length[r]]",
  pychk="float(len(r))")
P("matrixpower", "linalg", "MatrixPower[A300, 4]",
  "np.linalg.matrix_power(A300, 4)")
P("lu", "linalg", "LUDecomposition[A500]", None, chk="N[Length[r]]",
  pychk="float(len(r))")
P("qr", "linalg", "QRDecomposition[A500]", "np.linalg.qr(A500)",
  chk="N[Length[r]]", pychk="float(len(r))")
P("cholesky", "linalg", "CholeskyDecomposition[A500]",
  "np.linalg.cholesky(A500)")
P("svd", "linalg", "SingularValueDecomposition[A300]",
  "np.linalg.svd(A300)", chk="N[Length[r]]", pychk="float(len(r))")
P("singularvaluelist", "linalg", "SingularValueList[A300]",
  "np.linalg.svd(A300, compute_uv=False)")
P("eigenvalues", "linalg", "Eigenvalues[A300]",
  "np.linalg.eigvalsh(A300)", chk="N[Length[r]]", pychk="float(len(r))")
P("eigenvectors", "linalg", "Eigenvectors[A300]", None,
  chk="N[Length[r]]", pychk="float(len(r))")
P("tr", "linalg", "Tr[A1000]", "np.trace(A1000)")
P("norm_matrix", "linalg", "Norm[A300]", "np.linalg.norm(A300, 2)")
P("norm_frobenius", "linalg", "Norm[A300, \"Frobenius\"]",
  "np.linalg.norm(A300, 'fro')")
P("pseudoinverse", "linalg", "PseudoInverse[A300]", "np.linalg.pinv(A300)")
P("matrixexp", "linalg", "MatrixExp[A6]", "sla.expm(A6)")
P("kroneckerproduct", "linalg", "KroneckerProduct[A200, IdentityMatrix[2]]",
  "np.kron(A200, np.eye(2))")
P("cross", "linalg", "Cross[{1., 2., 3.}, {4., 5., 6.}]",
  "np.cross([1.,2.,3.],[4.,5.,6.])")
P("conjugatetranspose", "linalg", "ConjugateTranspose[A1000]",
  "A1000.conj().T.copy()")
P("orthogonalize", "linalg", "Orthogonalize[A200]", None)
P("symmetricmatrixq", "linalg", "SymmetricMatrixQ[A1000]",
  "bool(np.allclose(A1000, A1000.T))", chk="ckt[r]", pychk="float(r)")
P("positivedefiniteq", "linalg", "PositiveDefiniteMatrixQ[A300]",
  "bool(np.all(np.linalg.eigvalsh(A300) > 0))", chk="ckt[r]", pychk="float(r)")

# --- SMALL matrices: the packing threshold is a floor here (experiment 18) ---
P("det6", "linalg-small", "Table[Det[A6], {200}]",
  "np.array([np.linalg.det(A6) for _ in range(200)])")
P("inverse6", "linalg-small", "Table[Inverse[A6], {200}]",
  "np.array([np.linalg.inv(A6) for _ in range(200)])")
P("linearsolve6", "linalg-small", "Table[LinearSolve[A6, b6], {200}]",
  "np.array([np.linalg.solve(A6, b6) for _ in range(200)])")
P("dot6", "linalg-small", "Table[A6 . A6, {200}]",
  "np.array([A6 @ A6 for _ in range(200)])")

# --- transforms ---------------------------------------------------------------
P("fourier", "transform", "Fourier[v]", "np.fft.fft(v)/np.sqrt(len(v))",
  chk="N[Abs[Total[Flatten[{r}]]]]", pychk="float(np.abs(np.sum(r)))")
P("inversefourier", "transform", "InverseFourier[v]",
  "np.fft.ifft(v)*np.sqrt(len(v))",
  chk="N[Abs[Total[Flatten[{r}]]]]", pychk="float(np.abs(np.sum(r)))")
P("fourierdct", "transform", "FourierDCT[bv]", None)
P("fourierdst", "transform", "FourierDST[bv]", None)
P("listconvolve_1d", "transform", "ListConvolve[k5v, v]",
  "np.convolve(v, k5v, mode='valid')")
P("listcorrelate_1d", "transform", "ListCorrelate[k5v, v]",
  "np.correlate(v, k5v, mode='valid')")
P("listcorrelate_2d", "transform", "ListCorrelate[k5, im]",
  "ssig.correlate2d(im, k5, mode='valid')")

# --- special functions --------------------------------------------------------
P("gamma", "special", "Gamma[p]", "sp.gamma(p)")
P("loggamma", "special", "LogGamma[p]", "sp.gammaln(p)")
P("polygamma", "special", "PolyGamma[p]", "sp.psi(p)")
P("erf", "special", "Erf[v]", "sp.erf(v)")
P("erfc", "special", "Erfc[v]", "sp.erfc(v)")
P("erfi", "special", "Erfi[v]", "sp.erfi(v)")
P("inverseerf", "special", "InverseErf[u*0.99]", "sp.erfinv(u*0.99)")
P("expintegralei", "special", "ExpIntegralEi[p]", "sp.expi(p)")
P("logintegral", "special", "LogIntegral[1. + p]",
  "sp.expi(np.log(1.0+p))")
P("sinintegral", "special", "SinIntegral[p]", "sp.sici(p)[0]")
P("cosintegral", "special", "CosIntegral[p]", "sp.sici(p)[1]")
P("sinhintegral", "special", "SinhIntegral[p]", "sp.shichi(p)[0]")
P("coshintegral", "special", "CoshIntegral[p]", "sp.shichi(p)[1]")
P("fresnelc", "special", "FresnelC[p]", "sp.fresnel(p)[1]")
P("fresnels", "special", "FresnelS[p]", "sp.fresnel(p)[0]")
P("productlog", "special", "ProductLog[p]", "np.real(sp.lambertw(p))")
P("besselj", "special", "BesselJ[0, p]", "sp.jv(0, p)")
P("bessely", "special", "BesselY[0, p]", "sp.yv(0, p)")
P("besseli", "special", "BesselI[0, p]", "sp.iv(0, p)")
P("besselk", "special", "BesselK[0, p]", "sp.kv(0, p)")
P("airyai", "special", "AiryAi[u]", "sp.airy(u)[0]")
P("airybi", "special", "AiryBi[u]", "sp.airy(u)[2]")
P("zeta", "special", "Zeta[1. + p]", "sp.zeta(1.0+p)")
P("beta_fn", "special", "Beta[p, p]", "sp.beta(p, p)")
P("harmonicnumber", "special", "HarmonicNumber[p]", "sp.psi(p+1)+np.euler_gamma")
P("legendrep", "special", "LegendreP[2, u]", "sp.eval_legendre(2, u)")
P("polylog", "special", "PolyLog[2, v*0.9]", "sp.spence(1.0 - v*0.9)")
P("hyper2f1", "special", "Hypergeometric2F1[0.5, 0.5, 1.5, v*0.5]",
  "sp.hyp2f1(0.5, 0.5, 1.5, v*0.5)")
P("hyper1f1", "special", "Hypergeometric1F1[0.5, 1.5, v]",
  "sp.hyp1f1(0.5, 1.5, v)")
P("gammaregularized", "special", "GammaRegularized[2., p]",
  "sp.gammaincc(2.0, p)")
P("betaregularized", "special", "BetaRegularized[v, 2., 3.]",
  "sp.betainc(2.0, 3.0, v)")
P("factorial_real", "special", "Factorial[10. v]", "sp.gamma(10.0*v + 1.0)")
P("binomial_real", "special", "Binomial[10., 10. v]",
  "sp.binom(10.0, 10.0*v)")
P("pochhammer", "special", "Pochhammer[p, 2.]", "sp.poch(p, 2.0)")

# --- integer / number theory over int64 arrays -------------------------------
P("gcd_arr", "integer", "GCD[cv, 1234]", "np.gcd(cv, 1234)")
P("lcm_arr", "integer", "LCM[cv, 12]", "np.lcm(cv, 12)")
P("mod_int", "integer", "Mod[iv, 1000]", "np.mod(iv, 1000)")
P("quotient_int", "integer", "Quotient[iv, 1000]", "iv // 1000")
P("bitand", "integer", "BitAnd[iv, 255]", "np.bitwise_and(iv, 255)")
P("bitor", "integer", "BitOr[iv, 255]", "np.bitwise_or(iv, 255)")
P("bitxor", "integer", "BitXor[iv, 255]", "np.bitwise_xor(iv, 255)")
P("bitshiftleft", "integer", "BitShiftLeft[iv, 3]", "np.left_shift(iv, 3)")
P("integerlength", "integer", "IntegerLength[cv]",
  "np.floor(np.log10(cv)).astype(np.int64)+1")
P("powermod", "integer", "PowerMod[cv, 3, 1009]", "pow(3, 1, 1009)*0 + "
  "np.array([pow(int(x), 3, 1009) for x in cv])")
P("evenq_int", "integer", "Map[EvenQ, cv]", "cv % 2 == 0", chk="ckt[r]")
P("divisible", "integer", "Map[Divisible[#, 3] &, cv]", "cv % 3 == 0", chk="ckt[r]")
P("primeq", "integer", "Map[PrimeQ, cv]", None, chk="ckt[r]")
P("eulerphi", "integer", "EulerPhi[Take[cv, 20000]]", None)
P("moebiusmu", "integer", "MoebiusMu[Take[cv, 20000]]", None)
P("divisorsigma", "integer", "DivisorSigma[1, Take[cv, 20000]]", None)
P("fibonacci_arr", "integer", "Fibonacci[Take[cv, 2000]]", None)
P("factorial_int", "integer", "Factorial[Take[cv, 2000]]", None)
P("prime_arr", "integer", "Prime[Take[cv, 20000]]", None)
P("integerdigits", "integer", "IntegerDigits[Take[cv, 20000]]", None,
  chk="N[Length[r]]", pychk="float(len(r))")

# --- numerics ------------------------------------------------------------------
P("n_of_list", "numerics", "N[iv]", "iv.astype(np.float64)")
P("interpolation_build", "numerics",
  "Interpolation[Transpose[{Range[1., 10000.], Take[v, 10000]}]]", None,
  chk="0.", pychk="0.0")
P("fit_linear", "numerics",
  "Fit[Transpose[{Range[1., 10000.], Take[v, 10000]}], {1, x}, x]", None,
  chk="0.", pychk="0.0")

GROUPS = sorted({p["group"] for p in PROBES})


# ---------------------------------------------------------------------------
# Mathilda side
# ---------------------------------------------------------------------------

def m_source(probes):
    """One timing block per probe.

    The checksum is bounded before it is printed. ck[r] is Total[Flatten[...]],
    which reduces to a number for every result the probe is supposed to
    produce -- but a head that comes back UNEVALUATED leaves Total unevaluated
    too, and Print then dumps the whole 10^6-element expression. That is not a
    hypothetical: the first run of this sweep wrote 103 MB before the guard
    existed. NumberQ collapses it to a marker, which is also the signal the
    probe wants: "this head did not answer".
    """
    out = [MPRE]
    for pr in probes:
        chk = pr["chk"] if pr["chk"] != "0." else "0."
        out.append(f"""
Clear[r];
tt = Table[First[AbsoluteTiming[r = {pr['m']};]], {{{REPS}}}];
cc = {chk};
Print["@", "{pr['id']}", "\\t", Min[tt], "\\t", Max[tt], "\\t",
      If[NumberQ[cc], cc, "UNEVAL"]];
Clear[r, tt, cc];
""")
    return "\n".join(out)


def parse_m(text):
    out = {}
    for line in text.splitlines():
        if not line.startswith("@"):
            continue
        parts = line[1:].split("\t")
        if len(parts) < 4:
            continue
        pid = parts[0].strip()
        try:
            out[pid] = {
                "min": float(parts[1]),
                "max": float(parts[2]),
                "chk": parts[3].strip(),
            }
        except ValueError:
            continue
    return out


def run_mathilda(probes, env_extra=None, tmpdir=None):
    """Batched, with a per-probe retry so one crash costs one row, not twelve."""
    results = {}
    env = dict(os.environ)
    if env_extra:
        env.update(env_extra)

    def run_group(group):
        src = m_source(group)
        path = os.path.join(tmpdir, "probe.m")
        with open(path, "w") as fh:
            fh.write(src)
        try:
            cp = subprocess.run([MATHILDA, "-file", path], capture_output=True,
                                text=True, timeout=TIMEOUT, env=env)
            return parse_m(cp.stdout)
        except subprocess.TimeoutExpired:
            return {}

    for i in range(0, len(probes), BATCH):
        chunk = probes[i:i + BATCH]
        got = run_group(chunk)
        missing = [p for p in chunk if p["id"] not in got]
        results.update(got)
        if missing and len(chunk) > 1:
            # A crash or a timeout took the rest of the batch with it. Re-run
            # the ones that did not report, one process each, so the row that
            # actually failed is the only one lost.
            for p in missing:
                results.update(run_group([p]))
                sys.stderr.write(f"  mathilda retry {p['id']}\n")
        if CHECKPOINT:
            with open(CHECKPOINT, "w") as fh:
                json.dump(results, fh)
        sys.stderr.write(f"  mathilda {min(i+BATCH, len(probes))}/{len(probes)}\n")
    return results


# ---------------------------------------------------------------------------
# Python side
# ---------------------------------------------------------------------------

PY_RUNNER = r"""
def _time(fn, reps):
    ts = []
    try:
        fn()
    except Exception as e:
        return None, None, "ERR:" + type(e).__name__
    for _ in range(reps):
        t0 = time.perf_counter()
        r = fn()
        ts.append(time.perf_counter() - t0)
    return min(ts), max(ts), r

OUT = {}
"""


def py_source(probes):
    out = [PYPRE, PY_RUNNER]
    for pr in probes:
        if not pr["py"]:
            continue
        chk = pr["pychk"]
        out.append(f"""
try:
    _lo, _hi, r = _time(lambda: {pr['py']}, {REPS})
    if _lo is None:
        OUT[{pr['id']!r}] = {{"err": r}}
    else:
        OUT[{pr['id']!r}] = {{"min": _lo, "max": _hi, "chk": {chk}}}
except Exception as _e:
    OUT[{pr['id']!r}] = {{"err": type(_e).__name__ + ": " + str(_e)[:60]}}
""")
    out.append("print('@@JSON@@' + json.dumps(OUT))")
    return "\n".join(out)


def run_python(probes, tmpdir):
    path = os.path.join(tmpdir, "probe.py")
    with open(path, "w") as fh:
        fh.write(py_source(probes))
    try:
        cp = subprocess.run([PYTHON, path], capture_output=True, text=True,
                            timeout=TIMEOUT * 3)
    except subprocess.TimeoutExpired:
        return {}
    for line in cp.stdout.splitlines():
        if line.startswith("@@JSON@@"):
            return json.loads(line[len("@@JSON@@"):])
    sys.stderr.write(cp.stderr[-2000:] + "\n")
    return {}


# ---------------------------------------------------------------------------

def agree(a, b, tol=1e-5):
    """Do the two systems' checksums describe the same computation?

    THE TOLERANCE IS SET BY THE PRINTING, NOT BY THE ARITHMETIC. Mathilda
    prints six significant figures; NumPy's repr gives seventeen. So the same
    number arrives here as "1.1752e+06" and "1175201.4651842169", which differ
    by 1.2e-6 relative -- and at tol=1e-6 the first run of this sweep reported
    eleven such rows as value mismatches while every one of them agreed.

    performance.md §8 records the identical bug in hpc_bench.py ("the harness
    disagreed with itself about agreement"), which is why it is worth a comment
    rather than a quiet constant: any comparison against a printed value must be
    no tighter than the shorter printing, and six figures is 1e-5.

    A checksum is a sum over up to 10^6 elements, so this is a check that the two
    systems ran the same ALGORITHM, not that they agree to the last ulp.
    """
    if a is None or b is None:
        return None
    if a == "UNEVAL":
        return False          # the head did not answer; never "agreement"
    try:
        x, y = float(a), float(b)
    except (TypeError, ValueError):
        return str(a).strip() == str(b).strip()
    if x == y:
        return True
    import math
    if math.isinf(x) and math.isinf(y):
        return (x > 0) == (y > 0)
    if math.isnan(x) and math.isnan(y):
        return True
    d = abs(x - y)
    s = max(abs(x), abs(y))
    return d <= tol * s if s else d <= tol


def fmt_t(t):
    if t is None:
        return "—"
    if t < 1e-3:
        return f"{t*1e6:.1f} µs"
    if t < 1.0:
        return f"{t*1e3:.2f} ms"
    return f"{t:.3f} s"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="comma-separated probe ids")
    ap.add_argument("--group", help="comma-separated groups: " + ",".join(GROUPS))
    ap.add_argument("--json", metavar="PATH")
    ap.add_argument("--diag", action="store_true",
                    help="add a MATHILDA_NO_PACK=1 column")
    ap.add_argument("--system", default="mathilda,python")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--emit", metavar="DIR",
                    help="write the .m and .py the sweep would run, so the "
                         "probe set is readable and runnable on its own")
    args = ap.parse_args()

    probes = PROBES
    if args.group:
        want = set(args.group.split(","))
        probes = [p for p in probes if p["group"] in want]
    if args.only:
        want = set(args.only.split(","))
        probes = [p for p in probes if p["id"] in want]

    if args.list:
        for p in probes:
            print(f"{p['id']:<26} {p['group']:<14} {p['m']}")
        print(f"\n{len(probes)} probes in {len(GROUPS)} groups")
        return 0

    if args.emit:
        os.makedirs(args.emit, exist_ok=True)
        mpath = os.path.join(args.emit, "numeric_coverage_sweep.m")
        ppath = os.path.join(args.emit, "numeric_coverage_sweep.py")
        with open(mpath, "w") as fh:
            fh.write("(* Generated by tools/numeric_sweep.py --emit.\n"
                     "   Runs unmodified in Mathilda and in Mathematica, which\n"
                     "   is the point: the two CAS columns cannot silently be\n"
                     "   timing different programs. Each line prints\n"
                     "   @id <tab> min <tab> max <tab> checksum. *)\n")
            fh.write(m_source(probes))
        with open(ppath, "w") as fh:
            fh.write('"""Generated by tools/numeric_sweep.py --emit.\n\n'
                     "The same algorithms in NumPy/SciPy, on the same data, in\n"
                     "the same order. Prints one JSON object of\n"
                     "{id: {min, max, chk}}.\n"
                     '"""\n')
            fh.write(py_source(probes))
        print(f"wrote {mpath}\nwrote {ppath}")
        return 0

    systems = set(args.system.split(","))
    mres, pres, dres = {}, {}, {}

    with tempfile.TemporaryDirectory() as td:
        if "mathilda" in systems:
            sys.stderr.write(f"mathilda: {len(probes)} probes\n")
            mres = run_mathilda(probes, tmpdir=td)
        if "python" in systems:
            sys.stderr.write("python: running\n")
            pres = run_python(probes, td)
        if args.diag:
            sys.stderr.write("mathilda (no pack): running\n")
            dres = run_mathilda(probes, env_extra={"MATHILDA_NO_PACK": "1"},
                                tmpdir=td)

    rows = []
    for p in probes:
        m = mres.get(p["id"], {})
        q = pres.get(p["id"], {})
        d = dres.get(p["id"], {})
        mt = m.get("min")
        pt = q.get("min")
        row = {
            "id": p["id"], "group": p["group"], "expr": p["m"],
            "mathilda": mt, "mathilda_max": m.get("max"),
            "python": pt, "nopack": d.get("min"),
            "mathilda_chk": m.get("chk"), "python_chk": q.get("chk"),
            "agree": agree(m.get("chk"), q.get("chk")),
            "ratio": (mt / pt) if (mt and pt) else None,
            "packgain": (d.get("min") / mt) if (d.get("min") and mt) else None,
            "py_err": q.get("err"),
        }
        rows.append(row)

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(rows, fh, indent=1)

    hdr = f"{'probe':<26} {'group':<13} {'Mathilda':>11} {'NumPy':>11} {'ratio':>9}"
    if args.diag:
        hdr += f" {'no-pack':>11} {'gain':>8}"
    hdr += "  ok"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        rat = "—"
        if r["ratio"]:
            rat = f"{r['ratio']:.2f}x" if r["ratio"] >= 1 else f"1/{1/r['ratio']:.2f}x"
        line = (f"{r['id']:<26} {r['group']:<13} {fmt_t(r['mathilda']):>11} "
                f"{fmt_t(r['python']):>11} {rat:>9}")
        if args.diag:
            g = f"{r['packgain']:.1f}x" if r["packgain"] else "—"
            line += f" {fmt_t(r['nopack']):>11} {g:>8}"
        ok = {True: "  ", False: " ✗", None: " ?"}[r["agree"]]
        line += ok
        print(line)

    missing = [r["id"] for r in rows if r["mathilda"] is None]
    bad = [r for r in rows if r["agree"] is False]
    slow = [r for r in rows if r["ratio"] and r["ratio"] > 3.0]
    print()
    print(f"{len(rows)} probes; {len(missing)} no Mathilda result; "
          f"{len(bad)} value mismatches; {len(slow)} more than 3x behind NumPy")
    if missing:
        print("  no result: " + ", ".join(missing))
    if bad:
        for r in bad:
            m = str(r["mathilda_chk"])[:40]
            p = str(r["python_chk"])[:40]
            print(f"  MISMATCH {r['id']}: mathilda={m} numpy={p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

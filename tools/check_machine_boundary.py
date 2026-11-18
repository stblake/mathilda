#!/usr/bin/env python3
"""
check_machine_boundary.py -- do all five surfaces agree at the edges of the
machine range?

WHAT IT GUARDS. The numeric-tower review (docs/design/numeric_tower_review.md)
fixed one policy at one set of choke points: a machine number has a 53-bit
mantissa and a normal exponent, and anything an operation produces outside
[DBL_MIN, DBL_MAX] (or as a subnormal, or as zero from nonzero operands) is
recomputed as a 53-bit MPFR number. Several of its findings were the SAME
bug on a different surface:

    C3  packed {1.*^308} * 10 gave {inf.0} while the scalar promoted
    E5  packed {1., 2.} * 10^400 gave {inf.0, inf.0}
    E9  Compile[{{x, _Real}}, x x][1.*^-200] gave 0. (the interpreter 1.*^-400)

Each surface has its own fast path -- a scalar builtin, the boxed List walk, a
packed buffer behind the transparency gate, a visible NDArray, and the Compile
VM -- so a policy honoured on one can be broken on another, and none of the
existing audits looks at boundary VALUES: check-nd-surfaces runs ordinary data.

WHAT IT DOES. For every head with an element-wise NDArray kernel (the REG_U /
REG_B registry, read from the source as numeric_coverage.py does) plus the core
arithmetic heads that ndarray.c dispatches directly, and for every boundary
input, it evaluates

    S    the scalar           H[a]                    (and H[filler])
    L    a boxed List          H[{a, filler}]
    P    a packed List         H[ToPackedArray[{a, filler}]]
    N    a visible NDArray     Normal[H[NDArray[{a, filler}]]]
    CS   Compile, scalar       Compile[{{x, _Real}}, H[x]][a]
    CV   Compile, rank 1       Compile[{{x, _Real, 1}}, H[x]][{a, filler}]

and requires every lane to equal the scalar answer for that element, as an
InputForm string -- bit for bit, as CLAUDE.md requires of the surfaces. A
binary head is probed with the array in either position. Inputs include the
overflow and underflow edges (+-1.*^308, 800., 1.*^-200, DBL_MIN, -0.) for
real heads and the int64 edges (2^63-1, -2^63, 3037000500) for integer ones,
and a huge exact operand (10^400) against a machine array for binary heads.

A head whose Compile[] lowering does not exist falls back to the interpreter
and so agrees trivially; the gate is about values, and check-compile-coverage
owns lowering. N-ary kernels (the hypergeometric and elliptic families) are
not probed: their buffer paths vary a parameter vector, not a boundary value.

RATCHET. A disagreement is a finding keyed `Head|probe`. Findings already in
BASELINE are reported but do not fail; a new one fails; a BASELINE entry that
no longer reproduces is reported so the line can be deleted.

Usage:
    python3 tools/check_machine_boundary.py               # the gate
    python3 tools/check_machine_boundary.py --only Exp,Times -v
    python3 tools/check_machine_boundary.py --list        # every finding
"""
import argparse
import concurrent.futures as cf
import os
import re
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numeric_coverage as nc  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "Mathilda")
TIMEOUT = 900
# Per-expression budget. A scalar answer that is itself out of reach (an exact
# Fibonacci[2^63-1]) is not a surface question, so that probe is skipped.
PROBE_SECONDS = 10
# Per-process resident-memory cap. TimeConstrained cannot interrupt a single
# GMP call (an exact Power[3037000500, 3037000500] is ~12 GB and one mpz_pow),
# and Mathilda has no working MemoryConstrained, so the driver polls each
# child's RSS and kills it past this cap. A first full run without it reached
# 47 GB. The head is then reported as MEMCAP, a finding like any other.
MEM_CAP_MB = int(os.environ.get("MB_MEM_CAP_MB", "2048"))
POLL_SECONDS = 0.2

# Core heads that ndarray.c / the Compile VM handle without a registered kernel.
CORE_UNARY = ["Minus", "Sqrt", "Abs", "Sign", "Exp", "Log"]
# Max/Min are not here: they are not element-wise (Max[{a, f}, b] is one scalar),
# so the lane comparison does not apply; the reductions own them.
CORE_BINARY = ["Plus", "Times", "Power", "Subtract", "Divide",
               "Mod", "Quotient", "ArcTan", "Log"]

# Registry entries that are not heads (macro parameters picked up by the scan).
NOT_HEADS = {"NAME"}

# Real boundary inputs, as literal text. The filler keeps the second lane in
# every head's ordinary domain-free range.
REAL_INPUTS = [
    ("big", "1.*^308"), ("-big", "-1.*^308"), ("exp", "800."), ("-exp", "-800."),
    ("tiny", "1.*^-200"), ("-tiny", "-1.*^-200"),
    ("dblmin", "2.2250738585072014*^-308"), ("negzero", "-0."),
]
REAL_FILLER = "0.5"
INT_INPUTS = [("imax", "9223372036854775807"), ("imin", "-9223372036854775808"),
              ("isqrt", "3037000500"), ("i62", "4611686018427387904")]
INT_FILLER = "3"

# Binary probes: (id, a, b). The array goes in one position, the other operand
# stays a scalar.
REAL_PAIRS = [
    ("big*10", "1.*^308", "10."), ("-big*10", "-1.*^308", "10."),
    ("big,big", "1.*^308", "1.*^308"), ("tiny,tiny", "1.*^-200", "1.*^-200"),
    ("dblmin,tiny", "2.2250738585072014*^-308", "1.*^-10"),
    ("1,huge", "1.", "10^400"), ("huge,2", "10^400", "2."),
    ("2,-big", "2.", "-1.*^308"),
]
INT_PAIRS = [("imax,2", "9223372036854775807", "2"),
             ("imin,-1", "-9223372036854775808", "-1"),
             ("isqrt,isqrt", "3037000500", "3037000500")]
# Heads whose exact answer at an int64-edge pair is astronomically large
# (Power[3037000500, 3037000500] has ~10^11 bits). That is an exact-arithmetic
# size question, not a machine-boundary one, so these heads get the real pairs
# and the bounded int pair "imax,2" only.
INT_GROWTH = {"Power", "Binomial", "Pochhammer", "Factorial2", "Beta",
              "Multinomial", "StirlingS1", "StirlingS2", "BellB", "Fibonacci",
              "LucasL", "Factorial", "Gamma", "Subfactorial", "BitSet",
              "BitShiftLeft", "PartitionsP", "PartitionsQ", "HurwitzZeta"}

# ---------------------------------------------------------------------------
# Known disagreements. Each line is `Head|probe|surface` with a reason; delete
# a line when the check reports it fixed.
# Ids G1-G10 are the gate's findings, recorded with repros in
# docs/design/numeric_tower_review.md section G.
# ---------------------------------------------------------------------------
BASELINE = {
    "ArcCos|u:-tiny|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCos|u:dblmin|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCos|u:negzero|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCos|u:tiny|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:-big|N1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:-big|P1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:-exp|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:-exp|N1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:-exp|P1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:-tiny|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:-tiny|N1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:-tiny|P1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:big|N1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:big|P1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:dblmin|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:dblmin|N1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:dblmin|P1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:exp|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:exp|N1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:exp|P1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:negzero|CS0": "G7: signed zero: -0. read as 0 by the scalar, as -0. by the kernels",
    "ArcCot|u:negzero|CV0": "G7: signed zero: -0. read as 0 by the scalar, as -0. by the kernels",
    "ArcCot|u:negzero|CV1": "G7: signed zero: -0. read as 0 by the scalar, as -0. by the kernels",
    "ArcCot|u:negzero|N0": "G7: signed zero: -0. read as 0 by the scalar, as -0. by the kernels",
    "ArcCot|u:negzero|N1": "G7: signed zero: -0. read as 0 by the scalar, as -0. by the kernels",
    "ArcCot|u:negzero|P0": "G7: signed zero: -0. read as 0 by the scalar, as -0. by the kernels",
    "ArcCot|u:negzero|P1": "G7: signed zero: -0. read as 0 by the scalar, as -0. by the kernels",
    "ArcCot|u:tiny|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:tiny|N1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcCot|u:tiny|P1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcSin|u:-tiny|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcSin|u:dblmin|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcSin|u:negzero|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcSin|u:tiny|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTanh|u:-big|N0": "G5: buffer surface leaks a subnormal machine real",
    "ArcTanh|u:-big|P0": "G5: buffer surface leaks a subnormal machine real",
    "ArcTanh|u:big|N0": "G5: buffer surface leaks a subnormal machine real",
    "ArcTanh|u:big|P0": "G5: buffer surface leaks a subnormal machine real",
    "ArcTan|u:-big|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:-big|N1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:-big|P1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:-exp|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:-exp|N1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:-exp|P1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:-tiny|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:-tiny|N1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:-tiny|P1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:big|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:big|N1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:big|P1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:dblmin|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:dblmin|N1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:dblmin|P1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:exp|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:exp|N1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:exp|P1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:negzero|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:negzero|N1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:negzero|P1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:tiny|CV1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:tiny|N1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "ArcTan|u:tiny|P1": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "Arg|u:negzero|CS0": "G7: signed zero: -0. read as 0 by the scalar, as -0. by the kernels",
    "Arg|u:negzero|CV0": "G7: signed zero: -0. read as 0 by the scalar, as -0. by the kernels",
    "Arg|u:negzero|N0": "G7: signed zero: -0. read as 0 by the scalar, as -0. by the kernels",
    "Arg|u:negzero|P0": "G7: signed zero: -0. read as 0 by the scalar, as -0. by the kernels",
    "Binomial|b:2,-big|CS0": "G8: Binomial[2., -1.*^308]: scalar Indeterminate, kernels 0.",
    "Binomial|b:2,-big|CV0": "G8: Binomial[2., -1.*^308]: scalar Indeterminate, kernels 0.",
    "Binomial|b:2,-big|CV1": "G8: Binomial[2., -1.*^308]: scalar Indeterminate, kernels 0.",
    "Binomial|b:2,-big|N0": "G8: Binomial[2., -1.*^308]: scalar Indeterminate, kernels 0.",
    "Binomial|b:2,-big|N1": "G8: Binomial[2., -1.*^308]: scalar Indeterminate, kernels 0.",
    "Binomial|b:2,-big|P0": "G8: Binomial[2., -1.*^308]: scalar Indeterminate, kernels 0.",
    "Binomial|b:2,-big|P1": "G8: Binomial[2., -1.*^308]: scalar Indeterminate, kernels 0.",
    "Binomial|b:2,-big|rN0": "G8: Binomial[2., -1.*^308]: scalar Indeterminate, kernels 0.",
    "Binomial|b:2,-big|rP0": "G8: Binomial[2., -1.*^308]: scalar Indeterminate, kernels 0.",
    "Binomial|b:dblmin,tiny|CV1": "G8: Binomial[2., -1.*^308]: scalar Indeterminate, kernels 0.",
    "Binomial|b:dblmin,tiny|N1": "G8: Binomial[2., -1.*^308]: scalar Indeterminate, kernels 0.",
    "Binomial|b:dblmin,tiny|P1": "G8: Binomial[2., -1.*^308]: scalar Indeterminate, kernels 0.",
    "Binomial|b:dblmin,tiny|rN1": "G8: Binomial[2., -1.*^308]: scalar Indeterminate, kernels 0.",
    "Binomial|b:dblmin,tiny|rP1": "G8: Binomial[2., -1.*^308]: scalar Indeterminate, kernels 0.",
    "Binomial|b:tiny,tiny|rN1": "G8: Binomial[2., -1.*^308]: scalar Indeterminate, kernels 0.",
    "Binomial|b:tiny,tiny|rP1": "G8: Binomial[2., -1.*^308]: scalar Indeterminate, kernels 0.",
    "Csch|u:-big|CS0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Csch|u:-big|CV0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Csch|u:-exp|CS0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Csch|u:-exp|CV0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Csch|u:big|CS0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Csch|u:big|CV0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Csch|u:exp|CS0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Csch|u:exp|CV0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Erfc|u:big|CS0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Erfc|u:big|CV0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Erfc|u:exp|CS0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Erfc|u:exp|CV0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "FresnelS|u:-tiny|N0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "FresnelS|u:-tiny|P0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "FresnelS|u:dblmin|N0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "FresnelS|u:dblmin|P0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "FresnelS|u:tiny|N0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "FresnelS|u:tiny|P0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "HurwitzZeta|b:big,big|N0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "HurwitzZeta|b:big,big|P0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "InverseErf|u:dblmin|N0": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "InverseErf|u:dblmin|P0": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "InverseJacobiDN|b:dblmin,tiny|N0": "G2: packed/NDArray InverseJacobiDN drops the imaginary part (ordinary data)",
    "InverseJacobiDN|b:dblmin,tiny|N1": "G2: packed/NDArray InverseJacobiDN drops the imaginary part (ordinary data)",
    "InverseJacobiDN|b:dblmin,tiny|P0": "G2: packed/NDArray InverseJacobiDN drops the imaginary part (ordinary data)",
    "InverseJacobiDN|b:dblmin,tiny|P1": "G2: packed/NDArray InverseJacobiDN drops the imaginary part (ordinary data)",
    "InverseJacobiDN|b:dblmin,tiny|rN0": "G2: packed/NDArray InverseJacobiDN drops the imaginary part (ordinary data)",
    "InverseJacobiDN|b:dblmin,tiny|rN1": "G2: packed/NDArray InverseJacobiDN drops the imaginary part (ordinary data)",
    "InverseJacobiDN|b:dblmin,tiny|rP0": "G2: packed/NDArray InverseJacobiDN drops the imaginary part (ordinary data)",
    "InverseJacobiDN|b:dblmin,tiny|rP1": "G2: packed/NDArray InverseJacobiDN drops the imaginary part (ordinary data)",
    "JacobiAmplitude|b:dblmin,tiny|rN1": "G5: buffer surface leaks a subnormal machine real",
    "JacobiAmplitude|b:dblmin,tiny|rP1": "G5: buffer surface leaks a subnormal machine real",
    "Log10|u:-big|N0": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "Log10|u:-big|P0": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "Log10|u:-exp|N0": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "Log10|u:-exp|P0": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "Log10|u:-tiny|N0": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "Log10|u:-tiny|P0": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "Log2|u:-big|N0": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "Log2|u:-big|P0": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "Log2|u:-exp|N0": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "Log2|u:-exp|P0": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "Log2|u:-tiny|N0": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "Log2|u:-tiny|P0": "G4: buffer/Compile kernel 1 ulp off the correctly rounded scalar",
    "Log|b:big,big|N0": "G9: interpreter Log[x, x] is exact 1 for machine x",
    "Log|b:big,big|P0": "G9: interpreter Log[x, x] is exact 1 for machine x",
    "Log|b:big,big|rN0": "G9: interpreter Log[x, x] is exact 1 for machine x",
    "Log|b:big,big|rP0": "G9: interpreter Log[x, x] is exact 1 for machine x",
    "Log|b:tiny,tiny|N0": "G9: interpreter Log[x, x] is exact 1 for machine x",
    "Log|b:tiny,tiny|P0": "G9: interpreter Log[x, x] is exact 1 for machine x",
    "Log|b:tiny,tiny|rN0": "G9: interpreter Log[x, x] is exact 1 for machine x",
    "Log|b:tiny,tiny|rP0": "G9: interpreter Log[x, x] is exact 1 for machine x",
    "Mod|b:dblmin,tiny|CV1": "G3: packed/NDArray/Compile Mod[x, 1.*^-10] is 0. (scalar exact fmod)",
    "Mod|b:dblmin,tiny|N1": "G3: packed/NDArray/Compile Mod[x, 1.*^-10] is 0. (scalar exact fmod)",
    "Mod|b:dblmin,tiny|P1": "G3: packed/NDArray/Compile Mod[x, 1.*^-10] is 0. (scalar exact fmod)",
    "Sech|u:-big|CS0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Sech|u:-big|CV0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Sech|u:-exp|CS0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Sech|u:-exp|CV0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Sech|u:big|CS0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Sech|u:big|CV0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Sech|u:exp|CS0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Sech|u:exp|CV0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "SphericalBesselJ|b:big*10|CS0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "SphericalBesselJ|b:big*10|CV0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "SphericalBesselJ|b:big*10|N0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "SphericalBesselJ|b:big*10|P0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "SphericalBesselJ|b:big*10|rN0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "SphericalBesselJ|b:big*10|rN1": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "SphericalBesselJ|b:big*10|rP0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "SphericalBesselJ|b:big*10|rP1": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "StruveH|b:-big*10|CS0": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveH|b:-big*10|CV0": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveH|b:-big*10|CV1": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveH|b:-big*10|N0": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveH|b:-big*10|N1": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveH|b:-big*10|P0": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveH|b:-big*10|P1": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveH|b:-big*10|rL0": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveH|b:-big*10|rL1": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveH|b:-big*10|rN0": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveH|b:-big*10|rN1": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveH|b:-big*10|rP0": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveH|b:-big*10|rP1": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveL|b:-big*10|CS0": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveL|b:-big*10|CV0": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveL|b:-big*10|CV1": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveL|b:-big*10|N0": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveL|b:-big*10|N1": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveL|b:-big*10|P0": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveL|b:-big*10|P1": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveL|b:-big*10|rL0": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveL|b:-big*10|rL1": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveL|b:-big*10|rN0": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveL|b:-big*10|rN1": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveL|b:-big*10|rP0": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "StruveL|b:-big*10|rP1": "G10: order -1.*^308: array/Compile surfaces exceed the probe budget",
    "Subtract|b:imax,2|CV0": "G1: Compile integer rank-1 Subtract[v, s] is garbage (ordinary data)",
    "Subtract|b:imax,2|CV1": "G1: Compile integer rank-1 Subtract[v, s] is garbage (ordinary data)",
    "Subtract|b:isqrt,isqrt|CV0": "G1: Compile integer rank-1 Subtract[v, s] is garbage (ordinary data)",
    "Subtract|b:isqrt,isqrt|CV1": "G1: Compile integer rank-1 Subtract[v, s] is garbage (ordinary data)",
    "Zeta|b:big,big|N0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
    "Zeta|b:big,big|P0": "G6: underflow flushed to 0. instead of Underflow[] / an MPFR value",
}

PRELUDE = r"""
sh[e_] := ToString[e, InputForm];
lanes[r_] := If[Head[r] === List && Length[r] == 2, sh /@ r,
  {"SHAPE:" <> sh[r], "SHAPE:" <> sh[r]}];
emit[id_, surf_, s_] := Print["@@", id, "\t", surf, "\t", s];
emit2[id_, surf_, r_] := (emit[id, surf <> "0", lanes[r][[1]]];
                          emit[id, surf <> "1", lanes[r][[2]]]);
pk[l_] := Developer`ToPackedArray[l];
nd[l_] := NDArray[l, DataType -> If[VectorQ[l, IntegerQ], "int64", "float64"]];
nrm[r_] := If[Head[r] === NDArray, Normal[r], r];
(* An array lane NDArray cannot hold (an exact 10^400) leaves the constructor
   unevaluated; that surface is then not applicable, not a disagreement. *)
SetAttributes[ndrun, HoldFirst];
ndrun[e_, l_] := If[NDArrayQ[nd[l]], nrm[e], {"SKIP", "SKIP"}];
"""


def unary_block(h, pid, a, f, ty):
    q = lambda e: "Quiet[TimeConstrained[" + e + ", " + str(PROBE_SECONDS) + ", \"TC\"]]"
    return "\n".join(([f'emit["{pid}", "INT", "1"];'] if ty == "_Integer" else []) + [
        f'emit["{pid}", "S0", sh[{q(f"{h}[{a}]")}]];',
        f'emit["{pid}", "S1", sh[{q(f"{h}[{f}]")}]];',
        f'emit2["{pid}", "L", {q(f"{h}[{{{a}, {f}}}]")}];',
        f'emit2["{pid}", "P", {q(f"{h}[pk[{{{a}, {f}}}]]")}];',
        f'emit2["{pid}", "N", {q(f"ndrun[{h}[nd[{{{a}, {f}}}]], {{{a}, {f}}}]")}];',
        f'emit["{pid}", "CS0", sh[{q(f"Compile[{{{{x, {ty}}}}}, {h}[x]][{a}]")}]];',
        f'emit2["{pid}", "CV", {q(f"Compile[{{{{x, {ty}, 1}}}}, {h}[x]][{{{a}, {f}}}]")}];',
    ])


def binary_block(h, pid, a, b, f, ty):
    q = lambda e: "Quiet[TimeConstrained[" + e + ", " + str(PROBE_SECONDS) + ", \"TC\"]]"
    out = ([f'emit["{pid}", "INT", "1"];'] if ty == "_Integer" else []) + [
        f'emit["{pid}", "S0", sh[{q(f"{h}[{a}, {b}]")}]];',
        f'emit["{pid}", "S1", sh[{q(f"{h}[{f}, {b}]")}]];',
        f'emit["{pid}", "T1", sh[{q(f"{h}[{a}, {f}]")}]];',
        f'emit2["{pid}", "L", {q(f"{h}[{{{a}, {f}}}, {b}]")}];',
        f'emit2["{pid}", "P", {q(f"{h}[pk[{{{a}, {f}}}], {b}]")}];',
        f'emit2["{pid}", "N", {q(f"ndrun[{h}[nd[{{{a}, {f}}}], {b}], {{{a}, {f}}}]")}];',
        f'emit2["{pid}", "rL", {q(f"{h}[{a}, {{{b}, {f}}}]")}];',
        f'emit2["{pid}", "rP", {q(f"{h}[{a}, pk[{{{b}, {f}}}]]")}];',
        f'emit2["{pid}", "rN", {q(f"ndrun[{h}[{a}, nd[{{{b}, {f}}}]], {{{b}, {f}}}]")}];',
        f'emit["{pid}", "rS1", sh[{q(f"{h}[{a}, {f}]")}]];',
        f'emit["{pid}", "CS0", sh[{q(f"Compile[{{{{x, {ty}}}, {{y, {ty}}}}}, {h}[x, y]][{a}, {b}]")}]];',
        f'emit2["{pid}", "CV", {q(f"Compile[{{{{x, {ty}, 1}}, {{y, {ty}}}}}, {h}[x, y]][{{{a}, {f}}}, {b}]")}];',
    ]
    return "\n".join(out)


def script_for(h, arities):
    parts = [PRELUDE]
    if "unary" in arities:
        for pid, a in REAL_INPUTS:
            parts.append(unary_block(h, "u:" + pid, a, REAL_FILLER, "_Real"))
        for pid, a in INT_INPUTS:
            if h in INT_GROWTH:
                continue
            parts.append(unary_block(h, "u:" + pid, a, INT_FILLER, "_Integer"))
    if "binary" in arities:
        for pid, a, b in REAL_PAIRS:
            parts.append(binary_block(h, "b:" + pid, a, b, REAL_FILLER, "_Real"))
        for pid, a, b in INT_PAIRS:
            if h in INT_GROWTH and pid != "imax,2":
                continue
            parts.append(binary_block(h, "b:" + pid, a, b, INT_FILLER, "_Integer"))
    return "\n".join(parts) + "\n"


# An exact answer the scalar leaves symbolic (Sin[2^62]) is a machine number in
# Compile[] by design -- the compiled body is typed _Real throughout -- so on an
# _Integer probe the Compile surfaces are compared only where the scalar answer
# is itself a number.
NUMBER = re.compile(r"-?[0-9.]+(\*\^-?[0-9]+)?(`[0-9.]*)?")
EXACT_INT = re.compile(r"-?[0-9]+")
MACHINE = re.compile(r"-?[0-9]+\.[0-9]*(\*\^-?[0-9]+)?")


def machine_typed_equal(got, want):
    """A Compile lane is machine-typed by design: on an _Integer probe whose
    scalar answer is an exact integer, the compiled machine real that is that
    integer correctly rounded (Divide[6, 3] -> 2., Sqrt[2^62] -> 2147483648.)
    agrees. Anything else -- a wrong value, an exact integer that differs -- is
    still compared as a string."""
    if not (want and got and EXACT_INT.fullmatch(want) and MACHINE.fullmatch(got)):
        return False
    return float(got.replace("*^", "e")) == float(int(want))


def judge(rows):
    """rows: {probe: {surface: string}} -> list of (probe, surface, got, want)."""
    bad = []
    for pid, s in rows.items():
        want = {"0": s.get("S0"), "1": s.get("S1")}
        rwant = {"0": s.get("S0"), "1": s.get("rS1")}
        for surf, got in s.items():
            if surf in ("S0", "S1", "T1", "rS1", "INT"):
                continue
            lane = surf[-1]
            w = rwant[lane] if surf.startswith("r") else want[lane]
            if got == '"SKIP"' or w == '"TC"':
                continue
            # A list holds BOTH lanes: when the boundary lane's own scalar ran
            # out of budget, the whole list does too, and lane 1 has nothing to
            # say. (HermiteH[-1.*^308, 10.] alone takes 13 s.)
            if got == 'SHAPE:"TC"' and '"TC"' in (s.get("S0"), s.get("T1")):
                continue
            if surf.startswith("C") and s.get("INT") and not NUMBER.fullmatch(w or ""):
                continue
            if surf.startswith("C") and machine_typed_equal(got, w):
                continue
            if got != w:
                bad.append((pid, surf, got, w))
    return bad


def rss_mb(pid):
    """Resident set size of pid in MB (ps works on macOS and Linux)."""
    try:
        r = subprocess.run(["ps", "-o", "rss=", "-p", str(pid)],
                           capture_output=True, text=True)
        return int(r.stdout.strip() or 0) // 1024
    except (ValueError, OSError):
        return 0


def watch(p):
    """Wait for p, killing it past TIMEOUT or MEM_CAP_MB. Returns an error tag."""
    t = 0.0
    while p.poll() is None:
        if rss_mb(p.pid) > MEM_CAP_MB:
            p.kill(); p.wait()
            return f"MEMCAP>{MEM_CAP_MB}MB"
        if t > TIMEOUT:
            p.kill(); p.wait()
            return "TIMEOUT"
        time.sleep(POLL_SECONDS)
        t += POLL_SECONDS
    return None


def run_head(h, arities):
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as fh:
        fh.write(script_for(h, arities))
        path = fh.name
    outf = tempfile.TemporaryFile("w+")
    try:
        p = subprocess.Popen([BIN, "-file", path], stdout=outf,
                             stderr=subprocess.DEVNULL, text=True, cwd=ROOT)
        err = watch(p)
        outf.seek(0)
        out, rc = outf.read(), p.returncode
    finally:
        os.unlink(path)
        outf.close()
    if err:
        return h, None, err
    rows = {}
    for line in out.splitlines():
        if not line.startswith("@@"):
            continue
        parts = line[2:].split("\t", 2)
        if len(parts) == 3:
            rows.setdefault(parts[0], {})[parts[1]] = parts[2]
    if rc != 0 and not rows:
        return h, None, f"CRASH rc={rc}"
    return h, rows, None


def heads():
    k = nc.kernel_heads()
    out = {}
    for h, ar in k.items():
        if h in NOT_HEADS:
            continue
        a = {x for x in ar if x in ("unary", "binary")}
        if a:
            out[h] = set(a)
    for h in CORE_UNARY:
        out.setdefault(h, set()).add("unary")
    for h in CORE_BINARY:
        out.setdefault(h, set()).add("binary")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="comma-separated heads")
    ap.add_argument("--list", action="store_true", help="print every finding")
    ap.add_argument("-v", action="store_true", help="print got/want per finding")
    ap.add_argument("-j", type=int, default=max(1, min(4, (os.cpu_count() or 2) // 4)))
    args = ap.parse_args()
    if not os.path.exists(BIN):
        print("check-machine-boundary: SKIP (no ./Mathilda)")
        return 0
    hs = heads()
    if args.only:
        want = set(args.only.split(","))
        hs = {h: a for h, a in hs.items() if h in want}
    findings = {}
    with cf.ThreadPoolExecutor(args.j) as ex:
        for h, rows, err in ex.map(lambda kv: run_head(*kv), sorted(hs.items())):
            if err:
                findings[f"{h}|*|{err}"] = (err, "")
                continue
            for pid, surf, got, w in judge(rows):
                findings[f"{h}|{pid}|{surf}"] = (got, w)
    new = sorted(k for k in findings if k not in BASELINE)
    fixed = sorted(k for k in BASELINE if k not in findings
                   and (not args.only or k.split("|")[0] in hs))
    shown = sorted(findings) if args.list else new
    for k in shown:
        tag = "" if k in new else " (baseline)"
        print(f"DISAGREE {k}{tag}")
        if args.v:
            got, w = findings[k]
            print(f"    got  {got}\n    want {w}")
    for k in fixed:
        print(f"FIXED {k} -- delete it from BASELINE")
    print(f"check-machine-boundary: {len(hs)} heads, {len(findings)} disagreements "
          f"({len(new)} new, {len(findings) - len(new)} baseline), {len(fixed)} fixed")
    return 1 if new else 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Kahan-style numerics stress corpus (generates tests/test_numeric_kahan.c).

The tests William Kahan would aim at a numeric system: IEEE 754 machine
properties in the manner of Paranoia (radix, precision, ties-to-even, Sterbenz,
exactness theorems), correct rounding of + - * / sqrt against the EXACT result,
faithful elementary functions at hard arguments (huge-argument reduction,
log1p-type and near-pole points), correctly rounded decimal<->binary conversion
and round trips, catastrophic-cancellation classics (Rump, Muller, the needle
triangle, the Ramanujan near-integers, compound interest) through N's adaptive
precision, branch-cut continuity and overflow-safe complex arithmetic, and
numerically stable statistics (two-pass variance).

Every case is a boolean whose truth follows from exact mathematics. Equal and
SameQ on approximate numbers are tolerant by design, so exactness is asserted
through exact binary values (SetPrecision[x, Infinity]) and an exact ulp.

Usage: python3 tools/numeric_kahan_cases.py OUT.m [tests/test_numeric_kahan.c]
xfail=<id> marks a known open bug (docs/design/numeric_tower_review.md).
"""
import sys

PRELUDE = [
    # exact binary exponent of a nonzero dyadic s: 2^kx <= |s| < 2^(kx+1)
    "kx[s_] := IntegerLength[Floor[Abs[s] 2^1200], 2] - 1201",
    # ulp of a machine-precision value (53-bit significand), exactly
    "kulp[x_] := 2^(kx[SetPrecision[x, Infinity]] - 52)",
    # exact equality of two approximate numbers (== and === are tolerant)
    "keq[a_, b_] := SetPrecision[a, Infinity] === SetPrecision[b, Infinity]",
    # x is the correctly rounded (to nearest) value of the exact q
    "kcr[x_, q_] := Abs[SetPrecision[x, Infinity] - q] <= kulp[x]/2",
    # x is within k ulps of the exact q
    "kfa[x_, q_, k_] := Abs[SetPrecision[x, Infinity] - q] <= k kulp[x]",
    # the exact value (to 60 digits) of f at the exact binary value of x
    "kref[f_, x_] := SetPrecision[N[f[SetPrecision[x, Infinity]], 60], Infinity]",
    # f[x] at machine precision is within k ulps of the true f(x)
    "kacc[f_, x_, k_] := kfa[f[x], kref[f, x], k]",
    # Sqrt[x] is correctly rounded: (s - u/2)^2 <= x <= (s + u/2)^2, exactly
    "kcrsqrt[x_] := With[{s = SetPrecision[Sqrt[x], Infinity], X = SetPrecision[x, Infinity], u = kulp[Sqrt[x]]}, (s - u/2)^2 <= X <= (s + u/2)^2]",
    # N certifies its digits: the 20- and 40-digit values agree to 18 digits
    "kagree[e_] := Abs[N[e, 20]/N[e, 40] - 1] < 10^-18",
    # InputForm round trip
    "krt[x_] := ToExpression[ToString[InputForm[x]]] === x",
]

CASES = []
def sec(name): CASES.append(("#", name))
def t(expr, xfail=None): CASES.append(("T", expr) if xfail is None else ("X", expr, xfail))

# ---------------------------------------------------------------- P. Paranoia
sec("P. machine arithmetic in the manner of Paranoia")
t("Module[{a = 1., b = 1.}, While[keq[((a + 1.) - a) - 1., 0.], a *= 2.]; While[keq[(a + b) - a, 0.], b += 1.]; keq[(a + b) - a, 2.]]")  # radix 2
t("Module[{p = 0, e = 1.}, While[!keq[1. + e, 1.], e /= 2.; p++]; p == 53]")                      # 53-bit significand
t("keq[$MachineEpsilon, 2.^-52]")
t("keq[1. + 2.^-53, 1.]")                                         # tie to even (down)
t("keq[1. + 3. 2.^-53, 1. + 2.^-51]")                             # tie to even (up)
t("With[{t = 2.^-53 + 2.^-60}, keq[1. + t, 1. + 2.^-52]]")        # above the tie: up
t("With[{t = 2.^-53 - 2.^-60}, keq[1. + t, 1.]]")                 # below the tie: down
t("keq[(1. - 2.^-54) + 2.^-54, 1.]")
t("keq[3. (1./3.), 1.]")                                          # 1 - 2^-54 exactly: tie, to 1
for x, y in [("1.7", "1.1"), ("0.75", "0.5"), ("1.*^300", "0.6*^300"), ("3.", "2."), ("1.0000000000000002", "0.9999999999999999")]:
    t(f"With[{{x = {x}, y = {y}}}, SetPrecision[x - y, Infinity] === SetPrecision[x, Infinity] - SetPrecision[y, Infinity]]")  # Sterbenz
for x in ["0.1", "1./3.", "1.*^-300", "1.*^300", "-7.25"]:
    # x/x = 1 exactly is an IEEE guarantee; a/b as a * (1/b) breaks it (K1)
    t(f"With[{{x = {x}}}, keq[x/x, 1.] && keq[x - x, 0.] && keq[x 1., x] && keq[x + 0., x]]",
      xfail="K1" if x == "1.*^-300" else None)
for x in ["0.1", "1./7.", "1.*^307"]:
    t(f"SetPrecision[2. {x}, Infinity] === 2 SetPrecision[{x}, Infinity]")
t("keq[(2.^53 - 1.) + 1., 2.^53] && keq[(2.^53 + 2.) - 2., 2.^53]")
t("keq[3. 7., 21.] && keq[2.^26 2.^26, 2.^52] && keq[(2.^26 + 1.)^2, 2.^52 + 2.^27 + 1.]")
t("MachineNumberQ[2.^1023 (2. - 2.^-52)] && !MachineNumberQ[2.^1023 2.]")
t("MachineNumberQ[2.^-1022] && !MachineNumberQ[2.^-1022/2.]")
t("keq[2.^-1022 2.^1022, 1.] && keq[2.^1023 2.^-1023, 1.]")

# ---------------------------------------------------------------- C. correctly rounded + - * /
sec("C. correctly rounded + - * / against the exact result")
OPS = [("+", "+"), ("-", "-"), ("*", "*"), ("/", "/")]
PAIRS = [("0.1", "0.2"), ("0.1", "0.3"), ("1./3.", "2./3."), ("1./7.", "1./3."),
         ("1.*^16 + 1.", "3."), ("123456.789", "0.001"), ("9.87654321*^200", "1.23456789*^-100"),
         ("1.0000000000000002", "0.9999999999999999"), ("2.^52 + 1.", "0.5"), ("3.141592653589793", "2.718281828459045"),
         ("1.*^-10", "7."), ("1.*^200", "1.*^200"), ("1.*^-200", "3.*^-200"), ("6.02214076*^23", "1.602176634*^-19"),
         ("-0.3", "0.1")]
K1_DIV = {("1.*^16 + 1.", "3."), ("1.*^-10", "7."), ("-0.3", "0.1")}   # found by this corpus
for a, b in PAIRS:
    for op, xop in OPS:
        xf = "K1" if (op == "/" and (a, b) in K1_DIV) else None
        t(f"With[{{a = {a}, b = {b}}}, kcr[a {op} b, SetPrecision[a, Infinity] {xop} SetPrecision[b, Infinity]]]", xfail=xf)
for a, b in PAIRS:   # Divide is one correctly rounded division
    t(f"With[{{a = {a}, b = {b}}}, kcr[Divide[a, b], SetPrecision[a, Infinity]/SetPrecision[b, Infinity]]]")
t("keq[-0.3/0.1, Divide[-0.3, 0.1]] && keq[(1.*^16 + 1.)/3., Divide[1.*^16 + 1., 3.]]", xfail="K1")

# ---------------------------------------------------------------- Q. correctly rounded sqrt
sec("Q. correctly rounded Sqrt")
for x in ["2.", "3.", "0.5", "10.", "1.*^-300", "1.*^300", "2.^52 + 1.", "0.1", "1./3.", "1.0000000000000002", "123456789.", "2.^1023"]:
    t(f"kcrsqrt[{x}]")
t("keq[Sqrt[9.], 3.] && keq[Sqrt[2.^100], 2.^50] && keq[Sqrt[0.25], 0.5]")
for x in ["0.1", "3.", "1.*^150", "1.*^-150", "1./3.", "123.456"]:
    t(f"keq[Sqrt[{x}^2], Abs[{x}]]")                                   # sqrt(x*x) = |x| exactly
for x, y in [("1.", "1.*^-20"), ("3.", "4."), ("1.*^-200", "1.*^-200"), ("0.1", "0.")]:
    t(f"SetPrecision[{x}/Sqrt[{x}^2 + {y}^2], Infinity] <= 1")

# ---------------------------------------------------------------- T. faithful elementary functions
sec("T. elementary functions within two ulps at hard arguments (the platform libm's bound)")
FUN = [
    ("Sin", ["1.*^22", "2.^1000", "1.*^300", "3.141592653589793", "1.*^-10", "100000.", "6.283185307179586", "1.5707963267948966"]),
    ("Cos", ["1.*^22", "2.^1000", "1.5707963267948966", "1.*^-8", "355.", "1.*^300"]),
    ("Tan", ["1.5707963267948966", "1.*^22", "0.7853981633974483", "1.*^-10", "1.*^300"]),
    ("Exp", ["709.", "-745.", "1.*^-20", "0.5", "100.", "-1.*^-20", "700.5"]),
    ("Log", ["1.0000000000000002", "0.9999999999999999", "1.*^300", "1.*^-300", "2.", "0.5", "1.*^-320"]),
    ("ArcTan", ["1.", "1.*^20", "1.*^-20", "0.4142135623730951"]),
    ("Sinh", ["1.*^-10", "20.", "700."]),
    ("Tanh", ["1.*^-10", "0.5", "19."]),
    ("ArcSin", ["0.5", "0.9999999", "1.*^-10"]),
    ("ArcCos", ["0.5", "-0.9999999"]),
]
for f, xs in FUN:
    for x in xs:
        t(f"kacc[{f}, {x}, 2]")

# ---------------------------------------------------------------- N. cancellation classics through N
sec("N. cancellation classics: N certifies its digits")
for e in ["Exp[Pi Sqrt[163]] - 640320^3 - 744", "Exp[Pi Sqrt[67]] - 5280^3 - 744", "Sqrt[10^20 + 1] - 10^10",
          "Cos[10^-20] - 1", "Log[1 + 10^-25]", "Exp[10^-30] - 1", "ArcTan[10^30] - Pi/2",
          "(1 + 10^-15)^(10^15) - E", "Sin[10^30]", "Sin[2^1000]", "Tan[Pi/2 - 10^-10] - 10^10",
          "Sinh[10^-20] - 10^-20", "Log[10^1000 + 1] - 1000 Log[10]", "1 - Cos[1/10^10]",
          "Sqrt[1 + 10^-40] - 1"]:
    t(f"kagree[{e}]")
t("N[Sin[Pi/7]^2 + Cos[Pi/7]^2 - 1, 20] == 0")
t("N[E^(I Pi/4) - (-1)^(1/4), 20] == 0")
t("N[Sqrt[2] Sqrt[3] - Sqrt[6], 30] == 0")
t("With[{a = 77617, b = 33096}, 1335/4 b^6 + a^2 (11 a^2 b^2 - b^6 - 121 b^4 - 2) + 11/2 b^8 + a/(2 b) === -54767/66192]")  # Rump
t("With[{a = 77617, b = 33096}, Abs[N[1335/4 b^6 + a^2 (11 a^2 b^2 - b^6 - 121 b^4 - 2) + 11/2 b^8 + a/(2 b), 20] + 54767/66192] < 10^-18]")
t("Module[{x = {4, 17/4}}, Do[x = {x[[2]], 108 - (815 - 1500/x[[1]])/x[[2]]}, {19}]; x[[2]] === (5^21 + 3^21)/(5^20 + 3^20)]")  # Muller, exact
t("Module[{x = {4., 4.25}}, Do[x = {x[[2]], 108. - (815. - 1500./x[[1]])/x[[2]]}, {40}]; Abs[x[[2]] - 100.] < 10^-6]")  # IEEE drifts to 100
t("With[{a = 100000, b = 9999999979/100000, c = 29/100000, s = (100000 + 9999999979/100000 + 29/100000)/2}, s (s - a) (s - b) (s - c) === (a + (b + c)) (c - (a - b)) (c + (a - b)) (a + (b - c))/16]")  # needle triangle
t("With[{a = 100000, b = 9999999979/100000, c = 29/100000}, kagree[Sqrt[(a + (b + c)) (c - (a - b)) (c + (a - b)) (a + (b - c))]/4]]")
t("Abs[N[(1 + 10^-15)^(10^15), 30] - N[E (1 - 1/(2 10^15)), 30]] < 10^-29")
t("(Expand[(x - 1)^7] /. x -> 10001/10000) === 10^-28")
t("With[{r = Min[x /. NSolve[x^2 - 10^8 x + 1 == 0, x]]}, Abs[r/N[(10^8 - Sqrt[10^16 - 4])/2, 30] - 1] < 10^-12]")
t("With[{r = Min[x /. NSolve[x^2 - 10.^8 x + 1. == 0, x]]}, Abs[r/N[(10^8 - Sqrt[10^16 - 4])/2, 30] - 1] < 10^-12]")
t("Abs[N[Exp[Pi Sqrt[163]] - 640320^3 - 744, 30] + 75/10^14] < 10^-14")

# ---------------------------------------------------------------- D. decimal -> binary
sec("D. correctly rounded decimal input")
DEC = [("0.1", "1/10"), ("0.2", "2/10"), ("0.3", "3/10"), ("1.7976931348623157*^308", "17976931348623157 10^292"),
       ("2.2250738585072014*^-308", "22250738585072014/10^324"), ("2.2250738585072011*^-308", "22250738585072011/10^324"),
       ("4.9406564584124654*^-324", "49406564584124654/10^340"), ("9007199254740993.", "9007199254740993"),
       ("1.0000000000000002", "10000000000000002/10^16"), ("8.98846567431158*^307", "898846567431158 10^293"),
       ("1.*^23", "10^23"), ("0.30000000000000004", "30000000000000004/10^17"), ("3.0000000000000004", "30000000000000004/10^16"),
       ("123456789.12345678", "12345678912345678/10^8"), ("7.0e-10", "7/10^10"), ("6.02214076*^23", "602214076 10^15"),
       ("1.602176634*^-19", "1602176634/10^28"), ("9.999999999999999*^22", "9999999999999999 10^7")]
for lit, q in DEC:
    t(f"kcr[{lit}, {q}]")
t("SetPrecision[9007199254740995., Infinity] === 2^53 + 4")         # tie to even
t("SetPrecision[9007199254740993., Infinity] === 2^53")

# ---------------------------------------------------------------- R. binary -> decimal round trip
sec("R. InputForm round trips")
for x in ["1./3.", "2./3.", "N[Pi]", "0.1 + 0.2", "1.*^23", "2.^-1022", "$MaxMachineNumber", "1.*^-320",
          "Exp[1000.]", "N[Pi, 30]", "N[1/3, 100]", "1.0000000000000002", "5.*^-324", "-0.5", "123.456", "2.^1023 3."]:
    t(f"krt[{x}]")
t("SeedRandom[1977]; And @@ Table[krt[RandomReal[] 10.^RandomInteger[{-300, 300}]], {300}]")

# ---------------------------------------------------------------- X. complex: branch cuts, robust division
sec("X. complex arithmetic: branch cuts and robustness")
t("Abs[Sqrt[-1. + 1.*^-300 I] - I] < 10^-15")                        # approached from above the cut
t("Abs[Sqrt[-1. - 1.*^-300 I] + I] < 10^-15")                        # from below
t("Abs[Log[-1. + 1.*^-300 I] - Pi I] < 10^-14")
t("Abs[Log[-1. - 1.*^-300 I] + Pi I] < 10^-14")
t("Abs[Arg[-1. - 1.*^-300 I] + Pi] < 10^-14 && Abs[Arg[-1. + 1.*^-300 I] - Pi] < 10^-14")
for z in ["2. + 3. I", "-4. + 0.5 I", "-4. - 0.5 I", "1.*^-200 + 1.*^-200 I", "1.*^200 - 3.*^200 I", "-1. + 1.*^-300 I"]:
    # |Log[z]| ~ 460 at |z| ~ 10^+-200 amplifies its rounding: 10^-12, not 10^-14
    # complex Sqrt as exp(log(z)/2) loses |log z| ulps at extreme |z| (K2)
    t(f"With[{{z = {z}}}, Abs[Exp[Log[z]] - z] < 10^-12 Abs[z] && Abs[Sqrt[z]^2 - z] < 10^-14 Abs[z] && Re[Sqrt[z]] >= 0]",
      xfail="K2" if z in ("1.*^-200 + 1.*^-200 I", "1.*^200 - 3.*^200 I") else None)
t("Abs[(1.*^300 + 1.*^300 I)/(1.*^300 + 1.*^300 I) - 1] < 10^-15")    # Smith-style robust division
t("Abs[(1.*^-300 + 1.*^-300 I)/(1.*^-300 + 1.*^-300 I) - 1] < 10^-15")
t("Abs[(1.*^300 + 1. I)/(1.*^300 + 1. I) - 1] < 10^-15")
t("Abs[(1. + 2. I)/(3. + 4. I) - (11/25 + 2/25 I)] < 10^-16")
t("Abs[Complex[1.*^308, 1.*^308]] == Sqrt[2] 10^308")
t("Abs[(1.*^-200 I)^2 + 10^-400] < 10^-415")
t("Abs[Sin[ArcSin[2.]] - 2] < 10^-15 && Abs[Cos[ArcCos[3.]] - 3] < 10^-14")
t("Abs[Sin[ArcSin[0.5 + 3. I]] - (0.5 + 3. I)] < 10^-14")
t("Abs[Tan[ArcTan[100.]] - 100.] < 10^-11 && Abs[ArcTan[Tan[1.5]] - 1.5] < 10^-15")   # sec^2 = 10^4 amplifies one ulp
t("Abs[1/(1/(1.*^300 + 1.*^300 I)) - (1.*^300 + 1.*^300 I)] < 10^286")

# ---------------------------------------------------------------- I. rounding to integers
sec("I. rounding to integers")
t("Round[0.5] === 0 && Round[1.5] === 2 && Round[2.5] === 2 && Round[-1.5] === -2 && Round[-2.5] === -2")
t("Floor[-0.5] === -1 && Ceiling[-0.5] === 0 && Floor[0.5] === 0 && Ceiling[1.5] === 2")
t("Round[2.^52 + 0.5] === 2^52")                                    # 2^52 + 1/2 is itself a tie, rounded to even
t("Floor[2.^62 + 1024.] === 2^62 + 1024")
t("IntegerPart[-2.^60 - 0.5] === -2^60")
t("Floor[N[10^30, 50] + 1/3] === 10^30 && Ceiling[N[10^30, 50] + 1/3] === 10^30 + 1")
t("Round[N[10^20, 40] + 1/2] === 10^20 && Round[N[10^20, 40] + 3/2] === 10^20 + 2")
t("Round[1.*^16 + 1.] === 10^16 && Round[1.*^16 + 2.] === 10^16 + 2")
t("Floor[1.9999999999999998] === 1 && Ceiling[1.0000000000000002] === 2")
t("Floor[-1.*^-300] === -1 && Ceiling[1.*^-300] === 1")
t("Round[0.49999999999999994] === 0")                               # classic round(x) = floor(x + 0.5) bug
t("Round[-0.49999999999999994] === 0")
t("IntegerPart[9007199254740993.] === 2^53")

# ---------------------------------------------------------------- M. arbitrary precision at fixed width
sec("M. MPFR: values and operations to the stated precision")
for q in ["1/3", "2/7", "10^100/3", "1/(3 10^100)", "22/7", "355/113"]:
    for p in ["20", "50", "200"]:
        t(f"Abs[SetPrecision[N[{q}, {p}], Infinity] - {q}] <= 10^-{p} Abs[{q}]")
for op in ["+", "-", "*", "/"]:
    t(f"With[{{x = N[1/3, 40], y = N[1/7, 40]}}, Abs[SetPrecision[x {op} y, Infinity] - (SetPrecision[x, Infinity] {op} SetPrecision[y, Infinity])] <= 10^-40 Abs[SetPrecision[x {op} y, Infinity]]]")
t("Abs[SetPrecision[Sqrt[N[2, 50]], Infinity]^2 - 2] < 10^-49")
t("Abs[SetPrecision[N[Pi, 100], Infinity] - N[Pi, 120]] < 10^-99")

# ---------------------------------------------------------------- S. statistics and summation
sec("S. numerically stable statistics and summation")
t("SetPrecision[Fold[Plus, 0., ConstantArray[0.1, 10]], Infinity] === 1 - 2^-53")
t("Abs[SetPrecision[Total[ConstantArray[0.1, 10]], Infinity] - 1] <= 2^-52")
t("keq[Variance[1.*^8 + {1., 2., 3.}], 1.]")                     # two-pass, not sum of squares
t("keq[Variance[1.*^9 + {4., 7., 13., 16.}], 30.]")
t("keq[Variance[1.*^15 + {4., 7., 13., 16.}], 30.]")
t("keq[Mean[1.*^15 + {4., 7., 13., 16.}], 1.*^15 + 10.]")
t("Abs[Total[1./Range[10^6]]/N[HarmonicNumber[10^6], 30] - 1] < 10^-13")
t("Mean[{1.*^308, 1.*^308, 1.*^308}] == 10^308")
t("N[Total[{1, 10^100, 1, -10^100}]] === 2.")

# ------------------------------------------------------------------------- emit
def m_escape(s): return s.replace("\\", "\\\\").replace('"', '\\"')
c_escape = m_escape

def emit_m(path):
    with open(path, "w") as f:
        for p in PRELUDE: f.write(p + ";\n")
        f.write('kN = 0; kFail = 0;\n')
        f.write('kt[s_String] := (kN++; If[!TrueQ[Quiet[ToExpression[s]]], kFail++; Print["FAIL: ", s]]);\n')
        for c in CASES:
            if c[0] == "#": f.write(f'(* {c[1]} *)\n')
            elif c[0] == "T": f.write(f'kt["{m_escape(c[1])}"];\n')
            else: f.write(f'(* xfail {c[2]} *) kt["{m_escape(c[1])}"];\n')
        f.write('Print["ran ", kN, " failed ", kFail];\n')

C_HEAD = r'''/* Mathilda -- Kahan-style numerics stress corpus (GENERATED by
 * tools/numeric_kahan_cases.py; add cases there and regenerate:
 *   python3 tools/numeric_kahan_cases.py /tmp/k.m tests/test_numeric_kahan.c).
 *
 * IEEE 754 machine properties in the manner of Paranoia, correct rounding of
 * + - * / sqrt against the exact result, faithful elementary functions at hard
 * arguments, decimal<->binary conversion and round trips, the cancellation
 * classics (Rump, Muller, the needle triangle, Ramanujan's near-integers)
 * through N, branch-cut continuity and robust complex arithmetic, and stable
 * statistics. Every case is a boolean whose truth follows from exact
 * mathematics; exactness is asserted through exact binary values
 * (SetPrecision[x, Infinity]), since Equal / SameQ are tolerant by design.
 * XFAIL cases are known open bugs (docs/design/numeric_tower_review.md); one
 * that starts to pass prints XPASS.
 */
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef struct { const char* expr; const char* xfail; } KCase;

static const char* PRELUDE[] = {
'''
C_MID = r'''};

static const KCase CASES[] = {
'''
C_TAIL = r'''};

static Expr* eval_src(const char* src) {
    Expr* e = parse_expression(src);
    if (!e) return NULL;
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

static int run_case(const char* src) {
    size_t n = strlen(src) + 16;
    char* wrapped = malloc(n);
    snprintf(wrapped, n, "TrueQ[Quiet[%s]]", src);
    Expr* r = eval_src(wrapped);
    free(wrapped);
    if (!r) return -1;
    int ok = r->type == EXPR_SYMBOL && strcmp(r->data.symbol.name, "True") == 0;
    expr_free(r);
    return ok;
}

int main(void) {
    symtab_init();
    core_init();
    for (size_t i = 0; i < sizeof PRELUDE / sizeof PRELUDE[0]; i++) {
        Expr* r = eval_src(PRELUDE[i]);
        if (r) expr_free(r);
    }
    size_t n = sizeof CASES / sizeof CASES[0];
    int fail = 0, xfail = 0, xpass = 0;
    for (size_t i = 0; i < n; i++) {
        int ok = run_case(CASES[i].expr);
        if (CASES[i].xfail) {
            if (ok == 1) { xpass++; printf("XPASS [%s]: %s\n", CASES[i].xfail, CASES[i].expr); }
            else xfail++;
        } else if (ok != 1) {
            fail++;
            printf("FAIL%s: %s\n", ok < 0 ? " (parse)" : "", CASES[i].expr);
        }
    }
    printf("Kahan numerics stress: %zu cases, %d failed, %d expected failures, %d unexpected passes\n",
           n, fail, xfail, xpass);
    return fail ? 1 : 0;
}
'''

def emit_c(path):
    with open(path, "w") as f:
        f.write(C_HEAD)
        for p in PRELUDE: f.write(f'    "{c_escape(p)}",\n')
        f.write(C_MID)
        for c in CASES:
            if c[0] == "#": f.write(f"    /* ---- {c[1]} */\n")
            elif c[0] == "T": f.write(f'    {{ "{c_escape(c[1])}", NULL }},\n')
            else: f.write(f'    {{ "{c_escape(c[1])}", "{c[2]}" }},\n')
        f.write(C_TAIL)

if __name__ == "__main__":
    emit_m(sys.argv[1])
    if len(sys.argv) > 2: emit_c(sys.argv[2])
    n = sum(1 for c in CASES if c[0] != "#")
    print(n, "cases;", sum(1 for c in CASES if c[0] == "X"), "xfail")

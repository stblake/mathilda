#!/usr/bin/env python3
"""Numeric-tower stress corpus (generates tests/test_numeric_tower_stress.c).

Usage: python3 tools/numeric_tower_stress_cases.py OUT.m [tests/test_numeric_tower_stress.c]

Every case is a boolean expression; xfail=<id> marks a known open bug
(docs/design/numeric_tower_review.md, section S).

Original summary: every case is a boolean Mathilda expression whose
truth follows from exact mathematics (identities, exact binary values, known
digits, cross-surface agreement) -- no expected value is taken from any other
system. Emits a .m runner (fast triage against ./Mathilda) and the ctest C file.
"""
import sys

P = "Developer`ToPackedArray"
CASES = []
def sec(name): CASES.append(("#", name))
def t(expr, xfail=None): CASES.append(("T", expr) if xfail is None else ("X", expr, xfail))

# ---------------------------------------------------------------- A. int64 / BigInt boundary
sec("A. int64 / BigInt boundary")
t("(2^63 - 1) + 1 === 2^63")
t("(-2^63) - 1 + 1 === -2^63")
t("(2^63) (2^63) === 2^126")
t("(2^63 - 1) (2^63 - 1) === 2^126 - 2^64 + 1")
t("(-2^63) (-1) === 2^63")
t("(-2^63)/(-1) === 2^63")
t("Quotient[-2^63, -1] === 2^63")
t("Mod[-2^63, -1] === 0")
t("Quotient[2^63, 2] === 2^62")
t("Mod[2^64 + 5, 2^32] === 5")
t("Abs[-2^63] === 2^63")
t("Sign[-2^63] === -1")
t("2^63 - 2^63 === 0")
t("Sqrt[2^126] === 2^63")
t("Sqrt[-2^63] === 2147483648 I Sqrt[2]")
t("(-2^63)^(1/3) === 2097152 (-1)^(1/3)")
t("(2^63)^(1/3) === 2097152")
t("(-2^63)^2 === 2^126")
t("(-2^63)^-1 === -1/2^63")
t("Numerator[-1/2^63] === -1 && Denominator[-1/2^63] === 2^63")
t("GCD[2^63, 2^62] === 2^62")
t("GCD[-2^63, -2^63] === 2^63")
t("LCM[2^62, 3] === 3 2^62")
t("LCM[-2^63, -2^63] === 2^63")
t("GCD[2^63 - 1, 2^62 - 1] === 1")                       # gcd(2^a-1,2^b-1) = 2^gcd(a,b)-1
t("GCD[2^64 - 1, 2^32 - 1] === 2^32 - 1")
t("With[{g = ExtendedGCD[2^64 + 1, 2^32 + 3]}, g[[1]] == g[[2, 1]] (2^64 + 1) + g[[2, 2]] (2^32 + 3)]")
t("PowerMod[3, 10^20, 10^9 + 7] === PowerMod[PowerMod[3, 10^10, 10^9 + 7], 10^10, 10^9 + 7]")
t("Mod[2 PowerMod[2, -1, 2^63 + 1], 2^63 + 1] === 1")
t("FromDigits[IntegerDigits[2^63 + 12345, 7], 7] === 2^63 + 12345")
t("FromDigits[IntegerDigits[-2^63]] === 2^63")
t("IntegerLength[2^63] === 19 && IntegerLength[10^19] === 20")
t("BitAnd[2^63 + 5, 2^63 - 1] === 5")
t("BitOr[2^63, 1] === 2^63 + 1")
t("BitXor[2^64 - 1, 2^63] === 2^63 - 1")
t("BitNot[-2^63] === 2^63 - 1")
t("BitShiftLeft[1, 64] === 2^64")
t("BitShiftRight[2^64 + 2^63, 63] === 3")
t("BitLength[2^63] === 64 && BitLength[2^63 - 1] === 63")
t("Floor[(2^63 + 1)/2] === 2^62 && Ceiling[(2^63 + 1)/2] === 2^62 + 1")
t("Round[(2^63 + 1)/2] === 2^62")                         # 2^62 + 1/2: half to even
t("Round[5/2] === 2 && Round[7/2] === 4 && Round[-5/2] === -2")
t("IntegerPart[-7/2] === -3 && FractionalPart[-7/2] === -1/2")
t("Mod[-7, 3] === 2 && Mod[7, -3] === -2")
t("Quotient[-7, 3] === -3 && Quotient[7, -3] === -3")
t("With[{a = -10^30 - 7, b = 2^63 + 9}, a == b Quotient[a, b] + Mod[a, b] && 0 <= Mod[a, b] < b]")
t("With[{a = 10^40 + 3, b = -(2^70 + 1)}, a == b Quotient[a, b] + Mod[a, b] && b < Mod[a, b] <= 0]")
t("With[{a = -2^63, b = 7}, a == b Quotient[a, b] + Mod[a, b] && 0 <= Mod[a, b] < b]")
t("Factorial[25] === 15511210043330985984000000")
t("Binomial[2^63, 2] === 2^62 (2^63 - 1)")
t("Binomial[100, 50] === 100891344545564193334812497256")
t("Fibonacci[100] === 354224848179261915075")
t("(2^63 + 1)/(2^63 + 1) === 1")
t("3^40/3^38 === 9")
t("Divide[2^100, 2^99] === 2")
t("(2^64 - 1)/(2^32 + 1) === 2^32 - 1")
t("Rational[2^70, -2^68] === -4")
t("Rational[-6, 2^70] === -3/2^69")
t("Total[Range[2^62, 2^62 + 9]] === 10 2^62 + 45")
t("Total[ConstantArray[2^62, 8]] === 2^65")
t("Accumulate[{2^62, 2^62, 2^62}] === {2^62, 2^63, 3 2^62}")
t("Times @@ ConstantArray[2^31, 3] === 2^93")
t("{2^40, 3} . {2^40, 5} === 2^80 + 15")
t("Plus @@ Table[(-1)^k 2^60, {k, 1, 20}] === 0")
t("Max[2^63, 2^63 - 1] === 2^63 && Min[-2^63, -2^63 + 1] === -2^63")
t("Sort[{2^63, -2^63, 2^63 - 1, 0}] === {-2^63, 0, 2^63 - 1, 2^63}")
t("2^63 > 2^63 - 1 && -2^63 < -2^63 + 1")
t("Hash[2 2^62] === Hash[2^63]")
t("Length[Union[{2^63, 2 2^62, 4 2^61}]] === 1")
t("<|2^63 -> a|>[2 2^62] === a")
t("Count[{1, 000000000000000000001}, 1] === 2")
t("(-1)^(2^64) === 1 && (-1)^(2^64 + 1) === -1")
t("I^(2^64 + 1) === I", xfail="S1")
t("(-1)^(-2^64) === 1")
t("2^(-64) === 1/2^64 && (1/2)^64 === 1/2^64")
t("(2^64 + 2^32)/(2^32) === 2^32 + 1")
t("Mod[(2^64 + 1)^2, 2^64] === 1")
t("Quotient[2^127, 2^63 + 1] === Floor[2^127/(2^63 + 1)]")

# ---------------------------------------------------------------- B. exact rationals
sec("B. exact rationals")
t("1/3 + 1/6 === 1/2")
t("(1/3)^-2 === 9")
t("(-2/3)^3 === -8/27 && (-2/3)^-3 === -27/8")
t("Denominator[(-3)^-5] === 243 && Numerator[(-3)^-5] === -1")
t("2^64/3 - (2^64 - 1)/3 === 1/3")
t("1/2^63 + 1/2^63 === 1/2^62")
t("(2^70 + 2)/(2^70 + 4) === (2^69 + 1)/(2^69 + 2)")
t("Floor[-(2^70 + 1)/2] === -2^69 - 1")
t("Ceiling[(2^70 + 1)/2] === 2^69 + 1")
t("Round[(2^70 + 1)/2] === 2^69")
t("Round[-(2^70 + 3)/2] === -(2^69 + 2)")
t("Mod[1/3, 1/7] === 1/21")
t("Quotient[10/3, 1/2] === 6")
t("Mod[-1/3, 1/2] === 1/6")
t("Abs[-2^70/3] === 2^70/3")
t("Max[1/3, 0.333] === 1/3")
t("1/3 > 0.3333333")
t("Sort[{1/3, 0.3333333333333333, 1/4}] === {1/4, 0.3333333333333333, 1/3}")
t("Rationalize[0.125] === 1/8")
t("Rationalize[0.1, 0] === 1/10")
t("SetPrecision[0.5, Infinity] === 1/2")
t("SetPrecision[0.1, Infinity] === 3602879701896397/36028797018963968")
t("SetPrecision[N[1/7], Infinity] === 2573485501354569/18014398509481984")   # ulp 2^-55: round(2^55/7) = 5146971002709138
t("1/(1/3 - 1/3 + 1/2) === 2")
t("Together[1/(2^63) - 1/(2^63 + 1)] === 1/(2^63 (2^63 + 1))")
t("Numerator[(2^65)/(2^64 + 2^63)] === 4 && Denominator[(2^65)/(2^64 + 2^63)] === 3")
t("(10^30 + 1/3) - 10^30 === 1/3")
t("Sum[1/2^k, {k, 1, 100}] === 1 - 1/2^100")
t("Product[(k + 1)/k, {k, 1, 1000}] === 1001")

# ---------------------------------------------------------------- C. IEEE machine edge cases
sec("C. IEEE machine edge cases")
t("SetPrecision[0.1 + 0.2, Infinity] - SetPrecision[0.3, Infinity] === 2^-54")
t("Head[1.0 - 1.0] === Real && 1.0 - 1.0 == 0")
t("-0. == 0. && Sign[-0.] === 0")
t("1/(1. + 1.*^-16) == 1.")
t("SetPrecision[2.^53 + 1., Infinity] === 2^53")
t("SetPrecision[2.^53 + 2., Infinity] === 2^53 + 2")
t("SetPrecision[1. + $MachineEpsilon, Infinity] === 1 + 2^-52")
t("SetPrecision[$MachineEpsilon, Infinity] === 2^-52")
t("SetPrecision[$MaxMachineNumber, Infinity] === 2^1024 - 2^971")
t("SetPrecision[$MinMachineNumber, Infinity] === 2^-1022")
t("MachineNumberQ[$MaxMachineNumber] && MachineNumberQ[$MinMachineNumber]")
t("Round[2.5] === 2 && Round[3.5] === 4 && Round[-2.5] === -2 && Round[0.5] === 0")
t("Floor[-0.5] === -1 && Ceiling[-0.5] === 0")
t("Round[1.*^20] === 100000000000000000000")
t("Floor[-1.*^19] === -10000000000000000000")
t("IntegerPart[3. 2.^62] === 3 2^62")
t("Quotient[1.*^300, 1.] === Floor[1.*^300]")
t("Quotient[-7.5, 2.] === -4 && Mod[-7.5, 2.] == 0.5")
t("Mod[5.5, -2] == -0.5")
t("Max[1., 2^70] === 2^70 && Max[1.*^30, 2^70] === 1.*^30")
t("Sort[{2^53 + 1, 2.^53, 2^53}] === {2^53, 2.^53, 2^53 + 1}")
t("Order[2.^53, 2^53 + 1] === 1")
t("Max[0.1, 1/10] == 0.1 && 0.1 == 1/10 && !(0.1 === 1/10)")
t("SetPrecision[N[2/3], Infinity] === 6004799503160661/9007199254740992")
t("SetPrecision[2.^-1022 2^-52, Infinity] === 2^-1074")
t("SetPrecision[N[10^22], Infinity] === 10^22")                # 10^22 is exactly a double
t("SetPrecision[N[10^23], Infinity] === 99999999999999991611392")
t("Floor[1.9999999999999998] === 1 && Ceiling[1.0000000000000002] === 2")

# ---------------------------------------------------------------- D. overflow / underflow policy
sec("D. machine/MPFR boundary policy")
t("$MaxMachineNumber 2 == 2^1025 - 2^972")
t("$MaxMachineNumber + $MaxMachineNumber == 2^1025 - 2^972")
t("!MachineNumberQ[$MaxMachineNumber 2]")
t("$MinMachineNumber/2 == 2^-1023 && !MachineNumberQ[$MinMachineNumber/2]")
t("$MinMachineNumber $MinMachineNumber == 2^-2044")
t("1.*^-300 1.*^-300 == 10^-600")
t("1.*^308 1.*^308/1.*^308 == 10^308")
t("Exp[1000.]/Exp[999.] == E")
t("Log[Exp[800.]] == 800")
t("Sqrt[1.*^-400] == 10^-200")
t("(1.*^200)^2 == 10^400 && (1.*^-200)^2 == 10^-400")
t("Abs[-1.*^400] == 10^400")
t("Abs[Floor[1.*^400] - 10^400] < 10^385 && IntegerQ[Floor[1.*^400]]")   # 53-bit 1.*^400 lies just below 10^400
t("1.5 < 10^400 && -1.5 > -10^400 && 1.*^-300 > 10^-400")
t("Sort[{N[10^500], 1.5, N[10^400], -N[10^400]}] == {-10^400, 1.5, 10^400, 10^500}")
t("Exp[-1000.] > 0 && Exp[-1000.] == Exp[-1000]")
t("Exp[-1000.] Exp[1000.] == 1")
t("Gamma[200.] == 199!")
t("Factorial[200.] == 200!")
t("Binomial[1000., 500.] == Binomial[1000, 500]")
t("Sinh[1000.] == N[Sinh[1000], 30]")
t("Tanh[1000.] === 1.")
t("ArcTan[1.*^400] == Pi/2")
t("Erfc[30.] == N[Erfc[30], 30]")
t("Exp[1000. + 0.5 I] == N[Exp[1000 + I/2], 30]", xfail="S3")
t("1.*^308 10 > 1.*^308")
t("N[10^-320] == 10^-320 && !MachineNumberQ[N[10^-320]]")
t("N[3/2 2^-1074] == 3/2 2^-1074")
t("2.^-1074 == 2^-1074 && !MachineNumberQ[2.^-1074]")
t("2.^1024 == 2^1024")
t("10.^-308/10^10 == 10^-318")
t("Exp[1000.] 0. === 0.")
t("1.*^400 == 10^400 && 1.*^-400 == 10^-400")
t("1.*^-320 == 10^-320")
t("FreeQ[Hold[1.*^400], _DirectedInfinity]")
t("1.*^308 + 1.*^308 - 1.*^308 == 10^308")

# ---------------------------------------------------------------- E. complex machine arithmetic
sec("E. complex machine arithmetic")
t("Sqrt[-1.] === Complex[0., 1.] && Sqrt[-4.] === Complex[0., 2.]")
t("(-8.)^(1/3) == 1. + Sqrt[3.] I")
t("(-1.)^0.5 === Complex[0., 1.]")
t("(-1.)^1.5 === Complex[0., -1.]")
t("(-2.)^2.5 == 4 Sqrt[2.] I")
t("Log[-1.] == Pi I")
t("Exp[I Pi 1.] == -1", xfail="S2")
t("Complex[1.0, 2] === Complex[1.0, 2.0]")
t("Head[Complex[2^70, 0.]] === Real")
t("Abs[3. + 4. I] === 5.")
t("Abs[Complex[3.*^200, 4.*^200]] == 5*^200")
t("Abs[Complex[3.*^-200, 4.*^-200]] == 5*^-200")
t("Re[(1.*^200 + 1.*^200 I)^2] == 0 && Im[(1.*^200 + 1.*^200 I)^2] == 2 10^400", xfail="S3")
t("(1. + 2. I)(3. - 4. I) === 11. + 2. I")
t("(2. I)^2 === -4. + 0. I || (2. I)^2 === -4.")
t("Precision[Complex[1.5`30, 2.]] === MachinePrecision")
t("Precision[Abs[Complex[N[Pi, 30], 1.5]]] === MachinePrecision")
t("Arg[-1.] == Pi && Arg[-1. I] == -Pi/2")
t("Conjugate[1.*^400 + 2. I] == 10^400 - 2 I")

# ---------------------------------------------------------------- F. cross-surface: packed / NDArray / Compile
sec("F. cross-surface agreement")
for f in ["Sin", "Cos", "Exp", "Log", "Sqrt", "Abs", "Tanh", "ArcTan", "Floor", "Round"]:
    t(f"With[{{v = {P}[{{0.5, 2., 100., 1.*^-5}}]}}, {f}[v] == Map[{f}, Developer`FromPackedArray[v]]]")
for f in ["Sin", "Exp", "Log", "Sqrt", "ArcTan"]:
    t(f"With[{{c = Compile[{{{{x, _Real}}}}, {f}[x]]}}, Map[c, {{0.5, 2., 100., 1.*^-5}}] == Map[{f}, {{0.5, 2., 100., 1.*^-5}}]]")
t(f"{P}[{{1.*^308, 1.}}] 10 == {{10^309, 10}}")
t(f"Total[{P}[{{1.*^308, 1.*^308, -1.*^308}}]] == 10^308")
t(f"Mean[{P}[{{1.*^308, 1.*^308}}]] == 10^308")
t("Norm[{3.*^200, 4.*^200}] == 5*^200")
t("Norm[{3.*^-200, 4.*^-200}] == 5*^-200")
t("{1.*^200, 1.*^200} . {1.*^200, 1.*^200} == 2 10^400")
t("Variance[{1.*^200, 3.*^200}] == 2 10^400")
t("StandardDeviation[{1.*^200, 3.*^200}] == Sqrt[2] 10^200")
t(f"Accumulate[{P}[{{1.*^-200, 1.*^-200}}]] == {{10^-200, 2 10^-200}}")
t(f"{P}[{{1.*^-200}}]^2 == {{10^-400}}", xfail="S4")
t(f"{P}[{{1.*^200, 2.}}]^2 == {{10^400, 4}}")
t(f"Exp[{P}[{{-800., 0.}}]] == {{Exp[-800], 1}}", xfail="S5")
t(f"Exp[{P}[{{800., 0.}}]] == {{Exp[800], 1}}")
t(f"{P}[{{2.}}]^-1074 == {{2^-1074}}")
t("Compile[{{x, _Real}}, x^2][1.*^-200] == 10^-400")
t("Compile[{{x, _Real}}, Exp[x]][-800.] == Exp[-800]")
t("Compile[{{x, _Real}}, Exp[x]][800.] == Exp[800]")
t("Compile[{{x, _Real}}, x x][1.*^200] == 10^400")
t("Compile[{{v, _Real, 1}}, v v][{1.*^-200, 2.}] == {10^-400, 4}")
t("Compile[{{v, _Real, 1}}, v + v][{1.*^308, 2.}] == {2 10^308, 4}")
t("NDArray[{1.*^308, 2.}] NDArray[{10., 1.}] == {10^309, 2}")
t("Sqrt[NDArray[{-1., 4.}]] == {I, 2}")
t(f"Sqrt[{P}[{{-1., 4.}}]] == {{I, 2}}")
t(f"With[{{v = {P}[Range[-3., 3.]]}}, v^3 == Map[#^3 &, Developer`FromPackedArray[v]]]")
t(f"Developer`PackedArrayQ[{P}[{{1., 2.}}] 10.]")

# ---------------------------------------------------------------- G. iterators and reductions
sec("G. iterators and reductions")
t("Table[2.^k, {k, -1076, -1074}] == {2^-1076, 2^-1075, 2^-1074}")
t("Table[10.^k, {k, 307, 309}] == {10^307, 10^308, 10^309}")
t("Range[1.*^-300, 3.*^-300, 1.*^-300] == {10^-300, 2 10^-300, 3 10^-300}")
t("Length[Range[0, 1, 1/10]] === 11")
t("Length[Range[0., 1., 1./3]] === 4")
t("Length[Range[0., 1., 0.1]] === 11")
t("Length[Table[x, {x, 0., 10., 0.01}]] === 1001")
t("Sum[1.*^-300, {i, 1, 10}] == 10^-299")
t("Sum[x, {x, 0., 1., 0.25}] == 5/2")
t("Product[1.*^100, {i, 1, 5}] == 10^500")
t("Product[1.*^-100, {i, 1, 5}] == 10^-500", xfail="S6")
t("Length[Range[1.*^-200, 5.*^-200, 1.*^-200]] === 5")
t("Length[Range[1.*^20, 2.*^20, 1.*^19]] === 11")
t("Total[Table[1/k, {k, 1, 100}]] === HarmonicNumber[100]")
t("Total[Range[10^6]] === 500000500000")
t("Total[Range[10^6]^2] === 333333833333500000")
t("Total[Range[2^31, 2^31 + 99]^2] === Sum[k^2, {k, 2^31, 2^31 + 99}]")
t("Total[Table[2^62, {10}]] === 10 2^62")
t("Max[Range[-2^62, -2^62 + 5]] === -2^62 + 5")
t("Mean[{2^62, 2^62 + 2}] === 2^62 + 1")
t("Mean[{1/3, 2/3}] === 1/2")
t("Variance[{2^62, 2^62 + 2}] === 2")
t("Median[{2^63, 1, 2^64}] === 2^63")
t("Fold[Plus, 0, Table[1/2^k, {k, 1, 60}]] === 1 - 1/2^60")

# ---------------------------------------------------------------- H. MPFR precision / accuracy
sec("H. arbitrary precision")
t("Precision[N[Pi, 50]] == 50")
t("Precision[N[Pi, 50] + N[E, 30]] == 30")
t("Precision[N[Pi, 50] 2] == 50")
t("Precision[N[Pi, 50] + 1.] === MachinePrecision")
t("Precision[N[Pi, 30]^1.5] === MachinePrecision")
t("Precision[N[Pi, 30]^2] == 30")
t("Precision[1.5`20 + 1/3] == 20")
t("Accuracy[N[10^-400, 20]] == 420")
t("Accuracy[N[10^500, 20]] == -480")
t("Precision[SetAccuracy[Pi 10^50, 10]] > 60")
t("Accuracy[SetAccuracy[Pi, 10]] >= 10")
t("Accuracy[SetAccuracy[Pi/10^50, 10]] == 10")
t("N[Pi, 50] - Pi == 0")
t("Abs[N[Pi, 60] - 314159265358979323846264338327950288419716939937510582097494/10^59] < 10^-58")
t("Abs[N[E, 40] - 2718281828459045235360287471352662497757/10^39] < 10^-38")
t("Abs[N[Sqrt[2], 50]^2 - 2] < 10^-48")
t("Abs[N[Sin[10^30], 30]^2 + N[Cos[10^30], 30]^2 - 1] < 10^-25")
t("N[Exp[1000], 20]/N[Exp[999], 20] == E")
t("Abs[N[Log[10^400], 30] - 400 N[Log[10], 30]] < 10^-25")
t("N[1/3, 30] + N[2/3, 30] == 1")
t("With[{x = N[1/3, 30]}, SetPrecision[x, Infinity] =!= 1/3 && N[SetPrecision[x, Infinity], 30] == x]")
t("Precision[N[Pi, 10^3]] == 1000")
t("Precision[2.`20 3.`30] == 20")
t("N[2^(1/2), 40] == N[Sqrt[2], 40]")
t("Abs[N[Pi, 30] - N[Pi, 40]] < 10^-29")
t("N[10^400 + 1, 20] == 10^400")
t("Abs[N[(10^400 + 1)/10^400, 500] - 1 - 10^-400] < 10^-490")   # fixed precision: ~100 good digits survive the cancellation
t("Precision[N[1/7, 50] - 1/7] == 50 || N[1/7, 50] - 1/7 == 0")
t("Round[N[10^30, 40] + 1/2] - 10^30 === 0 || Round[N[10^30, 40] + 1/2] - 10^30 === 1")
t("Floor[N[10^30, 40] + 1/2] - 10^30 === 0")
t("Quotient[N[10^30, 40] + 1/2, 1] - 10^30 === 0")
t("IntegerPart[N[-7/2, 30]] === -3")
t("Mod[N[10, 30], 3] == 1")
t("SetPrecision[2^100, 10] == 2^100")
t("Precision[SetPrecision[1.5, 40]] == 40")

# ---------------------------------------------------------------- I. round trips and bounds
sec("I. round trips, $MinPrecision / $MaxPrecision")
t("ToExpression[ToString[InputForm[N[Pi, 40]]]] === N[Pi, 40]")
t("ToExpression[ToString[InputForm[Exp[1000.]]]] == Exp[1000.]")
t("ToExpression[ToString[InputForm[1.*^-320]]] == 10^-320")
t("ToExpression[ToString[InputForm[0.1]]] === 0.1")
t("ToExpression[ToString[InputForm[$MaxMachineNumber]]] === $MaxMachineNumber")
t("ToExpression[ToString[InputForm[2.^-1022]]] === 2.^-1022")
t("ToExpression[ToString[InputForm[1/3 + 0.``20]]] == 1/3")
t("ToExpression[ToString[InputForm[-2^63 - 1/7]]] === -2^63 - 1/7")
t("Block[{$MaxPrecision = 20}, Precision[Sqrt[2.`30]]] == 20")
t("Quiet[Block[{$MinPrecision = 30}, Precision[N[Pi, 10]]]] == 30")
t("Quiet[Block[{$MinPrecision = 25, $MaxPrecision = 25}, Precision[N[E, 50] N[Pi, 10]]]] == 25")
t("Block[{$MinPrecision = 30}, Precision[1.5 2.]] === MachinePrecision")
t("Quiet[Block[{$MaxPrecision = 20}, N[Pi, 50]]]; Precision[N[Pi, 50]] == 50")
t("{$MinPrecision, $MaxPrecision} === {0, Infinity}")

# ---------------------------------------------------------------- J. parser literals
sec("J. parser literals")
t("-2.`30^2 < 0 && -2`10^2 < 0 && -2.5*^3^2 == -6.25*^6")
t("1.5*^3 === 1500.")
t("Hold[-2^2] === Hold[Times[-1, Power[2, 2]]]")
t("Precision[2.5`10 + 0] == 10")
t("000000000000000000001 === 1")
t("Head[123456789012345678901234567890] === Integer")
t("-9223372036854775808 === -2^63")
t("9223372036854775807 + 1 === 2^63")
t("-9223372036854775809 === -2^63 - 1")
t("Hold[1.5`+nan] === Hold[1.5 + nan]")
t("1.`30 == 1 && Precision[1.`30] == 30")
t("2.*^-3 === 0.002")
t("Accuracy[1.5``20] == 20 || Abs[Accuracy[1.5``20] - 20] < 1")
t("16^^FF === 255 && 2^^1010 === 10", xfail="S7")
t("16^^FFFFFFFFFFFFFFFF === 2^64 - 1", xfail="S7")

# ------------------------------------------------------------------------- emit
def m_escape(s): return s.replace("\\", "\\\\").replace('"', '\\"')
def c_escape(s): return s.replace("\\", "\\\\").replace('"', '\\"')

def emit_m(path):
    with open(path, "w") as f:
        f.write('stressN = 0; stressFail = 0;\n')
        f.write('st[s_String] := (stressN++; If[!TrueQ[Quiet[ToExpression[s]]], stressFail++; Print["FAIL: ", s]]);\n')
        for c in CASES:
            if c[0] == "#": f.write(f'(* {c[1]} *)\n')
            elif c[0] == "T": f.write(f'st["{m_escape(c[1])}"];\n')
            else: f.write(f'(* xfail {c[2]} *) st["{m_escape(c[1])}"];\n')
        f.write('Print["ran ", stressN, " failed ", stressFail];\n')

def count(): return sum(1 for c in CASES if c[0] != "#")

C_HEAD = r'''/* Mathilda -- numeric-tower STRESS corpus (generated; do not hand-edit the
 * case list). Generator: tools/numeric_tower_stress_cases.py -- add cases
 * there and regenerate:
 *   python3 tools/numeric_tower_stress_cases.py /tmp/stress.m tests/test_numeric_tower_stress.c
 * (the .m runner triages new cases against ./Mathilda -file in seconds).
 *
 * Every case is a boolean expression whose truth follows from exact
 * mathematics: identities ((a + b) - b === a, a == b Quotient[a, b] + Mod[a, b]),
 * exact binary values (SetPrecision[x, Infinity]), digits of known constants,
 * and agreement across the scalar / List / packed / NDArray / Compile
 * surfaces. No expected value is taken from another system.
 *
 * XFAIL cases are known open bugs, keyed by id in
 * docs/design/numeric_tower_review.md (section S). One that starts to pass
 * prints XPASS so its marker can be removed; any other failure fails the
 * suite. That makes the corpus a ratchet.
 */
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef struct { const char* expr; const char* xfail; } StressCase;

static const StressCase CASES[] = {
'''
C_TAIL = r'''};

static int run_case(const char* src) {
    /* Quiet so expected diagnostics (preclg, ovfl, ...) do not flood the log. */
    size_t n = strlen(src) + 16;
    char* wrapped = malloc(n);
    snprintf(wrapped, n, "TrueQ[Quiet[%s]]", src);
    Expr* e = parse_expression(wrapped);
    free(wrapped);
    if (!e) return -1;
    Expr* r = evaluate(e);
    int ok = r && r->type == EXPR_SYMBOL && strcmp(r->data.symbol.name, "True") == 0;
    expr_free(e);
    expr_free(r);
    return ok;
}

int main(void) {
    symtab_init();
    core_init();
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
    printf("numeric-tower stress: %zu cases, %d failed, %d expected failures, %d unexpected passes\n",
           n, fail, xfail, xpass);
    return fail ? 1 : 0;
}
'''

def emit_c(path):
    with open(path, "w") as f:
        f.write(C_HEAD)
        for c in CASES:
            if c[0] == "#": f.write(f"    /* ---- {c[1]} */\n")
            elif c[0] == "T": f.write(f'    {{ "{c_escape(c[1])}", NULL }},\n')
            else: f.write(f'    {{ "{c_escape(c[1])}", "{c[2]}" }},\n')
        f.write(C_TAIL)

if __name__ == "__main__":
    emit_m(sys.argv[1])
    if len(sys.argv) > 2: emit_c(sys.argv[2])
    print(count(), "cases;", sum(1 for c in CASES if c[0] == "X"), "xfail")

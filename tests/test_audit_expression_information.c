/* Builtin-audit regressions: expression information (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* ByteCount#1 (high): "ByteCount segfaults on a deeply nested expression
 * (unbounded C recursion)".  The walk is iterative, and its total is the sum
 * over nodes, so f-nesting adds a fixed amount per level. */
static void test_bytecount_deep_nesting(void) {
    audit_run("auditDeep = Nest[f, x, 200000]");
    audit_true("ByteCount#1", "IntegerQ[ByteCount[auditDeep]]");
    audit_true("ByteCount#1",
               "ByteCount[auditDeep] - ByteCount[Nest[f, x, 199999]] == "
               "ByteCount[f[x]] - ByteCount[x]");
    audit_chk("ByteCount#1", "LeafCount[auditDeep]", "200001");
    audit_run("Clear[auditDeep]");
    /* Shallow values: a list's count is its own node plus its parts'. */
    audit_true("ByteCount#1",
               "ByteCount[{a, \"bc\"}] > ByteCount[a] + ByteCount[\"bc\"] + ByteCount[List]");
    audit_true("ByteCount#1", "ByteCount[2^1000] > ByteCount[2]");
}

/* Positive#1 (high): "Positive returns True for numeric expressions that are
 * exactly zero (sign taken from rounding noise)".  The sign of an exact
 * closed form is certified (ball arithmetic, else the exact zero test); a
 * hidden zero is zero. */
static void test_positive_exact_zero(void) {
    audit_chk("Positive#1", "Positive[Sqrt[2] + Sqrt[3] - Sqrt[5 + 2 Sqrt[6]]]", "False");
    audit_chk("Positive#1", "Positive[Sin[Pi/7]^2 + Cos[Pi/7]^2 - 1]", "False");
    audit_chk("Positive#1", "Negative[Log[6] - Log[2] - Log[3]]", "False");
    audit_chk("Positive#1", "Negative[Sqrt[5 + 2 Sqrt[6]] - Sqrt[2] - Sqrt[3]]", "False");
    audit_chk("Positive#1", "NonNegative[Sqrt[2] + Sqrt[3] - Sqrt[5 + 2 Sqrt[6]]]", "True");
    audit_chk("Positive#1", "NonPositive[Sin[Pi/7]^2 + Cos[Pi/7]^2 - 1]", "True");
    /* Genuine signs, complex values and approximate numbers still decide. */
    audit_chk("Positive#1", "Positive[Pi - 3]", "True");
    audit_chk("Positive#1", "Negative[Pi - 22/7]", "True");
    audit_chk("Positive#1", "Positive[Pi - 80143857/25510582]", "True");
    audit_chk("Positive#1", "Positive[Zeta[3] - 1]", "True");
    audit_chk("Positive#1", "Positive[E - Pi + 1/2]", "True");
    audit_chk("Positive#1", "Positive[Pi + I]", "False");
    audit_chk("Positive#1", "Positive[Sqrt[2] - 1.4142]", "True");
    audit_chk("Positive#1", "Positive[x]", "Positive[x]");
    audit_chk("Positive#1", "Positive[{1, -2., 0, 3/2}]", "{True, False, False, True}");
}

/* EndPackage#1 (high): "EndPackage with an unclosed Begin[\"`Private`\"]
 * leaves $Context and $ContextPath in package mode".  EndPackage unwinds the
 * open Begin frames and restores the state saved by BeginPackage. */
static void test_endpackage_unclosed_begin(void) {
    audit_run("auditCtx0 = $Context; auditPath0 = $ContextPath");
    audit_run("BeginPackage[\"AuditPkgA`\"]");
    audit_run("Begin[\"`Private`\"]");
    audit_run("EndPackage[]");
    audit_chk("EndPackage#1", "$Context === auditCtx0", "True");
    audit_chk("EndPackage#1", "$ContextPath === Prepend[auditPath0, \"AuditPkgA`\"]", "True");
    /* The ordinary BeginPackage / Begin / End / EndPackage sequence. */
    audit_run("BeginPackage[\"AuditPkgB`\"]");
    audit_run("Begin[\"`Private`\"]");
    audit_run("End[]");
    audit_run("EndPackage[]");
    audit_chk("EndPackage#1", "$Context === auditCtx0", "True");
    audit_chk("EndPackage#1", "First[$ContextPath]", "\"AuditPkgB`\"");
    /* EndPackage with no open package: EndPackage::noctx, and a lone Begin
     * frame is left for End[] to close. */
    audit_run("Begin[\"AuditCtxC`\"]");
    audit_chk("EndPackage#1", "Check[EndPackage[], fired, EndPackage::noctx]", "fired");
    audit_chk("EndPackage#1", "$Context", "\"AuditCtxC`\"");
    audit_run("End[]");
    audit_chk("EndPackage#1", "$Context === auditCtx0", "True");
}

/* MemberQ#1 (medium): -n is levels 1 through -n. */
static void test_s14_memberq_levels(void) {
    audit_chk("MemberQ#1", "MemberQ[{{a}}, {a}, -1]", "True");
    audit_chk("MemberQ#1", "MemberQ[{{a}}, {a}, {1, -1}]", "True");
    audit_chk("MemberQ#1", "MemberQ[{{a}}, {a}, {-1}]", "False");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("MemberQ#3", "Quiet[Check[MemberQ[{1}, 1, x], fired, MemberQ::level]]", "fired");
    audit_chk("MemberQ#3", "Quiet[Check[MemberQ[{1}, 1, 1, 1], fired, MemberQ::nonopt]]", "fired");
    audit_chk("MemberQ#3", "Quiet[Check[MemberQ[], fired, MemberQ::argb]]", "fired");
    audit_chk("FreeQ#1", "Quiet[Check[FreeQ[f[x], x, a], fired, FreeQ::level]]", "fired");
    audit_chk("FreeQ#1", "Quiet[Check[FreeQ[f[x], x, Foo -> 1], fired, FreeQ::optx]]", "fired");
    audit_chk("FreeQ#1", "Quiet[Check[FreeQ[f[x], x, 1, 2], fired, FreeQ::nonopt]]", "fired");
    audit_chk("FreeQ#1", "FreeQ[f[a], f, Heads -> False]", "True");
    audit_chk("SetAttributes#3", "Quiet[Check[SetAttributes[f], fired, SetAttributes::argr]]", "fired");
    audit_chk("SetAttributes#3", "Quiet[SetAttributes[f]]", "SetAttributes[f]");
    audit_chk("Activate#2", "Quiet[Check[Activate[a, b, c], fired, Activate::argt]]", "fired");
    audit_chk("Activate#2", "Quiet[Activate[a, b, c]]", "Activate[a, b, c]");
    audit_chk("AtomQ#2", "Quiet[Check[AtomQ[1, 2], fired, AtomQ::argx]]", "fired");
    audit_chk("AtomQ#2", "Quiet[AtomQ[1, 2]]", "AtomQ[1, 2]");
    audit_chk("ByteCount#3", "Quiet[Check[ByteCount[1, 2], fired, ByteCount::argx]]", "fired");
    audit_chk("ByteCount#3", "Quiet[ByteCount[1, 2]]", "ByteCount[1, 2]");
    audit_chk("ClearAttributes#5", "Quiet[Check[ClearAttributes[f, Listable, Flat], fired, ClearAttributes::argrx]]", "fired");
    audit_chk("ClearAttributes#5", "Quiet[ClearAttributes[f, Listable, Flat]]", "ClearAttributes[f, Listable, Flat]");
    audit_chk("EvenQ#1", "Quiet[Check[EvenQ[2, 3], fired, EvenQ::argx]]", "fired");
    audit_chk("EvenQ#1", "Quiet[EvenQ[2, 3]]", "EvenQ[2, 3]");
    audit_chk("ExactNumberQ#2", "Quiet[Check[ExactNumberQ[], fired, ExactNumberQ::argx]]", "fired");
    audit_chk("ExactNumberQ#2", "Quiet[ExactNumberQ[]]", "ExactNumberQ[]");
    audit_chk("Hash#2", "Quiet[Check[Hash[1, \"Expression\", \"Integer\", 4], fired, Hash::argb]]", "fired");
    audit_chk("Hash#2", "Quiet[Hash[1, \"Expression\", \"Integer\", 4]]", "Hash[1, \"Expression\", \"Integer\", 4]");
    audit_chk("InexactNumberQ#2", "Quiet[Check[InexactNumberQ[1., 2.], fired, InexactNumberQ::argx]]", "fired");
    audit_chk("IntegerQ#1", "Quiet[Check[IntegerQ[1, 2], fired, IntegerQ::argx]]", "fired");
    audit_chk("IntegerQ#1", "Quiet[IntegerQ[1, 2]]", "IntegerQ[1, 2]");
    audit_chk("ListQ#1", "Quiet[Check[ListQ[{1}, {2}], fired, ListQ::argx]]", "fired");
    audit_chk("ListQ#1", "Quiet[ListQ[{1}, {2}]]", "ListQ[{1}, {2}]");
    audit_chk("MachineNumberQ#2", "Quiet[Check[MachineNumberQ[], fired, MachineNumberQ::argx]]", "fired");
    audit_chk("MachineNumberQ#2", "Quiet[MachineNumberQ[]]", "MachineNumberQ[]");
    audit_chk("MatrixQ#1", "Quiet[Check[MatrixQ[{{1}}, NumberQ, 3], fired, MatrixQ::argt]]", "fired");
    audit_chk("MatrixQ#1", "Quiet[MatrixQ[{{1}}, NumberQ, 3]]", "MatrixQ[{{1}}, NumberQ, 3]");
    audit_chk("NumberQ#2", "Quiet[Check[NumberQ[1, 2], fired, NumberQ::argx]]", "fired");
    audit_chk("NumberQ#2", "Quiet[NumberQ[1, 2]]", "NumberQ[1, 2]");
    audit_chk("NumericQ#4", "Quiet[Check[NumericQ[1, 2], fired, NumericQ::argx]]", "fired");
    audit_chk("NumericQ#4", "Quiet[NumericQ[1, 2]]", "NumericQ[1, 2]");
    audit_chk("OddQ#1", "Quiet[Check[OddQ[], fired, OddQ::argx]]", "fired");
    audit_chk("OddQ#1", "Quiet[OddQ[]]", "OddQ[]");
    audit_chk("StringQ#1", "Quiet[Check[StringQ[Sequence[\"a\", \"b\"]], fired, StringQ::argx]]", "fired");
    audit_chk("StringQ#1", "Quiet[StringQ[Sequence[\"a\", \"b\"]]]", "StringQ[\"a\", \"b\"]");
}

/* Wave 4 G4: transparent display wrappers in TeX, Unevaluated through MatchQ,
 * attributes and options in Definition / Information, OutputForm and
 * StringForm, genuine-number tests for Rational/Complex, context-name checks,
 * composite Hold attributes, ComposeList evaluated step by step, Names
 * patterns, NumericQ of infinities, Row separators, SubValues[s] = list,
 * Symbol names and SyntaxQ on comments. */
static void test_wave4_g4_expression(void) {
    audit_run("ClearAll[f, x, a, b, gi4, f3g4, f4g4, gca4, hca4, k]");
    audit_chk("HoldForm#1", "ToString[TeXForm[HoldForm[1 + 1]]]", "\"1+1\"");
    audit_chk("HoldForm#1", "ToString[HoldForm[x^2], TeXForm]", "\"x^{2}\"");
    audit_chk("TraditionalForm#1", "ToString[TeXForm[TraditionalForm[x^2]]]", "\"x^{2}\"");
    audit_chk("HoldPattern#1", "MatchQ[Unevaluated[1 + 1], HoldPattern[1 + 1]]", "True");
    audit_run("SetAttributes[gi4, Listable]");
    audit_chk("Information#1", "Information[gi4]", "\"Attributes[gi4] = {Listable}\"");
    audit_chk("Print#1", "ToString[OutputForm[x^2]]", "\"x^2\"");
    audit_chk("Print#2", "ToString[StringForm[\"`` and ``\", a, b]]", "\"a and b\"");
    audit_chk("AtomQ#1", "{AtomQ[Rational[a, b]], AtomQ[Complex[a, b]], AtomQ[1/2], AtomQ[I]}", "{False, False, True, True}");
    audit_chk("Begin#1", "Check[Begin[\"1abc`\"], fired, Begin::cxt]", "fired");
    audit_chk("Begin#1", "Quiet[{Begin[\"a b`\"], Begin[\"a``\"], $Context}]", "{Begin[\"a b`\"], Begin[\"a``\"], \"Global`\"}");
    audit_run("SetAttributes[gca4, HoldAll]; ClearAttributes[gca4, HoldFirst]");
    audit_chk("ClearAttributes#1", "Attributes[gca4]", "{HoldAll}");
    audit_run("SetAttributes[hca4, HoldFirst]; ClearAttributes[hca4, HoldAll]");
    audit_chk("ClearAttributes#2", "Attributes[hca4]", "{HoldFirst}");
    audit_chk("ClearAttributes#2", "SetAttributes[hca4, HoldAll]; ClearAttributes[hca4, HoldAll]; Attributes[hca4]", "{}");
    audit_chk("ComposeList#1", "k = 0; ComposeList[{(k++; #1) &, (k++; #1) &, (k++; #1) &}, x]; k", "3");
    audit_chk("ComposeList#2", "Last[ComposeList[Table[#1 + 1 &, {2000}], 0]]", "2000");
    audit_chk("ComposeList#2", "ComposeList[{f, g}, x]", "{x, f[x], g[f[x]]}");
    audit_chk("Context#1", "Context[\"`x\"]", "\"Global`\"");
    audit_run("SetAttributes[f3g4, Listable]; f3g4[x_] := x^2");
    audit_chk("Definition#1", "ToString[Definition[f3g4]]", "\"Attributes[f3g4] = {Listable}\\n\\nf3g4[x_] := x^2\"");
    audit_run("Options[f4g4] = {opt -> 1}; f4g4[x_] := x");
    audit_chk("Definition#1", "ToString[Definition[f4g4]]", "\"f4g4[x_] := x\\n\\nOptions[f4g4] = {opt -> 1}\"");
    audit_chk("ExactNumberQ#1", "{ExactNumberQ[Rational[a, b]], ExactNumberQ[1/2]}", "{False, True}");
    audit_chk("InexactNumberQ#1", "{InexactNumberQ[Complex[a, 1.]], InexactNumberQ[1. + I]}", "{False, True}");
    audit_chk("Names#1", "Names[\"sin\", IgnoreCase -> True]", "{\"Sin\"}");
    audit_chk("Names#2", "Names[\"Sin\" | \"Cos\"]", "{\"Cos\", \"Sin\"}");
    audit_chk("Names#2", "Names[\"ArcSin\" ~~ __]", "{\"ArcSinh\"}");
    audit_chk("NumberQ#1", "{NumberQ[Complex[1, x]], NumberQ[Rational[1, x]], NumberQ[Complex[Pi, 1]], NumberQ[Rational[Pi, 2]], NumberQ[2 + 3 I]}", "{False, False, False, False, True}");
    audit_chk("NumericQ#1", "{NumericQ[Infinity], NumericQ[-Infinity], NumericQ[ComplexInfinity], NumericQ[DirectedInfinity[1]], NumericQ[Pi]}", "{False, False, False, False, True}");
    audit_chk("NumericQ#2", "{NumericQ[Rational[x, 2]], NumericQ[Complex[x, 1]], NumericQ[Complex[Pi, 1]]}", "{False, False, True}");
    audit_chk("Row#1", "ToString[Row[{a, b, c}, x]]", "\"axbxc\"");
    audit_chk("SubValues#1", "SubValues[usv4] = {HoldPattern[usv4[x_][y_]] :> x - y}; usv4[5][2]", "3");
    audit_chk("SubValues#1", "SubValues[usv4]", "{HoldPattern[usv4[x_][y_]] :> x - y}");
    audit_chk("Symbol#1", "Check[Symbol[\"`\"], fired, Symbol::symname]", "fired");
    audit_chk("Symbol#1", "Quiet[Symbol[\"a`\"]]", "Symbol[\"a`\"]");
    audit_chk("SyntaxQ#1", "{SyntaxQ[\"(* unterminated\"], SyntaxQ[\"(* a (* nested *) b\"], SyntaxQ[\"1 (* c *)\"], SyntaxQ[\"\\\"(*\\\"\"]}", "{False, False, True, True}");
    audit_chk("EndPackage#2", "Begin[\"XepG4`\"]; Quiet[EndPackage[]]; c0 = $Context; End[]; c0", "\"XepG4`\"");
    audit_chk("InputForm#1", "ToString[Hold[(a -> b) -> c], InputForm]", "\"Hold[(a -> b) -> c]\"");
    audit_chk("InputForm#1", "FullForm[ToExpression[ToString[Hold[(a -> b) -> c], InputForm]]]", "Hold[Rule[Rule[a, b], c]]");
    audit_chk("InputForm#1", "{ToString[Hold[(a = b) = c], InputForm], ToString[Hold[(a :> b) :> c], InputForm], ToString[Hold[a -> (b -> c)], InputForm]}", "{\"Hold[(a = b) = c]\", \"Hold[(a :> b) :> c]\", \"Hold[a -> b -> c]\"}");
    audit_chk("InputForm#2", "{ToString[Hold[(a == b) == c], InputForm], ToString[Hold[(a > b) > c], InputForm], ToString[Hold[(a === b) === c], InputForm]}", "{\"Hold[(a == b) == c]\", \"Hold[(a > b) > c]\", \"Hold[(a === b) === c]\"}");
    audit_chk("InputForm#2", "ToString[Hold[a == b == c], InputForm]", "\"Hold[a == b == c]\"");
    audit_chk("InputForm#3", "{ToString[Hold[Plus[]], InputForm], ToString[Hold[{Times[], 2}], InputForm], ToString[Hold[Plus[] + 1], InputForm], ToString[Hold[f[Plus[]]], InputForm]}", "{\"Hold[Plus[]]\", \"Hold[{Times[], 2}]\", \"Hold[Plus[] + 1]\", \"Hold[f[Plus[]]]\"}");
    audit_chk("FortranForm#1", "{ToString[FortranForm[1/10]], ToString[FortranForm[1/10^30]]}", "{\"0.1\", \"1.e-30\"}");
    audit_chk("FortranForm#3", "{ToString[FortranForm[N[1/3, 30]]], ToString[FortranForm[N[10^40, 20]]]}", "{\"0.333333333333333333333333333333\", \"1.e40\"}");
    audit_chk("NumberForm#1", "{ToString[NumberForm[0.1, 20]], ToString[NumberForm[N[Pi], 30]]}", "{\"0.1\", \"3.141592653589793\"}");
    audit_chk("PossibleZeroQ#3", "{PossibleZeroQ[a == b], PossibleZeroQ[x < y], PossibleZeroQ[a - a]}", "{False, False, True}");
}

int main(void) {
    audit_init();
    TEST(test_bytecount_deep_nesting);
    TEST(test_positive_exact_zero);
    TEST(test_endpackage_unclosed_begin);
    TEST(test_s14_memberq_levels);
    TEST(test_wave3_arity_options);
    TEST(test_wave4_g4_expression);
    return audit_finish("audit_expression_information");
}

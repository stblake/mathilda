/* Builtin-audit regressions: uncategorised heads (review/issues/<Head>.md;
 * tests/test_audit.h). */
#include "test_audit.h"

/* SlotSequence#1 (high): "## with no remaining arguments is left
 * unsubstituted instead of becoming an empty sequence". */
static void test_slotsequence_empty(void) {
    audit_run("ClearAll[f, a, b, c]");
    audit_chk("SlotSequence#1", "{##} &[]", "{}");
    audit_chk("SlotSequence#1", "{a, ##} &[]", "{a}");
    audit_chk("SlotSequence#1", "{##3} &[a, b]", "{}");
    audit_chk("SlotSequence#1", "f[## &[]]", "f[]");
    /* sanity */
    audit_chk("SlotSequence#1", "{##2} &[a, b, c]", "{b, c}");
    audit_chk("SlotSequence#1", "f[##] &[a, b]", "f[a, b]");
}

/* Decrement#1 (high): "Decrement evaluates a Part l-value's index twice:
 * returns one element and decrements another".  Shared by every
 * read-modify-write operator (increment_core, AppendTo/PrependTo). */
static void test_decrement_part_once(void) {
    audit_run("ClearAll[di, dl]; di = 1; dl = {10, 20, 30}");
    audit_chk("Decrement#1", "{dl[[di++]]--, di, dl}", "{10, 2, {9, 20, 30}}");
    audit_run("di = 1; dl = {10, 20, 30}");
    audit_chk("Decrement#1", "{dl[[di++]]++, di, dl}", "{10, 2, {11, 20, 30}}");
    audit_run("di = 1; dl = {10, 20, 30}");
    audit_chk("Decrement#1", "{++dl[[di++]], di, dl}", "{11, 2, {11, 20, 30}}");
    audit_run("di = 1; dl = {10, 20, 30}");
    audit_chk("Decrement#1", "{--dl[[di++]], di, dl}", "{9, 2, {9, 20, 30}}");
    audit_run("di = 1; dl = {{1}, {2}}");
    audit_chk("Decrement#1", "AppendTo[dl[[di++]], 9]; {dl, di}", "{{{1, 9}, {2}}, 2}");
    audit_run("di = 1; dl = {{1}, {2}}");
    audit_chk("Decrement#1", "PrependTo[dl[[di++]], 9]; {dl, di}", "{{{9, 1}, {2}}, 2}");
    audit_run("di = 1; dl = {{10, 20}, {30, 40}}");
    audit_chk("Decrement#1", "dl[[di++, di++]] -= 1; {dl, di}", "{{{10, 19}, {30, 40}}, 3}");
    /* sanity */
    audit_run("di = 5");
    audit_chk("Decrement#1", "{di--, di}", "{5, 4}");
    audit_run("dl = {1, 2, 3}");
    audit_chk("Decrement#1", "dl[[2]] += 10; dl", "{1, 12, 3}");
    audit_chk("Decrement#1", "Catch[dl[[Throw[7]]]++]", "7");
}

/* Divide#2 (high): "Divide[0., 0.] (and 0/0., 0./0) gives ComplexInfinity
 * instead of Indeterminate". */
static void test_divide_zero_zero_real(void) {
    audit_chk("Divide#2", "Quiet[Divide[0., 0.]]", "Indeterminate");
    audit_chk("Divide#2", "Quiet[Divide[0, 0.]]", "Indeterminate");
    audit_chk("Divide#2", "Quiet[Divide[0., 0]]", "Indeterminate");
    audit_chk("Divide#2", "Check[Divide[0., 0.], fired, Infinity::indet]", "fired");
    /* sanity: a nonzero dividend is still ComplexInfinity */
    audit_chk("Divide#2", "Quiet[Divide[1., 0.]]", "ComplexInfinity");
    audit_chk("Divide#2", "Quiet[Divide[0, 0]]", "Indeterminate");
    audit_chk("Divide#2", "Divide[1., 4]", "0.25");
}

/* TeXForm#1 (high): "2 3^x typesets as 23^x".  Math mode ignores spaces,
 * so a factor that starts with a digit gets an explicit "\ ". */
static void test_texform_digit_factor(void) {
    /* "\\ " in the language is one backslash and a space */
    audit_true("TeXForm#1", "ToString[TeXForm[2*3^x]] === \"2\\\\ 3^{x}\"");
    audit_true("TeXForm#1", "ToString[TeXForm[2*3^x/y]] === \"\\\\frac{2\\\\ 3^{x}}{y}\"");
    /* sanity: a symbol factor needs no explicit space */
    audit_true("TeXForm#1", "ToString[TeXForm[2 x]] === \"2 x\"");
}

/* ListConvolve#1 (medium): the level argument lev below the list's rank. */
static void test_s14_listconvolve_level(void) {
    audit_run("ClearAll[a, b, c, d, e, f, x, y]");
    audit_true("ListConvolve#1",
               "Expand[ListConvolve[{x, y}, {{a, b}, {c, d}, {e, f}}, {-1, 1}, 0, Times, Plus, 1]] === "
               "Expand[{{c x + a y, d x + b y}, {e x + c y, f x + d y}}]");
    audit_true("ListConvolve#1", "ListConvolve[{1, 2}, {3, 4, 5}, {-1, 1}, 0, Times, Plus, 1] === "
                                 "ListConvolve[{1, 2}, {3, 4, 5}]");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("Divide#5", "Quiet[Check[Divide[x], fired, Divide::argr]]", "fired");
    audit_chk("Divide#5", "Quiet[Divide[x]]", "Divide[x]");
    audit_chk("RuleDelayed#3", "Quiet[Check[RuleDelayed[a], fired, RuleDelayed::argr]]", "fired");
    audit_chk("RuleDelayed#3", "Quiet[RuleDelayed[a]]", "RuleDelayed[a]");
    audit_chk("Subtract#3", "Quiet[Check[Subtract[a, b, c], fired, Subtract::argrx]]", "fired");
    audit_chk("Subtract#3", "Quiet[Subtract[a, b, c]]", "Subtract[a, b, c]");
    audit_chk("BooleanQ#1", "Quiet[Check[BooleanQ[Sequence[]], fired, BooleanQ::argx]]", "fired");
    audit_chk("BooleanQ#1", "Quiet[BooleanQ[Sequence[]]]", "BooleanQ[]");
    audit_chk("Decrement#4", "Quiet[Check[Decrement[s12x, 1], fired, Decrement::argx]]", "fired");
    audit_chk("Decrement#4", "Quiet[Decrement[s12x, 1]]", "Decrement[s12x, 1]");
    audit_chk("Increment#2", "Quiet[Check[Increment[a, b], fired, Increment::argx]]", "fired");
    audit_chk("Increment#2", "Quiet[Increment[a, b]]", "Increment[a, b]");
    audit_chk("PreDecrement#4", "Quiet[Check[PreDecrement[s12x, 1], fired, PreDecrement::argx]]", "fired");
    audit_chk("PreDecrement#4", "Quiet[PreDecrement[s12x, 1]]", "PreDecrement[s12x, 1]");
    audit_chk("PreIncrement#3", "Quiet[Check[PreIncrement[], fired, PreIncrement::argx]]", "fired");
    audit_chk("PreIncrement#3", "Quiet[PreIncrement[]]", "PreIncrement[]");
    audit_chk("SplitBy#2", "Quiet[Check[SplitBy[{1}], fired, SplitBy::argr]]", "fired");
    audit_chk("SplitBy#2", "Quiet[SplitBy[{1}]]", "SplitBy[{1}]");
    audit_chk("SymbolName#2", "Quiet[Check[SymbolName[Sequence[a, b]], fired, SymbolName::argx]]", "fired");
    audit_chk("SymbolName#2", "Quiet[SymbolName[Sequence[a, b]]]", "SymbolName[a, b]");
    audit_chk("FromPackedArray#2", "Quiet[Check[FromPackedArray[1, 2], fired, FromPackedArray::argx]]", "fired");
    audit_chk("FromPackedArray#2", "Quiet[FromPackedArray[1, 2]]", "FromPackedArray[1, 2]");
    audit_chk("ToNDArray#4", "Quiet[Check[ToNDArray[{1, 2}, Foo -> 1], fired, ToNDArray::optx]]", "fired");
    audit_chk("IntegerPart#3", "Quiet[Check[IntegerPart[1, 2], fired, IntegerPart::argx]]", "fired");
    audit_chk("IntegerPart#3", "Quiet[IntegerPart[1, 2]]", "IntegerPart[1, 2]");
}

/* ================================================================== */
/* Wave 3, shared fix S11: Indeterminate / infinity / pole values     */
/* (src/infinity_rules.c and the per-head rules it does not cover).   */
/* ================================================================== */

static void test_s11_special_values(void) {
    /* Subtract#1: infinities in different complex directions have no sum */
    audit_true("Subtract#1", "Quiet[Subtract[DirectedInfinity[I], DirectedInfinity[I]]] === Indeterminate");
    audit_true("Subtract#1", "Check[DirectedInfinity[I] - DirectedInfinity[I], \"fired\", Infinity::indet] === \"fired\"");
    audit_true("Subtract#1", "Times[-1, DirectedInfinity[I]] === DirectedInfinity[-I]");
    audit_true("Subtract#1", "DirectedInfinity[I] + DirectedInfinity[I] === DirectedInfinity[I]");
    audit_true("Subtract#1", "Quiet[DirectedInfinity[1] + DirectedInfinity[I]] === Indeterminate");
    /* Divide#4: a complex number over 0 */
    audit_true("Divide#4", "Quiet[Divide[I, 0]] === ComplexInfinity");
    audit_true("Divide#4", "Quiet[Divide[2 + 3 I, 0]] === ComplexInfinity");
    audit_true("Divide#4", "Check[Divide[I, 0], \"fired\", Power::infy] === \"fired\"");
    /* SphericalHankelH2#2: a pole at 0 for every real order */
    audit_true("SphericalHankelH2#2", "SphericalHankelH2[-1, 0] === ComplexInfinity");
    audit_true("SphericalHankelH2#2", "SphericalHankelH2[1/2, 0] === ComplexInfinity");
    audit_true("SphericalHankelH2#2", "SphericalHankelH1[-1, 0] === ComplexInfinity");
    audit_true("SphericalHankelH2#2", "SphericalHankelH2[0, 0] === ComplexInfinity");
    /* SphericalHankelH2#3: an inexact zero argument */
    audit_true("SphericalHankelH2#3", "SphericalHankelH2[0, 0.] === ComplexInfinity");
}

/* Wave 4 G8: uncategorised findings. */
static void test_wave4_g8_uncategorised(void) {
    audit_chk("RuleDelayed#1", "ToString[InputForm[(a :> b) :> c]]", "\"(a :> b) :> c\"");
    audit_chk("RuleDelayed#1", "ToExpression[ToString[InputForm[(a :> b) :> c]]] === ((a :> b) :> c)", "True");
    /* Decompose: no linear component, a composite pure power split */
    audit_chk("Decompose#1", "Decompose[x^2 + 1, x]", "{1 + x^2}");
    audit_chk("Decompose#1", "Decompose[x^3 + 2, x]", "{2 + x^3}");
    audit_chk("Decompose#1", "Decompose[(x^2 + x)^3 + 1, x]", "{1 + x^3, x + x^2}");
    audit_chk("Decompose#1", "Map[Exponent[#, x] &, Decompose[Expand[(2 x^2 + x)^2], x]]", "{2, 2}");
    audit_chk("Decompose#1", "Expand[Fold[#2 /. x -> #1 &, x, Reverse[Decompose[Expand[(2 x^2 + x)^2], x]]] - (2 x^2 + x)^2]", "0");
    audit_chk("Decompose#2", "Decompose[x^4, x]", "{x^2, x^2}");
    audit_chk("Decompose#2", "Decompose[x^9, x]", "{x^3, x^3}");
    audit_chk("Decompose#2", "Decompose[x^8 + x^4, x]", "{x + x^2, x^2, x^2}");
    audit_chk("Decompose#2", "Decompose[x^6 + 1, x]", "{1 + x^2, x^3}");
    audit_chk("Decompose#2", "Decompose[x + 1, x]", "{1 + x}");
    /* the increment family on indexed l-values */
    audit_chk("Decrement#2", "w4f[1] = 5; {w4f[1]--, w4f[1]}", "{5, 4}");
    audit_chk("PreDecrement#1", "w4g[1] = 3; {--w4g[1], w4g[1]}", "{2, 2}");
    audit_chk("PreIncrement#1", "w4c[_] = 0; {++w4c[a], w4c[a]}", "{1, 1}");
    audit_chk("PreIncrement#1", "w4h[1] = 5; ++w4h[1]", "6");
    audit_chk("PreIncrement#1", "w4t = <|\"a\" -> 1|>; {++w4t[\"a\"], w4t}", "{2, <|\"a\" -> 2|>}");
    /* End[] restores $Context only */
    audit_chk("End#1", "Begin[\"W4A`\"]; AppendTo[$ContextPath, \"W4B`\"]; End[]; MemberQ[$ContextPath, \"W4B`\"]", "True");
    audit_run("$ContextPath = DeleteCases[$ContextPath, \"W4B`\"]");
    /* Series of a Jacobi quotient function at its pole */
    audit_chk("JacobiNS#1", "Head[Series[JacobiNS[u, m], {u, 0, 3}]]", "SeriesData");
    audit_chk("JacobiNS#1", "Expand[Normal[Series[JacobiNS[u, m], {u, 0, 3}]] - (1/u + (1 + m) u/6 + (7 - 22 m + 7 m^2) u^3/360)]", "0");
    audit_chk("JacobiNS#1", "Expand[Normal[Series[JacobiDS[u, m], {u, 0, 2}]] - (1/u + (1 - 2 m) u/6)]", "0");
    /* PolynomialSqrt with algebraic-number coefficients */
    audit_chk("PolynomialSqrt#1", "PolynomialSqrt[Expand[(x + Sqrt[2])^2]]", "Sqrt[2] + x");
    audit_chk("PolynomialSqrt#1", "PolynomialSqrt[x^4 + 2 Sqrt[2] x^2 + 2]", "Sqrt[2] + x^2");
    audit_chk("PolynomialSqrt#1", "PolynomialSqrt[x^2 + 1]", "$Failed");
    /* RatCanonPrototype's fresh names avoid the input's own symbols */
    audit_chk("RatCanonPrototype#1", "RatCanonPrototype[$rcp0$/Log[y]]", "$rcp0$/Log[y]");
    /* SphericalHankelH1 just off the real axis */
    audit_chk("SphericalHankelH1#1", "Abs[SphericalHankelH1[0, 1. + 1.*^-8 I] - N[-I Exp[I z]/z /. z -> 1 + I/10^8, 20]] < 10^-14", "True");
    audit_chk("SphericalHankelH1#1", "Abs[SphericalHankelH1[2, 1.5 + 0.3 I] - N[SphericalHankelH1[2, 3/2 + 3/10 I], 20]] < 10^-14", "True");
    audit_chk("SphericalHankelH1#1", "Abs[SphericalHankelH2[3, 0.5 - 2. I] - N[SphericalHankelH2[3, 1/2 - 2 I], 20]] < 10^-13", "True");
    /* Subsets by position */
    audit_chk("Subsets#1", "Subsets[{a, b, c, d}, All, {5}]", "{{d}}");
    audit_chk("Subsets#1", "Subsets[{a, b, c, d}, All, {3, 5}]", "{{b}, {c}, {d}}");
    audit_chk("Subsets#1", "Subsets[{a, b, c, d}, All, {2, 8, 3}]", "{{a}, {d}, {a, d}}");
    /* TeXForm */
    audit_chk("TeXForm#3", "ToString[TeXForm[\"a_b%\"]]", "\"\\\\text{a\\\\_b\\\\%}\"");
    audit_chk("TeXForm#3", "ToString[TeXForm[{x, \"a&b\"}]]", "\"\\\\{x,\\\\text{a\\\\&b}\\\\}\"");
    audit_chk("TeXForm#4", "ToString[TeXForm[\\[Alpha]]]", "\"\\\\alpha \"");
    audit_chk("TeXForm#4", "ToString[TeXForm[\\[Beta] + \\[Gamma]]]", "\"\\\\beta +\\\\gamma \"");
    audit_chk("TeXForm#5", "ToString[TeXForm[HoldForm[Integrate[f[x], {x, 0, 1}]]]]", "\"\\\\int_{0}^{1} f(x) \\\\, dx\"");
    audit_chk("TeXForm#5", "ToString[TeXForm[HoldForm[Sum[a[i], {i, 1, n}]]]]", "\"\\\\sum _{i=1}^{n} a(i)\"");
    audit_chk("TeXForm#5", "ToString[TeXForm[HoldForm[Limit[f[x], x -> 0]]]]", "\"\\\\lim_{x\\\\to 0} f(x)\"");
    audit_chk("TeXForm#5", "ToString[TeXForm[Integrate[f[x], x]]]", "\"\\\\int f(x) \\\\, dx\"");
}

int main(void) {
    audit_init();
    TEST(test_slotsequence_empty);
    TEST(test_decrement_part_once);
    TEST(test_divide_zero_zero_real);
    TEST(test_texform_digit_factor);
    TEST(test_s14_listconvolve_level);
    TEST(test_wave3_arity_options);
    TEST(test_s11_special_values);
    TEST(test_wave4_g8_uncategorised);
    return audit_finish("audit_uncategorised");
}

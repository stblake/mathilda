/* Builtin-audit regressions, shared root cause S7: attributes and Protected
 * (review/issues/<Head>.md; tests/test_audit.h).  One suite for the class
 * rather than one per category, because the findings span ten categories
 * and share three fixes: the attribute table (src/attr.c), the assignment
 * target of f[] (src/eval.c), and the Protected checks of the non-Set
 * writers (Clear, ClearAll, Unset, AppendTo/PrependTo, Options, Default). */
#include "test_audit.h"

/* List#1 Blank#1 BlankSequence#1 BlankNullSequence#1 Head#1 Append#1
 * Delete#5 First#1 Last#1 Most#1 Rest#1 Longest#1 Shortest#2 Through#1
 * Complex#5 Length#2 Optional#2 Insert#5 Prepend#2 Dimensions#3 UpTo#2:
 * System symbols with no Protected attribute. */
static void test_system_symbols_protected(void) {
    static const char* heads[] = {
        "List", "Blank", "BlankSequence", "BlankNullSequence", "Head", "Append",
        "Delete", "First", "Last", "Most", "Rest", "Longest", "Shortest",
        "Through", "Complex", "Length", "Optional", "Insert", "Prepend",
        "Dimensions", "UpTo", "Integer", "Real", "String", "True", "False",
        "Null", "Infinity", "All", "None", "Automatic", "Integers", "Reals",
        NULL };
    char buf[160];
    for (int i = 0; heads[i]; i++) {
        snprintf(buf, sizeof buf, "MemberQ[Attributes[%s], Protected]", heads[i]);
        audit_true("S7-protected", buf);
    }
    audit_chk("List#1", "Check[List = 5, fired, Set::wrsym]", "fired");
    audit_chk("List#1", "Quiet[List = 5]; {3, 1, 2}", "{3, 1, 2}");
    audit_chk("Blank#1", "Quiet[Blank[] = 5]; MatchQ[1, _]", "True");
    audit_chk("BlankSequence#1", "Quiet[BlankSequence[] = 5]; MatchQ[f[1, 2], f[__]]", "True");
    audit_chk("First#1", "Quiet[First = 5]; First[{7, 8}]", "7");
    /* Deliberately settable: the $-variables, MakeBoxes, Derivative. */
    audit_chk("S7-protected", "MemberQ[Attributes[$RecursionLimit], Protected]", "False");
    audit_chk("Derivative#5", "Derivative[1][hh] = cc; hh'[x]", "cc[x]");
}

/* Attribute sets that were wrong rather than missing. */
static void test_attribute_corrections(void) {
    /* Equivalent#1: not Flat -- the nested equivalence is a truth value. */
    audit_chk("Equivalent#1", "Attributes[Equivalent]", "{Orderless, Protected}");
    audit_chk("Equivalent#1", "Equivalent[a, Equivalent[b, c]]", "Equivalent[a, Equivalent[b, c]]");
    audit_chk("Equivalent#1", "Equivalent[False, Equivalent[False, False]]", "False");
    audit_chk("Equivalent#1", "Equivalent[True, Equivalent[False, False]]", "True");
    /* Join#2, Subtract#2, Attributes#4, Out#2.  Not changed: PossibleZeroQ#4
     * (deliberately not Listable -- it threads its first argument itself so an
     * Assumptions argument is never threaded; the list form already works) and
     * HeavisideTheta#4 (NumericFunction needs a Compile lowering first). */
    audit_true("Join#2", "SubsetQ[Attributes[Join], {Flat, OneIdentity}]");
    audit_chk("Join#2", "Join[{1}, Join[{2}, {3}]]", "{1, 2, 3}");
    audit_chk("Join#2", "Join[{{1}}, {{2}}, 2]", "{{1, 2}}");
    audit_true("Subtract#2", "MemberQ[Attributes[Subtract], OneIdentity]");
    audit_chk("Attributes#4", "Attributes[{Sin, Hold}]",
              "{{Listable, NumericFunction, Protected}, {HoldAll, Protected}}");
    audit_chk("S7-protected", "PossibleZeroQ[{0, 1, Sin[Pi]}]", "{True, False, True}");
    audit_true("Out#2", "SubsetQ[Attributes[Out], {Listable, NHoldFirst, Protected}]");
    audit_chk("Out#2", "Out[{1, 2}]", "{Out[1], Out[2]}");
    /* FirstCase#1: HoldRest, so the default runs only on a miss. */
    audit_chk("FirstCase#1", "n = 0; FirstCase[{2, 3}, _?EvenQ, n = 1; 0]; n", "0");
    audit_chk("FirstCase#1", "n = 0; FirstCase[{1, 3}, _?EvenQ, n = 1; -1]", "-1");
    audit_chk("FirstCase#1", "p = _?OddQ; FirstCase[{2, 3}, p]", "3");
}

/* Alphabet#1 Kernels#1 MemoryInUse#1 SessionTime#1 TimeUsed#1: f[] = v
 * found no assignment target and so bypassed Protected. */
static void test_zero_argument_lhs(void) {
    audit_chk("Alphabet#1", "Check[Alphabet[] = 3, fired, Set::write]", "fired");
    audit_chk("Alphabet#1", "Quiet[Alphabet[] = 3]; Length[Alphabet[]]", "26");
    audit_chk("Kernels#1", "Check[Kernels[] := 4, fired, SetDelayed::write]", "fired");
    audit_chk("MemoryInUse#1", "Quiet[MemoryInUse[] = -1]; MemoryInUse[] > 0", "True");
    audit_chk("SessionTime#1", "Quiet[SessionTime[] := -5]; SessionTime[] >= 0", "True");
    audit_chk("TimeUsed#1", "Quiet[TimeUsed[] = -5]; TimeUsed[] >= 0", "True");
    /* An unprotected f[] is still an ordinary definition. */
    audit_chk("Alphabet#1", "zz[] = 9; zz[]", "9");
}

/* Options#1 Default#2: the options and defaults of a Protected symbol. */
static void test_options_default_protected(void) {
    audit_chk("Options#1", "Check[Options[Sin] = {q -> 1}, fired, Set::write]", "fired");
    audit_chk("Options#1", "Quiet[Options[Sin] = {q -> 1}]; Options[Sin]", "{}");
    audit_chk("Options#1", "Options[oo] = {a -> 1}; Options[oo]", "{a -> 1}");
    audit_chk("Default#2", "Check[Default[Plus] = 1, fired, Set::write]", "fired");
    audit_chk("Default#2", "Quiet[Default[Plus] = 1]; Replace[y, x_. + y -> x]", "0");
    audit_chk("Default#2", "Default[gg] = 5; gg[x_, y_.] := {x, y}; gg[1]", "{1, 5}");
}

/* Clear#1 ClearAll#4 Unset#2 Unset#4 PrependTo#1: the other writers. */
static void test_other_writers(void) {
    audit_chk("Clear#1", "Check[g1[1] = 2; Protect[g1]; Clear[g1], fired, Clear::wrsym]", "fired");
    audit_chk("Clear#1", "g2[1] = 2; Protect[g2]; Quiet[Clear[g2]]; g2[1]", "2");
    audit_chk("Clear#1", "v1 = 1; Clear[v1]; v1", "v1");
    audit_chk("ClearAll#4", "Check[u1 = 1; Protect[u1]; ClearAll[u1], fired, ClearAll::wrsym]", "fired");
    audit_chk("ClearAll#4", "u2 = 1; Protect[u2]; Quiet[ClearAll[u2]]; u2", "1");
    audit_chk("ClearAll#4", "Check[ClearAll[Sin], fired, ClearAll::wrsym]", "fired");
    audit_chk("Unset#4", "Quiet[Sin =.]", "$Failed");
    audit_chk("Unset#4", "Check[Sin =., fired, Unset::wrsym]", "fired");
    audit_chk("Unset#2", "N[nn] = 3; N[nn] =.; N[nn]", "nn");
    audit_chk("Unset#2", "N[kk, 30] = 7; N[kk, 30] =.; N[kk, 30]", "kk");
    audit_chk("Unset#2", "Quiet[N[zz] =.]", "$Failed");
    audit_chk("PrependTo#1", "Check[pq = {1}; Protect[pq]; PrependTo[pq, 0], fired, Set::wrsym]", "fired");
    audit_chk("PrependTo#1", "pr = {1}; Protect[pr]; Quiet[PrependTo[pr, 0]]; pr", "{1}");
    audit_chk("PrependTo#1", "ps = {1}; Protect[ps]; Quiet[AppendTo[ps, 0]]; ps", "{1}");
    audit_chk("PrependTo#1", "pt = {1, 2}; AppendTo[pt, 3]; PrependTo[pt, 0]; pt", "{0, 1, 2, 3}");
    audit_chk("PrependTo#1", "pw = {{1}}; AppendTo[pw[[1]], 2]; pw", "{{1, 2}}");
}

int main(void) {
    audit_init();
    TEST(test_system_symbols_protected);
    TEST(test_attribute_corrections);
    TEST(test_zero_argument_lhs);
    TEST(test_options_default_protected);
    TEST(test_other_writers);
    return audit_finish("audit_protected");
}

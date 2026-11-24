/* Builtin-audit regressions: control flow (review/issues/<Head>.md;
 * tests/test_audit.h). */
#include "test_audit.h"

/* For#1 (high, hang): "Throw in For's start is swallowed, so
 * Catch[For[Throw[7], True, i++]] loops forever".  The start's value was
 * freed with no flow-control check.  The same holds for the test (For#2) and
 * for the evaluated iterator bounds of Do / Table / Sum / Product. */
static void test_for_start_throw(void) {
    audit_run("ClearAll[fi]");
    audit_chk("For#1", "Catch[For[Throw[7], True, fi++]]", "7");
    audit_chk("For#2", "Catch[For[fi = 0, Throw[5], fi++]]", "5");
    audit_chk("For#2", "Catch[For[fi = 0, If[fi == 2, Throw[fi], fi < 5], fi++]]", "2");
    audit_chk("For#1", "Quiet[CheckAbort[For[Abort[], True, fi++], aborted]]", "aborted");
    audit_chk("For#1", "Catch[Do[1, {fi, Throw[3]}]]", "3");
    audit_chk("For#1", "Catch[Table[1, {fi, 1, Throw[4]}]]", "4");
    audit_chk("For#1", "Catch[Sum[fi, {fi, Throw[6]}]]", "6");
    audit_chk("For#1", "Catch[Product[fi, {fi, Throw[8]}]]", "8");
    /* sanity: ordinary loops */
    audit_chk("For#1", "Module[{s = 0}, For[fi = 1, fi <= 4, fi++, s += fi]; s]", "10");
    audit_chk("For#1", "For[fi = 0, fi < 3, fi++]", "Null");
    audit_chk("For#1", "f0[] := For[fi = 0, True, fi++, If[fi == 3, Return[fi]]]; f0[]", "3");
    audit_chk("For#1", "Table[fi, {fi, 3}]", "{1, 2, 3}");
}

/* Break#1 / Continue#1 (high): "Break[] / Continue[] evaluated inside a
 * function argument does not exit / continue the enclosing loop". */
static void test_break_continue_in_argument(void) {
    audit_run("ClearAll[bcr, bci, bcg]");
    audit_chk("Break#1", "bcr = {}; Do[AppendTo[bcr, bci]; f[If[bci == 2, Break[]]], {bci, 3}]; bcr", "{1, 2}");
    audit_chk("Break#1", "Do[g[If[bci == 2, Break[]]], {bci, 5}]", "Null");
    audit_chk("Break#1", "bcr = 0; While[True, bcr++; h[If[bcr > 3, Break[]]]]; bcr", "4");
    audit_chk("Continue#1", "bcr = {}; Do[AppendTo[bcr, If[bci == 2, Continue[], bci]], {bci, 3}]; bcr", "{1, 3}");
    audit_chk("Continue#1", "bcr = {}; Do[Map[If[# == 2, Continue[]] &, {bci}]; AppendTo[bcr, bci], {bci, 3}]; bcr", "{1, 3}");
    audit_chk("Continue#1", "bcr = {}; For[bci = 0, bci < 4, bci++, AppendTo[bcr, If[OddQ[bci], Continue[], bci]]]; bcr", "{0, 2}");
    /* through a builtin that consumes the element without re-evaluating it:
     * inside a loop, Break / Continue are unwinds like Throw */
    audit_chk("Break#1", "bcr = 0; Do[Select[{1, 2}, (If[# == 2, Break[]]; True) &]; bcr = bci, {bci, 3}]; bcr", "0");
    audit_chk("Break#1", "bcr = 0; Do[Count[{1, 2}, _?(If[# == 2, Break[]; True] &)]; bcr = bci, {bci, 3}]; bcr", "0");
    audit_chk("Break#1", "bci = 0; While[True, bci++; Select[{bci}, (If[# > 2, Break[]]; True) &]]; bci", "3");
    audit_chk("Continue#1", "bcr = 0; Do[Sort[{2, 1}, (Continue[]; #1 < #2) &]; bcr++, {3}]; bcr", "0");
    audit_chk("Continue#1", "bcr = {}; For[bci = 1, bci <= 4, bci++, Scan[If[# == 3, Continue[]] &, {bci}]; AppendTo[bcr, bci]]; bcr", "{1, 2, 4}");
    /* an inner loop consumes its own Break; Catch does not trap one */
    audit_chk("Break#1", "Do[Do[If[bcj == 2, Break[]]; bcr = {bci, bcj}, {bcj, 3}], {bci, 2}]; bcr", "{2, 1}");
    audit_chk("Break#1", "Catch[Do[Catch[Break[]]; Throw[bad], {2}]; ok]", "ok");
    /* sanity: outside any loop the markers stay values; ordinary loops work */
    audit_chk("Break#1", "Quiet[f[Break[]]]", "f[Break[]]");
    audit_chk("Break#1", "Do[If[bci == 2, Break[]], {bci, 3}]", "Null");
    audit_chk("Break#1", "Catch[Do[Throw[bci], {bci, 3}]]", "1");
    audit_chk("Break#1", "Table[bci, {bci, 3}]", "{1, 2, 3}");
}

/* Abort#1 (high): "Abort[] inside an NIntegrate integrand is swallowed and
 * a wrong numeric value is returned". */
static void test_abort_in_nintegrate(void) {
    audit_chk("Abort#1", "CheckAbort[NIntegrate[If[x > 0.5, Abort[], x], {x, 0, 1}], c]", "c");
    audit_chk("Abort#1", "Catch[NIntegrate[If[x > 0.5, Throw[t], x], {x, 0, 1}]]", "t");
    /* sanity */
    audit_true("Abort#1", "Abs[NIntegrate[x, {x, 0, 1}] - 0.5] < 10^-10");
}

/* While#1 / For#3 (high): the compiled numeric loop path turned exact
 * variables into Reals.  It now engages only when every assigned variable is
 * provably inexact; the rest run in the interpreter. */
static void test_numloop_exact_types(void) {
    audit_run("ClearAll[nx, ny, nc, ni, nw]");
    audit_chk("While#1", "nx = 1.; ny = 0; While[nx < 3., nx = nx + 1.; ny = ny + 1]; {nx, ny, Head[ny]}",
              "{3.0, 2, Integer}");
    audit_chk("While#1", "nx = 1.; nc = 0; While[nx < 100., nx = nx*1.5; nc = nc + 1]; nc", "12");
    audit_chk("While#1", "nx = 0.; ny = 0; Do[nx = nx + 1.; ny = ny + 1, {2}]; {nx, ny}", "{2.0, 2}");
    audit_chk("While#1", "nx = 5; While[nx < 3., nx = nx + 1.]; nx", "5");
    audit_chk("For#3", "nx = 2.5; For[ni = 0, ni < 3, ni++, nx = ni]; {nx, Head[nx]}", "{2, Integer}");
    audit_chk("For#3", "ny = 0.5; For[ni = 0, ni < 2, ni++, ny = ni/2]; ny", "1/2");
    audit_chk("For#3", "nw = 1.; For[ni = 0, ni < 3, ni++, nw = Sqrt[ni]]; nw", "Sqrt[2]");
    audit_chk("For#3", "nx = 2.5; For[ni = 0, ni < 3, ni++, ny = nx; nx = ni]; {nx, ny}", "{2, 1}");
    /* sanity: real loops keep the fast path's value */
    audit_chk("For#3", "nx = 0.; Do[nx = nx + 0.5, {10}]; nx", "5.0");
    audit_chk("For#3", "nx = 0; Do[nx = nx + ni, {ni, 10}]; nx", "55");
    audit_chk("For#3", "nx = 0.; Do[nx = nx + ni, {ni, 10}]; nx", "55.0");
}

/* ConditionalExpression#1 (high): "not propagated through arithmetic
 * (Plus/Times/Power) or NumericFunction heads". */
static void test_conditionalexpression_propagates(void) {
    audit_run("ClearAll[x, y, c]");
    audit_true("ConditionalExpression#1",
               "1 + ConditionalExpression[x, x > 0] === ConditionalExpression[1 + x, x > 0]");
    audit_true("ConditionalExpression#1",
               "ConditionalExpression[x, x > 0] + ConditionalExpression[y, y > 0] === "
               "ConditionalExpression[x + y, x > 0 && y > 0]");
    audit_true("ConditionalExpression#1",
               "Sin[ConditionalExpression[x, x > 0]] === ConditionalExpression[Sin[x], x > 0]");
    audit_true("ConditionalExpression#1",
               "2 ConditionalExpression[x, x > 0] === ConditionalExpression[2 x, x > 0]");
    audit_true("ConditionalExpression#1",
               "ConditionalExpression[x, x > 0]^2 === ConditionalExpression[x^2, x > 0]");
    /* sanity: a decided condition still reduces */
    audit_chk("ConditionalExpression#1", "1 + ConditionalExpression[x, True]", "1 + x");
    audit_chk("ConditionalExpression#1", "ConditionalExpression[x, x > 0] /. x -> 2", "2");
}

/* Continue#2 (low): a Continue[] at top level is reported and held. */
static void test_s14_continue_toplevel(void) {
    mth_msg_suppress_push();
    audit_chk("Continue#2", "Continue[]", "Hold[Continue[]]");
    audit_chk("Continue#2", "Break[]", "Hold[Break[]]");
    mth_msg_suppress_pop();
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("Catch#1", "Quiet[Check[Catch[1, 2, 3, 4], fired, Catch::argb]]", "fired");
    audit_chk("Catch#1", "Quiet[Catch[1, 2, 3, 4]]", "Catch[1, 2, 3, 4]");
    audit_chk("For#5", "Quiet[Check[For[s12i = 0, s12i < 3], fired, For::argb]]", "fired");
    audit_chk("For#5", "Quiet[For[s12i = 0, s12i < 3]]", "For[s12i = 0, s12i < 3]");
    audit_chk("If#1", "Quiet[Check[If[True], fired, If::argbu]]", "fired");
    audit_chk("If#1", "Quiet[If[True]]", "If[True]");
    audit_chk("Not#1", "Quiet[Check[Not[a, b], fired, Not::argx]]", "fired");
    audit_chk("Not#1", "Quiet[Not[a, b]]", "Not[a, b]");
    audit_chk("While#2", "Quiet[Check[While[1, 2, 3], fired, While::argb]]", "fired");
    audit_chk("While#2", "Quiet[While[1, 2, 3]]", "While[1, 2, 3]");
    audit_chk("AbortProtect#2", "Quiet[Check[AbortProtect[1, 2], fired, AbortProtect::argx]]", "fired");
    audit_chk("AbortProtect#2", "Quiet[AbortProtect[1, 2]]", "AbortProtect[1, 2]");
    audit_chk("Boole#3", "Quiet[Check[Boole[True, False], fired, Boole::argx]]", "fired");
    audit_chk("Boole#3", "Quiet[Boole[True, False]]", "Boole[True, False]");
    audit_chk("CheckAbort#1", "Quiet[Check[CheckAbort[x], fired, CheckAbort::argr]]", "fired");
    audit_chk("CheckAbort#1", "Quiet[CheckAbort[x]]", "CheckAbort[x]");
    audit_chk("Goto#1", "Quiet[Check[Goto[1, 2], fired, Goto::argx]]", "fired");
    audit_chk("Goto#1", "Quiet[Goto[1, 2]]", "Goto[1, 2]");
    audit_chk("Implies#1", "Quiet[Check[Implies[p], fired, Implies::argr]]", "fired");
    audit_chk("Implies#1", "Quiet[Implies[p]]", "Implies[p]");
    audit_chk("Monitor#1", "Quiet[Check[Monitor[x], fired, Monitor::argr]]", "fired");
    audit_chk("Monitor#1", "Quiet[Monitor[x]]", "Monitor[x]");
    audit_chk("Quiet#1", "Quiet[Check[Quiet[1/0, All, None, x], fired, Quiet::argb]]", "fired");
    audit_chk("TrueQ#1", "Quiet[Check[TrueQ[True, False], fired, TrueQ::argx]]", "fired");
    audit_chk("TrueQ#1", "Quiet[TrueQ[True, False]]", "TrueQ[True, False]");
    audit_chk("Dialog#1", "Quiet[Check[Dialog[Sequence[1, 2]], fired, Dialog::nonopt]]", "fired");
    audit_chk("Dialog#1", "Quiet[Dialog[Sequence[1, 2]]]", "Dialog[1, 2]");
    audit_chk("Quiet#1", "Quiet[Quiet[1/0, All, None, x]]", "Quiet[1/0, All, None, x]");
}

int main(void) {
    audit_init();
    TEST(test_for_start_throw);
    TEST(test_break_continue_in_argument);
    TEST(test_abort_in_nintegrate);
    TEST(test_numloop_exact_types);
    TEST(test_conditionalexpression_propagates);
    TEST(test_s14_continue_toplevel);
    TEST(test_wave3_arity_options);
    return audit_finish("audit_control_flow");
}

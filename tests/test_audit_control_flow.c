/* Builtin-audit regressions: control flow (review/issues/<Head>.md;
 * tests/test_audit.h). */
#include "test_audit.h"
#include "compile/compiled_function.h"   /* compiled_function_disassemble */

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

/* Wave 4 G8: control flow.
 *  - Abort[] is sticky, as an interrupt is: a callee that drops the sentinel (a numeric
 *    solver that saw a non-numeric sample) no longer swallows it -- the next evaluator
 *    poll re-raises it, and CheckAbort / the top level consume it;
 *  - CheckAbort inside AbortProtect absorbs the abort it catches;
 *  - Goto leaves Do / While / For / Scan for a Label around the loop;
 *  - a rejected Quit["a"] is inert and the CompoundExpression goes on;
 *  - Do and Table run a symbolic range of known length ({i, a, a + 2}). */
static void test_wave4_g8_control(void) {
    audit_chk("Abort#2", "Quiet[CheckAbort[FindMinimum[Abort[], {x, 2}], c]]", "c");
    audit_chk("Abort#3", "Quiet[CheckAbort[FindRoot[Abort[], {x, 1}], c]]", "c");
    audit_chk("Abort#3", "Quiet[CheckAbort[NSum[If[k > 3, Abort[], 1/k^2], {k, 1, Infinity}], c]]", "c");
    audit_chk("Abort#3", "Quiet[CheckAbort[NIntegrate[Abort[], {x, 0, 1}], c]]", "c");
    audit_chk("Abort#3", "Quiet[Check[CheckAbort[FindRoot[Abort[], {x, 1}], c], fired]]", "c");
    audit_chk("Abort#4", "Quiet[CheckAbort[NDSolve[{y'[t] == If[t > 0.5, Abort[], y[t]], y[0] == 1}, y, {t, 0, 1}], c]]", "c");
    audit_chk("Abort#4", "1 + 1", "2");                      /* nothing left pending */
    audit_run("w4$log = {}");
    audit_chk("AbortProtect#1",
              "CheckAbort[AbortProtect[CheckAbort[Abort[], AppendTo[w4$log, \"inner\"]];"
              " AppendTo[w4$log, \"cont\"]; Null], \"outer\"]", "Null");
    audit_chk("AbortProtect#1", "w4$log", "{\"inner\", \"cont\"}");
    audit_chk("AbortProtect#1", "CheckAbort[AbortProtect[Abort[]; 7], d]", "d");   /* still deferred */
    audit_chk("Label#1", "(Do[If[i == 3, Goto[w4$out]], {i, 10}]; \"no\"; Label[w4$out]; \"jumped\")", "\"jumped\"");
    audit_chk("Label#1", "w4$n = 0; (While[w4$n < 1000000, w4$n++; If[w4$n > 3, Goto[w4$out]]]; w4$n = -1;"
              " Label[w4$out]; w4$n)", "4");
    audit_chk("Label#1", "(For[w4$i = 1, w4$i < 10, w4$i++, If[w4$i == 2, Goto[w4$out]]]; w4$i = -1;"
              " Label[w4$out]; w4$i)", "2");
    audit_chk("Label#1", "(Scan[If[# == 2, Goto[w4$out]] &, {1, 2, 3}]; \"no\"; Label[w4$out]; \"scan\")", "\"scan\"");
    audit_chk("Quit#1", "Quiet[Quit[\"a\"]; 5]", "5");
    audit_chk("Quit#1", "Quiet[Check[Quit[\"a\"], fired, Quit::intm]]", "fired");
    audit_chk("Do#3", "w4$acc = {}; Do[AppendTo[w4$acc, i], {i, a, a + 2}]; w4$acc", "{a, 1 + a, 2 + a}");
    audit_chk("Do#3", "Do[Null, {i, a, a + 2}]", "Null");
    audit_chk("Do#3", "Table[i, {i, a, a + 2}]", "{a, 1 + a, 2 + a}");
    audit_chk("Do#3", "Table[i, {i, a, b}]", "Table[i, {i, a, b}]");
}

/* Wave 4 G8: ConditionalExpression and Boole under assumptions, Normal and D. */
static void test_wave4_g8_conditions(void) {
    audit_chk("ConditionalExpression#2", "Normal[ConditionalExpression[x, x > 0]]", "x");
    audit_chk("ConditionalExpression#3", "Simplify[ConditionalExpression[x, x > 0], x > 0]", "x");
    audit_chk("ConditionalExpression#3", "Refine[ConditionalExpression[x, x > 0], x > 0]", "x");
    audit_chk("ConditionalExpression#3", "Refine[ConditionalExpression[x, x > 0], x < 0]", "Undefined");
    audit_chk("ConditionalExpression#4", "D[ConditionalExpression[x^2, x > 0], x]", "ConditionalExpression[2 x, x > 0]");
    audit_chk("Boole#1", "Integrate[Boole[x < 1], {x, 0, 2}]", "1");
    audit_chk("Boole#1", "Integrate[x Boole[x > 1], {x, 0, 3}]", "4");
    audit_chk("Boole#2", "Refine[Boole[x > 0], x > 0]", "1");
    audit_chk("Boole#2", "Simplify[Boole[x > 0], x > 0]", "1");
    audit_chk("Boole#2", "Refine[Boole[x > 0], x < -1]", "0");
}

/* Wave 4 G8: Compile -- 0.^0. is Indeterminate (the call goes to the interpreter), True|False
 * is the Boolean type spec, CompilePrint names the managed types and their uses. */
static void test_wave4_g8_compile(void) {
    audit_chk("Compile#1", "Quiet[Compile[{x}, 0.^x][0.]]", "Indeterminate");
    audit_chk("Compile#1", "Quiet[Compile[{x}, x^0.][0.]]", "Indeterminate");
    audit_chk("Compile#1", "Quiet[Check[Compile[{x}, x^0.][0.], fired, Power::indet]]", "fired");
    audit_chk("Compile#1", "Compile[{x}, x^0.5][4.]", "2.0");
    audit_chk("CompileDiagnostics#1", "\"Compiled\" /. CompileDiagnostics[{{b, True | False}}, If[b, 1.0, 2.0]]", "True");
    audit_chk("CompileDiagnostics#1", "Compile[{{b, True | False}}, If[b, 1.0, 2.0]][False]", "2.0");
    audit_chk("CompilePrint#1", "StringContainsQ[ToString[CompilePrint[Compile[{x}, x^2, WorkingPrecision -> 30]]], \"unused\"]", "False");
    audit_chk("CompilePrint#1", "StringContainsQ[ToString[CompilePrint[Compile[{x}, x^2, WorkingPrecision -> 30]]], \"Unknown\"]", "False");
}

/* Wave 4 G8: Message[f::m::lang, ...] (parsed as one MessageName with three arguments),
 * ParallelMap over a rule, the built-in slash commands reserved in every mode. */
static void test_wave4_g8_misc_control(void) {
    audit_chk("Message#1", "FullForm[Hold[f::m::French]]", "Hold[MessageName[f, \"m\", \"French\"]]");
    audit_run("w4m::m = \"bad `1`\"");
    audit_chk("Message#1", "Quiet[Check[Message[w4m::m::French, 1], fired, w4m::m]]", "fired");
    audit_chk("ParallelMap#1", "ParallelMap[f, a -> b]", "f[a] -> f[b]");
    audit_chk("ParallelMap#1", "ParallelMap[f, {1, 2}]", "{f[1], f[2]}");
    audit_chk("SlashCommand#1", "Quiet[SlashCommand[\"help\", f]]", "$Failed");
    audit_chk("SlashCommand#1", "Quiet[Check[SlashCommand[\"help\", f], fired, SlashCommand::builtin]]", "fired");
    audit_chk("SlashCommand#1", "First[SlashCommand[]]", "\"/help\" -> \"list the commands and keys\"");
}

/* Wave 5 G8: control-flow findings -- argument counts, missing messages,
 * General::stop, Which/Unevaluated, Compile signatures and values. */
static void test_wave5_g8_control(void) {
    audit_chk("Do#5", "Quiet[Check[Do[Print[i], {i, w5x}], fired, Do::iterb]]", "fired");
    audit_chk("Do#5", "Quiet[Do[Print[i], {i, w5x}]]", "Do[Print[i], {i, w5x}]");
    audit_chk("FixedPoint#2", "Quiet[Check[FixedPoint[f], fired, FixedPoint::argtu]]", "fired");
    audit_chk("FixedPoint#2", "Quiet[FixedPoint[f]]", "FixedPoint[f]");
    audit_chk("FixedPoint#2", "Quiet[Check[FixedPoint[f, x, 2, 3], fired, FixedPoint::nonopt]]", "fired");
    audit_chk("FixedPoint#2", "Quiet[Check[FixedPoint[f, x, -1], fired, FixedPoint::intnm]]", "fired");
    audit_chk("FixedPoint#2", "Quiet[FixedPoint[f, x, -1]]", "FixedPoint[f, x, -1]");
    audit_chk("FixedPoint#2", "FixedPoint[#/2 &, 1.0, SameTest -> (Abs[#1 - #2] < 0.1 &)]", "0.0625");
    audit_chk("Throw#1", "Quiet[Check[Throw[], fired, Throw::argb]]", "fired");
    audit_chk("Throw#1", "Quiet[Throw[1, 2, 3, 4]]", "Throw[1, 2, 3, 4]");
    audit_chk("Throw#1", "Catch[Throw[1, a]; 2, a]", "1");
    audit_chk("Which#1", "Which[Unevaluated[1 > 0], 1]", "1");
    audit_chk("Which#1", "Which[Unevaluated[1 < 0], 1, True, 2]", "2");
    audit_chk("Check#2", "Quiet[Check[Check[w5x], fired, Check::argmu]]", "fired");
    audit_chk("Check#2", "Quiet[Check[Check[], fired, Check::argm]]", "fired");
    audit_chk("Check#2", "Quiet[Check[w5x]]", "Check[w5x]");
    audit_chk("Exit#1", "Quiet[Check[Exit[w5x], fired, Exit::intm]]", "fired");
    audit_chk("Exit#1", "Quiet[Exit[w5x]]", "Exit[w5x]");
    audit_chk("Quit#2", "Quiet[Check[Quit[1, 2], fired, Quit::argt]]", "fired");
    audit_chk("Quit#2", "Quiet[Quit[1, 2]; 5]", "5");
    audit_chk("ParallelEvaluate#1", "Quiet[Check[ParallelEvaluate[1, 2, 3], fired, ParallelEvaluate::argt]]", "fired");
    audit_chk("ParallelEvaluate#1", "Quiet[ParallelEvaluate[1, 2, 3]]", "ParallelEvaluate[1, 2, 3]");
    audit_chk("Message#2", "Quiet[Check[Message[3], fired, Message::name]]", "fired");
    audit_chk("Message#2", "Quiet[Check[Message[{f::m}], fired, Message::name]]", "fired");
    audit_chk("Message#3", "Quiet[Check[Message[], fired, Message::argm]]", "fired");
    audit_chk("Message#3", "Quiet[Message[]]", "Message[]");
    audit_chk("Quiet#2", "Quiet[Check[Quiet[1/0, Power], fired, Quiet::anmlist]]", "fired");
    audit_chk("Quiet#2", "Quiet[Quiet[1/0, Power]]", "Quiet[1/0, Power]");
    audit_chk("Quiet#2", "Quiet[1/0, {Power::infy}]", "ComplexInfinity");
    /* General::stop: the fourth and later printings are held back, Check still sees them */
    audit_run("w5m::m = \"bad `1`\"");
    /* General::stop follows the third printing; the later two are generated
     * (in $MessageList, seen by Check) but not printed */
    audit_chk("Message#4", "Block[{$MessageList = {}}, Do[Message[w5m::m, i], {i, 5}]; $MessageList]",
              "{w5m::m, w5m::m, w5m::m, General::stop, w5m::m, w5m::m}");
    audit_chk("Off#1", "Block[{$MessageList = {}}, Do[1/0, {5}]; Position[$MessageList, HoldForm[General::stop]]]", "{{4}}");
    audit_chk("Off#1", "Check[Do[1/0, {5}], caught]", "caught");
    audit_chk("Off#1", "Block[{$MessageList = {}}, Off[General::stop]; Do[1/0, {5}]; On[General::stop]; "
              "MemberQ[$MessageList, HoldForm[General::stop]]]", "False");
    audit_chk("Off#1", "Block[{$MessageList = {}}, 1/0; Count[$MessageList, HoldForm[Power::infy]]]", "1");
}

/* CompilePrint writes to stdout and returns Null, so its text is checked at
 * the C level: `src` must build a CompiledFunction whose disassembly holds
 * `needle`. */
static void audit_disasm_has(const char* id, const char* src, const char* needle) {
    Expr* e = parse_expression(src);
    Expr* r = e ? evaluate(e) : NULL;
    char* d = (r && r->type == EXPR_COMPILED) ? compiled_function_disassemble(r->data.compiled) : NULL;
    bool ok = d && strstr(d, needle);
    if (!ok) fprintf(stderr, "  [%s] disassembly lacks \"%s\":\n%s\n", id, needle, d ? d : "(none)");
    free(d);
    if (r) expr_free(r);
    if (e) expr_free(e);
    audit_chk(id, ok ? "True" : "False", "True");
}

static void test_wave5_g8_compile(void) {
    audit_chk("Compile#2", "Compile[{}, 3.][]", "3.0");
    audit_chk("CompileDiagnostics#4", "\"Compiled\" /. CompileDiagnostics[{}, 1.0]", "True");
    audit_chk("Compile#3", "Compile[{x}, x][Pi] === N[Pi]", "True");
    audit_chk("Compile#3", "Compile[{x}, x + 1][Sqrt[2]] === N[1 + Sqrt[2]]", "True");
    audit_chk("Compile#3", "Head[Compile[{x}, x][Pi]]", "Real");
    audit_chk("Compile#3", "Abs[Compile[{x}, Sin[x]][Pi]] < 10^-15", "True");
    audit_chk("Compile#3", "Compile[{x}, x + 1][w5y]", "1 + w5y");
    audit_chk("CompiledFunction#1", "Quiet[Check[Compile[{x}, x^2 + 1][1, 2], fired, CompiledFunction::cfct]]", "fired");
    audit_chk("CompiledFunction#1", "Head[Quiet[Compile[{x}, x^2 + 1][1, 2]]]", "CompiledFunction[{x}, x^2 + 1]");
    audit_chk("CompiledFunction#2", "Compile[{x}, x^2] === Compile[{x}, x^2]", "True");
    audit_chk("CompiledFunction#2", "Compile[{x}, x^2] === Compile[{x}, x^3]", "False");
    audit_chk("CompiledFunction#2", "Compile[{{x, _Integer}}, x^2] === Compile[{x}, x^2]", "False");
    audit_chk("CompiledFunction#2", "Length[Union[{Compile[{x}, x], Compile[{x}, x]}]]", "1");
    audit_chk("CompiledFunction#3", "Compile[{x}, x, CompilationTarget -> \"C\"][1.]", "1.0");
    audit_chk("CompiledFunction#3", "Compile[{x}, x, Parallelization -> True][1.]", "1.0");
    audit_chk("CompiledFunction#3", "Compile[{x}, x, {{y, _Real}}][1.]", "1.0");
    audit_chk("CompileDiagnostics#3", "Block[{w5g = 3.0}, \"Reason\" /. CompileDiagnostics[{{x, _Real}}, x + w5g]]",
              "\"symbol is not a declared argument or a local variable of the compiled code\"");
    audit_disasm_has("CompilePrint#2", "Compile[{{x, _Real}}, x + 1, WorkingPrecision -> 30]", "T2 = R0 + T1");
    audit_disasm_has("CompilePrint#2", "Compile[{{x, _Real}}, x + 1, WorkingPrecision -> 30]", "T1 = 1.0000");
    audit_disasm_has("CompilePrint#2", "Compile[{{x, _Real}}, Sin[x] + 1, WorkingPrecision -> 30]", "= Sin[R0]");
}

static void test_wave5_g8_misc_control(void) {
    audit_chk("SlashCommand#2", "Quiet[SlashCommand[\"abcdefghijabcdefghijabcdefghijabc\", f]]", "$Failed");
    audit_chk("SlashCommand#2", "StringContainsQ[Information[SlashCommand], \"at most 32 characters\"]", "True");
    audit_chk("CompileDiagnostics#2", "\"Compiled\" /. CompileDiagnostics[{{z, _Complex}}, Zeta[z]]", "True");
    audit_chk("CompileDiagnostics#2", "Abs[Compile[{{z, _Complex}}, Zeta[z]][2. + I] - N[Zeta[2 + I]]] < 10^-12", "True");
}

/* Wave 6 G1: the spec page's Zeta example (a complex Zeta compiles). */
static void test_wave6_g1_control(void) {
    audit_chk("CompileDiagnostics#2", "\"Compiled\" /. CompileDiagnostics[{{z, _Complex}}, Zeta[z]]", "True");
    audit_chk("CompileDiagnostics#2", "Abs[Compile[{{z, _Complex}}, Zeta[z]][2. + I] - N[Zeta[2 + I]]] < 10^-12", "True");
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
    TEST(test_wave4_g8_control);
    TEST(test_wave4_g8_conditions);
    TEST(test_wave4_g8_compile);
    TEST(test_wave4_g8_misc_control);
    TEST(test_wave5_g8_control);
    TEST(test_wave5_g8_compile);
    TEST(test_wave5_g8_misc_control);
    TEST(test_wave6_g1_control);
    return audit_finish("audit_control_flow");
}

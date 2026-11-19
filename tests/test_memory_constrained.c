/*
 * test_memory_constrained.c -- MemoryConstrained[expr, b, failexpr].
 *
 * Expected values are the reference language's (probed in its kernel), except
 * where an answer depends on byte counts -- every system charges a different
 * number of bytes per object -- in which case the test states the bound it
 * relies on rather than the reference's exact figure.  After every group the
 * accounting must be back at rest: no frame open, no trip pending, the limit
 * back to "no constraint".
 *
 * Side effects are observed through a log list (AppendTo) rather than Echo, so
 * the order of events is part of the checked value.
 */
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "memacct.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

static void assert_at_rest(const char* where) {
    if (cf_depth() != 0 || g_cf_tripped != -1 || g_mth_mem_limit != INT64_MAX) {
        fprintf(stderr, "FAIL (%s): accounting not at rest: depth %d, tripped %d, limit %lld\n",
                where, cf_depth(), g_cf_tripped, (long long)g_mth_mem_limit);
        abort();
    }
}

/* The call is left unevaluated and the named message is issued. */
static void assert_message(const char* input, const char* tag) {
    char src[512];
    snprintf(src, sizeof src, "Quiet[Check[%s, \"msg\", MemoryConstrained::%s]]", input, tag);
    assert_eval_eq(src, "\"msg\"", 0);
}

static void test_attributes_and_arity(void) {
    assert_eval_eq("Attributes[MemoryConstrained]", "{HoldAll, Protected}", 0);
    assert_eval_eq("StringQ[Information[MemoryConstrained, \"Usage\"]] || "
                   "StringLength[ToString[Information[MemoryConstrained]]] > 0", "True", 0);
    assert_message("MemoryConstrained[]", "argt");
    assert_message("MemoryConstrained[1]", "argtu");
    assert_message("MemoryConstrained[1, 2, 3, 4]", "argt");
    assert_eval_eq("Quiet[MemoryConstrained[]]", "MemoryConstrained[]", 0);
    assert_eval_eq("Quiet[MemoryConstrained[1, 2, 3, 4]]", "MemoryConstrained[1, 2, 3, 4]", 0);
    assert_at_rest("arity");
}

static void test_budget_argument(void) {
    /* Positive machine integer or Infinity; anything else: ipnfm, unevaluated. */
    assert_message("MemoryConstrained[1, -1]", "ipnfm");
    assert_message("MemoryConstrained[1, 0, \"f\"]", "ipnfm");
    assert_message("MemoryConstrained[1, 1.5]", "ipnfm");
    assert_message("MemoryConstrained[1, x]", "ipnfm");
    assert_message("MemoryConstrained[1, 10^20]", "ipnfm");
    assert_eval_eq("Quiet[MemoryConstrained[1, -1]]", "MemoryConstrained[1, -1]", 0);
    /* An unevaluated call never ran its body or its failexpr. */
    assert_eval_eq("bud = 0; Quiet[MemoryConstrained[bud = 1, -1, bud = 2]]; bud", "0", 0);
    /* The budget is evaluated (HoldAll holds it only until the call runs). */
    assert_eval_eq("bb = 10^7; MemoryConstrained[Length[Range[10]], bb]", "10", 0);
    assert_eval_eq("MemoryConstrained[Length[Range[10^6]], 2*bb]", "1000000", 0);
    /* Infinity imposes no constraint. */
    assert_eval_eq("MemoryConstrained[Length[Range[10^6]], Infinity]", "1000000", 0);
    assert_eval_eq("MemoryConstrained[Length[Range[10^6]], Infinity, \"f\"]", "1000000", 0);
    assert_eval_eq("MemoryConstrained[Length[Range[10^6]], DirectedInfinity[1]]", "1000000", 0);
    assert_at_rest("budget");
}

/* The documented examples. */
static void test_reference_examples(void) {
    assert_eval_eq("MemoryConstrained[Range[10^6], 10000]", "$Aborted", 0);
    assert_eval_eq("MemoryConstrained[Range[10], 10, \"too much memory\"]", "\"too much memory\"", 0);
    /* The loop runs until the memory for Range's result passes the limit: k
     * stops part-way (the reference, which charges more per element, stops
     * at 324; the exact value depends on bytes per element). */
    assert_eval_eq("Clear[k]; MemoryConstrained[Table[k = i^2; Range[k], {i, 100}], 10000]", "$Aborted", 0);
    assert_eval_eq("IntegerQ[k] && 25 <= k <= 1600", "True", 0);
    /* failexpr is evaluated only if the limit is reached. */
    assert_eval_eq("fl = {}; {MemoryConstrained[Sin[Pi], 100, AppendTo[fl, \"failure\"]], fl}", "{0, {}}", 0);
    assert_eval_eq("fl = {}; {MemoryConstrained[Range[100], 100, AppendTo[fl, \"failure\"]; \"failure\"], fl}",
                   "{\"failure\", {\"failure\"}}", 0);
    assert_eval_eq("MemoryConstrained[Range[10], 10^6]", "{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}", 0);
    assert_at_rest("examples");
}

static void test_net_accounting(void) {
    /* Memory allocated and freed again inside the body is not charged. */
    assert_eval_eq("MemoryConstrained[Do[Range[10^4], {1000}], 10^6]", "Null", 0);
    assert_eval_eq("MemoryConstrained[Do[RandomReal[1, 10^4], {1000}], 10^6]", "Null", 0);
    assert_eval_eq("MemoryConstrained[Do[x = Range[10^4]; Clear[x], {1000}], 10^6]", "Null", 0);
    assert_eval_eq("MemoryConstrained[Do[2^10000, {10^4}], 10^5]", "Null", 0);
    assert_eval_eq("MemoryConstrained[Do[N[Pi, 1000], {100}], 10^5]", "Null", 0);
    assert_eval_eq("MemoryConstrained[Do[StringJoin[Table[\"ab\", {1000}]], {1000}], 10^5]", "Null", 0);
    assert_eval_eq("MemoryConstrained[Do[Module[{t = Table[i, {i, 1000}]}, Total[t]], {200}], 10^6]", "Null", 0);
    /* ... but memory that is kept is. */
    assert_eval_eq("MemoryConstrained[Table[RandomReal[1, 10^3], {100}]; 1, 10^5]", "$Aborted", 0);
    assert_eval_eq("MemoryConstrained[Length[Table[f[i], {i, 10^5}]], 10^5]", "$Aborted", 0);
    assert_eval_eq("MemoryConstrained[Length[Table[f[i], {i, 10^5}]], 10^8]", "100000", 0);
    /* A block that is requested counts even if it is freed before any poll. */
    assert_eval_eq("MemoryConstrained[StringJoin[Table[\"ab\", {10^5}]]; 1, 10^5]", "$Aborted", 0);
    assert_at_rest("net");
}

static void test_kinds_of_memory(void) {
    /* Big-number limbs (through GMP's allocator). */
    assert_eval_eq("MemoryConstrained[IntegerLength[2^(10^7)], 10^4, \"big\"]", "\"big\"", 0);
    assert_eval_eq("MemoryConstrained[IntegerLength[2^(10^5)], 10^6]", "30103", 0);
    /* MPFR numbers. */
    assert_eval_eq("MemoryConstrained[N[Pi, 10^5]; 1, 10^4, \"mpfr\"]", "\"mpfr\"", 0);
    /* Strings. */
    assert_eval_eq("MemoryConstrained[StringLength[StringJoin[Table[\"abc\", {10^5}]]], 10^4, \"str\"]",
                   "\"str\"", 0);
    /* Boxed expressions, no packing. */
    assert_eval_eq("MemoryConstrained[Length[Table[{i, x}, {i, 10^5}]], 10^5, \"boxed\"]", "\"boxed\"", 0);
    /* Packed arrays of each kind. */
    assert_eval_eq("MemoryConstrained[RandomInteger[10, 10^6], 10^4, \"ri\"]", "\"ri\"", 0);
    assert_eval_eq("MemoryConstrained[RandomComplex[1 + I, 10^6], 10^4, \"rc\"]", "\"rc\"", 0);
    assert_eval_eq("MemoryConstrained[ConstantArray[0., {1000, 1000}], 10^4, \"ca\"]", "\"ca\"", 0);
    assert_eval_eq("MemoryConstrained[Range[0., 1., 10^-6], 10^4, \"rr\"]", "\"rr\"", 0);
    assert_eval_eq("MemoryConstrained[Table[Sin[x], {x, 0., 1., 10^-6}], 10^4, \"tbl\"]", "\"tbl\"", 0);
    /* A request far beyond the machine is refused before it is made. */
    assert_eval_eq("MemoryConstrained[RandomReal[1, 10^10], 10^3, \"huge\"]", "\"huge\"", 0);
    assert_eval_eq("MemoryConstrained[ConstantArray[1, {10^5, 10^5}], 10^3, \"huge\"]", "\"huge\"", 0);
    assert_eval_eq("MemoryConstrained[RandomInteger[1, {10^5, 10^5}], 10^3, \"huge\"]", "\"huge\"", 0);
    assert_at_rest("kinds");
}

static void test_abort_protect(void) {
    /* A nested AbortProtect blocks the abort until it finishes. */
    assert_eval_eq("lg = {}; r = MemoryConstrained[AbortProtect[x = Range[10^4]; AppendTo[lg, \"eval\"]], "
                   "10^2, AppendTo[lg, \"constraint\"]; \"constraint\"]; {r, lg}",
                   "{\"constraint\", {\"eval\", \"constraint\"}}", 0);
    /* An enclosing AbortProtect is overridden ... */
    assert_eval_eq("lg = {}; r = AbortProtect[MemoryConstrained[RandomReal[1, 10^5]; AppendTo[lg, \"eval\"], "
                   "10^4, AppendTo[lg, \"constraint\"]; \"constraint\"]]; {r, lg}",
                   "{\"constraint\", {\"constraint\"}}", 0);
    /* ... while the code around MemoryConstrained runs to completion. */
    assert_eval_eq("lg = {}; r = AbortProtect[MemoryConstrained[RandomReal[1, 10^5]; AppendTo[lg, \"eval\"], "
                   "10^4, AppendTo[lg, \"constraint\"]]; AppendTo[lg, \"past\"]; \"past\"]; {r, lg}",
                   "{\"past\", {\"constraint\", \"past\"}}", 0);
    /* Nested AbortProtects inside: deferred to the outermost of them. */
    assert_eval_eq("lg = {}; MemoryConstrained[AbortProtect[AbortProtect[x = Range[10^4]; AppendTo[lg, 1]]; "
                   "AppendTo[lg, 2]]; AppendTo[lg, 3], 10^2, lg]", "{1, 2}", 0);
    assert_at_rest("abortprotect");
}

static void test_check_abort(void) {
    /* CheckAbort catches the abort and the body continues. */
    assert_eval_eq("lg = {}; MemoryConstrained[CheckAbort[Range[10^6]; AppendTo[lg, \"first\"], "
                   "AppendTo[lg, \"second\"], PropagateAborts -> False]; AppendTo[lg, \"finished\"]; lg, 10^4]",
                   "{\"second\", \"finished\"}", 0);
    /* The default, Automatic, runs failexpr and then hands a constraint's
     * abort on to its owner (the reference); a plain Abort[] would stop. */
    assert_eval_eq("lg = {}; r = MemoryConstrained[CheckAbort[Range[10^6]; AppendTo[lg, \"first\"], "
                   "AppendTo[lg, \"second\"]]; AppendTo[lg, \"finished\"]; lg, 10^4, \"mc\"]; {r, lg}",
                   "{\"mc\", {\"second\"}}", 0);
    assert_eval_eq("lg = {}; r = MemoryConstrained[CheckAbort[Range[10^6], AppendTo[lg, \"second\"], "
                   "PropagateAborts -> Automatic]; \"fin\", 10^4, \"mc\"]; {r, lg}",
                   "{\"mc\", {\"second\"}}", 0);
    /* PropagateAborts -> True passes it on after failexpr. */
    assert_eval_eq("lg = {}; r = MemoryConstrained[CheckAbort[Range[10^6], AppendTo[lg, \"second\"], "
                   "PropagateAborts -> True]; AppendTo[lg, \"finished\"], 10^4, \"mc\"]; {r, lg}",
                   "{\"mc\", {\"second\"}}", 0);
    /* The first argument of CheckAbort is constrained, the second is not: it
     * runs to completion, and then the constraint takes effect. */
    assert_eval_eq("lg = {}; r = MemoryConstrained[CheckAbort[RandomReal[1, 10^10]; AppendTo[lg, \"eval\"], "
                   "RandomReal[1, 10^5]; AppendTo[lg, \"constraint\"]], 10^3]; {r, lg}",
                   "{$Aborted, {\"constraint\"}}", 0);
    /* A CheckAbort outside the MemoryConstrained sees no abort. */
    assert_eval_eq("CheckAbort[MemoryConstrained[Range[10^6], 10^4, \"mc\"], \"ca\"]", "\"mc\"", 0);
    assert_eval_eq("CheckAbort[MemoryConstrained[Range[10^6], 10^4], \"ca\"]", "$Aborted", 0);
    assert_at_rest("checkabort");
}

static void test_aborts_that_are_not_ours(void) {
    /* A plain Abort[] is not caught by MemoryConstrained. */
    assert_eval_eq("CheckAbort[MemoryConstrained[Abort[], 10^6, \"f\"], \"outer\"]", "\"outer\"", 0);
    assert_eval_eq("lg = {}; CheckAbort[MemoryConstrained[Abort[], 10^6, AppendTo[lg, 1]], lg]", "{}", 0);
    /* Throw passes through. */
    assert_eval_eq("Catch[MemoryConstrained[Throw[3], 10^6, \"f\"]]", "3", 0);
    assert_eval_eq("Catch[MemoryConstrained[Throw[3, t], 10^6, \"f\"], t]", "3", 0);
    assert_at_rest("not ours");
}

static void test_nesting(void) {
    /* The frame whose budget broke answers. */
    assert_eval_eq("MemoryConstrained[MemoryConstrained[Range[10^5], 10^9, \"inner\"], 10^4, \"outer\"]",
                   "\"outer\"", 0);
    assert_eval_eq("MemoryConstrained[MemoryConstrained[Range[10^5], 10^4, \"inner\"], 10^9, \"outer\"]",
                   "\"inner\"", 0);
    assert_eval_eq("MemoryConstrained[{MemoryConstrained[Range[10^5], 10^4, \"inner\"], 5}, 10^9, \"outer\"]",
                   "{\"inner\", 5}", 0);
    /* failexpr runs under the enclosing budget. */
    assert_eval_eq("MemoryConstrained[MemoryConstrained[Range[10^5], 10^4, Range[10^6]], 10^5, \"outer\"]",
                   "\"outer\"", 0);
    /* Three deep, middle frame broken. */
    assert_eval_eq("MemoryConstrained[MemoryConstrained[MemoryConstrained[Range[10^5], 10^9, \"a\"], "
                   "10^4, \"b\"], 10^9, \"c\"]", "\"b\"", 0);
    /* Sequential constraints are independent. */
    assert_eval_eq("{MemoryConstrained[Range[10^6], 10^4, \"a\"], Length[MemoryConstrained[Range[10^5], 10^7, \"b\"]]}",
                   "{\"a\", 100000}", 0);
    assert_at_rest("nesting");
}

static void test_with_time_constrained(void) {
    assert_eval_eq("TimeConstrained[MemoryConstrained[While[True], 10^6, \"mem\"], 0.3, \"time\"]", "\"time\"", 0);
    assert_at_rest("tc outside");
    assert_eval_eq("MemoryConstrained[TimeConstrained[Range[10^6], 5, \"time\"], 10^4, \"mem\"]", "\"mem\"", 0);
    assert_eval_eq("MemoryConstrained[TimeConstrained[Length[Range[10^5]], 5, \"time\"], 10^7, \"mem\"]",
                   "100000", 0);
    /* A timeout unwinding a MemoryConstrained frame leaves no limit behind. */
    assert_eval_eq("TimeConstrained[MemoryConstrained[MemoryConstrained[While[True], 10^6], 10^7], 0.2, \"t\"]; "
                   "Length[Range[10^6]]", "1000000", 0);
    assert_at_rest("tc inside");
}

static void test_state_restored(void) {
    assert_eval_eq("MemoryConstrained[Range[10^6], 10^4, \"x\"]; Length[Range[10^7]]", "10000000", 0);
    assert_eval_eq("Quiet[MemoryConstrained[1, -1]]; Length[Range[10^6]]", "1000000", 0);
    assert_eval_eq("CheckAbort[MemoryConstrained[Abort[], 10^4], 0]; Length[Range[10^6]]", "1000000", 0);
    /* Ordinary aborts still behave as before after a constraint fired. */
    assert_eval_eq("MemoryConstrained[Range[10^6], 10^4]; CheckAbort[Abort[], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("MemoryConstrained[Range[10^6], 10^4]; AbortProtect[1 + 1]", "2", 0);
    assert_at_rest("restored");
}

static void test_compile(void) {
    assert_eval_eq("mcf = Compile[{{n, _Integer}}, MemoryConstrained[Total[Range[n]], 10^4, -1]]; "
                   "{mcf[10], mcf[10^6], mcf[100]}", "{55, -1, 5050}", 0);
    assert_eval_eq("Lookup[CompileDiagnostics[{{n, _Integer}}, MemoryConstrained[Total[Range[n]], 10^4, -1]], "
                   "\"Compiled\"]", "True", 0);
    assert_eval_eq("Lookup[CompileDiagnostics[{{x, _Real}}, MemoryConstrained[Sin[x], 10^4, 0.]], \"Compiled\"]",
                   "True", 0);
    assert_eval_eq("Lookup[CompileDiagnostics[{{v, _Real, 1}}, MemoryConstrained[Total[v^2], 10^4, 0.]], "
                   "\"Compiled\"]", "True", 0);
    assert_eval_eq("Lookup[CompileDiagnostics[{{v, _Real, 1}}, MemoryConstrained[v^2, 10^4]], \"Compiled\"]",
                   "True", 0);
    /* Two-argument form: an overrun answers $Aborted (through the interpreter). */
    assert_eval_eq("mcg = Compile[{{n, _Integer}}, MemoryConstrained[Total[Range[n]], 10^4]]; {mcg[10], mcg[10^6]}",
                   "{55, $Aborted}", 0);
    /* Nested compiled frames: the broken budget's failexpr. */
    assert_eval_eq("mch = Compile[{{n, _Integer}}, MemoryConstrained[MemoryConstrained[Total[Range[n]], 10^9, -1], "
                   "10^4, -2]]; {mch[10], mch[10^6]}", "{55, -2}", 0);
    assert_eval_eq("mci = Compile[{{n, _Integer}}, MemoryConstrained[MemoryConstrained[Total[Range[n]], 10^4, -1], "
                   "10^9, -2]]; {mci[10], mci[10^6]}", "{55, -1}", 0);
    /* An interpreted constraint around a compiled call. */
    assert_eval_eq("MemoryConstrained[mcf[10^6], 10^9]", "-1", 0);
    assert_eval_eq("mcj = Compile[{{n, _Integer}}, Total[Range[n]]]; MemoryConstrained[mcj[10^6], 10^4, \"outer\"]",
                   "\"outer\"", 0);
    assert_eval_eq("MemoryConstrained[mcj[10^6], 10^8, \"outer\"]", "500000500000", 0);
    /* Infinity compiles to the bare body; a runtime budget comes from an argument. */
    assert_eval_eq("mck = Compile[{{n, _Integer}}, MemoryConstrained[Total[Range[n]], Infinity, -1]]; mck[10^6]",
                   "500000500000", 0);
    assert_eval_eq("mcb = Compile[{{n, _Integer}, {b, _Integer}}, MemoryConstrained[Total[Range[n]], b, -1]]; "
                   "{mcb[10^5, 10^4], mcb[10^5, 10^7]}", "{-1, 5000050000}", 0);
    /* A non-positive runtime budget goes to the interpreter, which complains. */
    assert_eval_eq("Quiet[mcb[10, 0]]", "MemoryConstrained[Total[Range[10]], 0, -1]", 0);
    /* Array growth (Table with an unknown length) is charged too. */
    assert_eval_eq("mcw = Compile[{{n, _Integer}}, MemoryConstrained[Total[Table[i, {i, n}]], 10^4, -1]]; "
                   "{mcw[10], mcw[10^5]}", "{55, -1}", 0);
    assert_at_rest("compile");
}

static void test_autocompile(void) {
    /* Table / Sum auto-compile their bodies; the constraint holds through them. */
    assert_eval_eq("Table[MemoryConstrained[Sin[x], 10^6, 0.], {x, 0., 1., 0.5}]",
                   "{0.0, 0.479426, 0.841471}", 0);
    assert_eval_eq("MemoryConstrained[Total[Table[Sin[x], {x, 0., 1., 10^-3}]] > 0, 10^5, \"tbl\"]", "True", 0);
    assert_eval_eq("MemoryConstrained[Table[Sin[x] + x^2, {x, 0., 1., 10^-6}], 10^5, \"tbl\"]", "\"tbl\"", 0);
    assert_eval_eq("MemoryConstrained[Abs[Sum[Sin[1. k], {k, 1, 10^5}]] < 2, 10^5, \"sum\"]", "True", 0);
    assert_at_rest("autocompile");
}

int main(void) {
    symtab_init();
    core_init();
    assert_at_rest("start");
    TEST(test_attributes_and_arity);
    TEST(test_budget_argument);
    TEST(test_reference_examples);
    TEST(test_net_accounting);
    TEST(test_kinds_of_memory);
    TEST(test_abort_protect);
    TEST(test_check_abort);
    TEST(test_aborts_that_are_not_ours);
    TEST(test_nesting);
    TEST(test_with_time_constrained);
    TEST(test_state_restored);
    TEST(test_compile);
    TEST(test_autocompile);
    printf("All MemoryConstrained tests passed.\n");
    return 0;
}

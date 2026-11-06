/*
 * test_tc_binding_unwind.c -- a TimeConstrained timeout must not leak the
 * temporary value of an iterator or a numerical-method variable.
 *
 * Table / Do / Sum / Product, the N-functions (NIntegrate, NSum, FindRoot,
 * NMinimize, NDSolve, ND, NLimit, NSeries, NResidue) and the plotters bind
 * their variable for the duration of the call.  A timeout siglongjmps past the
 * C frame that would put the variable back, so the binding has to live on the
 * dynamic-scope stack that the TimeConstrained unwind drains (modular.c,
 * mth_bind_push / mth_bind_pop).  Before that, `k` was left at whatever value
 * the loop had reached (k = 41006 after an interrupted Do) for the rest of the
 * session.
 *
 * Every case asserts both halves: that the timeout really fired (the failexpr
 * `ab` is returned) and that the variable is back -- unset, or its previous
 * value.
 */
#include "test_utils.h"
#include "symtab.h"
#include "core.h"

#include <stdio.h>
#include <string.h>

/* A body that is slow per step, so a 10^6-step loop is certain to time out. */
#define SLOW(v) "Expand[(1 + " v ")^200]"

static void test_iterators_unset(void) {
    assert_eval_eq("{TimeConstrained[Do[" SLOW("k") ", {k, 10^6}], 0.1, ab], k}",
                   "{ab, k}", 0);
    assert_eval_eq("{TimeConstrained[Table[" SLOW("j") ", {j, 10^6}], 0.1, ab], j}",
                   "{ab, j}", 0);
    assert_eval_eq("{TimeConstrained[Sum[" SLOW("s") ", {s, 10^6}], 0.1, ab], s}",
                   "{ab, s}", 0);
    assert_eval_eq("{TimeConstrained[Product[" SLOW("p") ", {p, 10^6}], 0.1, ab], p}",
                   "{ab, p}", 0);
    /* Several iterators, list iterators, real iterators. */
    assert_eval_eq("{TimeConstrained[Table[" SLOW("a1 b1") ", {a1, 1000}, {b1, 1000}], 0.1, ab], a1, b1}",
                   "{ab, a1, b1}", 0);
    assert_eval_eq("{TimeConstrained[Do[" SLOW("c1") ", {c1, Range[10^5]}], 0.1, ab], c1}",
                   "{ab, c1}", 0);
    assert_eval_eq("{TimeConstrained[Do[" SLOW("r1") ", {r1, 0.5, 10^6}], 0.1, ab], r1}",
                   "{ab, r1}", 0);
}

static void test_iterators_restore_previous(void) {
    assert_eval_eq("k2 = 7; {TimeConstrained[Do[" SLOW("k2") ", {k2, 10^6}], 0.1, ab], k2}",
                   "{ab, 7}", 0);
    assert_eval_eq("j2 = \"old\"; {TimeConstrained[Table[" SLOW("j2") ", {j2, 10^6}], 0.1, ab], j2}",
                   "{ab, \"old\"}", 0);
    /* An indexed iterator v[1] is a DownValue, restored the same way. */
    assert_eval_eq("vv[1] = 5; {TimeConstrained[Do[" SLOW("vv[1]") ", {vv[1], 10^6}], 0.1, ab], vv[1]}",
                   "{ab, 5}", 0);
    /* The session is healthy afterwards: the same symbol iterates again and
     * its outer value is still put back on a normal exit. */
    assert_eval_eq("Table[k2^2, {k2, 3}]", "{1, 4, 9}", 0);
    assert_eval_eq("Do[Null, {k2, 5}]; k2", "7", 0);
    assert_eval_eq("Sum[j2, {j2, 4}]", "10", 0);
    assert_eval_eq("j2", "\"old\"", 0);
}

static void test_block_interplay(void) {
    /* Iterator over a Block local, interrupted: the iterator frame is newer,
     * so it is drained first, then the Block's -- x gets its global value. */
    assert_eval_eq("x3 = 1; {TimeConstrained[Block[{x3 = 2}, Table[" SLOW("x3") ", {x3, 10^6}]], 0.1, ab], x3}",
                   "{ab, 1}", 0);
    /* Block inside the loop body, over the same symbol. */
    assert_eval_eq("y3 = 1; {TimeConstrained[Do[Block[{y3 = 0}, " SLOW("y3") "], {y3, 10^6}], 0.1, ab], y3}",
                   "{ab, 1}", 0);
    /* Interleaved: iterator, Block, iterator. */
    assert_eval_eq("z3 = 1; w3 = 2; {TimeConstrained[Do[Block[{w3}, Do[" SLOW("z3 w3") ", {w3, 10^6}]], {z3, 10^6}], 0.1, ab], z3, w3}",
                   "{ab, 1, 2}", 0);
    /* Normal (uninterrupted) nesting keeps working. */
    assert_eval_eq("Block[{x3 = 2}, {Table[x3, {x3, 3}], x3}]", "{{1, 2, 3}, 2}", 0);
    assert_eval_eq("Table[Table[q3, {q3, i3}], {i3, 3}]", "{{1}, {1, 2}, {1, 2, 3}}", 0);
    assert_eval_eq("Table[q3, {q3, 2}, {q3, 2}]", "{{1, 2}, {1, 2}}", 0);
}

static void test_nested_timeconstrained(void) {
    /* An inner timeout drains only what it jumped over: the outer Table's
     * binding of i4 stays in force and is put back normally afterwards. */
    assert_eval_eq("i4 = 9; {Table[{i4, TimeConstrained[Do[" SLOW("m4") ", {m4, 10^6}], 0.05, ab]}, {i4, 2}], i4, m4}",
                   "{{{1, ab}, {2, ab}}, 9, m4}", 0);
    assert_eval_eq("{TimeConstrained[Table[TimeConstrained[Do[" SLOW("n4") ", {n4, 10^6}], 0.02, in], {o4, 10^6}], 0.2, out], n4, o4}",
                   "{out, n4, o4}", 0);
}

static void test_numerical_methods(void) {
    assert_eval_eq("{TimeConstrained[NSum[1/(m^2 + Sin[m]), {m, 1, Infinity}, WorkingPrecision -> 60], 0.1, ab], m}",
                   "{ab, m}", 0);
    /* WorkingPrecision -> 400 keeps the integral well past the budget (at 40
     * digits it now finishes in ~0.2 s, inside it). */
    assert_eval_eq("{TimeConstrained[NIntegrate[Sin[x^2] Exp[-x/1000] BesselJ[3, x], {x, 0, Infinity}, WorkingPrecision -> 400], 0.2, ab], x}",
                   "{ab, x}", 0);
    assert_eval_eq("xr = 3; {TimeConstrained[NIntegrate[Sin[xr^2] Exp[-xr/1000] BesselJ[3, xr], {xr, 0, Infinity}, WorkingPrecision -> 400], 0.2, ab], xr}",
                   "{ab, 3}", 0);
}

/* A timeout that lands while MPFR recomputes a cached constant (pi, log 2,
 * Euler's gamma) at a higher precision must not leave the cache poisoned:
 * N[Pi, 30] was Indeterminate for the rest of the session after a timed-out
 * N[Pi, 2000000], and a later Levin rule built its nodes from that NaN pi
 * (BesselJ then aborted the kernel on a negative MPFR precision). */
static void test_mpfr_constant_cache_after_timeout(void) {
    assert_eval_eq("TimeConstrained[N[Pi, 2000000], 0.05, ab]", "ab", 0);
    assert_eval_eq("N[Pi, 30]", "3.14159265358979323846264338328", 0);
    assert_eval_eq("TimeConstrained[N[EulerGamma, 300000], 0.05, ab]", "ab", 0);
    assert_eval_eq("N[EulerGamma, 30]", "0.577215664901532860606512090082", 0);
    assert_eval_eq("TimeConstrained[N[Log[2], 3000000], 0.05, ab]", "ab", 0);
    assert_eval_eq("N[Log[2], 30]", "0.693147180559945309417232121458", 0);
    /* a NaN argument of the arbitrary-precision BesselJ is a NaN, not a crash */
    assert_eval_eq("TimeConstrained[N[Pi, 2000000], 0.05, ab]", "ab", 0);
    assert_eval_eq("Abs[N[BesselJ[3, Pi/7], 30] - 185968813646768863267810859058/10^32] < 10^-31", "True", 0);
}

/* Each N-function, driven slowly enough to time out (a body that is slow to
 * evaluate at every sample), must leave its variable untouched. */
static void test_numerical_methods_slow_body(void) {
    const char* rows[][2] = {
        { "FindRoot[slowf[u] - 2, {u, 1}]", "u" },
        { "NMinimize[slowf[u]^2, u]", "u" },
        { "FindMinimum[slowf[u]^2, {u, 1}]", "u" },
        { "NDSolve[{yy'[u] == slowf[u], yy[0] == 0}, yy, {u, 0, 1}]", "u" },
        { "ND[slowf[u], u, 1]", "u" },
        { "NLimit[slowf[u], u -> 0]", "u" },
        { "NSeries[slowf[u], {u, 0, 4}]", "u" },
        { "NResidue[slowf[u]/u, {u, 0}]", "u" },
        { "NSum[slowf[u], {u, 1, Infinity}]", "u" },
        { "NIntegrate[slowf[u], {u, 0, 1}]", "u" },
        { "Plot[slowf[u], {u, 0, 1}]", "u" },
    };
    /* slowf is numeric-only, so the methods sample it; every sample costs
     * hundreds of Expands of a polynomial in the sample value (no memo hit). */
    assert_eval_eq("slowf[t_?NumericQ] := (Do[Expand[(t + i5 + tt)^40], {i5, 300}]; Sin[t]); Clear[u]; Null", "Null", 0);
    char src[512];
    for (size_t i = 0; i < sizeof rows / sizeof rows[0]; i++) {
        snprintf(src, sizeof src, "{TimeConstrained[%s, 0.1, ab], %s}", rows[i][0], rows[i][1]);
        printf("  %s\n", rows[i][0]);
        assert_eval_eq(src, "{ab, u}", 0);
    }
    assert_eval_eq("{Table[u, {u, 2}], i5}", "{{1, 2}, i5}", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_iterators_unset);
    TEST(test_iterators_restore_previous);
    TEST(test_block_interplay);
    TEST(test_nested_timeconstrained);
    TEST(test_numerical_methods);
    TEST(test_numerical_methods_slow_body);
    TEST(test_mpfr_constant_cache_after_timeout);
    printf("All tc_binding_unwind tests passed.\n");
    return 0;
}

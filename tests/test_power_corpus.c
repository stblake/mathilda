/*
 * test_power_corpus.c -- Regression tests for Power algorithmic fixes.
 *
 * Bug 3 -- (I*Pi)^(1/2) -> (-1)^(1/4) Sqrt[Pi]
 *
 *   builtin_power now distributes a rational exponent over a Times-base
 *   when the base contains a pure-imaginary factor (Complex[0, k]) and
 *   every other factor is known strictly positive (Pi, E, positive
 *   integer/rational). The branch falls through unchanged when the
 *   factor list is mixed (Sqrt[I x] etc.) so we don't cross a branch
 *   cut for unknown-sign inputs.
 *
 * Bug 4 -- Power[2, -5/3] -> 1/(2 * 2^(2/3))
 *
 *   The integer-part extraction in builtin_power had a guard `p >= q`
 *   which excluded negative `p`. We relax to `|p| >= q` so the C99
 *   truncating divide produces a_int = -1, b_rem = -2 for p=-5, q=3,
 *   yielding the canonical Times[Rational[1,2], Power[2, -2/3]].
 *
 *   The Times-side canonicalizer (src/times.c radical-canonical pass)
 *   was simultaneously updated to NOT pull factors of `b` from den
 *   into the exponent past the proper-fraction boundary; otherwise the
 *   1/2 * 2^(-2/3) form would round-trip back into Power[2, -5/3] and
 *   loop. This test pins the new fixed-point.
 *
 * Bug 5 -- Power[0, 0] -> Indeterminate
 *
 *   Previously returned 1 because Power[_, 0] -> 1 fired before the
 *   0^anything dispatch. Now 0^0 is intercepted first and emits the
 *   `Power::indet` warning.
 */

#include "test_utils.h"
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"

#include "mth_machine_pow.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* Bug 3: pure-imaginary times positive base, rational exponent. */
/* A machine complex power that overflows a double is Overflow[], not the
 * Indeterminate an inf*0 inside binary powering would give (NIntegrate of a
 * 10^1426-sized contour integrand relies on it). */
static void test_complex_power_overflow(void) {
    assert_eval_eq("(1. + 0.9 I)^20000", "Overflow[]", 0);
    assert_eval_eq("Head[(1. + 0.9 I)^20]", "Complex", 0);
}

static void test_bug3_imag_times_positive(void) {
    assert_eval_eq("(I*Pi)^(1/2)",
                   "(-1)^(1/4) Sqrt[Pi]", 0);
    /* Sqrt[I*Pi] = (I*Pi)^(1/2). */
    assert_eval_eq("Sqrt[I*Pi]",
                   "(-1)^(1/4) Sqrt[Pi]", 0);
    /* Multiple positive factors. */
    assert_eval_eq("(I*Pi*2)^(1/2)",
                   "(-1)^(1/4) Sqrt[2] Sqrt[Pi]", 0);
    /* -I (negative imaginary). Power[-1, -1/4] canonicalises (floor of the
     * exponent into [0,1)) to -(-1)^(3/4), matching Mathematica:
     * Sqrt[Pi]/(-1)^(1/4) == (-1)^(-1/4) Sqrt[Pi] == -(-1)^(3/4) Sqrt[Pi]. */
    assert_eval_eq("(-I*Pi)^(1/2)",
                   "-(-1)^(3/4) Sqrt[Pi]", 0);
    /* Cube root case. */
    assert_eval_eq("(I*Pi)^(1/3)",
                   "(-1)^(1/6) Pi^(1/3)", 0);
    /* Negative branch: x of unknown sign -- must NOT distribute. */
    assert_eval_eq("Sqrt[I*x]",
                   "Sqrt[I x]", 0);
    /* No imaginary factor: leave as-is. */
    assert_eval_eq("Sqrt[2*Pi]",
                   "Sqrt[2 Pi]", 0);
}

/* Bug 4: integer-part extraction for negative rational exponents. */
static void test_bug4_negative_integer_part(void) {
    /* Print form: 1/2/2^(2/3) (Times-print of Rational[1,2] * Power[2,-2/3]). */
    assert_eval_eq("Power[2, -5/3]",       "1/(2 2^(2/3))", 0);
    /* FullForm pinpoints the Times[Rational[1,2], Power[2, Rational[-2,3]]]. */
    assert_eval_eq("FullForm[Power[2, -5/3]]",
                   "Times[Rational[1, 2], Power[2, Rational[-2, 3]]]", 0);
    /* Negative-half exponent: 2^(-3/2) = 1/2 * 2^(-1/2). */
    assert_eval_eq("Power[2, -3/2]", "1/(2 Sqrt[2])", 0);
    /* Sqrt[2]/4 must canonicalise to the same proper-fraction form. */
    assert_eval_eq("Sqrt[2]/4", "1/(2 Sqrt[2])", 0);
    /* Positive side still works (regression). */
    assert_eval_eq("Power[2, 5/3]", "2 2^(2/3)", 0);
    /* Proper-fraction stays as-is. */
    assert_eval_eq("Power[2, -1/3]", "1/2^(1/3)", 0);
    /* Pulling factors from numerator coefficient still works. */
    assert_eval_eq("8 * Power[2, -1/3]", "4 2^(2/3)", 0);
}

/* Bug 5: 0^0 is Indeterminate. */
static void test_bug5_zero_to_zero(void) {
    /* 0^0 must produce the symbol Indeterminate. */
    assert_eval_eq("Power[0, 0]", "Indeterminate", 0);
    assert_eval_eq("Power[0., 0]", "Indeterminate", 0);
    assert_eval_eq("Power[0, 0.]", "Indeterminate", 0);
    /* 0^positive integer still 0. */
    assert_eval_eq("Power[0, 3]", "0", 0);
    /* Positive base ^ 0 still 1 (regression). */
    assert_eval_eq("Power[5, 0]", "1", 0);
    /* 0^negative still ComplexInfinity (regression). */
    assert_eval_eq("Power[0, -1]", "ComplexInfinity", 0);
}

/* Rational-base + rational-exponent canonicalisation (2026-05-12).
 * Powers of Rational[n, d] with non-integer exponent now distribute
 * to Power[n, p/q] * Power[d, -p/q] when at least one piece has a
 * perfect q-th-power extraction (or |n| == 1). Cases without any
 * extraction (Sqrt[2/3], (4/9)^(2/3)) stay unevaluated. */
static void test_rational_base_rational_exp(void) {
    /* Numerator-only trigger (|n| == 1). */
    assert_eval_eq("Power[1/54, 2/3]", "1/(9 2^(2/3))", 0);
    assert_eval_eq("Power[1/8, 2/3]",  "1/4", 0);
    assert_eval_eq("Power[1/4, 1/2]",  "1/2", 0);

    /* Both numerator and denominator extract. */
    assert_eval_eq("Power[8/27, 2/3]", "4/9", 0);
    assert_eval_eq("Power[4/9, 1/2]",  "2/3", 0);

    /* Numerator-only perfect extraction with non-perfect denominator. */
    assert_eval_eq("Power[8/3, 1/2]", "2 Sqrt[2/3]", 0);

    /* No extraction possible -- left unevaluated. */
    assert_eval_eq("Power[2/3, 1/2]",  "Sqrt[2/3]", 0);
    assert_eval_eq("Power[4/9, 2/3]",  "(4/9)^(2/3)", 0);

    /* Negative numerator routes through the integer negative-base path. */
    assert_eval_eq("Power[-1/8, 1/3]", "1/2 (-1)^(1/3)", 0);
    assert_eval_eq("Power[-1/4, 1/2]", "1/2*I", 0);
}

/* Power[Times[positive_factors], p/q] distribution (2026-05-12).
 * Distributes when every factor is known positive AND at least one
 * factor cleanly reduces to a rational coefficient under Power[_, p/q]. */
static void test_power_of_times_positives(void) {
    /* Integer perfect-square factor extracts cleanly. */
    assert_eval_eq("Sqrt[4 Pi]",   "2 Sqrt[Pi]", 0);
    assert_eval_eq("(Pi/8)^(2/3)", "1/4 Pi^(2/3)", 0);

    /* Conservative when no factor cleanly reduces. */
    assert_eval_eq("Sqrt[2 Pi]",       "Sqrt[2 Pi]", 0);
    assert_eval_eq("Power[4 Pi, 2/3]", "(4 Pi)^(2/3)", 0);
    assert_eval_eq("Sqrt[2 Sqrt[3]]",  "Sqrt[2 Sqrt[3]]", 0);

    /* Chains with rational-base distribution: Sqrt[(1/54)^(2/3)] cascades
     * (1/54)^(2/3) -> 1/9 * 2^(-2/3), then Sqrt of that distributes. */
    assert_eval_eq("Sqrt[Power[1/54, 2/3]]", "1/(3 2^(1/3))", 0);

    /* Nested clean extraction. */
    assert_eval_eq("Sqrt[4 Sqrt[3]]", "2 3^(1/4)", 0);
}

/* Power[Times[positives], p/q] distribution when the trigger is a
 * Power[positive, p'/q'] factor that composes cleanly (composed
 * reduced denominator <= q'), not a perfect-power rational factor
 * (2026-05-13).  Closes Sqrt[k/2^(2/3)] -> Sqrt[k]/2^(1/3) for
 * arbitrary positive k. */
static void test_power_of_times_power_factor_composes(void) {
    /* Power factor triggers distribution; rational coefficient need
     * not factor cleanly under the outer p/q. */
    assert_eval_eq("Sqrt[3/2^(2/3)]",      "Sqrt[3]/2^(1/3)",       0);
    assert_eval_eq("Sqrt[1/3/2^(2/3)]",    "1/(2^(1/3) Sqrt[3])",   0);
    assert_eval_eq("Sqrt[Pi/2^(2/3)]",     "Sqrt[Pi]/2^(1/3)",      0);

    /* Pre-existing perfect-square cases still work. */
    assert_eval_eq("Sqrt[1/2^(2/3)]",      "1/2^(1/3)",             0);
    assert_eval_eq("Sqrt[9/2^(2/3)]",      "3/2^(1/3)",             0);

    /* Adjacency cases must remain unchanged: composition would
     * introduce a new (larger-denominator) radical. */
    assert_eval_eq("Sqrt[2 Sqrt[3]]",      "Sqrt[2 Sqrt[3]]",       0);
    assert_eval_eq("Sqrt[2 Pi]",           "Sqrt[2 Pi]",            0);
    assert_eval_eq("Power[4 Pi, 2/3]",     "(4 Pi)^(2/3)",          0);
}

/* Nested rational powers with |inner exponent| < 1 compose on the principal
 * branch for ANY base (no positivity assumption), and a positive numeric
 * coefficient that fully reduces under the outer denominator splits out of a
 * mixed Times base. Together these reduce Sqrt[(1/27/a)^(2/3)] to
 * 1/3 (1/a)^(1/3) for symbolic a (2026-06-15). */
static void test_power_nested_symbolic_compose(void) {
    /* Main case: compose (z^(2/3))^(1/2) -> z^(1/3), then split out 1/27. */
    assert_eval_eq("Sqrt[(1/27/a)^(2/3)]", "1/3 (1/a)^(1/3)", 0);
    /* Composition alone, symbolic base (|2/3| < 1). */
    assert_eval_eq("Sqrt[a^(2/3)]",        "a^(1/3)",         0);
    /* Repeated Sqrt folds: |1/2| < 1 at each level. */
    assert_eval_eq("Sqrt[Sqrt[Sqrt[a]]]",  "a^(1/8)",         0);
    /* Inner exponent |p| >= q must NOT compose for unknown-sign base. */
    assert_eval_eq("Sqrt[a^(3/2)]",        "Sqrt[a^(3/2)]",   0);
    /* Symbolic inner exponent must NOT compose (branch-cut unsafe). */
    assert_eval_eq("Sqrt[Power[2, a]]",    "Sqrt[2^a]",       0);
}

/* Power-of-Power composition with positive base + rational outer (2026-05-12). */
static void test_power_of_power_positive_base(void) {
    /* Positive base with a *real* (rational) inner exponent merges:
     * (2^(2/3))^(1/2) = 2^(1/3), since 2^(2/3) is a positive real. */
    assert_eval_eq("Sqrt[Power[2, 2/3]]", "2^(1/3)", 0);

    /* Positive base but *symbolic* inner exponent must NOT merge: a may be
     * complex, and Sqrt[2^a] != 2^(a/2) on the principal branch (e.g.
     * a = 3 Pi I / Log[2] gives Sqrt[2^a] = I but 2^(a/2) = -I). Matches
     * Mathematica, which leaves Sqrt[2^a] unevaluated -- PowerExpand is the
     * documented route to the merged form. */
    assert_eval_eq("Sqrt[Power[2, a]]", "Sqrt[2^a]", 0);

    /* Pi positive -- Sqrt[Pi^2] reduces to Pi (not |Pi|). */
    assert_eval_eq("Sqrt[Power[Pi, 2]]", "Pi", 0);

    /* Unknown-sign symbol -- must NOT reduce (Sqrt[x^4] != x^2 for x<0). */
    assert_eval_eq("Sqrt[Power[x, 4]]", "Sqrt[x^4]", 0);

    /* Integer outer exponent on negative base still works (pre-existing rule). */
    assert_eval_eq("Power[Power[-1, 2], 3]", "1", 0);

    /* (E^x)^y with symbolic exponents must stay nested -- x may be complex,
     * so the core evaluator must not fold it to E^(x y). PowerExpand handles
     * the merge when the user asserts it is safe. */
    assert_eval_eq("(E^x)^y", "(E^x)^y", 0);
    assert_eval_eq("(E^x)^y", "Power[Power[E, x], y]", 1);
}

/* Perfect-power extraction on BigInt bases: previously skipped because
 * the gate insisted on EXPR_INTEGER. */
static void test_bigint_perfect_root(void) {
    assert_eval_eq("Sqrt[10^20]",  "10000000000", 0);
    assert_eval_eq("Sqrt[10^100]", "100000000000000000000000000000000000000000000000000", 0);
    assert_eval_eq("Sqrt[2^64]",   "4294967296", 0);
    assert_eval_eq("Sqrt[2^100]",  "1125899906842624", 0);
    assert_eval_eq("Power[10^60, 1/3]",  "100000000000000000000", 0);
    assert_eval_eq("Power[10^100, 1/2]", "100000000000000000000000000000000000000000000000000", 0);
    assert_eval_eq("Power[2^100, 1/4]",  "33554432", 0);
    assert_eval_eq("Power[10^100, 1/5]", "100000000000000000000", 0);
    assert_eval_eq("Power[10^21, 1/3]",  "10000000", 0);
    /* Mixed exponent: 10^20 = 10^(20), so (10^20)^(2/3) = 10^(40/3) which
     * factors as 10^13 * 10^(1/3). */
    assert_eval_eq("Power[10^20, 2/3]",  "10000000000000 10^(1/3)", 0);

    /* Negative BigInt base with q == 2: Sqrt[-bigint] = I * Sqrt[|bigint|]
     * via the I^p * |n|^(p/2) identity. Previously left unreduced because
     * the negative-base path gated on EXPR_INTEGER only. */
    assert_eval_eq("Sqrt[-(10^100)]", "100000000000000000000000000000000000000000000000000*I", 0);
    assert_eval_eq("Sqrt[-(10^50)]",  "10000000000000000000000000*I", 0);
    assert_eval_eq("Sqrt[-(2^64)]",   "4294967296*I", 0);
    /* (-bigint)^(3/2): I^3 = -I. */
    assert_eval_eq("Power[-(10^100), 3/2]",
        "-1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000*I", 0);
}

/* Constant complex-power fold in Simplify: Power[c1, c2] -> Exp[c2 Log[c1]]
 * for a variable-free non-positive-real base, kept only on a whole-expression
 * complexity win. */
void test_const_cpow_fold(void) {
    /* Targets: principal-value identities that collapse to 0 / 1. */
    assert_eval_eq("Simplify[I^I - Exp[-Pi/2]]", "0", 0);
    assert_eval_eq("Simplify[2 I^I - 2 Exp[-Pi/2]]", "0", 0);
    assert_eval_eq("Simplify[(-1)^I Exp[Pi]]", "1", 0);
    assert_eval_eq("Simplify[I^(2 I) Exp[Pi]]", "1", 0);

    /* Soundness / non-regression: must stay unchanged (the Exp form is not
     * strictly simpler, or the base is a positive real, or a surd). */
    assert_eval_eq("Simplify[2^I]", "2^I", 0);
    assert_eval_eq("Simplify[I^I]", "I^I", 0);
    assert_eval_eq("Simplify[(-1)^(1/3)]", "(-1)^(1/3)", 0);
    assert_eval_eq("Simplify[(-8)^(1/3)]", "2 (-1)^(1/3)", 0);
}

/* ---- machine Power: one helper on every surface ------------------------- *
 * A machine real (or complex) raised to an integer must be the SAME double on
 * every surface: the scalar builtin, Power with an integer-valued Real
 * exponent, a plain List, a packed List, a visible NDArray, Compile[] at rank 1
 * and at scalar shape, and an auto-compiled Table body. Before
 * src/mth_machine_pow.h the interpreter took libm pow, the VM binary powering
 * and the ND kernel x*x*x / pow, so x^3 differed in the last bit for ~26% of
 * random reals (x^5 ~51%; x^2 ~0.1%, Apple's pow(x, 2) not being x*x). The
 * reference is the helper itself; its accuracy is checked against MPFR. */
#define PWS_N 640
static double pws_x[PWS_N];

static uint64_t pws_lcg(uint64_t* s) {
    *s = *s * 6364136223846793005ull + 1442695040888963407ull;
    return *s;
}
static double pws_unit(uint64_t* s) {             /* uniform in [-2, 2) */
    return ((double)(pws_lcg(s) >> 11) * 0x1p-53) * 4.0 - 2.0;
}

static void pws_inputs(void) {
    uint64_t s = 0x9E3779B97F4A7C15ull;
    int k = 0;
    for (; k < 520; k++) pws_x[k] = pws_unit(&s);
    for (int i = 0; i < 40; i++) pws_x[k++] = ldexp(1.0 + i / 41.0, 4 * i - 80);  /* 2^-80 .. 2^76 */
    for (int i = 0; i < 20; i++) pws_x[k++] = -ldexp(1.5 + i / 50.0, 15 * i - 150);
    /* straddle the double-double range bounds of |n| = 3, 5, 32 */
    const double b[3] = { 0x1p313, 0x1p188, 0x1p29 };
    for (int j = 0; j < 3; j++) {
        pws_x[k++] = b[j];        pws_x[k++] = b[j] * (1 + 0x1p-52);
        pws_x[k++] = 1.0 / b[j];  pws_x[k++] = (1.0 / b[j]) * (1 - 0x1p-53);
    }
    while (k < PWS_N) { pws_x[k] = 1.0 + 1e-3 * k; k++; }
}

static int pws_same_bits(double a, double b) { return memcmp(&a, &b, sizeof a) == 0; }

static Expr* pws_eval(const char* src) {
    Expr* p = parse_expression(src);
    Expr* e = evaluate(p);
    expr_free(p);
    return e;
}

/* Bind the symbol `name` to List[elems...] (the elements are consumed). */
static void pws_bind(const char* name, Expr** elems, int n) {
    Expr* l = expr_new_function(expr_new_symbol("List"), elems, (size_t)n);
    Expr* sa[2] = { expr_new_symbol(name), l };
    Expr* set = expr_new_function(expr_new_symbol("Set"), sa, 2);
    expr_free(evaluate(set));
    expr_free(set);
}

/* `src` must be a List of PWS_N machine reals equal bit for bit to want[]. */
static void pws_check_real(const char* src, const double* want) {
    Expr* e = pws_eval(src);
    ASSERT_MSG(e->type == EXPR_FUNCTION && e->data.function.arg_count == PWS_N,
               "%s: not a list of %d", src, PWS_N);
    for (int i = 0; i < PWS_N; i++) {
        if (!isnormal(want[i])) continue;          /* overflow/underflow: see below */
        Expr* r = e->data.function.args[i];
        ASSERT_MSG(r->type == EXPR_REAL && pws_same_bits(r->data.real, want[i]),
                   "%s: element %d (x = %.17g) is %.17g, helper gives %.17g",
                   src, i + 1, pws_x[i], r->type == EXPR_REAL ? r->data.real : NAN, want[i]);
    }
    expr_free(e);
}

static void test_power_machine_surfaces_bitexact(void) {
    static const int64_t ns[] = { 2, 3, 4, 5, 7, 10, 16, 32, 33, -1, -2, -3, -5, -32, -40 };
    static const char* surf[] = {
        "Developer`FromPackedArray[pwx]^%lld",
        "Developer`FromPackedArray[Developer`ToPackedArray[pwx]^%lld]",
        "Developer`FromPackedArray[Developer`ToPackedArray[pwx]^(%lld.)]",
        "Normal[NDArray[pwx]^%lld]",
        "Developer`FromPackedArray[Compile[{{v, _Real, 1}}, v^%lld][pwx]]",
        "Developer`FromPackedArray[Compile[{{v, _Real, 1}}, 0. + v^%lld][pwx]]",
        "Developer`FromPackedArray[Map[Compile[{{u, _Real}}, u^%lld], pwx]]",
        "Developer`FromPackedArray[Map[Compile[{{u, _Real}, {e, _Real}}, u^e][#, %lld.] &, pwx]]",
        "Developer`FromPackedArray[Table[pwx[[i]]^%lld, {i, Length[pwx]}]]",
    };
    double want[PWS_N];
    char src[512];
    pws_inputs();
    Expr* el[PWS_N];
    for (int i = 0; i < PWS_N; i++) el[i] = expr_new_real(pws_x[i]);
    pws_bind("pwx", el, PWS_N);
    for (size_t j = 0; j < sizeof ns / sizeof *ns; j++) {
        int64_t n = ns[j];
        for (int i = 0; i < PWS_N; i++) want[i] = mth_powi_machine(pws_x[i], n);
        /* the scalar builtin, one call per element: Integer and Real exponent.
         * A power that overflows or underflows a double is no machine number
         * on the scalar surface (it promotes to an extended-exponent real), so
         * those cells are skipped here and in pws_check_real. */
        for (int i = 0; i < PWS_N; i++) {
            if (!isnormal(want[i])) continue;
            for (int realexp = 0; realexp < 2; realexp++) {
                Expr* a[2] = { expr_new_real(pws_x[i]),
                               realexp ? expr_new_real((double)n) : expr_new_integer(n) };
                Expr* call = expr_new_function(expr_new_symbol("Power"), a, 2);
                Expr* r = evaluate(call);
                ASSERT_MSG(r->type == EXPR_REAL && pws_same_bits(r->data.real, want[i]),
                           "Power[%.17g, %lld%s] = %.17g, helper gives %.17g", pws_x[i],
                           (long long)n, realexp ? "." : "",
                           r->type == EXPR_REAL ? r->data.real : NAN, want[i]);
                expr_free(call); expr_free(r);
            }
        }
        for (size_t q = 0; q < sizeof surf / sizeof *surf; q++) {
            snprintf(src, sizeof src, surf[q], (long long)n);
            pws_check_real(src, want);
        }
    }
    /* constant-folded inside Compile[] (the optimiser's OP_POWI_R fold) */
    for (int i = 0; i < 40; i++) {
        snprintf(src, sizeof src, "Compile[{{u, _Real}}, u + %.17g^5][0.]", pws_x[i]);
        Expr* r = pws_eval(src);
        double w = mth_powi_machine(pws_x[i], 5);
        ASSERT_MSG(r->type == EXPR_REAL && pws_same_bits(r->data.real, w),
                   "%s = %.17g, helper %.17g", src, r->type == EXPR_REAL ? r->data.real : NAN, w);
        expr_free(r);
    }
}

/* Complex base, integer exponent: Re and Im bit for bit on every surface. The
 * scalar builtin answers a Real when the imaginary part is exactly 0. */
#define PWC_N 160
static void test_power_machine_complex_surfaces_bitexact(void) {
    static const int64_t ns[] = { 2, 3, 5, 7, 16, -1, -2, -3, -6 };
    static const char* surf[] = {
        "Developer`FromPackedArray[Table[Apply[Power, {pwz[[i]], %lld}], {i, Length[pwz]}]]",
        "Developer`FromPackedArray[pwz]^%lld",
        "Developer`FromPackedArray[Developer`ToPackedArray[pwz]^%lld]",
        "Normal[NDArray[pwz]^%lld]",
        "Developer`FromPackedArray[Compile[{{v, _Complex, 1}}, v^%lld][pwz]]",
        "Developer`FromPackedArray[Map[Compile[{{u, _Complex}}, u^%lld], pwz]]",
    };
    double zr[PWC_N], zi[PWC_N], wr[PWC_N], wi[PWC_N];
    uint64_t s = 12345;
    for (int k = 0; k < PWC_N; k++) { zr[k] = pws_unit(&s); zi[k] = pws_unit(&s); }
    zr[0] = 0.0;   zi[0] = 1.0;                    /* I^n exactly */
    zr[1] = 1e200; zi[1] = 3e199;                  /* Smith's inversion: |z|^2 overflows */
    Expr* el[PWC_N];
    for (int k = 0; k < PWC_N; k++) {
        Expr* c[2] = { expr_new_real(zr[k]), expr_new_real(zi[k]) };
        el[k] = expr_new_function(expr_new_symbol("Complex"), c, 2);
    }
    pws_bind("pwz", el, PWC_N);
    for (size_t j = 0; j < sizeof ns / sizeof *ns; j++) {
        for (int k = 0; k < PWC_N; k++) mth_powi_machine_c(zr[k], zi[k], ns[j], &wr[k], &wi[k]);
        for (size_t q = 0; q < sizeof surf / sizeof *surf; q++) {
            char src[256];
            snprintf(src, sizeof src, surf[q], (long long)ns[j]);
            Expr* e = pws_eval(src);
            ASSERT_MSG(e->type == EXPR_FUNCTION && e->data.function.arg_count == PWC_N,
                       "%s: not a list of %d", src, PWC_N);
            for (int k = 0; k < PWC_N; k++) {
                /* z = 1e200 overflows every power but z^-1 (which is the point
                 * of it: Smith's 1/z never forms |z|^2) */
                if (!isfinite(wr[k]) || !isfinite(wi[k])) continue;
                Expr* r = e->data.function.args[k];
                double re = NAN, im = NAN;
                if (r->type == EXPR_REAL) { re = r->data.real; im = 0.0; }
                else if (r->type == EXPR_FUNCTION && r->data.function.arg_count == 2 &&
                         r->data.function.args[0]->type == EXPR_REAL &&
                         r->data.function.args[1]->type == EXPR_REAL) {
                    re = r->data.function.args[0]->data.real;
                    im = r->data.function.args[1]->data.real;
                }
                ASSERT_MSG(pws_same_bits(re, wr[k]) &&
                           (pws_same_bits(im, wi[k]) || (r->type == EXPR_REAL && wi[k] == 0.0)),
                           "%s: element %d (z = %.17g + %.17g I) is %.17g + %.17g I, helper %.17g + %.17g I",
                           src, k + 1, zr[k], zi[k], re, im, wr[k], wi[k]);
            }
            expr_free(e);
        }
    }
}

/* The helper's accuracy, which is why it was chosen over libm pow and over
 * plain binary powering: correctly rounded against MPFR on random inputs
 * wherever its double-double arm applies. */
static void test_power_machine_helper_accuracy(void) {
#ifdef USE_MPFR
    static const int64_t ns[] = { 2, 3, 4, 5, 6, 7, 8, 11, 16, 25, 32, -1, -2, -3, -5, -9, -32 };
    mpfr_t x, r;
    mpfr_init2(x, 53); mpfr_init2(r, 53);
    uint64_t s = 777;
    for (size_t j = 0; j < sizeof ns / sizeof *ns; j++) {
        int64_t n = ns[j];
        unsigned m = (unsigned)(n < 0 ? -n : n);
        int miss = 0, tried = 0;
        for (int i = 0; i < 4000; i++) {
            double xv = pws_unit(&s);
            if (m >= 2 && n != 2 && (fabs(xv) > MTH_POWI_BOUND[m] || fabs(xv) * MTH_POWI_BOUND[m] < 1.0))
                continue;                            /* libm pow's range */
            if (xv == 0.0) continue;
            mpfr_set_d(x, xv, MPFR_RNDN);
            mpfr_pow_si(r, x, (long)n, MPFR_RNDN);  /* correctly rounded */
            tried++;
            if (!pws_same_bits(mth_powi_machine(xv, n), mpfr_get_d(r, MPFR_RNDN))) miss++;
        }
        ASSERT_MSG(tried > 3900 && miss == 0, "x^%lld: %d of %d not correctly rounded",
                   (long long)n, miss, tried);
    }
    mpfr_clear(x); mpfr_clear(r);
#endif
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_bug3_imag_times_positive);
    TEST(test_bug4_negative_integer_part);
    TEST(test_bug5_zero_to_zero);
    TEST(test_rational_base_rational_exp);
    TEST(test_power_of_times_positives);
    TEST(test_power_of_power_positive_base);
    TEST(test_power_of_times_power_factor_composes);
    TEST(test_power_nested_symbolic_compose);
    TEST(test_bigint_perfect_root);
    TEST(test_const_cpow_fold);
    TEST(test_power_machine_surfaces_bitexact);
    TEST(test_power_machine_complex_surfaces_bitexact);
    TEST(test_power_machine_helper_accuracy);
    TEST(test_complex_power_overflow);

    printf("All power_corpus tests passed!\n");
    return 0;
}

/* test_ffield.c -- the ffield library (src/poly/ffield/, Phase 4.5): the
 * generic ring Q(p_1..p_r)[t]/(mu) (ff_fnx.c), its Expr conversions
 * (ff_expr.c), and the FLINT` polynomial-matrix kernels over it
 * ("Parameters" -> {p...}, with or without Extension -> {t, mu}).
 *
 *  - ring laws on seeded random elements in C: associativity and
 *    distributivity of the product modulo mu, x * inv(x) = 1, (x / y) y = x,
 *    and a gcd over K[x] that recovers a planted common factor;
 *  - Expr round trips: a polynomial read into K[x] and printed back equals
 *    the input by value (denominators in the parameters, t reduced mod mu);
 *  - the kernels against the Trager package's own routines (tragerOption
 *    "NativeFn" False, the reference) and against their defining invariants;
 *  - the out-of-scope contract: a denominator with t or x, a reducible mu,
 *    an unknown symbol leave the call unevaluated.
 * The independent oracles (PARI/GP and the reference CAS) are run by
 * tools/ffield_stress/check.py, a developer tool.
 */
#include "trager_test_utils.h"

#ifdef USE_FLINT
#include "poly/ffield/ff_fnx.h"
#include "poly/ffield/ff_expr.h"
#include <flint/gr_poly.h>

static void ff_check(int cond, const char* what, int line) {
    tr_assertions++;
    if (!cond) { printf("FAIL %s:%d  %s\n", __FILE__, line, what); exit(1); }
}
#define FF_CHECK(c) ff_check((c), #c, __LINE__)

/* a random element: coefficients (small rational) / (1 + small linear form) */
static void ff_rand(gr_ptr e, gr_ctx_t K, flint_rand_t st) {
    ff_fnx_data* D = FF_FNX_DATA(K);
    fmpz_mpoly_q_t g, s; fmpz_mpoly_q_init(g, D->mctx); fmpz_mpoly_q_init(s, D->mctx);
    for (slong j = 0; j < D->n; j++) {
        fmpz_mpoly_q_struct* c = FF_FNX_COEFF(e, j);
        fmpz_mpoly_q_set_si(c, (slong)n_randint(st, 9) - 4, D->mctx);
        for (slong v = 0; v < D->nvars; v++) {
            fmpz_mpoly_q_gen(g, v, D->mctx);
            fmpz_mpoly_q_mul_si(g, g, (slong)n_randint(st, 5) - 2, D->mctx);
            fmpz_mpoly_q_add(c, c, g, D->mctx);
        }
        if (n_randint(st, 2)) {
            fmpz_mpoly_q_one(s, D->mctx);
            if (D->nvars > 0) { fmpz_mpoly_q_gen(g, 0, D->mctx); fmpz_mpoly_q_add(s, s, g, D->mctx); }
            fmpz_mpoly_q_div(c, c, s, D->mctx);
        }
    }
    fmpz_mpoly_q_clear(g, D->mctx); fmpz_mpoly_q_clear(s, D->mctx);
}

static void ring_laws(slong nvars, const char* mu_str) {
    fmpq_poly_t mu; fmpq_poly_init(mu);
    if (mu_str) fmpq_poly_set_str(mu, mu_str);
    gr_ctx_t K;
    FF_CHECK(ff_fnx_ctx_init(K, nvars, mu_str ? mu : NULL) == GR_SUCCESS);
    FF_CHECK(gr_ctx_is_field(K) == T_TRUE);
    flint_rand_t st; flint_rand_init(st);
    gr_ptr a = gr_heap_init(K), b = gr_heap_init(K), c = gr_heap_init(K);
    gr_ptr l = gr_heap_init(K), r = gr_heap_init(K), u = gr_heap_init(K);
    for (int it = 0; it < 25; it++) {
        ff_rand(a, K, st); ff_rand(b, K, st); ff_rand(c, K, st);
        int s = GR_SUCCESS;
        /* (a b) c = a (b c) */
        s |= gr_mul(l, a, b, K); s |= gr_mul(l, l, c, K);
        s |= gr_mul(r, b, c, K); s |= gr_mul(r, a, r, K);
        FF_CHECK(s == GR_SUCCESS && gr_equal(l, r, K) == T_TRUE);
        /* a (b + c) = a b + a c */
        s |= gr_add(u, b, c, K); s |= gr_mul(l, a, u, K);
        s |= gr_mul(r, a, b, K); s |= gr_mul(u, a, c, K); s |= gr_add(r, r, u, K);
        FF_CHECK(s == GR_SUCCESS && gr_equal(l, r, K) == T_TRUE);
        if (gr_is_zero(a, K) != T_TRUE) {
            /* a inv(a) = 1 and (b / a) a = b */
            FF_CHECK(gr_inv(u, a, K) == GR_SUCCESS);
            FF_CHECK(gr_mul(l, u, a, K) == GR_SUCCESS && gr_is_one(l, K) == T_TRUE);
            FF_CHECK(gr_div(l, b, a, K) == GR_SUCCESS && gr_mul(l, l, a, K) == GR_SUCCESS);
            FF_CHECK(gr_equal(l, b, K) == T_TRUE);
        }
    }
    /* inverse of zero is a domain error, not a crash */
    FF_CHECK(gr_zero(a, K) == GR_SUCCESS && gr_inv(u, a, K) == GR_DOMAIN);
    /* gcd over K[x] of (x - a)(x - b) and (x - a)(x + c) is x - a (for b != -c) */
    gr_poly_t P, Q, G, L1, L2, L3; gr_poly_init(P, K); gr_poly_init(Q, K); gr_poly_init(G, K);
    gr_poly_init(L1, K); gr_poly_init(L2, K); gr_poly_init(L3, K);
    for (int it = 0; it < 6; it++) {
        int s = GR_SUCCESS;
        ff_rand(a, K, st); ff_rand(b, K, st); ff_rand(c, K, st);
        s |= gr_add(u, b, c, K);
        if (gr_is_zero(u, K) == T_TRUE) continue;
        s |= gr_one(l, K);
        s |= gr_neg(r, a, K); s |= gr_poly_set_coeff_scalar(L1, 1, l, K); s |= gr_poly_set_coeff_scalar(L1, 0, r, K);
        s |= gr_neg(r, b, K); s |= gr_poly_set_coeff_scalar(L2, 1, l, K); s |= gr_poly_set_coeff_scalar(L2, 0, r, K);
        s |= gr_poly_set_coeff_scalar(L3, 1, l, K); s |= gr_poly_set_coeff_scalar(L3, 0, c, K);
        s |= gr_poly_mul(P, L1, L2, K); s |= gr_poly_mul(Q, L1, L3, K);
        s |= gr_poly_gcd(G, P, Q, K); s |= gr_poly_make_monic(G, G, K);
        FF_CHECK(s == GR_SUCCESS && gr_poly_equal(G, L1, K) == T_TRUE);
    }
    gr_poly_clear(P, K); gr_poly_clear(Q, K); gr_poly_clear(G, K);
    gr_poly_clear(L1, K); gr_poly_clear(L2, K); gr_poly_clear(L3, K);
    gr_heap_clear(a, K); gr_heap_clear(b, K); gr_heap_clear(c, K);
    gr_heap_clear(l, K); gr_heap_clear(r, K); gr_heap_clear(u, K);
    flint_rand_clear(st);
    gr_ctx_clear(K);
    fmpq_poly_clear(mu);
}

static void test_ring(void) {
    ring_laws(1, NULL);            /* Q(a) */
    ring_laws(2, NULL);            /* Q(a, b) */
    ring_laws(1, "3  -2 0 1");     /* Q(Sqrt[2])(a) */
    ring_laws(2, "3  1 0 1");      /* Q(i)(a, b) */
    ring_laws(1, "4  -2 0 0 1");   /* Q(2^(1/3))(a) */
    /* a non-monic mu is refused */
    fmpq_poly_t mu; fmpq_poly_init(mu); fmpq_poly_set_str(mu, "3  -2 0 2");
    gr_ctx_t K;
    FF_CHECK(ff_fnx_ctx_init(K, 1, mu) == GR_DOMAIN);
    fmpq_poly_clear(mu);
}
#endif

static void test_round_trip(void) {
    /* QuotientRemainder by 1 returns the input itself, through K[x] and back */
    TR_TRUE("Module[{e = (3 a x^2)/(1 + a) - x/(2 a) + 5}, Together[First[Last[{"
            "FLINT`PolynomialQuotientRemainder[e, 1, x, \"Parameters\" -> {a}]}]] - e] === 0]");
    TR_TRUE("Module[{e = t x^2/(a - b) + (1 + t) x/(a b) - 3/(1 + b)}, Together[First["
            "FLINT`PolynomialQuotientRemainder[e, 1, x, \"Parameters\" -> {a, b}, Extension -> {t, t^2 - 2}]] - e] === 0]");
    /* t^3 reduced modulo t^2 - 2 */
    TR_TRUE("Together[First[FLINT`PolynomialQuotientRemainder[t^3 x, 1, x, \"Parameters\" -> {a}, Extension -> {t, t^2 - 2}]] - 2 t x] === 0");
}

static void test_kernels_vs_package(void) {
    /* Hermite forms over Q(a) and Q(Sqrt[2])(a) against the package's hnf */
    TR_TRUE("Module[{M = {{x^2 - a, x + 1/a}, {(x - 1/(a + 1)) x, x^2 + a x}}, n, p},"
            " n = FLINT`PolynomialHermiteForm[M, x, \"Parameters\" -> {a}];"
            " p = hnf[fieldFn[{a}], M, x];"
            " ListQ[n] && AllTrue[Flatten[n - p], Together[#] === 0 &]]");
    TR_TRUE("Module[{M = {{x^2 - a t, x + t/a}, {(x - (a - t)/(a^2 - 2)) x, x^2 + a x}}, n, p},"
            " n = FLINT`PolynomialHermiteForm[M, x, Extension -> {t, t^2 - 2}, \"Parameters\" -> {a}];"
            " p = hnf[fieldFn[{a}, {t, t^2 - 2, Sqrt[2]}], M, x];"
            " ListQ[n] && AllTrue[Flatten[n - p], Expand[PolynomialRemainder[Numerator[Together[#]], t^2 - 2, t]] === 0 &]]");
    /* extended gcd: g monic, s A + u B = g, g | A */
    TR_TRUE("Module[{F = (x - a)/(1 + a), A, B, r}, A = Expand[F (x + 1)]; B = Expand[F (x^2 - b)];"
            " r = FLINT`PolynomialExtendedGCD[A, B, x, \"Parameters\" -> {a, b}];"
            " Together[r[[1]] - (x - a)] === 0 && Together[r[[2, 1]] A + r[[2, 2]] B - r[[1]]] === 0]");
    /* determinant and adjugate */
    TR_TRUE("Module[{M = {{x + a, t}, {1/a, x - t}}, r},"
            " r = FLINT`PolynomialDetAdjugate[M, x, \"Parameters\" -> {a}, Extension -> {t, t^2 + 1}];"
            " AllTrue[Flatten[{r[[1]] - Det[M], r[[2]] . M - r[[1]] IdentityMatrix[2]}],"
            "   Expand[PolynomialRemainder[Numerator[Together[#]], t^2 + 1, t]] === 0 &]]");
    /* the ideal operations through the native route agree with the package by value */
    TR_TRUE("Module[{F = fieldFn[{a}, {t, t^2 - 2, Sqrt[2]}], M = {{x^2 - a t, x + t/a}, {0, x^2 - 1}}, r1, r2},"
            " tragerOption[\"NativeFn\"] = False; r1 = idealReduce[F, <|\"m\" -> 2, \"x\" -> x|>, <|\"M\" -> M, \"den\" -> x - 1|>];"
            " tragerOption[\"NativeFn\"] = True;  r2 = idealReduce[F, <|\"m\" -> 2, \"x\" -> x|>, <|\"M\" -> M, \"den\" -> x - 1|>];"
            " tragerOption[\"NativeFn\"] = False;"
            " Together[r1[\"den\"] - r2[\"den\"]] === 0 && AllTrue[Flatten[r1[\"M\"] - r2[\"M\"]], Together[#] === 0 &]]");
}

static void test_out_of_scope(void) {
    /* t in a denominator: unevaluated */
    TR_TRUE("Head[FLINT`PolynomialHermiteForm[{{x, 1/(a + t)}, {0, x}}, x, Extension -> {t, t^2 - 2}, \"Parameters\" -> {a}]] === FLINT`PolynomialHermiteForm");
    /* x in a denominator */
    TR_TRUE("Head[FLINT`PolynomialQuotientRemainder[1/x, 1, x, \"Parameters\" -> {a}]] === FLINT`PolynomialQuotientRemainder");
    /* a symbol that is not a parameter */
    TR_TRUE("Head[FLINT`PolynomialQuotientRemainder[c x, 1, x, \"Parameters\" -> {a}]] === FLINT`PolynomialQuotientRemainder");
    /* a reducible mu */
    TR_TRUE("Head[FLINT`PolynomialQuotientRemainder[t x, 1, x, \"Parameters\" -> {a}, Extension -> {t, t^2 - 1}]] === FLINT`PolynomialQuotientRemainder");
    /* parameters with a modulus */
    TR_TRUE("Head[FLINT`PolynomialQuotientRemainder[a x, 1, x, \"Parameters\" -> {a}, Modulus -> 7]] === FLINT`PolynomialQuotientRemainder");
}

int main(void) {
    tr_init();
#ifdef USE_FLINT
    test_ring();
#endif
    test_round_trip();
    test_kernels_vs_package();
    test_out_of_scope();
    printf("ffield_tests: all %d assertions passed\n", tr_assertions);
    return 0;
}

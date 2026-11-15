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
 *    an unknown symbol leave the call unevaluated;
 *  - the certified multi-modular Hermite form (ff_hnf.c) against the direct
 *    elimination on seeded random matrices over Q and number fields whose
 *    minimal polynomials split mod p in different ways, its declines (prime
 *    budget, non-monic modulus), and the "Modular" switch of the builtin;
 *  - the certified multi-modular extended gcd (ff_xgcd.c) against gr_poly_xgcd
 *    on planted common factors, and its declines (a zero input, cofactors
 *    that are not unique).
 * The independent oracles (PARI/GP and the reference CAS) are run by
 * tools/ffield_stress/check.py, a developer tool.
 */
#include "trager_test_utils.h"

#ifdef USE_FLINT
#include "poly/ffield/ff_fnx.h"
#include "poly/ffield/ff_expr.h"
#include "poly/ffield/ff_hnf.h"
#include "poly/ffield/ff_xgcd.h"
#include <flint/gr_poly.h>
#include <flint/nf.h>
#include <flint/nf_elem.h>

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

#ifdef USE_FLINT
/* ---- the multi-modular Hermite form (ff_hnf.c) ---- */

/* a random t-polynomial of degree < n with coefficients up to B / (1..3) */
static void hm_rand_tp(fmpq_poly_t a, slong n, slong B, flint_rand_t st) {
    fmpq_poly_zero(a);
    fmpq_t q; fmpq_init(q);
    for (slong j = 0; j < n; j++) {
        if (n_randint(st, 3) == 0) continue;
        fmpq_set_si(q, (slong)n_randint(st, (ulong)(2 * B + 1)) - B, 1 + (ulong)n_randint(st, 3));
        fmpq_poly_set_coeff_fmpq(a, j, q);
    }
    fmpq_clear(q);
}

static void hm_rand_poly(ff_qtpoly_struct* a, slong maxlen, slong n, slong B, flint_rand_t st) {
    ff_qtpoly_set_length(a, (slong)n_randint(st, (ulong)maxlen + 1));
    for (slong d = 0; d < a->len; d++) hm_rand_tp(a->c + d, n, B, st);
    ff_qtpoly_normalise(a);
}

/* K = Q (nf NULL) or Q[t]/(mu) as a gr field, and qtpoly <-> gr_poly */
static void hm_to_gr(gr_poly_t P, const ff_qtpoly_struct* a, gr_ctx_t K, nf_t nf, int isnf) {
    slong sz = gr_ctx_sizeof_elem(K);
    FF_CHECK(gr_poly_zero(P, K) == GR_SUCCESS);
    if (a->len == 0) return;
    gr_poly_fit_length(P, a->len, K);
    fmpq_t q; fmpq_init(q);
    for (slong d = 0; d < a->len; d++) {
        gr_ptr c = (char*)P->coeffs + d * sz;
        if (isnf) nf_elem_set_fmpq_poly((nf_elem_struct*)c, a->c + d, nf);
        else { fmpq_poly_get_coeff_fmpq(q, a->c + d, 0); FF_CHECK(gr_set_fmpq(c, q, K) == GR_SUCCESS); }
    }
    fmpq_clear(q);
    _gr_poly_set_length(P, a->len, K);
    _gr_poly_normalise(P, K);
}

/* one random case: the modular form equals the direct form (both are the
 * Hermite form of span(rows) + D K[u]^nc); returns the modular verdict */
static int hm_case(const char* mu_str, slong nr, slong nc, slong maxlen, slong B, flint_rand_t st,
                   ff_mm_stats* stats) {
    fmpq_poly_t mu; fmpq_poly_init(mu);
    int isnf = mu_str != NULL;
    slong n = 1;
    if (isnf) { FF_CHECK(fmpq_poly_set_str(mu, mu_str) == 0); n = fmpq_poly_degree(mu); }
    gr_ctx_t K; nf_t nf;
    if (isnf) { nf_init(nf, mu); gr_ctx_init_nf(K, mu); } else gr_ctx_init_fmpq(K);
    ff_qtpoly_struct* rows = malloc((size_t)(nr * nc) * sizeof(ff_qtpoly_struct));
    ff_qtpoly_struct* H = malloc((size_t)(nc * nc) * sizeof(ff_qtpoly_struct));
    ff_qtpoly_struct D;
    for (slong i = 0; i < nr * nc; i++) { ff_qtpoly_init(rows + i); hm_rand_poly(rows + i, maxlen, n, B, st); }
    for (slong i = 0; i < nc * nc; i++) ff_qtpoly_init(H + i);
    ff_qtpoly_init(&D);
    /* D monic of degree 2..5 */
    slong dd = 2 + (slong)n_randint(st, 4);
    ff_qtpoly_set_length(&D, dd + 1);
    for (slong d = 0; d < dd; d++) hm_rand_tp(D.c + d, n, B, st);
    fmpq_poly_one(D.c + dd);
    int ok = ff_hnf_modular(H, rows, nr, nc, &D, isnf ? mu : NULL, stats);
    if (ok) {
        /* the direct elimination over K */
        slong R = nr + nc;
        gr_poly_struct* e = malloc((size_t)(R * nc) * sizeof(gr_poly_struct));
        for (slong i = 0; i < R * nc; i++) gr_poly_init(e + i, K);
        gr_poly_t Dg, h; gr_poly_init(Dg, K); gr_poly_init(h, K);
        hm_to_gr(Dg, &D, K, nf, isnf);
        for (slong i = 0; i < nr * nc; i++) {
            hm_to_gr(e + i, rows + i, K, nf, isnf);
            FF_CHECK(gr_poly_rem(e + i, e + i, Dg, K) == GR_SUCCESS);
        }
        for (slong k = 0; k < nc; k++) FF_CHECK(gr_poly_set(e + (nr + k) * nc + k, Dg, K) == GR_SUCCESS);
        FF_CHECK(ff_hnf_gr(e, R, nc, Dg, K));
        for (slong i = 0; i < nc * nc; i++) {
            hm_to_gr(h, H + i, K, nf, isnf);
            FF_CHECK(gr_poly_equal(h, e + i, K) == T_TRUE);
        }
        for (slong i = 0; i < R * nc; i++) gr_poly_clear(e + i, K);
        free(e);
        gr_poly_clear(Dg, K); gr_poly_clear(h, K);
    }
    for (slong i = 0; i < nr * nc; i++) ff_qtpoly_clear(rows + i);
    for (slong i = 0; i < nc * nc; i++) ff_qtpoly_clear(H + i);
    ff_qtpoly_clear(&D);
    free(rows); free(H);
    gr_ctx_clear(K);
    if (isnf) nf_clear(nf);
    fmpq_poly_clear(mu);
    return ok;
}

static void test_hnf_modular(void) {
    flint_rand_t st; flint_rand_init(st);
    /* Q, Q(Sqrt[2]), Q(omega), and a cubic with Galois group S3 (mu mod p
     * irreducible, linear times quadratic, or split, by p) */
    const char* fields[] = { NULL, "3  -2 0 1", "3  1 1 1", "4  -1 -1 0 1" };
    for (int f = 0; f < 4; f++)
        for (int it = 0; it < 25; it++) {
            ff_mm_stats S;
            slong nc = 1 + (slong)n_randint(st, 3), nr = 1 + (slong)n_randint(st, 4);
            /* the kernel accepts every case of this size within its budget */
            FF_CHECK(hm_case(fields[f], nr, nc, 6, 40, st, &S));
            FF_CHECK(S.primes >= 1 && S.certificates >= 1);
        }
    /* large heights need several primes */
    {
        ff_mm_stats S;
        FF_CHECK(hm_case("3  1 1 1", 4, 3, 9, 1000000, st, &S));
        FF_CHECK(S.primes >= 2);
    }
    /* a budget of one prime cannot reconstruct a form with large heights: it declines */
    ff_mm_set_max_primes(1);
    {
        ff_mm_stats S;
        FF_CHECK(!hm_case("3  1 1 1", 4, 3, 9, 1000000, st, &S));
    }
    ff_mm_set_max_primes(0);
    /* a modulus that is not monic is declined */
    {
        ff_qtpoly_struct r, H, D;
        ff_qtpoly_init(&r); ff_qtpoly_init(&H); ff_qtpoly_init(&D);
        ff_qtpoly_set_length(&r, 1); fmpq_poly_set_si(r.c + 0, 1);
        ff_qtpoly_set_length(&D, 2); fmpq_poly_set_si(D.c + 1, 2);
        FF_CHECK(!ff_hnf_modular(&H, &r, 1, 1, &D, NULL, NULL));
        ff_qtpoly_clear(&r); ff_qtpoly_clear(&H); ff_qtpoly_clear(&D);
    }
    flint_rand_clear(st);
}

/* the modular extended gcd against gr_poly_xgcd, normalised (g monic) */
static void hm_mul(ff_qtpoly_struct* r, const ff_qtpoly_struct* a, const ff_qtpoly_struct* b,
                   const fmpq_poly_t mu, int isnf) {
    ff_qtpoly_set_length(r, a->len && b->len ? a->len + b->len - 1 : 0);
    fmpq_poly_t m; fmpq_poly_init(m);
    for (slong i = 0; i < a->len; i++)
        for (slong j = 0; j < b->len; j++) {
            fmpq_poly_mul(m, a->c + i, b->c + j);
            if (isnf) fmpq_poly_rem(m, m, mu);
            fmpq_poly_add(r->c + i + j, r->c + i + j, m);
        }
    fmpq_poly_clear(m);
    ff_qtpoly_normalise(r);
}

static void test_xgcd_modular(void) {
    flint_rand_t st; flint_rand_init(st);
    const char* fields[] = { NULL, "3  -2 0 1", "3  1 1 1", "4  -1 -1 0 1" };
    for (int f = 0; f < 4; f++)
        for (int it = 0; it < 25; it++) {
            fmpq_poly_t mu; fmpq_poly_init(mu);
            int isnf = fields[f] != NULL;
            slong n = 1;
            if (isnf) { FF_CHECK(fmpq_poly_set_str(mu, fields[f]) == 0); n = fmpq_poly_degree(mu); }
            gr_ctx_t K; nf_t nf;
            if (isnf) { nf_init(nf, mu); gr_ctx_init_nf(K, mu); } else gr_ctx_init_fmpq(K);
            ff_qtpoly_struct c, a0, b0, a, b, g, s, t;
            ff_qtpoly_init(&c); ff_qtpoly_init(&a0); ff_qtpoly_init(&b0); ff_qtpoly_init(&a);
            ff_qtpoly_init(&b); ff_qtpoly_init(&g); ff_qtpoly_init(&s); ff_qtpoly_init(&t);
            /* a = c a0, b = c b0: a planted common factor of degree 0..3 */
            hm_rand_poly(&c, 4, n, 30, st);
            hm_rand_poly(&a0, 6, n, 30, st);
            hm_rand_poly(&b0, 6, n, 30, st);
            if (c.len == 0) { ff_qtpoly_set_length(&c, 1); fmpq_poly_set_si(c.c, 1); }
            hm_mul(&a, &c, &a0, mu, isnf);
            hm_mul(&b, &c, &b0, mu, isnf);
            ff_mm_stats S;
            /* the gcd alone: always certified for nonzero inputs, and equal
             * to the monic gcd of gr_poly_gcd */
            if (a.len > 0 && b.len > 0) {
                ff_qtpoly_struct g1; ff_qtpoly_init(&g1);
                FF_CHECK(ff_gcd_modular(&g1, &a, &b, isnf ? mu : NULL, &S));
                gr_poly_t A1, B1, G1, h1;
                gr_poly_init(A1, K); gr_poly_init(B1, K); gr_poly_init(G1, K); gr_poly_init(h1, K);
                hm_to_gr(A1, &a, K, nf, isnf); hm_to_gr(B1, &b, K, nf, isnf);
                FF_CHECK(gr_poly_gcd(G1, A1, B1, K) == GR_SUCCESS && gr_poly_make_monic(G1, G1, K) == GR_SUCCESS);
                hm_to_gr(h1, &g1, K, nf, isnf);
                FF_CHECK(gr_poly_equal(h1, G1, K) == T_TRUE);
                gr_poly_clear(A1, K); gr_poly_clear(B1, K); gr_poly_clear(G1, K); gr_poly_clear(h1, K);
                ff_qtpoly_clear(&g1);
            }
            /* exact division: a / c = a0 (c nonzero), and a / c not exact
             * after adding 1 to a when deg c >= 1 */
            if (a.len > 0) {
                ff_qtpoly_struct q; ff_qtpoly_init(&q);
                FF_CHECK(ff_divexact(&q, &a, &c, isnf ? mu : NULL));
                gr_poly_t Q1, A0; gr_poly_init(Q1, K); gr_poly_init(A0, K);
                hm_to_gr(Q1, &q, K, nf, isnf); hm_to_gr(A0, &a0, K, nf, isnf);
                FF_CHECK(gr_poly_equal(Q1, A0, K) == T_TRUE);
                gr_poly_clear(Q1, K); gr_poly_clear(A0, K);
                if (c.len >= 2) {
                    ff_qtpoly_struct a1; ff_qtpoly_init(&a1);
                    ff_qtpoly_set_length(&a1, a.len);
                    for (slong d = 0; d < a.len; d++) fmpq_poly_set(a1.c + d, a.c + d);
                    fmpq_poly_add_si(a1.c + 0, a1.c + 0, 1);
                    FF_CHECK(!ff_divexact(&q, &a1, &c, isnf ? mu : NULL));
                    ff_qtpoly_clear(&a1);
                }
                ff_qtpoly_clear(&q);
            }
            int ok = ff_xgcd_modular(&g, &s, &t, &a, &b, isnf ? mu : NULL, &S);
            if (a.len == 0 || b.len == 0) FF_CHECK(!ok);
            else {
                gr_poly_t A, B, G, S1, T1, h;
                gr_poly_init(A, K); gr_poly_init(B, K); gr_poly_init(G, K);
                gr_poly_init(S1, K); gr_poly_init(T1, K); gr_poly_init(h, K);
                hm_to_gr(A, &a, K, nf, isnf); hm_to_gr(B, &b, K, nf, isnf);
                FF_CHECK(gr_poly_xgcd(G, S1, T1, A, B, K) == GR_SUCCESS);
                slong dg = G->length - 1;
                if (dg == a.len - 1 && dg == b.len - 1) FF_CHECK(!ok);   /* cofactors not unique */
                else {
                    FF_CHECK(ok && S.certificates >= 1);
                    gr_ptr lc = gr_heap_init(K);
                    FF_CHECK(gr_inv(lc, (char*)G->coeffs + dg * gr_ctx_sizeof_elem(K), K) == GR_SUCCESS);
                    FF_CHECK(gr_poly_mul_scalar(G, G, lc, K) == GR_SUCCESS);
                    FF_CHECK(gr_poly_mul_scalar(S1, S1, lc, K) == GR_SUCCESS);
                    FF_CHECK(gr_poly_mul_scalar(T1, T1, lc, K) == GR_SUCCESS);
                    gr_heap_clear(lc, K);
                    hm_to_gr(h, &g, K, nf, isnf); FF_CHECK(gr_poly_equal(h, G, K) == T_TRUE);
                    hm_to_gr(h, &s, K, nf, isnf); FF_CHECK(gr_poly_equal(h, S1, K) == T_TRUE);
                    hm_to_gr(h, &t, K, nf, isnf); FF_CHECK(gr_poly_equal(h, T1, K) == T_TRUE);
                }
                gr_poly_clear(A, K); gr_poly_clear(B, K); gr_poly_clear(G, K);
                gr_poly_clear(S1, K); gr_poly_clear(T1, K); gr_poly_clear(h, K);
            }
            ff_qtpoly_clear(&c); ff_qtpoly_clear(&a0); ff_qtpoly_clear(&b0); ff_qtpoly_clear(&a);
            ff_qtpoly_clear(&b); ff_qtpoly_clear(&g); ff_qtpoly_clear(&s); ff_qtpoly_clear(&t);
            gr_ctx_clear(K);
            if (isnf) nf_clear(nf);
            fmpq_poly_clear(mu);
        }
    flint_rand_clear(st);
}
#endif /* USE_FLINT */

static void test_hnf_modular_builtin(void) {
    /* the builtin: modular and direct forms agree over Q(omega) with large heights */
    TR_TRUE("Module[{M = {{12345 x^3 - 678/7 t x + 9, x^2 - 3 t}, {x - 1111/13, 98765 t x^2 + 1},"
            " {x^4 + t, 77/3 x}}, D = x^5 + 2 x + 3, a, b},"
            " a = FLINT`PolynomialHermiteForm[M, x, Extension -> {t, t^2 + t + 1}, \"DeterminantMultiple\" -> D];"
            " b = FLINT`PolynomialHermiteForm[M, x, Extension -> {t, t^2 + t + 1}, \"DeterminantMultiple\" -> D, \"Modular\" -> False];"
            " ListQ[a] && a === b]");
    /* over Q */
    TR_TRUE("Module[{M = {{3 x^2 - 22/7, x + 5}, {x^3 - 1, 1/9 x}}, D = x^4 - 2, a, b},"
            " a = FLINT`PolynomialHermiteForm[M, x, \"DeterminantMultiple\" -> D];"
            " b = FLINT`PolynomialHermiteForm[M, x, \"DeterminantMultiple\" -> D, \"Modular\" -> False];"
            " ListQ[a] && a === b]");
    /* over Q(theta)(p) with entries free of the parameter: the modular form, by value */
    TR_TRUE("Module[{M = {{x^2 - 7/3 t, x + 1}, {x^3 + t, 2 x}}, D = x^5 - 1, a, b},"
            " a = FLINT`PolynomialHermiteForm[M, x, Extension -> {t, t^2 - 2}, \"Parameters\" -> {c}, \"DeterminantMultiple\" -> D];"
            " b = FLINT`PolynomialHermiteForm[M, x, Extension -> {t, t^2 - 2}, \"Parameters\" -> {c}, \"DeterminantMultiple\" -> D, \"Modular\" -> False];"
            " ListQ[a] && AllTrue[Flatten[Together[a - b]], # === 0 &]]");
    /* with a parameter the modular kernel declines and the direct form is returned */
    TR_TRUE("Module[{M = {{x^2 - c t, x + 1/c}, {x^3 + t, 2 x}}, D = x^5 - 1, a, b},"
            " a = FLINT`PolynomialHermiteForm[M, x, Extension -> {t, t^2 - 2}, \"Parameters\" -> {c}, \"DeterminantMultiple\" -> D];"
            " b = FLINT`PolynomialHermiteForm[M, x, Extension -> {t, t^2 - 2}, \"Parameters\" -> {c}, \"DeterminantMultiple\" -> D, \"Modular\" -> False];"
            " ListQ[a] && a === b]");
    /* the extended gcd: modular and direct agree (a planted factor, large heights) */
    TR_TRUE("Module[{c = x^2 - 1234567/89 t x + 5, a, b, r1, r2}, a = Expand[c (x^3 + 7/3 t)]; b = Expand[c (x^2 - 99999 x + t)];"
            " r1 = FLINT`PolynomialExtendedGCD[a, b, x, Extension -> {t, t^2 + t + 1}];"
            " r2 = FLINT`PolynomialExtendedGCD[a, b, x, Extension -> {t, t^2 + t + 1}, \"Modular\" -> False];"
            " ListQ[r1] && r1 === r2 && Exponent[r1[[1]], x] === 2]");
    TR_TRUE("Head[FLINT`PolynomialExtendedGCD[x, x^2, x, \"DeterminantMultiple\" -> x]] === FLINT`PolynomialExtendedGCD");
    /* the field gcd through ff_gcd_modular: a planted factor over Q(Sqrt[2] + Sqrt[3]), and coprime inputs */
    TR_TRUE("Module[{mu = t^4 - 10 t^2 + 1, c, a, b, red}, red[e_] := Expand[PolynomialRemainder[Expand[e], mu, t]];"
            " c = x^2 + 3 t x - 7 t^3; a = red[c (x^3 + t x + 5)]; b = red[c (x^2 - 11 t^2)];"
            " FLINT`PolynomialFieldGCD[a, b, x, Extension -> {t, mu}] === red[c]]");
    TR_TRUE("FLINT`PolynomialFieldGCD[x^3 + t, x^2 - 3, x, Extension -> {t, t^2 + t + 1}] === 1");
    /* a bad option value leaves the call unevaluated */
    TR_TRUE("Head[FLINT`PolynomialHermiteForm[{{x}}, x, \"DeterminantMultiple\" -> x^2, \"Modular\" -> 3]] === FLINT`PolynomialHermiteForm");
}

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
    test_hnf_modular();
    test_xgcd_modular();
#endif
    test_hnf_modular_builtin();
    test_round_trip();
    test_kernels_vs_package();
    test_out_of_scope();
    printf("ffield_tests: all %d assertions passed\n", tr_assertions);
    return 0;
}

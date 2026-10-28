/*
 * ratfun_gcd.inc -- PolynomialGCD and PolynomialLCM on the ratfun context.
 *
 * Included by ratfun.c (inside USE_FLINT) so it shares the context, the
 * conversion and the certified output of Together/Cancel/Apart:
 *
 *  - ONE context over all the arguments (rf_ctx_build_n): the same atoms,
 *    generators and canonical variable order Together uses, so E^(2x) - 1 and
 *    E^x - 1 share the generator E^x and Log[x] is a variable;
 *  - each argument converted ONCE (rf_from_expr, which keeps a product's
 *    factors) and expanded to c * p, c rational, p primitive with a positive
 *    leading coefficient;
 *  - the gcd folded over the p's with fmpz_mpoly_gcd (fmpz_poly_gcd when there
 *    is one variable and the polynomials are dense), stopping at 1, and the
 *    rational contents folded with fmpq_gcd;
 *  - the answer rendered once, canonical and certified (rf_direct_ok).
 *
 * The output contract is the classical path's (poly.c), which the tests pin:
 *
 *   PolynomialGCD: the expanded gcd times the gcd of the contents,
 *       PolynomialGCD[x/2 + 1/2, x/3 + 1/3] -> 1/6 + x/6.
 *   PolynomialLCM: the classical fold L <- (L/g) a_i or L (a_i/g), whichever
 *       quotient has fewer terms, g the gcd including its content, kept as a
 *       product:  PolynomialLCM[2x + 2, 4x - 4] -> (1 + x)(-4 + 4x).
 *
 * Under Modulus -> p (p a word-size prime) the same runs over GF(p) on
 * nmod_mpoly: the gcd is monic, every coefficient a residue 0..p-1.
 */

#include <flint/nmod.h>
#include <flint/nmod_mpoly.h>
#include <flint/ulong_extras.h>

static int rg_enabled(void) {
    static int on = -1;
    if (on < 0) {
        const char* e = getenv("MATHILDA_POLYGCD_FLINT");
        on = !(e && strcmp(e, "off") == 0);
    }
    return on;
}

/* ------------------------------------------------------------------ */
/*  Over Q: c * p                                                      */
/* ------------------------------------------------------------------ */

/* c * p, p primitive with a positive leading coefficient; c = 0 is the zero
 * polynomial (p is then left 1). */
typedef struct { fmpq_t c; fmpz_mpoly_t p; } Rg;

static void rg_init(Rg* g, const RfCtx* C) { fmpq_init(g->c); fmpz_mpoly_init(g->p, C->ctx); }
static void rg_clear(Rg* g, const RfCtx* C) { fmpq_clear(g->c); fmpz_mpoly_clear(g->p, C->ctx); }

/* g = e; 0 when e is not a polynomial in the context.  With want_poly = 0
 * only the content is set (the gcd is already 1). */
static int rg_from_expr(const Expr* e, Rg* g, RfCtx* C, int want_poly) {
    Rf R; rf_init(&R);
    int ok = rf_from_expr(e, &R, C) && !C->fail && rf_is_poly(&R);
    if (ok) {
        fmpq_set(g->c, R.c);
        if (fmpq_is_zero(g->c) || !want_poly) fmpz_mpoly_one(g->p, C->ctx);
        else rf_expand_side(g->p, &R, 0, C);
    }
    rf_clear(&R, C);
    return ok;
}

/* G = gcd(A, B), positive leading coefficient.  One dense variable goes
 * through fmpz_poly, whose gcd is the faster one there; a sparse one (x^n + 1,
 * large n) stays in fmpz_mpoly, which does not densify it. */
static int rg_gcd(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const RfCtx* C) {
    const fmpz_mpoly_ctx_struct* ctx = C->ctx;
    if (C->n == 1 && A->length > 1 && B->length > 1) {
        slong da = fmpz_mpoly_degree_si(A, 0, ctx), db = fmpz_mpoly_degree_si(B, 0, ctx);
        if (da < 16 * A->length + 64 && db < 16 * B->length + 64) {
            fmpz_poly_t a, b, g;
            fmpz_poly_init(a); fmpz_poly_init(b); fmpz_poly_init(g);
            int ok = fmpz_mpoly_get_fmpz_poly(a, A, 0, ctx) && fmpz_mpoly_get_fmpz_poly(b, B, 0, ctx);
            if (ok) {
                fmpz_poly_gcd(g, a, b);
                fmpz_mpoly_set_fmpz_poly(G, g, 0, ctx);
            }
            fmpz_poly_clear(a); fmpz_poly_clear(b); fmpz_poly_clear(g);
            if (ok) return 1;
        }
    }
    return fmpz_mpoly_gcd(G, A, B, ctx);
}

/* The canonical order of distinct monomials in SYMBOL atoms, read off the
 * exponents: ascending lex with the last canonical atom most significant,
 *     1 + 2 x + x^2 + 2 y + 2 x y + y^2,
 * which is what expr_orderless_sort gives a Plus of them -- without an
 * expr_compare per pair, which was two thirds of a 500-term gcd's time. */
typedef struct { const ulong* ex; slong nv; Expr* t; } RgTerm;

static int rg_cmp_terms(const void* pa, const void* pb) {
    const RgTerm* a = pa; const RgTerm* b = pb;
    for (slong v = a->nv - 1; v >= 0; v--)
        if (a->ex[v] != b->ex[v]) return a->ex[v] < b->ex[v] ? -1 : 1;
    return 0;
}

static int rg_symbol_atoms(const RfCtx* C) {
    for (size_t v = 0; v < C->n; v++)
        if (C->atom[v]->type != EXPR_SYMBOL) return 0;
    return 1;
}

/* terms[0..len) of p (FLINT's descending lex, atom 0 most significant) put in
 * canonical order: reversed for one variable, sorted on the exponents for
 * several symbols, expr_orderless_sort otherwise. */
static void rg_sort_terms(Expr** terms, const fmpz_mpoly_t p, const RfCtx* C) {
    slong len = p->length, nv = (slong)C->n;
    if (nv == 1) {
        for (slong i = 0, j = len - 1; i < j; i++, j--) { Expr* t = terms[i]; terms[i] = terms[j]; terms[j] = t; }
        return;
    }
    ulong* ex = rg_symbol_atoms(C) ? malloc((size_t)(len * nv) * sizeof(ulong)) : NULL;
    RgTerm* T = ex ? malloc((size_t)len * sizeof(RgTerm)) : NULL;
    if (!T) { free(ex); expr_orderless_sort(terms, (size_t)len); return; }
    for (slong i = 0; i < len; i++) {
        fmpz_mpoly_get_term_exp_ui(ex + i * nv, p, i, C->ctx);
        T[i].ex = ex + i * nv; T[i].nv = nv; T[i].t = terms[i];
    }
    qsort(T, (size_t)len, sizeof(RgTerm), rg_cmp_terms);
    for (slong i = 0; i < len; i++) terms[i] = T[i].t;
    free(T); free(ex);
}

/* Q = A/B when B divides A exactly; dense in one variable as rg_gcd. */
static int rg_divides(fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const RfCtx* C) {
    const fmpz_mpoly_ctx_struct* ctx = C->ctx;
    if (C->n == 1 && A->length > 1 && B->length > 1 &&
        fmpz_mpoly_degree_si(A, 0, ctx) < 16 * A->length + 64) {
        fmpz_poly_t a, b, q;
        fmpz_poly_init(a); fmpz_poly_init(b); fmpz_poly_init(q);
        int ok = fmpz_mpoly_get_fmpz_poly(a, A, 0, ctx) && fmpz_mpoly_get_fmpz_poly(b, B, 0, ctx) &&
                 fmpz_poly_divides(q, a, b);
        if (ok) fmpz_mpoly_set_fmpz_poly(Q, q, 0, ctx);
        fmpz_poly_clear(a); fmpz_poly_clear(b); fmpz_poly_clear(q);
        return ok;
    }
    return fmpz_mpoly_divides(Q, A, B, ctx);
}

/* c * p expanded: a sum of monomials with rational coefficients, canonical
 * and certified when the context allows it (rf_direct_ok), else built and
 * evaluated once.  MATHILDA_RATFUN_DIRECT_CHECK=1 compares the certified sum
 * with the evaluated one. */
static Expr* rg_render_raw(const fmpq_t c, const fmpz_mpoly_t p, const RfCtx* C, int direct);

static Expr* rg_render(const fmpq_t c, const fmpz_mpoly_t p, const RfCtx* C) {
    static int check = -1;
    if (check < 0) check = getenv("MATHILDA_RATFUN_DIRECT_CHECK") ? 1 : 0;
    int direct = rf_direct_ok(C);
    Expr* r = rg_render_raw(c, p, C, direct);
    if (r && direct && check) {
        Expr* e = rg_render_raw(c, p, C, 0);
        if (e && !expr_eq(r, e)) {
            char* sd = expr_to_string(r); char* se = expr_to_string(e);
            fprintf(stderr, "RATFUN_DIRECT_MISMATCH direct=%s evaluated=%s\n", sd, se);
            free(sd); free(se);
        }
        if (e) expr_free(e);
    }
    return r;
}

static Expr* rg_render_raw(const fmpq_t c, const fmpz_mpoly_t p, const RfCtx* C, int direct) {
    if (fmpq_is_zero(c) || p->length == 0) return expr_new_integer(0);
    slong len = p->length, nv = (slong)C->n;
    ulong* ex = malloc((size_t)nv * sizeof(ulong));
    Expr** terms = malloc((size_t)len * sizeof(Expr*));
    Expr** fs = malloc((size_t)(nv + 1) * sizeof(Expr*));
    if (!ex || !terms || !fs) { free(ex); free(terms); free(fs); return NULL; }
    fmpq_t q; fmpq_init(q);
    for (slong i = 0; i < len; i++) {
        fmpz_mpoly_get_term_exp_ui(ex, p, i, C->ctx);
        size_t nf = 0;
        int lead;
        if (fmpq_is_one(c)) {                     /* an integer coefficient: no fmpq */
            const fmpz* z = p->coeffs + i;
            lead = !fmpz_is_one(z);
            if (lead) fs[nf++] = rf_fmpz_to_expr(z);
        } else {
            fmpq_mul_fmpz(q, c, p->coeffs + i);
            lead = !fmpq_is_one(q);
            if (lead) fs[nf++] = direct ? rf_number_expr(q) : fb_fmpq_to_expr(q);
        }
        for (slong v = 0; v < nv; v++) {
            if (!ex[v]) continue;
            Expr* a = rf_atom_expr(C, v);
            if (ex[v] == 1) { fs[nf++] = a; continue; }
            Expr* pw = expr_new_function(expr_new_symbol(SYM_Power),
                           (Expr*[]){ a, expr_new_integer((int64_t)ex[v]) }, 2);
            fs[nf++] = direct ? rf_certified(pw) : pw;
        }
        if (direct) expr_orderless_sort(fs + lead, nf - (size_t)lead);
        if (nf == 0) terms[i] = expr_new_integer(1);
        else if (nf == 1) terms[i] = fs[0];
        else {
            Expr* t = expr_new_function(expr_new_symbol(SYM_Times), fs, nf);
            terms[i] = direct ? rf_certified(t) : t;
        }
    }
    fmpq_clear(q);
    Expr* r;
    if (len == 1) r = terms[0];
    else {
        if (direct) rg_sort_terms(terms, p, C);
        r = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)len);
        if (direct) rf_certified(r);
    }
    free(ex); free(terms); free(fs);
    return direct ? r : eval_and_free(r);
}

/* The GCD over Q. */
static Expr* rg_gcd_q(Expr* const* a, size_t n, RfCtx* C) {
    fmpq_t cg; fmpq_init(cg);
    fmpz_mpoly_t G; fmpz_mpoly_init(G, C->ctx);
    Rg t; rg_init(&t, C);
    int any = 0, ok = 1;
    for (size_t i = 0; i < n && ok; i++) {
        int one = any && fmpz_mpoly_is_one(G, C->ctx);
        if (!rg_from_expr(a[i], &t, C, !one)) { ok = 0; break; }
        if (fmpq_is_zero(t.c)) continue;
        if (!any) {
            fmpq_abs(cg, t.c);
            fmpz_mpoly_swap(G, t.p, C->ctx);
            any = 1;
            continue;
        }
        fmpq_gcd(cg, cg, t.c);
        if (!one && !rg_gcd(G, G, t.p, C)) ok = 0;
    }
    Expr* r = !ok ? NULL : !any ? expr_new_integer(0) : rg_render(cg, G, C);
    rg_clear(&t, C);
    fmpz_mpoly_clear(G, C->ctx);
    fmpq_clear(cg);
    return r;
}

/* The shown factors as one product, evaluated once; the factors are moved
 * into it (*ns = 0). */
static Expr* rg_product(Expr** show, size_t* ns) {
    size_t k = *ns;
    *ns = 0;
    if (k == 1) return show[0];
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), show, k));
}

/* The number of terms the classical path counts for a quotient: a sum's
 * length, 1 for anything else. */
static slong rg_terms(const fmpz_mpoly_t p) { return p->length > 1 ? p->length : 1; }

/* The LCM over Q: the classical fold, the shown factors kept as they come
 * (an argument as given, a quotient expanded).
 *
 * L is held as the shown factors' primitive parts F_1..F_k (and a content),
 * never expanded: gcd(L, a) is taken factor by factor,
 *     gcd(F_1 ... F_k, a) = h_1 h_2 ... h_k,  h_j = gcd(F_j, a / (h_1 ... h_(j-1))),
 * every gcd the size of one argument rather than of the growing LCM (the gcd
 * of a 10-argument LCM's running product was half its time), stopping once
 * the remaining cofactor is 1; that cofactor is a/g, the second quotient. */
static Expr* rg_lcm_q(Expr* const* a, size_t n, RfCtx* C) {
    const fmpz_mpoly_ctx_struct* ctx = C->ctx;
    fmpq_t Lc, gc, qc; fmpq_init(Lc); fmpq_init(gc); fmpq_init(qc);
    fmpz_mpoly_struct* F = malloc(n * sizeof *F);     /* primitive parts of show[] */
    fmpz_mpoly_struct* H = malloc(n * sizeof *H);     /* h_j */
    Expr** show = malloc(n * sizeof(Expr*));
    if (!F || !H || !show) { free(F); free(H); free(show); fmpq_clear(Lc); fmpq_clear(gc); fmpq_clear(qc); return NULL; }
    for (size_t k = 0; k < n; k++) { fmpz_mpoly_init(F + k, ctx); fmpz_mpoly_init(H + k, ctx); }
    fmpz_mpoly_t g, r, q1, u; fmpz_mpoly_init(g, ctx); fmpz_mpoly_init(r, ctx);
    fmpz_mpoly_init(q1, ctx); fmpz_mpoly_init(u, ctx);
    Rg t; rg_init(&t, C);
    size_t ns = 0;
    int zero = 0;
    int ok = rg_from_expr(a[0], &t, C, 1);
    if (ok) {
        zero = fmpq_is_zero(t.c);
        fmpq_set(Lc, t.c);
        fmpz_mpoly_swap(F + 0, t.p, ctx);
        show[ns++] = expr_copy(a[0]);
    }
    for (size_t i = 1; i < n && ok && !zero; i++) {
        if (!(ok = rg_from_expr(a[i], &t, C, 1))) break;
        if (fmpq_is_zero(t.c)) { zero = 1; break; }
        fmpq_gcd(gc, Lc, t.c);
        fmpz_mpoly_one(g, ctx);
        fmpz_mpoly_set(r, t.p, ctx);
        for (size_t j = 0; j < ns && ok; j++) {
            if (fmpz_mpoly_is_one(r, ctx)) { fmpz_mpoly_one(H + j, ctx); continue; }
            if (!(ok = rg_gcd(H + j, F + j, r, C))) break;
            if (fmpz_mpoly_is_one(H + j, ctx)) continue;
            fmpz_mpoly_mul(g, g, H + j, ctx);
            ok = rg_divides(r, r, H + j, C);
        }
        if (!ok) break;
        if (fmpq_is_one(gc) && fmpz_mpoly_is_one(g, ctx)) {
            fmpz_mpoly_swap(F + ns, t.p, ctx);
            show[ns++] = expr_copy(a[i]);
            fmpq_mul(Lc, Lc, t.c);
            continue;
        }
        /* Q1 = L/g = prod F_j/h_j, Q2 = a_i/g = r; the shorter quotient is
         * shown, the other side kept as it was */
        fmpz_mpoly_one(q1, ctx);
        for (size_t j = 0; j < ns && ok; j++) {
            if (fmpz_mpoly_is_one(H + j, ctx)) fmpz_mpoly_mul(q1, q1, F + j, ctx);
            else if ((ok = rg_divides(u, F + j, H + j, C))) fmpz_mpoly_mul(q1, q1, u, ctx);
        }
        if (!ok) break;
        if (rg_terms(q1) <= rg_terms(r)) {
            fmpq_div(qc, Lc, gc);
            Expr* e1 = rg_render(qc, q1, C);
            for (size_t k = 0; k < ns; k++) expr_free(show[k]);
            ns = 0;
            if (!(ok = e1 != NULL)) break;
            fmpz_mpoly_swap(F + 0, q1, ctx);
            fmpz_mpoly_swap(F + 1, t.p, ctx);
            show[ns++] = e1;
            show[ns++] = expr_copy(a[i]);
        } else {
            fmpq_div(qc, t.c, gc);
            Expr* e2 = rg_render(qc, r, C);
            if (!(ok = e2 != NULL)) break;
            fmpz_mpoly_swap(F + ns, r, ctx);
            show[ns++] = e2;
        }
        fmpq_mul(Lc, Lc, t.c);                 /* L = (L/g) a_i */
        fmpq_div(Lc, Lc, gc);
    }
    Expr* r_out = (ok && zero) ? expr_new_integer(0) : ok ? rg_product(show, &ns) : NULL;
    for (size_t k = 0; k < ns; k++) expr_free(show[k]);
    free(show);
    for (size_t k = 0; k < n; k++) { fmpz_mpoly_clear(F + k, ctx); fmpz_mpoly_clear(H + k, ctx); }
    free(F); free(H);
    fmpz_mpoly_clear(g, ctx); fmpz_mpoly_clear(r, ctx); fmpz_mpoly_clear(q1, ctx); fmpz_mpoly_clear(u, ctx);
    rg_clear(&t, C);
    fmpq_clear(Lc); fmpq_clear(gc); fmpq_clear(qc);
    return r_out;
}

/* ------------------------------------------------------------------ */
/*  Over GF(p)                                                         */
/* ------------------------------------------------------------------ */

/* A number's residue mod p; 0 on a denominator divisible by p. */
static int rgm_residue(const Expr* c, nmod_t mod, ulong* r) {
    if (!c) { *r = 1; return 1; }
    if (c->type == EXPR_INTEGER) {           /* p <= INT64_MAX (rg_entry) */
        int64_t m = c->data.integer % (int64_t)mod.n;
        *r = (ulong)(m < 0 ? m + (int64_t)mod.n : m);
        return 1;
    }
    if (c->type == EXPR_BIGINT) { *r = mpz_fdiv_ui(c->data.bigint, mod.n); return 1; }
    fmpq_t q; fmpq_init(q);
    int ok = fb_expr_to_fmpq(c, q);
    if (ok) {
        ulong d = fmpz_fdiv_ui(fmpq_denref(q), mod.n);
        ok = d != 0;
        if (ok) *r = nmod_mul(fmpz_fdiv_ui(fmpq_numref(q), mod.n), n_invmod(d, mod.n), mod);
    }
    fmpq_clear(q);
    return ok;
}

/* out = e mod p.  A sum of monomial terms is read straight into nmod_mpoly
 * (every coefficient reduced as it is met: no content, no big integers);
 * anything else goes through the rational form. */
static int rgm_from_expr(const Expr* e, nmod_mpoly_t out, const nmod_mpoly_ctx_t M,
                         RfCtx* C, nmod_t mod) {
    nmod_mpoly_zero(out, M);
    slong nv = (slong)C->n;
    ulong* ex = malloc((size_t)nv * sizeof(ulong));
    if (!ex) return 0;
    const Expr* one = e;
    size_t nt = 1;
    const Expr* const* ts = &one;
    if (rf_head_is(e, SYM_Plus)) { nt = e->data.function.arg_count; ts = (const Expr* const*)e->data.function.args; }
    int ok = 1, direct = 1;
    const Expr* cf;
    for (size_t k = 0; k < nt && direct; k++) {
        ulong r;
        if (!rf_monomial_term(ts[k], C, ex, &cf)) { direct = 0; break; }
        if (!rgm_residue(cf, mod, &r)) { ok = 0; break; }
        if (r) nmod_mpoly_push_term_ui_ui(out, r, ex, M);
    }
    if (ok && direct) {
        nmod_mpoly_sort_terms(out, M);
        nmod_mpoly_combine_like_terms(out, M);
    } else if (ok) {
        nmod_mpoly_zero(out, M);
        Rg t; rg_init(&t, C);
        ok = rg_from_expr(e, &t, C, 1);
        ulong c = 0;
        if (ok && !fmpq_is_zero(t.c)) {
            ulong d = fmpz_fdiv_ui(fmpq_denref(t.c), mod.n);
            ok = d != 0;
            if (ok) c = nmod_mul(fmpz_fdiv_ui(fmpq_numref(t.c), mod.n), n_invmod(d, mod.n), mod);
        }
        if (ok && c) {
            for (slong i = 0; i < t.p->length; i++) {
                ulong r = nmod_mul(fmpz_fdiv_ui(t.p->coeffs + i, mod.n), c, mod);
                if (!r) continue;
                fmpz_mpoly_get_term_exp_ui(ex, t.p, i, C->ctx);
                nmod_mpoly_push_term_ui_ui(out, r, ex, M);
            }
            nmod_mpoly_sort_terms(out, M);
            nmod_mpoly_combine_like_terms(out, M);
        }
        rg_clear(&t, C);
    }
    free(ex);
    return ok;
}

/* G = gcd(A, B), monic; nmod_poly in one variable. */
static int rgm_gcd(nmod_mpoly_t G, const nmod_mpoly_t A, const nmod_mpoly_t B,
                   const nmod_mpoly_ctx_t M, const RfCtx* C, ulong p) {
    if (C->n == 1) {
        nmod_poly_t a, b, g;
        nmod_poly_init(a, p); nmod_poly_init(b, p); nmod_poly_init(g, p);
        int ok = nmod_mpoly_get_nmod_poly(a, A, 0, M) && nmod_mpoly_get_nmod_poly(b, B, 0, M);
        if (ok) {
            nmod_poly_gcd(g, a, b);
            nmod_mpoly_set_nmod_poly(G, g, 0, M);
        }
        nmod_poly_clear(a); nmod_poly_clear(b); nmod_poly_clear(g);
        if (ok) return 1;
    }
    return nmod_mpoly_gcd(G, A, B, M);
}

/* The residues as an expanded sum (coefficients 0..p-1). */
static Expr* rgm_render(const nmod_mpoly_t P, const nmod_mpoly_ctx_t M, const RfCtx* C) {
    if (nmod_mpoly_is_zero(P, M)) return expr_new_integer(0);
    fmpz_mpoly_t Z; fmpz_mpoly_init(Z, C->ctx);
    ulong* ex = malloc(C->n * sizeof(ulong));
    if (!ex) { fmpz_mpoly_clear(Z, C->ctx); return NULL; }
    for (slong i = 0; i < P->length; i++) {
        nmod_mpoly_get_term_exp_ui(ex, P, i, M);
        fmpz_mpoly_push_term_ui_ui(Z, P->coeffs[i], ex, C->ctx);
    }
    free(ex);
    fmpz_mpoly_sort_terms(Z, C->ctx);
    fmpq_t one; fmpq_init(one); fmpq_one(one);
    Expr* r = rg_render(one, Z, C);
    fmpq_clear(one);
    fmpz_mpoly_clear(Z, C->ctx);
    return r;
}

static slong rgm_terms(const nmod_mpoly_t p) { return p->length > 1 ? p->length : 1; }

static Expr* rg_gcd_mod(Expr* const* a, size_t n, RfCtx* C, ulong p) {
    nmod_mpoly_ctx_t M; nmod_mpoly_ctx_init(M, (slong)C->n, ORD_LEX, p);
    nmod_t mod; nmod_init(&mod, p);
    nmod_mpoly_t G, T; nmod_mpoly_init(G, M); nmod_mpoly_init(T, M);
    int any = 0, ok = 1;
    for (size_t i = 0; i < n && ok; i++) {
        if (!(ok = rgm_from_expr(a[i], T, M, C, mod))) break;
        if (nmod_mpoly_is_zero(T, M)) continue;
        if (!any) { nmod_mpoly_make_monic(G, T, M); any = 1; continue; }
        if (!nmod_mpoly_is_one(G, M)) ok = rgm_gcd(G, G, T, M, C, p);
    }
    Expr* r = !ok ? NULL : !any ? expr_new_integer(0) : rgm_render(G, M, C);
    nmod_mpoly_clear(G, M); nmod_mpoly_clear(T, M);
    nmod_mpoly_ctx_clear(M);
    return r;
}

/* The LCM over GF(p): the fold of rg_lcm_q with a monic gcd; an argument is
 * shown reduced (its residues). */
static Expr* rg_lcm_mod(Expr* const* a, size_t n, RfCtx* C, ulong p) {
    nmod_mpoly_ctx_t M; nmod_mpoly_ctx_init(M, (slong)C->n, ORD_LEX, p);
    nmod_t mod; nmod_init(&mod, p);
    nmod_mpoly_t L, T, G, Q1, Q2;
    nmod_mpoly_init(L, M); nmod_mpoly_init(T, M); nmod_mpoly_init(G, M);
    nmod_mpoly_init(Q1, M); nmod_mpoly_init(Q2, M);
    Expr** show = malloc(n * sizeof(Expr*));
    size_t ns = 0;
    int ok = show != NULL, zero = 0;
    if (ok && (ok = rgm_from_expr(a[0], L, M, C, mod))) {
        zero = nmod_mpoly_is_zero(L, M);
        Expr* e = zero ? NULL : rgm_render(L, M, C);
        if (!zero && !e) ok = 0; else if (e) show[ns++] = e;
    }
    for (size_t i = 1; i < n && ok && !zero; i++) {
        if (!(ok = rgm_from_expr(a[i], T, M, C, mod))) break;
        if (nmod_mpoly_is_zero(T, M)) { zero = 1; break; }
        if (!(ok = rgm_gcd(G, L, T, M, C, p))) break;
        if (nmod_mpoly_is_one(G, M)) {
            Expr* e = rgm_render(T, M, C);
            if (!e) { ok = 0; break; }
            show[ns++] = e;
            nmod_mpoly_mul(L, L, T, M);
            continue;
        }
        if (!(ok = nmod_mpoly_divides(Q1, L, G, M) && nmod_mpoly_divides(Q2, T, G, M))) break;
        if (rgm_terms(Q1) <= rgm_terms(Q2)) {
            Expr* e1 = rgm_render(Q1, M, C);
            Expr* e2 = rgm_render(T, M, C);
            for (size_t k = 0; k < ns; k++) expr_free(show[k]);
            ns = 0;
            if (e1) show[ns++] = e1;
            if (e2) show[ns++] = e2;
            if (!e1 || !e2) ok = 0;
        } else {
            Expr* e2 = rgm_render(Q2, M, C);
            if (!e2) ok = 0; else show[ns++] = e2;
        }
        nmod_mpoly_mul(L, Q1, T, M);
    }
    Expr* r = (ok && zero) ? expr_new_integer(0) : ok ? rg_product(show, &ns) : NULL;
    for (size_t k = 0; k < ns; k++) expr_free(show[k]);
    free(show);
    nmod_mpoly_clear(L, M); nmod_mpoly_clear(T, M); nmod_mpoly_clear(G, M);
    nmod_mpoly_clear(Q1, M); nmod_mpoly_clear(Q2, M);
    nmod_mpoly_ctx_clear(M);
    return r;
}

/* ------------------------------------------------------------------ */
/*  Entry points                                                       */
/* ------------------------------------------------------------------ */

/* Number arguments have no variables, so no context: they are left to
 * the classical path (integer gcd/lcm), as is a non-prime or huge modulus. */
static Expr* rg_entry(Expr* const* a, size_t n, uint64_t p, int lcm) {
    if (n < 2 || !rg_enabled()) return NULL;
    size_t nz = 0;
    for (size_t i = 0; i < n; i++)
        nz += a[i]->type == EXPR_INTEGER && a[i]->data.integer == 0;
    if (nz == n) return expr_new_integer(0);          /* gcd(0, 0) = lcm(0, 0) = 0 */
    if (p && (p > (uint64_t)INT64_MAX || !n_is_prime((ulong)p))) return NULL;
    RfCtx C;
    if (!rf_ctx_build_n((const Expr* const*)a, n, &C, NULL)) return NULL;
    Expr* r = p ? (lcm ? rg_lcm_mod(a, n, &C, (ulong)p) : rg_gcd_mod(a, n, &C, (ulong)p))
                : (lcm ? rg_lcm_q(a, n, &C) : rg_gcd_q(a, n, &C));
    rf_ctx_clear(&C);
    return r;
}

Expr* ratfun_polynomial_gcd(Expr* const* a, size_t n, uint64_t p) { return rg_entry(a, n, p, 0); }
Expr* ratfun_polynomial_lcm(Expr* const* a, size_t n, uint64_t p) { return rg_entry(a, n, p, 1); }

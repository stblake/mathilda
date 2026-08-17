/* gbmod.c
 *
 * Buchberger Gröbner basis over the prime field GF(p).  See gbmod.h.
 *
 * The polynomial substrate mirrors groebner.c's GBPoly (sparse exponent
 * tuples, sorted descending by the active monomial order) but carries native
 * uint64 residue coefficients instead of GMP rationals.  GF(p) is a field, so
 * division is multiplication by the modular inverse and the ordinary
 * field-Buchberger algorithm (with the coprime-leading-term criterion)
 * applies unchanged.
 */

#include "gbmod.h"

#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "core.h"        /* tc_check_deadline() */

/* ------------------------------------------------------------------ */
/*  Modular arithmetic (p < 2^31, so products fit in uint64).          */
/* ------------------------------------------------------------------ */

static inline uint64_t gfp_add(uint64_t a, uint64_t b, uint64_t p) {
    uint64_t s = a + b;
    return s >= p ? s - p : s;
}
static inline uint64_t gfp_sub(uint64_t a, uint64_t b, uint64_t p) {
    return a >= b ? a - b : a + p - b;
}
static inline uint64_t gfp_mul(uint64_t a, uint64_t b, uint64_t p) {
    return (a * b) % p;
}

bool gfp_is_prime(uint64_t p) {
    if (p < 2) return false;
    if (p % 2 == 0) return p == 2;
    for (uint64_t d = 3; d * d <= p; d += 2)
        if (p % d == 0) return false;
    return true;
}

uint64_t gfp_inv(uint64_t a, uint64_t p) {
    /* Extended Euclid; assumes gcd(a, p) == 1 (p prime, a % p != 0). */
    int64_t t = 0, newt = 1;
    int64_t r = (int64_t)p, newr = (int64_t)(a % p);
    while (newr != 0) {
        int64_t q = r / newr;
        int64_t tmp;
        tmp = t - q * newt; t = newt; newt = tmp;
        tmp = r - q * newr; r = newr; newr = tmp;
    }
    if (t < 0) t += (int64_t)p;
    return (uint64_t)t;
}

/* ------------------------------------------------------------------ */
/*  Construction / storage                                             */
/* ------------------------------------------------------------------ */

static inline int* gfp_exp_at(const GFpPoly* a, size_t i) {
    return a->exps + i * (size_t)a->n_vars;
}

GFpPoly* gfp_poly_new(int n_vars, GBOrder order, uint64_t p) {
    GFpPoly* a = (GFpPoly*)calloc(1, sizeof(GFpPoly));
    a->n_vars = n_vars;
    a->order = order;
    a->p = p;
    return a;
}

static void gfp_poly_reserve(GFpPoly* a, size_t cap) {
    if (cap <= a->cap) return;
    a->exps  = (int*)realloc(a->exps, sizeof(int) * cap * (size_t)(a->n_vars > 0 ? a->n_vars : 1));
    a->coefs = (uint64_t*)realloc(a->coefs, sizeof(uint64_t) * cap);
    a->cap = cap;
}

GFpPoly* gfp_poly_copy(const GFpPoly* a) {
    GFpPoly* b = gfp_poly_new(a->n_vars, a->order, a->p);
    if (a->n_terms) {
        gfp_poly_reserve(b, a->n_terms);
        memcpy(b->exps, a->exps, sizeof(int) * a->n_terms * (size_t)(a->n_vars > 0 ? a->n_vars : 1));
        memcpy(b->coefs, a->coefs, sizeof(uint64_t) * a->n_terms);
        b->n_terms = a->n_terms;
    }
    return b;
}

void gfp_poly_free(GFpPoly* a) {
    if (!a) return;
    free(a->exps);
    free(a->coefs);
    free(a);
}

void gfp_poly_push_term(GFpPoly* a, const int* exps, uint64_t coef) {
    coef %= a->p;
    if (coef == 0) return;
    if (a->n_terms == a->cap) gfp_poly_reserve(a, a->cap ? a->cap * 2 : 8);
    if (a->n_vars > 0)
        memcpy(gfp_exp_at(a, a->n_terms), exps, sizeof(int) * (size_t)a->n_vars);
    a->coefs[a->n_terms] = coef;
    a->n_terms++;
}

/* ------------------------------------------------------------------ */
/*  Monomial order                                                     */
/* ------------------------------------------------------------------ */

/* Compare exponent vectors a, b under `order`: >0 if a is the larger
 * monomial, <0 if smaller, 0 if equal. */
static int gfp_cmp(const int* a, const int* b, int n, GBOrder order) {
    if (order == GB_ORDER_GREVLEX) {
        int da = 0, db = 0;
        for (int i = 0; i < n; i++) { da += a[i]; db += b[i]; }
        if (da != db) return da - db;
        /* reverse lex: last differing variable, smaller exp is larger */
        for (int i = n - 1; i >= 0; i--)
            if (a[i] != b[i]) return b[i] - a[i];
        return 0;
    }
    /* Lexicographic (default). */
    for (int i = 0; i < n; i++)
        if (a[i] != b[i]) return a[i] - b[i];
    return 0;
}

/* qsort-with-context via a file-static pointer (mirrors groebner.c's
 * sort_cmp_idx).  Used only inside gfp_poly_normalize, never re-entrant. */
static const GFpPoly* g_sort_ctx = NULL;
static int sort_cmp_idx(const void* pa, const void* pb) {
    size_t ia = *(const size_t*)pa, ib = *(const size_t*)pb;
    const GFpPoly* a = g_sort_ctx;
    /* descending: larger monomial first */
    return -gfp_cmp(gfp_exp_at(a, ia), gfp_exp_at(a, ib), a->n_vars, a->order);
}

void gfp_poly_normalize(GFpPoly* a) {
    size_t nt = a->n_terms;
    if (nt <= 1) {
        if (nt == 1 && (a->coefs[0] %= a->p) == 0) a->n_terms = 0;
        return;
    }
    size_t* idx = (size_t*)malloc(sizeof(size_t) * nt);
    for (size_t i = 0; i < nt; i++) idx[i] = i;
    g_sort_ctx = a;
    qsort(idx, nt, sizeof(size_t), sort_cmp_idx);
    g_sort_ctx = NULL;

    int nv = a->n_vars > 0 ? a->n_vars : 1;
    int* nexp = (int*)malloc(sizeof(int) * nt * (size_t)nv);
    uint64_t* ncoef = (uint64_t*)malloc(sizeof(uint64_t) * nt);
    size_t out = 0;
    for (size_t k = 0; k < nt; k++) {
        const int* e = gfp_exp_at(a, idx[k]);
        uint64_t c = a->coefs[idx[k]] % a->p;
        if (out > 0 && a->n_vars > 0
            && memcmp(nexp + (out - 1) * (size_t)a->n_vars, e,
                      sizeof(int) * (size_t)a->n_vars) == 0) {
            ncoef[out - 1] = gfp_add(ncoef[out - 1], c, a->p);
            if (ncoef[out - 1] == 0) out--;   /* terms cancelled */
        } else if (out > 0 && a->n_vars == 0) {
            ncoef[out - 1] = gfp_add(ncoef[out - 1], c, a->p);
            if (ncoef[out - 1] == 0) out--;
        } else {
            if (c == 0) continue;
            if (a->n_vars > 0)
                memcpy(nexp + out * (size_t)a->n_vars, e, sizeof(int) * (size_t)a->n_vars);
            ncoef[out] = c;
            out++;
        }
    }
    free(idx);
    free(a->exps);
    free(a->coefs);
    a->exps = nexp;
    a->coefs = ncoef;
    a->n_terms = out;
    a->cap = nt;
}

/* ------------------------------------------------------------------ */
/*  Queries                                                            */
/* ------------------------------------------------------------------ */

bool gfp_poly_is_zero(const GFpPoly* a) { return a->n_terms == 0; }

bool gfp_poly_is_constant(const GFpPoly* a) {
    if (a->n_terms == 0) return true;
    if (a->n_terms != 1) return false;
    const int* e = gfp_exp_at(a, 0);
    for (int i = 0; i < a->n_vars; i++) if (e[i] != 0) return false;
    return true;
}

const int* gfp_poly_lm(const GFpPoly* a) {
    return a->n_terms ? gfp_exp_at(a, 0) : NULL;
}
uint64_t gfp_poly_lc(const GFpPoly* a) {
    return a->n_terms ? a->coefs[0] : 0;
}

int gfp_degree_in(const GFpPoly* a, int var) {
    int d = 0;
    for (size_t t = 0; t < a->n_terms; t++) {
        int e = gfp_exp_at(a, t)[var];
        if (e > d) d = e;
    }
    return d;
}

bool gfp_contains_var(const GFpPoly* a, int var) {
    for (size_t t = 0; t < a->n_terms; t++)
        if (gfp_exp_at(a, t)[var] > 0) return true;
    return false;
}

/* ------------------------------------------------------------------ */
/*  Arithmetic                                                         */
/* ------------------------------------------------------------------ */

static void gfp_poly_make_monic(GFpPoly* a) {
    if (a->n_terms == 0) return;
    uint64_t inv = gfp_inv(a->coefs[0], a->p);
    if (inv == 1) return;
    for (size_t t = 0; t < a->n_terms; t++)
        a->coefs[t] = gfp_mul(a->coefs[t], inv, a->p);
}

/* result = a - x^add_exp * c * b  (leading-term-cancelling subtraction). */
static GFpPoly* gfp_poly_sub_mul(const GFpPoly* a, const int* add_exp,
                                 uint64_t c, const GFpPoly* b) {
    GFpPoly* r = gfp_poly_new(a->n_vars, a->order, a->p);
    gfp_poly_reserve(r, a->n_terms + b->n_terms);
    for (size_t t = 0; t < a->n_terms; t++)
        gfp_poly_push_term(r, gfp_exp_at(a, t), a->coefs[t]);
    int nv = a->n_vars;
    int* e = (int*)malloc(sizeof(int) * (size_t)(nv > 0 ? nv : 1));
    for (size_t t = 0; t < b->n_terms; t++) {
        const int* be = gfp_exp_at(b, t);
        for (int i = 0; i < nv; i++) e[i] = be[i] + add_exp[i];
        uint64_t nc = gfp_sub(0, gfp_mul(c, b->coefs[t], a->p), a->p);  /* negate */
        gfp_poly_push_term(r, e, nc);
    }
    free(e);
    gfp_poly_normalize(r);
    return r;
}

/* ------------------------------------------------------------------ */
/*  S-polynomial and reduction                                         */
/* ------------------------------------------------------------------ */

static GFpPoly* gfp_spoly(const GFpPoly* f, const GFpPoly* g) {
    int n = f->n_vars;
    const int* lmf = gfp_poly_lm(f);
    const int* lmg = gfp_poly_lm(g);
    int* lcm = (int*)malloc(sizeof(int) * (size_t)(n > 0 ? n : 1));
    int* af  = (int*)malloc(sizeof(int) * (size_t)(n > 0 ? n : 1));
    int* ag  = (int*)malloc(sizeof(int) * (size_t)(n > 0 ? n : 1));
    for (int i = 0; i < n; i++) {
        lcm[i] = lmf[i] > lmg[i] ? lmf[i] : lmg[i];
        af[i] = lcm[i] - lmf[i];
        ag[i] = lcm[i] - lmg[i];
    }
    /* S = (1/lc f) x^af f - (1/lc g) x^ag g.  Build (1/lc f) x^af f, then
     * subtract (1/lc g) x^ag g via gfp_poly_sub_mul. */
    uint64_t cf = gfp_inv(f->coefs[0], f->p);
    uint64_t cg = gfp_inv(g->coefs[0], g->p);
    GFpPoly* tf = gfp_poly_new(n, f->order, f->p);
    gfp_poly_reserve(tf, f->n_terms);
    for (size_t t = 0; t < f->n_terms; t++) {
        const int* fe = gfp_exp_at(f, t);
        for (int i = 0; i < n; i++) af[i] = (lcm[i] - lmf[i]) + fe[i];
        gfp_poly_push_term(tf, af, gfp_mul(cf, f->coefs[t], f->p));
    }
    gfp_poly_normalize(tf);
    /* ag holds the shift for g */
    for (int i = 0; i < n; i++) ag[i] = lcm[i] - lmg[i];
    GFpPoly* s = gfp_poly_sub_mul(tf, ag, cg, g);
    gfp_poly_free(tf);
    free(lcm); free(af); free(ag);
    return s;
}

static bool exp_divides(const int* b, const int* a, int* q, int n) {
    for (int i = 0; i < n; i++) {
        if (a[i] < b[i]) return false;
        q[i] = a[i] - b[i];
    }
    return true;
}

/* Full multivariate division: normal form of `p` modulo the basis. */
static GFpPoly* gfp_reduce(const GFpPoly* poly, GFpPoly* const* basis, size_t n) {
    GFpPoly* r = gfp_poly_copy(poly);
    if (r->n_terms == 0) return r;
    int nv = r->n_vars;
    int* q = (int*)malloc(sizeof(int) * (size_t)(nv > 0 ? nv : 1));
    bool reduced;
    do {
        reduced = false;
        for (size_t t = 0; t < r->n_terms; t++) {
            const int* re = gfp_exp_at(r, t);
            for (size_t bi = 0; bi < n; bi++) {
                const GFpPoly* g = basis[bi];
                if (g->n_terms == 0) continue;
                if (!exp_divides(gfp_poly_lm(g), re, q, nv)) continue;
                /* factor = r.coef[t] / lc(g) */
                uint64_t factor = gfp_mul(r->coefs[t], gfp_inv(g->coefs[0], g->p), g->p);
                GFpPoly* nr = gfp_poly_sub_mul(r, q, factor, g);
                gfp_poly_free(r);
                r = nr;
                reduced = true;
                goto restart;
            }
        }
restart: ;
    } while (reduced);
    free(q);
    return r;
}

/* ------------------------------------------------------------------ */
/*  Buchberger                                                         */
/* ------------------------------------------------------------------ */

static bool lm_coprime(const GFpPoly* f, const GFpPoly* g) {
    const int* lf = gfp_poly_lm(f);
    const int* lg = gfp_poly_lm(g);
    for (int i = 0; i < f->n_vars; i++)
        if (lf[i] > 0 && lg[i] > 0) return false;
    return true;
}

static void gfp_finalize_basis(GFpPoly** G, size_t* nG_io) {
    size_t nG = *nG_io;
    if (nG == 0) { *nG_io = 0; return; }
    int nv = G[0]->n_vars;

    /* Discard any g_i whose LM is divisible by another g_k's LM. */
    bool* keep = (bool*)malloc(sizeof(bool) * nG);
    for (size_t i = 0; i < nG; i++) keep[i] = true;
    for (size_t i = 0; i < nG; i++) {
        if (!keep[i]) continue;
        const int* lmi = gfp_poly_lm(G[i]);
        for (size_t k = 0; k < nG; k++) {
            if (k == i || !keep[k]) continue;
            const int* lmk = gfp_poly_lm(G[k]);
            bool div_ok = true;
            for (int v = 0; v < nv; v++) if (lmi[v] < lmk[v]) { div_ok = false; break; }
            if (div_ok) { keep[i] = false; break; }
        }
    }
    size_t out = 0;
    for (size_t i = 0; i < nG; i++) {
        if (keep[i]) G[out++] = G[i];
        else gfp_poly_free(G[i]);
    }
    nG = out;
    free(keep);

    /* Reduce each survivor by the others to a fixed point, keep monic. */
    bool changed;
    do {
        changed = false;
        for (size_t i = 0; i < nG; i++) {
            GFpPoly** other = (GFpPoly**)malloc(sizeof(GFpPoly*) * (nG ? nG : 1));
            size_t m = 0;
            for (size_t k = 0; k < nG; k++) if (k != i) other[m++] = G[k];
            GFpPoly* nr = gfp_reduce(G[i], other, m);
            gfp_poly_make_monic(nr);
            free(other);
            if (nr->n_terms == 0) {                 /* collapsed; drop */
                gfp_poly_free(nr);
                gfp_poly_free(G[i]);
                for (size_t k = i; k + 1 < nG; k++) G[k] = G[k + 1];
                nG--; i--; changed = true; continue;
            }
            bool same = (nr->n_terms == G[i]->n_terms);
            for (size_t t = 0; same && t < nr->n_terms; t++) {
                if (nr->coefs[t] != G[i]->coefs[t]) same = false;
                else if (nv > 0 && memcmp(gfp_exp_at(nr, t), gfp_exp_at(G[i], t),
                                          sizeof(int) * (size_t)nv) != 0) same = false;
            }
            if (!same) { gfp_poly_free(G[i]); G[i] = nr; changed = true; }
            else gfp_poly_free(nr);
        }
    } while (changed);

    /* Sort ascending by leading monomial (Mathematica convention). */
    for (size_t i = 0; i + 1 < nG; i++) {
        size_t pick = i;
        for (size_t j = i + 1; j < nG; j++)
            if (gfp_cmp(gfp_poly_lm(G[j]), gfp_poly_lm(G[pick]), nv, G[0]->order) < 0)
                pick = j;
        if (pick != i) { GFpPoly* tmp = G[i]; G[i] = G[pick]; G[pick] = tmp; }
    }
    *nG_io = nG;
}

GFpPoly** gfp_buchberger(GFpPoly* const* F, size_t n, size_t* out_n) {
    GFpPoly** G = NULL;
    size_t nG = 0, capG = 0;
    for (size_t i = 0; i < n; i++) {
        if (F[i]->n_terms == 0) continue;
        if (nG == capG) { capG = capG ? capG * 2 : 8; G = (GFpPoly**)realloc(G, sizeof(GFpPoly*) * capG); }
        G[nG] = gfp_poly_copy(F[i]);
        gfp_poly_make_monic(G[nG]);
        nG++;
    }
    if (nG == 0) { *out_n = 0; free(G); return NULL; }

    /* Pair queue: all (i, j), i < j.  Simple normal-strategy queue with the
     * coprime-leading-term (Buchberger product) criterion for pruning. */
    typedef struct { size_t i, j; } Pair;
    Pair* pairs = NULL; size_t nP = 0, capP = 0;
    for (size_t j = 1; j < nG; j++)
        for (size_t i = 0; i < j; i++) {
            if (nP == capP) { capP = capP ? capP * 2 : 16; pairs = (Pair*)realloc(pairs, sizeof(Pair) * capP); }
            pairs[nP].i = i; pairs[nP].j = j; nP++;
        }

    size_t head = 0;
    while (head < nP) {
        tc_check_deadline();
        Pair pr = pairs[head++];
        GFpPoly* gi = G[pr.i];
        GFpPoly* gj = G[pr.j];
        if (lm_coprime(gi, gj)) continue;            /* product criterion */
        GFpPoly* s = gfp_spoly(gi, gj);
        GFpPoly* r = gfp_reduce(s, G, nG);
        gfp_poly_free(s);
        if (r->n_terms == 0) { gfp_poly_free(r); continue; }
        gfp_poly_make_monic(r);
        if (nG == capG) { capG = capG ? capG * 2 : 8; G = (GFpPoly**)realloc(G, sizeof(GFpPoly*) * capG); }
        G[nG] = r;
        for (size_t i = 0; i < nG; i++) {
            if (nP == capP) { capP = capP ? capP * 2 : 16; pairs = (Pair*)realloc(pairs, sizeof(Pair) * capP); }
            pairs[nP].i = i; pairs[nP].j = nG; nP++;
        }
        nG++;
    }
    free(pairs);

    gfp_finalize_basis(G, &nG);
    *out_n = nG;
    return G;
}

void gfp_basis_free(GFpPoly** G, size_t n) {
    if (!G) return;
    for (size_t i = 0; i < n; i++) gfp_poly_free(G[i]);
    free(G);
}

/* ------------------------------------------------------------------ */
/*  Conversion from the rational engine                                */
/* ------------------------------------------------------------------ */

GFpPoly* gfp_from_gbpoly(const GBPoly* g, GBOrder order, uint64_t p) {
    GFpPoly* a = gfp_poly_new(g->n_vars, order, p);
    for (size_t t = 0; t < g->n_terms; t++) {
        unsigned long num = mpz_fdiv_ui(mpq_numref(g->coefs[t]), (unsigned long)p);
        unsigned long den = mpz_fdiv_ui(mpq_denref(g->coefs[t]), (unsigned long)p);
        if (den == 0) { gfp_poly_free(a); return NULL; }   /* pole mod p */
        uint64_t coef = gfp_mul((uint64_t)num, gfp_inv((uint64_t)den, p), p);
        gfp_poly_push_term(a, g->exps + t * (size_t)(g->n_vars > 0 ? g->n_vars : 1), coef);
    }
    gfp_poly_normalize(a);
    return a;
}

GBPoly* gbpoly_from_gfp(const GFpPoly* a) {
    GBPoly* g = gb_poly_new(a->n_vars, a->order, 0);
    for (size_t t = 0; t < a->n_terms; t++)
        gb_poly_push_term_si(g, gfp_exp_at(a, t), (int64_t)a->coefs[t], 1);
    gb_poly_normalize(g);
    return g;
}

/* ------------------------------------------------------------------ */
/*  Solving support                                                    */
/* ------------------------------------------------------------------ */

GFpPoly* gfp_poly_subst(const GFpPoly* a, int var, uint64_t value) {
    GFpPoly* r = gfp_poly_new(a->n_vars, a->order, a->p);
    gfp_poly_reserve(r, a->n_terms);
    int nv = a->n_vars;
    int* e = (int*)malloc(sizeof(int) * (size_t)(nv > 0 ? nv : 1));
    for (size_t t = 0; t < a->n_terms; t++) {
        const int* ae = gfp_exp_at(a, t);
        int deg = ae[var];
        uint64_t c = a->coefs[t];
        /* multiply coefficient by value^deg mod p */
        uint64_t vp = 1, base = value % a->p;
        for (int k = 0; k < deg; k++) vp = gfp_mul(vp, base, a->p);
        c = gfp_mul(c, vp, a->p);
        for (int i = 0; i < nv; i++) e[i] = ae[i];
        e[var] = 0;
        gfp_poly_push_term(r, e, c);
    }
    free(e);
    gfp_poly_normalize(r);
    return r;
}

bool gfp_univariate_in(const GFpPoly* a, int var) {
    for (size_t t = 0; t < a->n_terms; t++) {
        const int* e = gfp_exp_at(a, t);
        for (int i = 0; i < a->n_vars; i++)
            if (i != var && e[i] != 0) return false;
    }
    return true;
}

uint64_t gfp_eval_univariate(const GFpPoly* a, int var, uint64_t r) {
    uint64_t acc = 0, rr = r % a->p;
    for (size_t t = 0; t < a->n_terms; t++) {
        int deg = gfp_exp_at(a, t)[var];
        uint64_t vp = 1;
        for (int k = 0; k < deg; k++) vp = gfp_mul(vp, rr, a->p);
        acc = gfp_add(acc, gfp_mul(a->coefs[t], vp, a->p), a->p);
    }
    return acc;
}

uint64_t gfp_eval_full(const GFpPoly* a, const uint64_t* values) {
    uint64_t acc = 0;
    for (size_t t = 0; t < a->n_terms; t++) {
        const int* e = gfp_exp_at(a, t);
        uint64_t term = a->coefs[t];
        for (int i = 0; i < a->n_vars; i++) {
            uint64_t base = values[i] % a->p;
            for (int k = 0; k < e[i]; k++) term = gfp_mul(term, base, a->p);
        }
        acc = gfp_add(acc, term, a->p);
    }
    return acc;
}

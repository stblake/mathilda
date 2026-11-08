/*
 * gb_f4.c -- Faugere's F4 over GF(p), word-size p.
 *
 * Buchberger reduces one S-polynomial at a time by heap division; F4 reduces
 * all the S-polynomials of one degree together, as the rows of one matrix,
 * by linear algebra over GF(p).  The steps, per round:
 *
 *  1. SELECTION (normal strategy): every pair of the smallest degree --
 *     deg lcm in a degree order, the sugar degree in Lex.  Pairs are kept by
 *     the Gebauer-Moeller criteria (chain criterion and product criterion on
 *     the new pairs, criterion B on the old ones), as gb_flint.c does.
 *  2. SYMBOLIC PREPROCESSING: the two halves m_i f_i, m_j f_j of each pair
 *     are rows; then, until closure, every monomial of a row that some
 *     leading monomial divides gets a REDUCER row u g (the shortest such g;
 *     every element ever added may reduce, as in gb_flint.c).  A row's
 *     coefficients are its polynomial's: rows share them, only the column
 *     indices are built.
 *  3. The MATRIX: columns are the monomials in decreasing order; the reducer
 *     rows are the pivot rows (distinct leading columns, monic).  Each
 *     S-row is loaded into a dense accumulator and reduced, column by column
 *     from the left, by the pivot rows -- and by the S-rows reduced before it,
 *     which become pivots themselves (the GBLA structure: the known-pivot
 *     block first, then the echelon form of what is left).  The products are
 *     accumulated lazily in two words per column (unsigned 128-bit sums of
 *     products below p^2) and reduced when the column is reached, or when
 *     the sum nears 2^128 -- about one reduction per entry, not one per
 *     product.
 *  4. Every nonzero reduced S-row has a leading monomial no basis element
 *     divides: it is a new basis element (the pairs updated).
 *
 * The output is the reduced basis: the minimal one, each element's tail
 * reduced by the smaller ones (f4_interreduce).
 *
 * TimeConstrained: tc_deadline_passed is polled every round and every 256
 * rows; on expiry everything is freed and *expired returned.
 */
#include "gb_f4.h"

#ifdef USE_FLINT

#include <stdlib.h>
#include <string.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#include "core.h"           /* tc_deadline_passed */

static int f4_on = -1;
int gb_f4_enable(int on) {
    int was = f4_on < 0 ? 1 : f4_on;
    f4_on = on;
    return was;
}
/* GroebnerBasis[..., Method -> "Buchberger"]: > 0 while a call runs that must
 * not use F4.  A nesting counter (not a flag) registered with the
 * TimeConstrained unwind registry, so a timeout restores it. */
static int f4_off_depth = 0;
int* gb_f4_off_counter(void) { return &f4_off_depth; }
/* Bases F4 has produced: Method -> "F4" compares it across a call. */
static long f4_hits = 0;
long gb_f4_hits(void) { return f4_hits; }
static int f4_enabled(void) {
    if (f4_on < 0) {
        const char* e = getenv("MATHILDA_GB_F4");
        f4_on = !(e && strcmp(e, "off") == 0);
    }
    return f4_on && f4_off_depth == 0;
}

/* ------------------------------------------------------------------ */
/*  Monomials: a hash table of exponent vectors                        */
/* ------------------------------------------------------------------ */

typedef struct {
    slong nv;
    int ord;
    ulong* ex;          /* n * nv */
    ulong* deg;         /* total degree */
    ulong* mask;        /* divisibility mask: bit v set when exponent of v > 0 */
    ulong* hash;
    slong n, cap;
    slong* tab; slong tcap;     /* open addressing, id + 1 */
} F4Mons;

static ulong f4_hash(const ulong* e, slong nv) {
    ulong h = UWORD(1469598103934665603);
    for (slong v = 0; v < nv; v++) h = (h ^ e[v]) * UWORD(1099511628211);
    return h ^ (h >> 29);
}

static void fm_init(F4Mons* M, slong nv, int ord) {
    memset(M, 0, sizeof *M);
    M->nv = nv; M->ord = ord;
    M->tcap = 1024;
    M->tab = flint_calloc((size_t)M->tcap, sizeof(slong));
}
static void fm_clear(F4Mons* M) {
    flint_free(M->ex); flint_free(M->deg); flint_free(M->mask); flint_free(M->hash); flint_free(M->tab);
}

static ulong fm_mask(const ulong* e, slong nv) {
    ulong m = 0;
    for (slong v = 0; v < nv; v++) if (e[v]) m |= UWORD(1) << (v % FLINT_BITS);
    return m;
}

/* The id of exponent vector e, inserted when new. */
static slong fm_id(F4Mons* M, const ulong* e) {
    slong nv = M->nv;
    ulong h = f4_hash(e, nv);
    ulong at = h & (ulong)(M->tcap - 1);
    while (M->tab[at]) {
        slong k = M->tab[at] - 1;
        if (M->hash[k] == h && memcmp(M->ex + k * nv, e, (size_t)nv * sizeof(ulong)) == 0) return k;
        at = (at + 1) & (ulong)(M->tcap - 1);
    }
    if (M->n == M->cap) {
        M->cap = M->cap ? 2 * M->cap : 1024;
        M->ex = flint_realloc(M->ex, (size_t)(M->cap * nv) * sizeof(ulong));
        M->deg = flint_realloc(M->deg, (size_t)M->cap * sizeof(ulong));
        M->mask = flint_realloc(M->mask, (size_t)M->cap * sizeof(ulong));
        M->hash = flint_realloc(M->hash, (size_t)M->cap * sizeof(ulong));
    }
    slong k = M->n++;
    memcpy(M->ex + k * nv, e, (size_t)nv * sizeof(ulong));
    ulong d = 0;
    for (slong v = 0; v < nv; v++) d += e[v];
    M->deg[k] = d; M->mask[k] = fm_mask(e, nv); M->hash[k] = h;
    M->tab[at] = k + 1;
    if (2 * M->n > M->tcap) {                         /* rehash */
        flint_free(M->tab);
        M->tcap *= 2;
        M->tab = flint_calloc((size_t)M->tcap, sizeof(slong));
        for (slong i = 0; i < M->n; i++) {
            ulong a = M->hash[i] & (ulong)(M->tcap - 1);
            while (M->tab[a]) a = (a + 1) & (ulong)(M->tcap - 1);
            M->tab[a] = i + 1;
        }
    }
    return k;
}

/* The id of the product of monomials a (an id) and the exponent vector u. */
static slong fm_mul(F4Mons* M, slong a, const ulong* u, ulong* tmp) {
    const ulong* ea = M->ex + a * M->nv;
    for (slong v = 0; v < M->nv; v++) tmp[v] = ea[v] + u[v];
    return fm_id(M, tmp);
}

/* < 0 when monomial a is the larger (FLINT's descending term order). */
static int fm_cmp(const F4Mons* M, slong a, slong b) {
    if (a == b) return 0;
    slong nv = M->nv;
    if (M->ord != ORD_LEX && M->deg[a] != M->deg[b]) return M->deg[a] > M->deg[b] ? -1 : 1;
    const ulong* ea = M->ex + a * nv; const ulong* eb = M->ex + b * nv;
    if (M->ord == ORD_DEGREVLEX) {
        for (slong v = nv - 1; v >= 0; v--)
            if (ea[v] != eb[v]) return ea[v] < eb[v] ? -1 : 1;
        return 0;
    }
    for (slong v = 0; v < nv; v++)
        if (ea[v] != eb[v]) return ea[v] > eb[v] ? -1 : 1;
    return 0;
}

static int fm_divides(const F4Mons* M, slong a, slong b) {      /* a | b */
    if (M->mask[a] & ~M->mask[b]) return 0;
    const ulong* ea = M->ex + a * M->nv; const ulong* eb = M->ex + b * M->nv;
    for (slong v = 0; v < M->nv; v++) if (ea[v] > eb[v]) return 0;
    return 1;
}

/* Sort ids descending in the order (merge sort: qsort has no context). */
static void fm_sort(const F4Mons* M, slong* a, slong n, slong* tmp) {
    if (n < 2) return;
    if (n < 16) {
        for (slong i = 1; i < n; i++) {
            slong x = a[i], j = i;
            while (j > 0 && fm_cmp(M, a[j - 1], x) > 0) { a[j] = a[j - 1]; j--; }
            a[j] = x;
        }
        return;
    }
    slong h = n / 2;
    fm_sort(M, a, h, tmp);
    fm_sort(M, a + h, n - h, tmp);
    slong i = 0, j = h, k = 0;
    while (i < h && j < n) tmp[k++] = fm_cmp(M, a[i], a[j]) <= 0 ? a[i++] : a[j++];
    while (i < h) tmp[k++] = a[i++];
    while (j < n) tmp[k++] = a[j++];
    memcpy(a, tmp, (size_t)n * sizeof(slong));
}

/* ------------------------------------------------------------------ */
/*  The basis and the pairs                                            */
/* ------------------------------------------------------------------ */

typedef struct { slong len; slong* mon; ulong* cf; } F4Poly;     /* monic, descending */

typedef struct { slong i, j, lcm; ulong deg; } F4Pair;

typedef struct {
    F4Mons M;
    nmod_t mod;
    F4Poly* B; slong nb, cb;
    ulong* sugar;
    char* active;           /* in the minimal set the pairs are formed from */
    F4Pair* P; slong np, cp;
    ulong* tmp;             /* nv scratch */
} F4;

static void f4_poly_clear(F4Poly* f) { flint_free(f->mon); flint_free(f->cf); }

static ulong f4_pair_deg(const F4* S, slong i, slong j, slong l) {
    if (S->M.ord != ORD_LEX) return S->M.deg[l];
    /* sugar */
    ulong dl = S->M.deg[l];
    ulong si = S->sugar[i] + dl - S->M.deg[S->B[i].mon[0]];
    ulong sj = S->sugar[j] + dl - S->M.deg[S->B[j].mon[0]];
    return si > sj ? si : sj;
}

static slong f4_lcm(F4* S, slong a, slong b) {
    const ulong* ea = S->M.ex + a * S->M.nv; const ulong* eb = S->M.ex + b * S->M.nv;
    for (slong v = 0; v < S->M.nv; v++) S->tmp[v] = ea[v] > eb[v] ? ea[v] : eb[v];
    return fm_id(&S->M, S->tmp);
}

static int f4_coprime(const F4* S, slong a, slong b) {
    if (!(S->M.mask[a] & S->M.mask[b])) return 1;
    const ulong* ea = S->M.ex + a * S->M.nv; const ulong* eb = S->M.ex + b * S->M.nv;
    for (slong v = 0; v < S->M.nv; v++) if (ea[v] && eb[v]) return 0;
    return 1;
}

/* Register the new element h (Gebauer-Moeller). */
static void f4_update(F4* S, slong h) {
    slong t = S->B[h].mon[0];
    slong nc = 0;
    slong* cg = flint_malloc((size_t)(h + 1) * sizeof(slong));
    slong* cl = flint_malloc((size_t)(h + 1) * sizeof(slong));
    char* keep = flint_malloc((size_t)(h + 1));
    for (slong g = 0; g < h; g++) {
        if (!S->active[g]) continue;
        cg[nc] = g; cl[nc] = f4_lcm(S, S->B[g].mon[0], t); keep[nc] = 1; nc++;
    }
    /* the chain criterion on the new pairs (of equal lcms one is kept) */
    for (slong a = 0; a < nc; a++)
        for (slong b = 0; b < nc && keep[a]; b++) {
            if (a == b || !keep[b]) continue;
            if (!fm_divides(&S->M, cl[b], cl[a])) continue;
            if (cl[a] != cl[b] || b < a) keep[a] = 0;
        }
    /* the product criterion: an equal-lcm group with a coprime member goes */
    for (slong a = 0; a < nc; a++) {
        if (!keep[a]) continue;
        int cop = f4_coprime(S, S->B[cg[a]].mon[0], t);
        for (slong b = 0; b < nc && !cop; b++)
            if (b != a && cl[b] == cl[a] && f4_coprime(S, S->B[cg[b]].mon[0], t)) cop = 1;
        if (cop) keep[a] = 0;
    }
    /* criterion B on the old pairs */
    slong w = 0;
    for (slong k = 0; k < S->np; k++) {
        F4Pair* q = S->P + k;
        int dead = 0;
        if (fm_divides(&S->M, t, q->lcm)) {
            slong l1 = f4_lcm(S, S->B[q->i].mon[0], t);
            slong l2 = f4_lcm(S, S->B[q->j].mon[0], t);
            dead = l1 != q->lcm && l2 != q->lcm;
        }
        if (!dead) S->P[w++] = *q;
    }
    S->np = w;
    for (slong a = 0; a < nc; a++) {
        if (!keep[a]) continue;
        if (S->np == S->cp) {
            S->cp = S->cp ? 2 * S->cp : 64;
            S->P = flint_realloc(S->P, (size_t)S->cp * sizeof(F4Pair));
        }
        F4Pair* q = S->P + S->np++;
        q->i = cg[a]; q->j = h; q->lcm = cl[a]; q->deg = f4_pair_deg(S, cg[a], h, cl[a]);
    }
    flint_free(cg); flint_free(cl); flint_free(keep);
    for (slong g = 0; g < h; g++)
        if (S->active[g] && fm_divides(&S->M, t, S->B[g].mon[0])) S->active[g] = 0;
    S->active[h] = 1;
}

/* Append f (taken) with the given sugar and update the pairs. */
static void f4_add(F4* S, F4Poly* f, ulong sugar) {
    if (S->nb == S->cb) {
        S->cb = S->cb ? 2 * S->cb : 32;
        S->B = flint_realloc(S->B, (size_t)S->cb * sizeof(F4Poly));
        S->sugar = flint_realloc(S->sugar, (size_t)S->cb * sizeof(ulong));
        S->active = flint_realloc(S->active, (size_t)S->cb);
    }
    slong h = S->nb++;
    S->B[h] = *f;
    S->sugar[h] = sugar;
    S->active[h] = 0;
    f4_update(S, h);
}

/* ------------------------------------------------------------------ */
/*  One round: rows, symbolic preprocessing, the matrix                */
/* ------------------------------------------------------------------ */

typedef struct {
    slong poly;         /* basis element */
    slong mult;         /* monomial id of the multiplier */
    slong* col;         /* column of each term (after numbering) */
    slong* mon;         /* monomial id of each term */
    ulong sugar;
} F4Row;

typedef struct {
    F4Row* r; slong n, c;
} F4Rows;

static void rows_push(F4Rows* R, slong poly, slong mult, ulong sugar) {
    if (R->n == R->c) {
        R->c = R->c ? 2 * R->c : 64;
        R->r = flint_realloc(R->r, (size_t)R->c * sizeof(F4Row));
    }
    F4Row* w = R->r + R->n++;
    w->poly = poly; w->mult = mult; w->col = NULL; w->mon = NULL; w->sugar = sugar;
}

/* The monomials of row w (multiplier times its polynomial). */
static void row_monomials(F4* S, F4Row* w) {
    const F4Poly* f = S->B + w->poly;
    w->mon = flint_malloc((size_t)f->len * sizeof(slong));
    /* a copy: fm_mul may grow (move) the exponent table */
    ulong* u = flint_malloc((size_t)S->M.nv * sizeof(ulong));
    memcpy(u, S->M.ex + w->mult * S->M.nv, (size_t)S->M.nv * sizeof(ulong));
    ulong* tmp = flint_malloc((size_t)S->M.nv * sizeof(ulong));
    for (slong t = 0; t < f->len; t++) w->mon[t] = fm_mul(&S->M, f->mon[t], u, tmp);
    flint_free(tmp); flint_free(u);
}

/* The reducer for monomial m: the shortest basis element whose leading
 * monomial divides m, -1 when none. */
static slong f4_reducer(const F4* S, slong m) {
    slong best = -1;
    for (slong k = 0; k < S->nb; k++) {
        if (!fm_divides(&S->M, S->B[k].mon[0], m)) continue;
        if (best < 0 || S->B[k].len < S->B[best].len) best = k;
    }
    return best;
}

/* One F4 round over the pairs sel[0..ns).  New elements are added to S.
 * Returns 1, 0 on the deadline. */
static int f4_round(F4* S, const F4Pair* sel, slong ns) {
    slong nv = S->M.nv;
    F4Rows piv = {0}, red = {0};
    /* the monomial -> pivot row table, by monomial id (grown as ids appear) */
    slong pcap = S->M.n + 1024;
    slong* pivot_of = flint_malloc((size_t)pcap * sizeof(slong));
    for (slong k = 0; k < pcap; k++) pivot_of[k] = -2;      /* -2: unseen */
    ulong* u = flint_malloc((size_t)nv * sizeof(ulong));
    int ok = 1;
#define SEEN_GROW(id) do { \
        if ((id) >= pcap) { \
            slong nc = 2 * pcap + (id); \
            pivot_of = flint_realloc(pivot_of, (size_t)nc * sizeof(slong)); \
            for (slong z = pcap; z < nc; z++) pivot_of[z] = -2; \
            pcap = nc; \
        } } while (0)
    /* every monomial of the matrix is a column: todo lists them as seen */
    slong* todo = NULL; slong ntodo = 0, ctodo = 0;
#define TODO_PUSH(id) do { \
        if (ntodo == ctodo) { ctodo = ctodo ? 2 * ctodo : 1024; \
            todo = flint_realloc(todo, (size_t)ctodo * sizeof(slong)); } \
        todo[ntodo++] = (id); } while (0)
    /* the halves of the pairs: the first to reach an lcm is its pivot row */
    for (slong k = 0; k < ns; k++) {
        const F4Pair* q = sel + k;
        slong ij[2] = {q->i, q->j};
        for (int h = 0; h < 2; h++) {
            slong g = ij[h];
            const ulong* el = S->M.ex + q->lcm * nv;
            const ulong* eg = S->M.ex + S->B[g].mon[0] * nv;
            for (slong v = 0; v < nv; v++) u[v] = el[v] - eg[v];
            slong mult = fm_id(&S->M, u);
            ulong sug = S->sugar[g] + S->M.deg[mult];
            SEEN_GROW(q->lcm);
            if (pivot_of[q->lcm] == -2) TODO_PUSH(q->lcm);
            if (pivot_of[q->lcm] < 0) {
                rows_push(&piv, g, mult, sug);
                pivot_of[q->lcm] = piv.n - 1;
            } else {
                /* a duplicate of a pivot or reducible row is skipped */
                int dup = 0;
                F4Row* pr = piv.r + pivot_of[q->lcm];
                if (pr->poly == g && pr->mult == mult) dup = 1;
                for (slong t = 0; t < red.n && !dup; t++)
                    dup = red.r[t].poly == g && red.r[t].mult == mult;
                if (!dup) rows_push(&red, g, mult, sug);
            }
        }
    }
    /* symbolic preprocessing: a reducer for every reducible monomial */
    for (int pass = 0; pass < 2; pass++) {
        F4Rows* R = pass ? &red : &piv;
        for (slong k = 0; k < R->n; k++) {
            row_monomials(S, R->r + k);
            const F4Poly* f = S->B + R->r[k].poly;
            for (slong t = 0; t < f->len; t++) {
                slong m = R->r[k].mon[t];
                SEEN_GROW(m);
                if (pivot_of[m] == -2) { pivot_of[m] = -1; TODO_PUSH(m); }
            }
        }
    }
    for (slong a = 0; a < ntodo; a++) {
        if ((a & 255) == 255 && tc_deadline_passed()) { ok = 0; break; }
        slong m = todo[a];
        if (pivot_of[m] >= 0) continue;
        slong g = f4_reducer(S, m);
        if (g < 0) continue;
        const ulong* em = S->M.ex + m * nv;
        const ulong* eg = S->M.ex + S->B[g].mon[0] * nv;
        for (slong v = 0; v < nv; v++) u[v] = em[v] - eg[v];
        slong mult = fm_id(&S->M, u);
        rows_push(&piv, g, mult, 0);
        pivot_of[m] = piv.n - 1;
        F4Row* w = piv.r + piv.n - 1;
        row_monomials(S, w);
        const F4Poly* f = S->B + g;
        for (slong t = 1; t < f->len; t++) {
            slong mm = w->mon[t];
            SEEN_GROW(mm);
            if (pivot_of[mm] == -2) { pivot_of[mm] = -1; TODO_PUSH(mm); }
        }
    }
    slong ncols = ntodo;
    slong* colmon = todo; todo = NULL;          /* the column monomials */
    slong* colof = NULL;                        /* monomial id -> column */
    ulong* acc = NULL;                          /* 2 words per column: lazy sums */
    slong* pivrow_of_col = NULL; ulong** pivcf = NULL; slong** pivcol = NULL; slong* pivlen = NULL;
    if (!ok) goto out;
    {
        slong* tmp = flint_malloc((size_t)(ncols ? ncols : 1) * sizeof(slong));
        fm_sort(&S->M, colmon, ncols, tmp);
        flint_free(tmp);
    }
    colof = flint_malloc((size_t)pcap * sizeof(slong));
    for (slong c = 0; c < ncols; c++) colof[colmon[c]] = c;
    /* pivots by column: the reducer rows, then reduced S-rows join */
    pivrow_of_col = flint_malloc((size_t)(ncols ? ncols : 1) * sizeof(slong));
    for (slong c = 0; c < ncols; c++) pivrow_of_col[c] = -1;
    slong npv = 0, cpv = piv.n + red.n + 1;
    pivcf = flint_malloc((size_t)cpv * sizeof(ulong*));
    pivcol = flint_malloc((size_t)cpv * sizeof(slong*));
    pivlen = flint_malloc((size_t)cpv * sizeof(slong));
    char* owned = flint_calloc((size_t)cpv, 1);     /* pivcf/pivcol allocated here */
    for (slong k = 0; k < piv.n; k++) {
        F4Row* w = piv.r + k;
        const F4Poly* f = S->B + w->poly;
        w->col = flint_malloc((size_t)f->len * sizeof(slong));
        for (slong t = 0; t < f->len; t++) w->col[t] = colof[w->mon[t]];
        pivcf[npv] = f->cf; pivcol[npv] = w->col; pivlen[npv] = f->len;
        pivrow_of_col[w->col[0]] = npv++;
    }
    /* reduce the S-rows, in order of their leading column */
    acc = flint_calloc((size_t)(2 * (ncols ? ncols : 1)), sizeof(ulong));
    ulong p = S->mod.n;
    int small = p < (UWORD(1) << 32);
    ulong p2 = small ? p * p : 0;
    for (slong k = 0; k < red.n && ok; k++) {
        if ((k & 255) == 255 && tc_deadline_passed()) { ok = 0; break; }
        F4Row* w = red.r + k;
        const F4Poly* f = S->B + w->poly;
        slong first = ncols, last = -1;
        for (slong t = 0; t < f->len; t++) {
            slong c = colof[w->mon[t]];
            if (!small) { acc[2 * c] = f->cf[t]; acc[2 * c + 1] = 0; }
            if (c < first) first = c;
            if (c > last) last = c;
        }
        /* column by column: reduce the lazy sum, eliminate by its pivot */
        slong lead = -1;
        slong nnz = 0;
        if (small) {
            /* p < 2^32: one word per column, kept in [0, p^2): a product
             * x c < p^2 is subtracted and p^2 added back on a borrow --
             * branch-free, no reduction until the column is reached */
            ulong* a1 = acc;
            for (slong t = 0; t < f->len; t++) {
                slong c = colof[w->mon[t]];
                a1[c] = f->cf[t];
            }
            for (slong c = first; c <= last; c++) {
                ulong v = a1[c];
                if (!v) continue;
                ulong x = v % p;
                a1[c] = x;
                if (!x) continue;
                slong pr = pivrow_of_col[c];
                if (pr < 0) { if (lead < 0) lead = c; nnz++; continue; }
                const ulong* pc = pivcf[pr]; const slong* pcol = pivcol[pr]; slong pl = pivlen[pr];
                a1[c] = 0;
                for (slong t = 1; t < pl; t++) {
                    slong cc = pcol[t];
                    ulong o = a1[cc], d = o - x * pc[t];
                    a1[cc] = d + (d > o ? p2 : 0);
                }
                if (pcol[pl - 1] > last) last = pcol[pl - 1];
            }
        } else {
        for (slong c = first; c <= last; c++) {
            ulong hi = acc[2 * c + 1], lo = acc[2 * c];
            if (!hi && !lo) continue;
            ulong x;
            NMOD2_RED2(x, hi, lo, S->mod);       /* the true entry mod p */
            acc[2 * c] = x; acc[2 * c + 1] = 0;
            if (!x) continue;
            slong pr = pivrow_of_col[c];
            if (pr < 0) { if (lead < 0) lead = c; nnz++; continue; }
            /* acc -= x * pivot: add (p - x) * pivot as unsigned 128-bit sums,
             * reduced only when the high word nears overflow */
            ulong a = p - x;
            const ulong* pc = pivcf[pr]; const slong* pcol = pivcol[pr]; slong pl = pivlen[pr];
            acc[2 * c] = 0;
            for (slong t = 1; t < pl; t++) {
                slong cc = pcol[t];
                ulong ph, plo;
                umul_ppmm(ph, plo, a, pc[t]);
                add_ssaaaa(acc[2 * cc + 1], acc[2 * cc], acc[2 * cc + 1], acc[2 * cc], ph, plo);
                if (acc[2 * cc + 1] >= (UWORD(1) << (FLINT_BITS - 2))) {
                    /* past 2^126: reduced now, so the next product (< 2^126)
                     * cannot carry out of the high word */
                    ulong r;
                    NMOD2_RED2(r, acc[2 * cc + 1], acc[2 * cc], S->mod);
                    acc[2 * cc] = r; acc[2 * cc + 1] = 0;
                }
            }
            if (pcol[pl - 1] > last) last = pcol[pl - 1];
        }
        }
        /* the entries: acc[c] (small p) or acc[2 c] (two words), reduced */
        slong st = small ? 1 : 2;
        if (lead < 0) {                       /* reduced to zero */
            for (slong c = first; c <= last; c++) { acc[st * c] = 0; acc[st * c + st - 1] = 0; }
            continue;
        }
        /* a new pivot: monic, sparse */
        ulong inv = n_invmod(acc[st * lead], p);
        ulong* ncf = flint_malloc((size_t)nnz * sizeof(ulong));
        slong* ncl = flint_malloc((size_t)nnz * sizeof(slong));
        slong q = 0;
        for (slong c = lead; c <= last; c++) {
            ulong x = acc[st * c];
            acc[st * c] = 0; acc[st * c + st - 1] = 0;
            if (!x) continue;
            ncf[q] = nmod_mul(x, inv, S->mod); ncl[q] = c; q++;
        }
        for (slong c = first; c < lead; c++) { acc[st * c] = 0; acc[st * c + st - 1] = 0; }
        pivcf[npv] = ncf; pivcol[npv] = ncl; pivlen[npv] = q; owned[npv] = 1;
        pivrow_of_col[lead] = npv++;
        /* the new basis element: its monomials back from the columns */
        F4Poly nf; nf.len = q;
        nf.mon = flint_malloc((size_t)q * sizeof(slong));
        nf.cf = flint_malloc((size_t)q * sizeof(ulong));
        for (slong t = 0; t < q; t++) { nf.mon[t] = colmon[ncl[t]]; nf.cf[t] = ncf[t]; }
        ulong sug = w->sugar;
        if (S->M.deg[nf.mon[0]] > sug) sug = S->M.deg[nf.mon[0]];
        f4_add(S, &nf, sug);
    }
    for (slong k = 0; k < npv; k++) if (owned[k]) { flint_free(pivcf[k]); flint_free(pivcol[k]); }
    flint_free(owned);
out:
    for (slong k = 0; k < piv.n; k++) { flint_free(piv.r[k].col); flint_free(piv.r[k].mon); }
    for (slong k = 0; k < red.n; k++) { flint_free(red.r[k].col); flint_free(red.r[k].mon); }
    flint_free(piv.r); flint_free(red.r);
    flint_free(pivot_of); flint_free(u); flint_free(todo); flint_free(colmon); flint_free(colof);
    flint_free(acc); flint_free(pivrow_of_col); flint_free(pivcf); flint_free(pivcol); flint_free(pivlen);
#undef SEEN_GROW
#undef TODO_PUSH
    return ok;
}

/* ------------------------------------------------------------------ */
/*  The reduced basis                                                  */
/* ------------------------------------------------------------------ */

/* A max-heap of monomial ids in the order. */
typedef struct { slong* a; slong n, c; } F4Heap;
static void hp_push(F4Heap* H, const F4Mons* M, slong x) {
    if (H->n == H->c) { H->c = H->c ? 2 * H->c : 256; H->a = flint_realloc(H->a, (size_t)H->c * sizeof(slong)); }
    slong i = H->n++;
    while (i > 0) {
        slong up = (i - 1) / 2;
        if (fm_cmp(M, H->a[up], x) <= 0) break;       /* parent larger: stop */
        H->a[i] = H->a[up]; i = up;
    }
    H->a[i] = x;
}
static slong hp_pop(F4Heap* H, const F4Mons* M) {
    slong top = H->a[0], x = H->a[--H->n], i = 0;
    for (;;) {
        slong l = 2 * i + 1, r = l + 1, b = i;
        slong bv = x;
        if (l < H->n && fm_cmp(M, H->a[l], bv) < 0) { b = l; bv = H->a[l]; }
        if (r < H->n && fm_cmp(M, H->a[r], bv) < 0) { b = r; bv = H->a[r]; }
        if (b == i) break;
        H->a[i] = H->a[b]; i = b;
    }
    if (H->n) H->a[i] = x;
    return top;
}

/* Tail-reduce the minimal basis sel[0..ns) (ascending by leading monomial)
 * in place: each element by the ones before it, which are already reduced,
 * so every reducer is final.  A heap of monomial ids and a coefficient per
 * id: no FLINT division (whose stack use does not suit worker threads). */
static int f4_interreduce(F4* S, const slong* sel, slong ns) {
    slong nv = S->M.nv;
    ulong* coef = NULL; char* inh = NULL; slong ccap = 0;
    F4Heap H = {0};
    ulong* u = flint_malloc((size_t)nv * sizeof(ulong));
    ulong* tmp = flint_malloc((size_t)nv * sizeof(ulong));
    int ok = 1;
#define GROWC(id) do { if ((id) >= ccap) { slong nc = 2 * (id) + 1024; \
        coef = flint_realloc(coef, (size_t)nc * sizeof(ulong)); inh = flint_realloc(inh, (size_t)nc); \
        memset(coef + ccap, 0, (size_t)(nc - ccap) * sizeof(ulong)); memset(inh + ccap, 0, (size_t)(nc - ccap)); \
        ccap = nc; } } while (0)
    for (slong a = 0; a < ns && ok; a++) {
        if ((a & 15) == 15 && tc_deadline_passed()) { ok = 0; break; }
        F4Poly* f = S->B + sel[a];
        slong rlen = 0, rcap = f->len;
        slong* rm = flint_malloc((size_t)rcap * sizeof(slong));
        ulong* rc = flint_malloc((size_t)rcap * sizeof(ulong));
        rm[rlen] = f->mon[0]; rc[rlen] = f->cf[0]; rlen++;
        for (slong t = 1; t < f->len; t++) {
            GROWC(f->mon[t]);
            coef[f->mon[t]] = f->cf[t];
            inh[f->mon[t]] = 1;
            hp_push(&H, &S->M, f->mon[t]);
        }
        while (H.n) {
            slong m = hp_pop(&H, &S->M);
            inh[m] = 0;
            ulong c = coef[m]; coef[m] = 0;
            if (!c) continue;
            /* a reducer among the earlier (reduced) elements: the shortest */
            slong g = -1;
            for (slong b = 0; b < a; b++) {
                const F4Poly* h = S->B + sel[b];
                if (fm_divides(&S->M, h->mon[0], m) && (g < 0 || h->len < S->B[g].len)) g = sel[b];
            }
            if (g < 0) {                                /* irreducible: in the result */
                if (rlen == rcap) {
                    rcap *= 2;
                    rm = flint_realloc(rm, (size_t)rcap * sizeof(slong));
                    rc = flint_realloc(rc, (size_t)rcap * sizeof(ulong));
                }
                rm[rlen] = m; rc[rlen] = c; rlen++;
                continue;
            }
            const F4Poly* h = S->B + g;
            const ulong* em = S->M.ex + m * nv; const ulong* eh = S->M.ex + h->mon[0] * nv;
            for (slong v = 0; v < nv; v++) u[v] = em[v] - eh[v];
            ulong nc = nmod_neg(c, S->mod);
            for (slong t = 1; t < h->len; t++) {
                slong mm = fm_mul(&S->M, h->mon[t], u, tmp);
                GROWC(mm);
                coef[mm] = nmod_add(coef[mm], nmod_mul(nc, h->cf[t], S->mod), S->mod);
                if (!inh[mm]) { inh[mm] = 1; hp_push(&H, &S->M, mm); }
            }
        }
        f4_poly_clear(f);
        f->len = rlen; f->mon = rm; f->cf = rc;
    }
    if (!ok) {                                          /* leave nothing behind */
        while (H.n) { slong m = hp_pop(&H, &S->M); coef[m] = 0; inh[m] = 0; }
    }
#undef GROWC
    flint_free(coef); flint_free(inh); flint_free(H.a); flint_free(u); flint_free(tmp);
    return ok;
}

/* ------------------------------------------------------------------ */
/*  Driver                                                             */
/* ------------------------------------------------------------------ */

int gb_f4_nmod(const nmod_mpoly_struct* F, slong n, const nmod_mpoly_ctx_t ctx,
               nmod_mpoly_struct** G, slong* ng, int* expired) {
    *expired = 0; *G = NULL; *ng = 0;
    if (!f4_enabled()) return 0;
    slong nv = ctx->minfo->nvars;
    if (nv < 1) return 0;
    int ord = (int)ctx->minfo->ord;
    for (slong i = 0; i < n; i++)
        if (F[i].length && nmod_mpoly_degrees_fit_si(F + i, ctx)) {
            slong* d = flint_malloc((size_t)nv * sizeof(slong));
            nmod_mpoly_degrees_si(d, F + i, ctx);
            int big = 0;
            for (slong v = 0; v < nv; v++) if (d[v] > (WORD(1) << 20)) big = 1;
            flint_free(d);
            if (big) return 0;
        } else if (F[i].length) return 0;
    F4 S; memset(&S, 0, sizeof S);
    fm_init(&S.M, nv, ord);
    S.mod = ctx->mod;
    S.tmp = flint_malloc((size_t)nv * sizeof(ulong));
    ulong* ex = flint_malloc((size_t)nv * sizeof(ulong));
    /* the inputs, monic; each is its own row set of degree 0 pairs: they are
     * added directly (the first rounds reduce them against each other) */
    for (slong i = 0; i < n; i++) {
        if (!F[i].length) continue;
        F4Poly f; f.len = F[i].length;
        f.mon = flint_malloc((size_t)f.len * sizeof(slong));
        f.cf = flint_malloc((size_t)f.len * sizeof(ulong));
        ulong inv = n_invmod(F[i].coeffs[0], ctx->mod.n);
        for (slong t = 0; t < f.len; t++) {
            nmod_mpoly_get_term_exp_ui(ex, F + i, t, ctx);
            f.mon[t] = fm_id(&S.M, ex);
            f.cf[t] = nmod_mul(F[i].coeffs[t], inv, ctx->mod);
        }
        /* a duplicate leading monomial is reduced by the matrix through a
         * pair (lcm = the monomial itself); nothing special here */
        ulong sug = 0;
        for (slong t = 0; t < f.len; t++) if (S.M.deg[f.mon[t]] > sug) sug = S.M.deg[f.mon[t]];
        f4_add(&S, &f, sug);
    }
    int ok = 1;
    /* a unit input: the ideal is (1) */
    while (S.np > 0) {
        if (tc_deadline_passed()) { ok = 0; break; }
        int unit = 0;
        for (slong k = 0; k < S.nb && !unit; k++) unit = S.active[k] && S.M.deg[S.B[k].mon[0]] == 0;
        if (unit) break;
        ulong dmin = S.P[0].deg;
        for (slong k = 1; k < S.np; k++) if (S.P[k].deg < dmin) dmin = S.P[k].deg;
        slong ns = 0;
        F4Pair* sel = flint_malloc((size_t)S.np * sizeof(F4Pair));
        slong w = 0;
        for (slong k = 0; k < S.np; k++) {
            if (S.P[k].deg == dmin) sel[ns++] = S.P[k];
            else S.P[w++] = S.P[k];
        }
        S.np = w;
        ok = f4_round(&S, sel, ns);
        flint_free(sel);
        if (!ok) break;
    }
    if (ok) {
        /* the minimal basis: active elements whose leading monomial no other
         * divides (the inputs were not interreduced; of equal ones the first
         * stays), a unit alone */
        for (slong k = 0; k < S.nb; k++) {
            if (!S.active[k]) continue;
            for (slong j = 0; j < S.nb && S.active[k]; j++) {
                if (j == k || !S.active[j]) continue;
                slong lj = S.B[j].mon[0], lk = S.B[k].mon[0];
                if (fm_divides(&S.M, lj, lk) && (lj != lk || j < k)) S.active[k] = 0;
            }
        }
        slong unit = -1;
        for (slong k = 0; k < S.nb && unit < 0; k++)
            if (S.active[k] && S.M.deg[S.B[k].mon[0]] == 0) unit = k;
        slong cnt = 0;
        slong* sel = flint_malloc((size_t)(S.nb ? S.nb : 1) * sizeof(slong));
        if (unit >= 0) {                                /* (1): the constant, monic */
            S.B[unit].len = 1;
            sel[cnt++] = unit;
        } else
            for (slong k = 0; k < S.nb; k++) if (S.active[k]) sel[cnt++] = k;
        /* ascending by leading monomial */
        for (slong i = 1; i < cnt; i++)
            for (slong j = i; j > 0 && fm_cmp(&S.M, S.B[sel[j - 1]].mon[0], S.B[sel[j]].mon[0]) < 0; j--) {
                slong t = sel[j]; sel[j] = sel[j - 1]; sel[j - 1] = t;
            }
        ok = f4_interreduce(&S, sel, cnt);
        slong q = 0;
        if (ok) {
            *G = flint_malloc((size_t)(cnt ? cnt : 1) * sizeof(nmod_mpoly_struct));
            for (slong a = 0; a < cnt; a++) {
                nmod_mpoly_struct* g = *G + q++;
                nmod_mpoly_init(g, ctx);
                const F4Poly* f = S.B + sel[a];
                for (slong t = 0; t < f->len; t++)
                    nmod_mpoly_push_term_ui_ui(g, f->cf[t], S.M.ex + f->mon[t] * nv, ctx);
                /* pushed in descending order: no sort needed */
            }
        } else *expired = 1;
        flint_free(sel);
        *ng = q;
    } else *expired = 1;
    if (!ok) *expired = 1;
    for (slong k = 0; k < S.nb; k++) f4_poly_clear(S.B + k);
    flint_free(S.B); flint_free(S.sugar); flint_free(S.active); flint_free(S.P);
    flint_free(S.tmp); flint_free(ex);
    fm_clear(&S.M);
    if (ok) f4_hits++;
    return ok;
}

#endif /* USE_FLINT */

#ifndef USE_FLINT
static int f4_off_depth_stub = 0;
int* gb_f4_off_counter(void) { return &f4_off_depth_stub; }
long gb_f4_hits(void) { return 0; }
#endif

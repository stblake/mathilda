/* partitions.c — IntegerPartitions
 *
 * A faithful, efficient recreation of the Wolfram-Language IntegerPartitions.
 * The whole surface collapses onto a single count-vector enumerator over an
 * ordered set of allowed parts, run with exact GMP rational arithmetic so that
 * integers, big integers, rationals and negative values are all handled by the
 * same code path.
 *
 * Forms:
 *   IntegerPartitions[n]                     all partitions of n
 *   IntegerPartitions[n, k]                  into at most k parts
 *   IntegerPartitions[n, {k}]                into exactly k parts
 *   IntegerPartitions[n, {kmin, kmax}]       between kmin and kmax parts
 *   IntegerPartitions[n, {kmin, kmax, dk}]   kmin, kmin+dk, ... parts (dk < 0
 *                                            counts down from kmin to kmax)
 *   IntegerPartitions[n, kspec, sspec]       parts drawn only from sspec (a list,
 *                                            or a number s for Range[s])
 *   IntegerPartitions[n, kspec, sspec, m]    first m (m>0), last |m| (m<0), None,
 *                                            or a Take sequence {a}, {a, b}, {a, b, s}
 *
 * n and the s_i may be rational and/or negative. Results are in reverse
 * lexicographic order; within a partition the parts appear in the order of the
 * reversed sspec (descending for the default Range[n]).
 *
 * Ownership: this builtin only *reads* `res`. On every NULL return (bad
 * arguments, ::undef, symbolic input) the evaluator keeps `res` unevaluated.
 */

#include "partitions.h"
#include "symtab.h"
#include "sym_names.h"
#include "attr.h"
#include "print.h"
#include "eval.h"          /* eval_stamp_list_literal */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>     /* PRId64 */
#include <math.h>
#include <gmp.h>
#include "message.h"   /* mth_message: Quiet/Check funnel */
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* M_PI is POSIX, not C99; glibc hides it under -std=c99. */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static Expr* partitions_scalar(bool distinct, Expr* res);

/* ----- small numeric helpers ------------------------------------------- */

/* Read an Integer, BigInt or Rational[p,q] into a freshly-initialised mpq_t.
 * Returns false (leaving `out` untouched) for any other expression. */
static bool ip_to_mpq(const Expr* e, mpq_t out) {
    if (expr_is_integer_like(e)) {
        mpz_t z;
        expr_to_mpz(e, z);            /* inits z */
        mpq_init(out);
        mpq_set_z(out, z);
        mpz_clear(z);
        return true;
    }
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational &&
        e->data.function.arg_count == 2) {
        Expr* pe = e->data.function.args[0];
        Expr* qe = e->data.function.args[1];
        if (expr_is_integer_like(pe) && expr_is_integer_like(qe)) {
            mpz_t pz, qz;
            expr_to_mpz(pe, pz);
            expr_to_mpz(qe, qz);
            mpq_init(out);
            mpq_set_num(out, pz);
            mpq_set_den(out, qz);
            mpq_canonicalize(out);
            mpz_clear(pz);
            mpz_clear(qz);
            return true;
        }
    }
    return false;
}

/* Build an Integer/BigInt or Rational[n,d] Expr from an mpq_t. */
static Expr* ip_mpq_to_expr(const mpq_t q) {
    mpz_t num, den;
    mpz_init(num);
    mpz_init(den);
    mpq_get_num(num, q);
    mpq_get_den(den, q);
    if (mpz_cmp_ui(den, 1) == 0) {
        Expr* r = expr_bigint_normalize(expr_new_bigint_from_mpz(num));
        mpz_clear(num);
        mpz_clear(den);
        return r;
    }
    Expr* ne = expr_bigint_normalize(expr_new_bigint_from_mpz(num));
    Expr* de = expr_bigint_normalize(expr_new_bigint_from_mpz(den));
    Expr* args[2] = { ne, de };
    Expr* r = expr_new_function(expr_new_symbol(SYM_Rational), args, 2);
    mpz_clear(num);
    mpz_clear(den);
    return r;
}

/* An mpz that fits in int64 into *out; false otherwise. (Not mpz_get_si:
 * `long` is 32 bits on LLP64 Windows.) */
static bool ip_mpz_i64(const mpz_t z, int64_t* out) {
    if (mpz_sizeinbase(z, 2) > 63) return false;
    uint64_t u = 0;
    size_t cnt = 0;
    mpz_export(&u, &cnt, -1, sizeof u, 0, 0, z);
    *out = mpz_sgn(z) < 0 ? -(int64_t)u : (int64_t)u;
    return true;
}

/* int64 into a (caller-initialised) mpz. */
static void ip_i64_mpz(mpz_t z, int64_t v) {
    uint64_t u = v < 0 ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
    mpz_import(z, 1, -1, sizeof u, 0, 0, &u);
    if (v < 0) mpz_neg(z, z);
}

/* floor(a / b) clamped to [0, INT64_MAX]; b must be non-zero. */
static int64_t ip_floor_div(const mpq_t a, const mpq_t b) {
    mpq_t q;
    mpq_init(q);
    mpq_div(q, a, b);
    mpz_t fl;
    mpz_init(fl);
    mpz_fdiv_q(fl, mpq_numref(q), mpq_denref(q));
    int64_t r;
    if (mpz_sgn(fl) < 0)              r = 0;
    else if (!ip_mpz_i64(fl, &r))     r = INT64_MAX;
    mpz_clear(fl);
    mpq_clear(q);
    return r;
}

static bool ip_is_sym(const Expr* e, const char* s) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == s;
}

static bool ip_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && ip_is_sym(e->data.function.head, SYM_List);
}

/* ----- enumeration context --------------------------------------------- */

typedef struct {
    mpq_t*   parts;          /* explicit sspec, reversed (display order)    */
    bool*    suffix_nonneg;  /* all parts[j>=i] have sign >= 0              */
    bool*    suffix_nonpos;  /* all parts[j>=i] have sign <= 0              */
    int64_t* counts;         /* current count chosen per part               */
    size_t   nparts;

    /* The admissible lengths: lo <= len <= hi (no upper bound when hi_inf)
     * with |len - anchor| a multiple of step. {kmin, kmax, dk} with dk < 0
     * is the same set read downwards: lo = kmax, hi = kmin, anchor = kmin. */
    int64_t  lo, hi, anchor, step;
    bool     hi_inf;

    Expr**   out;            /* collected partitions, in order              */
    size_t   out_n, out_cap;
    int64_t  limit;          /* stop after this many; 0 = unlimited         */
    bool     stop;
    bool     error;          /* allocation failure                          */
} PartCtx;

static bool ip_len_ok(const PartCtx* c, int64_t len) {
    if (len < c->lo) return false;
    if (!c->hi_inf && len > c->hi) return false;
    int64_t d = len >= c->anchor ? len - c->anchor : c->anchor - len;
    return c->step <= 1 || d % c->step == 0;
}

/* Append a finished partition, honouring the first-m limit. Consumes it. */
static void ip_push(PartCtx* c, Expr* part) {
    if (!part) { c->error = true; return; }
    if (c->out_n == c->out_cap) {
        size_t cap = c->out_cap ? c->out_cap * 2 : 16;
        Expr** grown = realloc(c->out, cap * sizeof(Expr*));
        if (!grown) { expr_free(part); c->error = true; return; }
        c->out = grown;
        c->out_cap = cap;
    }
    c->out[c->out_n++] = part;
    if (c->limit > 0 && c->out_n >= (size_t)c->limit) c->stop = true;
}

/* ----- explicit sspec: exact count-vector enumerator --------------------- */

/* Materialise the current count vector as List[...] and append it. */
static void ip_emit(PartCtx* c, int64_t len) {
    Expr* part = expr_new_function(expr_new_symbol(SYM_List), NULL, (size_t)len);
    if (!part) { c->error = true; return; }
    size_t idx = 0;
    for (size_t i = 0; i < c->nparts; i++)
        for (int64_t k = 0; k < c->counts[i]; k++)
            part->data.function.args[idx++] = ip_mpq_to_expr(c->parts[i]);
    ip_push(c, part);
}

/* Assign a count to parts[i], recursing in descending-count order. This is
 * the general engine: rational and negative parts, unsorted or repeated
 * sspec entries (each entry is its own part, so {1, 1, 2} repeats results). */
static void ip_recurse(PartCtx* c, size_t i, const mpq_t remaining, int64_t len) {
    if (c->stop || c->error) return;

    if (i == c->nparts) {
        if (mpq_sgn(remaining) == 0 && ip_len_ok(c, len)) ip_emit(c, len);
        return;
    }

    const mpq_t* p = (const mpq_t*)&c->parts[i];
    int psgn = mpq_sgn(c->parts[i]);
    int rsgn = mpq_sgn(remaining);

    /* Length budget from a finite upper bound. */
    int64_t budget = -1;                    /* -1 == unlimited */
    if (!c->hi_inf) {
        budget = c->hi - len;
        if (budget < 0) budget = 0;
    }

    /* Tight numeric bound: only valid when the remaining parts are single
     * signed in p's direction, so contributions cannot overshoot and recover. */
    int64_t nb = -1;                        /* -1 == unbounded */
    if (psgn > 0 && c->suffix_nonneg[i]) {
        nb = (rsgn <= 0) ? 0 : ip_floor_div(remaining, *p);
    } else if (psgn < 0 && c->suffix_nonpos[i]) {
        nb = (rsgn >= 0) ? 0 : ip_floor_div(remaining, *p);
    }

    int64_t cmax;
    if (nb >= 0 && budget >= 0)      cmax = (nb < budget) ? nb : budget;
    else if (nb >= 0)                cmax = nb;
    else if (budget >= 0)            cmax = budget;
    else { c->error = true; return; }       /* unreachable: undef gate guards */

    mpq_t step, rem2;
    mpq_init(step);
    mpq_init(rem2);
    mpz_t cz;
    mpz_init(cz);
    for (int64_t cc = cmax; cc >= 0; cc--) {
        ip_i64_mpz(cz, cc);
        mpq_set_z(step, cz);
        mpq_mul(step, step, *p);            /* cc * p */
        mpq_sub(rem2, remaining, step);     /* remaining - cc*p */
        c->counts[i] = cc;
        ip_recurse(c, i + 1, rem2, len + cc);
        if (c->stop || c->error) break;
    }
    mpz_clear(cz);
    mpq_clear(step);
    mpq_clear(rem2);
}

/* ----- explicit sspec, machine integers ---------------------------------- */

/* The general enumerator pays exact rational arithmetic per step. When n
 * and every part are positive machine integers and the parts are strictly
 * descending (an ascending sspec), the output is plain reverse-lexicographic
 * order and a word-sized DFS over non-increasing part sequences produces it:
 *   - the first candidate at each depth is found by binary search;
 *   - a finite upper length prunes by ceil(r / p) parts still needed, a lower
 *     one by the most parts the remainder can still be cut into;
 *   - a reachability table (r a sum of parts S[j..]) prunes dead branches
 *     when the smallest part exceeds 1 and the table fits in 16 MB.
 * Each partition is one List node over refcounted copies of one shared
 * Integer per part value. Returns false (nothing emitted) when the inputs do
 * not qualify. */
static bool ip_fast(PartCtx* c, const mpq_t N) {
    if (mpz_cmp_ui(mpq_denref(N), 1) != 0) return false;
    int64_t n;
    if (!ip_mpz_i64(mpq_numref(N), &n) || n < 0) return false;
    size_t ns = c->nparts;
    int64_t* S = ns ? malloc(ns * sizeof(int64_t)) : NULL;
    if (ns && !S) return false;
    for (size_t j = 0; j < ns; j++) {
        if (mpz_cmp_ui(mpq_denref(c->parts[j]), 1) != 0 ||
            !ip_mpz_i64(mpq_numref(c->parts[j]), &S[j]) ||
            S[j] <= 0 || (j && S[j] >= S[j - 1])) { free(S); return false; }
    }

    /* n == 0: the single empty partition, subject to the length spec. */
    if (n == 0 || ns == 0) {
        if (n == 0 && ip_len_ok(c, 0))
            ip_push(c, expr_new_function(expr_new_symbol(SYM_List), NULL, 0));
        free(S);
        return true;
    }

    int64_t smin = S[ns - 1];
    int64_t maxdepth = n / smin;
    if (!c->hi_inf && c->hi < maxdepth) maxdepth = c->hi;
    if (maxdepth > ((int64_t)1 << 26)) { free(S); return false; }

    unsigned char* reach = NULL;
    if (smin > 1 && (uint64_t)ns * (uint64_t)(n + 1) <= (uint64_t)16 << 20) {
        reach = calloc(ns * (size_t)(n + 1), 1);
        if (reach) {
            for (size_t k = ns; k > 0; k--) {
                size_t j = k - 1;
                unsigned char* row = reach + j * (size_t)(n + 1);
                const unsigned char* below = (j + 1 < ns) ? row + (n + 1) : NULL;
                row[0] = 1;
                for (int64_t v = 1; v <= n; v++)
                    row[v] = (unsigned char)((below && below[v]) || (v >= S[j] && row[v - S[j]]));
            }
            if (!reach[n]) {                  /* no partition at all */
                free(reach); free(S);
                return true;
            }
        }
    }

    size_t* cur = malloc(((size_t)maxdepth + 1) * sizeof(size_t));
    size_t* nxt = malloc(((size_t)maxdepth + 1) * sizeof(size_t));
    Expr** node = calloc(ns, sizeof(Expr*));
    if (!cur || !nxt || !node) {
        free(cur); free(nxt); free(node); free(reach); free(S);
        return false;
    }

    Expr* lhead = expr_new_symbol(SYM_List);
    bool any_len = c->lo == 0 && c->hi_inf && c->step <= 1;
    int64_t d = 0, r = n;
    nxt[0] = 0;
    for (;;) {
        bool descend = false;
        if (r == 0) {
            if (any_len || ip_len_ok(c, d)) {
                Expr* part = expr_new_function(expr_copy(lhead), NULL, (size_t)d);
                if (part) {
                    for (int64_t k = 0; k < d; k++) {
                        size_t j = cur[k];
                        if (!node[j]) node[j] = expr_new_integer(S[j]);
                        part->data.function.args[k] = expr_copy(node[j]);
                    }
                }
                ip_push(c, part);
                if (c->stop || c->error) break;
            }
        } else if (d < maxdepth) {
            size_t lo = nxt[d];
            if (lo < ns && S[lo] > r) {          /* first index with S[j] <= r */
                size_t a = lo, b = ns;
                while (a < b) { size_t m = a + (b - a) / 2; if (S[m] > r) a = m + 1; else b = m; }
                lo = a;
            }
            for (size_t j = lo; j < ns; j++) {
                int64_t p = S[j], rest = r - p;
                /* upper length: r / p parts still needed; smaller p need more */
                if (!c->hi_inf && r / p + (r % p != 0) > c->hi - d) break;
                /* lower length: the most parts r can still become */
                if (d + 1 + rest / smin < c->lo) continue;
                if (reach && !reach[j * (size_t)(n + 1) + (size_t)rest]) continue;
                cur[d] = j;
                nxt[d] = j + 1;
                r = rest;
                d++;
                nxt[d] = j;                   /* parts are non-increasing */
                descend = true;
                break;
            }
        }
        if (descend) continue;
        if (d == 0) break;                    /* backtrack */
        d--;
        r += S[cur[d]];
    }

    for (size_t j = 0; j < ns; j++) if (node[j]) expr_free(node[j]);
    expr_free(lhead);
    free(node); free(cur); free(nxt); free(reach); free(S);
    return true;
}

/* ----- explicit sspec, machine integers in any order -------------------- */

/* ip_recurse's count-vector walk -- the count of priority part i chosen
 * largest first, then part i+1 -- in int64, for positive machine-integer
 * parts in ANY order (an unsorted or repeated sspec fixes the display order
 * and the duplicates, so ip_fast's reverse-lexicographic walk does not
 * apply). Pruned by a reachability table (rest a sum of parts i+1..) and by
 * the length bounds through the suffix extremes of the parts. */
typedef struct {
    PartCtx*        c;
    const int64_t*  P;          /* parts in priority order, all > 0        */
    size_t          ns;
    int64_t         n;
    int64_t*        cnt;
    const int64_t*  sufmin;     /* min P[i..], max P[i..]                   */
    const int64_t*  sufmax;
    const unsigned char* reach; /* reach[i * (n + 1) + v]: v from P[i..]    */
    Expr**          node;       /* shared Integer per entry                 */
    Expr*           lhead;
} IpCount;

static void ipc_emit(IpCount* k, int64_t len) {
    Expr* part = expr_new_function(expr_copy(k->lhead), NULL, (size_t)len);
    if (part) {
        size_t idx = 0;
        for (size_t i = 0; i < k->ns; i++) {
            if (!k->cnt[i]) continue;
            if (!k->node[i]) k->node[i] = expr_new_integer(k->P[i]);
            for (int64_t j = 0; j < k->cnt[i]; j++)
                part->data.function.args[idx++] = expr_copy(k->node[i]);
        }
    }
    ip_push(k->c, part);
}

static void ipc_rec(IpCount* k, size_t i, int64_t r, int64_t len) {
    PartCtx* c = k->c;
    if (c->stop || c->error) return;
    if (i == k->ns) {
        if (r == 0 && ip_len_ok(c, len)) ipc_emit(k, len);
        return;
    }
    int64_t p = k->P[i];
    int64_t cmax = r / p;
    if (!c->hi_inf && c->hi - len < cmax) cmax = c->hi - len;
    for (int64_t cc = cmax; cc >= 0; cc--) {
        int64_t rest = r - cc * p, l2 = len + cc;
        if (rest > 0) {
            if (i + 1 == k->ns) continue;           /* nothing left to fill it */
            if (k->reach && !k->reach[(i + 1) * (size_t)(k->n + 1) + (size_t)rest]) continue;
            int64_t mx = k->sufmax[i + 1];
            if (!c->hi_inf && l2 + rest / mx + (rest % mx != 0) > c->hi) continue;
            if (l2 + rest / k->sufmin[i + 1] < c->lo) continue;
        }
        k->cnt[i] = cc;
        ipc_rec(k, i + 1, rest, l2);
        if (c->stop || c->error) break;
    }
    k->cnt[i] = 0;
}

/* Returns false (nothing emitted) unless n and every part are positive
 * machine integers. */
static bool ip_counts_i64(PartCtx* c, const mpq_t N) {
    if (mpz_cmp_ui(mpq_denref(N), 1) != 0) return false;
    int64_t n;
    if (!ip_mpz_i64(mpq_numref(N), &n) || n < 0) return false;
    size_t ns = c->nparts;
    if (ns == 0) return false;
    int64_t* P = malloc(ns * sizeof(int64_t));
    int64_t* sufmin = malloc(ns * sizeof(int64_t));
    int64_t* sufmax = malloc(ns * sizeof(int64_t));
    int64_t* cnt = calloc(ns, sizeof(int64_t));
    Expr** node = calloc(ns, sizeof(Expr*));
    bool ok = P && sufmin && sufmax && cnt && node;
    for (size_t j = 0; ok && j < ns; j++)
        ok = mpz_cmp_ui(mpq_denref(c->parts[j]), 1) == 0 &&
             ip_mpz_i64(mpq_numref(c->parts[j]), &P[j]) && P[j] > 0;
    if (!ok) { free(P); free(sufmin); free(sufmax); free(cnt); free(node); return false; }
    for (size_t k = ns; k > 0; k--) {
        size_t j = k - 1;
        sufmin[j] = (j + 1 < ns && sufmin[j + 1] < P[j]) ? sufmin[j + 1] : P[j];
        sufmax[j] = (j + 1 < ns && sufmax[j + 1] > P[j]) ? sufmax[j + 1] : P[j];
    }
    unsigned char* reach = NULL;
    if ((uint64_t)(ns + 1) * (uint64_t)(n + 1) <= (uint64_t)16 << 20) {
        reach = calloc((ns + 1) * (size_t)(n + 1), 1);
        if (reach) {
            reach[ns * (size_t)(n + 1)] = 1;      /* only 0 from no parts */
            for (size_t k = ns; k > 0; k--) {
                size_t j = k - 1;
                unsigned char* row = reach + j * (size_t)(n + 1);
                const unsigned char* below = row + (n + 1);
                for (int64_t v = 0; v <= n; v++)
                    row[v] = (unsigned char)(below[v] || (v >= P[j] && row[v - P[j]]));
            }
        }
    }
    IpCount k = { c, P, ns, n, cnt, sufmin, sufmax, reach, node, expr_new_symbol(SYM_List) };
    if (!reach || reach[n]) ipc_rec(&k, 0, n, 0);
    for (size_t j = 0; j < ns; j++) if (node[j]) expr_free(node[j]);
    expr_free(k.lhead);
    free(P); free(sufmin); free(sufmax); free(cnt); free(node); free(reach);
    return true;
}

/* ----- range mode: the parts are 1..smax ---------------------------------- */

/* The default sspec (Range[n]) and a number s (Range[s]) need no part list:
 * at depth d with remainder r the candidates are the integers
 *     max(1, ceil(r / (hi - d)))  ..  min(previous part, r, smax, r - (lo - d - 1)),
 * largest first -- the length bounds prune exactly, so every node visited
 * leads to an admissible length (up to the step), and n need not be small:
 * IntegerPartitions[10^9, 2] is 5*10^8 + 1 rows, not a 10^9-entry part
 * table. The stack grows with the depth actually reached. */
static void ip_range_i64(PartCtx* c, int64_t n, int64_t smax) {
    if (n < 0) return;
    if (n == 0) {
        if (ip_len_ok(c, 0)) ip_push(c, expr_new_function(expr_new_symbol(SYM_List), NULL, 0));
        return;
    }
    if (smax > n) smax = n;
    if (smax < 1) return;

    /* One shared Integer per small part value. */
    size_t ncache = smax < 4096 ? (size_t)smax : 4096;
    Expr** node = calloc(ncache + 1, sizeof(Expr*));
    Expr* lhead = expr_new_symbol(SYM_List);
    size_t cap = 64;
    int64_t* cur = malloc(cap * sizeof(int64_t));  /* chosen part per depth   */
    int64_t* ub  = malloc(cap * sizeof(int64_t));  /* next candidate to try   */
    if (!node || !cur || !ub) {
        free(node); free(cur); free(ub); expr_free(lhead);
        c->error = true;
        return;
    }
    bool any_len = c->lo == 0 && c->hi_inf && c->step <= 1;

    int64_t d = 0, r = n;
    ub[0] = smax;
    for (;;) {
        bool descend = false;
        if (r == 0) {
            if (any_len || ip_len_ok(c, d)) {
                Expr* part = expr_new_function(expr_copy(lhead), NULL, (size_t)d);
                if (part) {
                    for (int64_t k = 0; k < d; k++) {
                        int64_t v = cur[k];
                        if ((size_t)v <= ncache) {
                            if (!node[v]) node[v] = expr_new_integer(v);
                            part->data.function.args[k] = expr_copy(node[v]);
                        } else {
                            part->data.function.args[k] = expr_new_integer(v);
                        }
                    }
                }
                ip_push(c, part);
                if (c->stop || c->error) break;
            }
        } else if (c->hi_inf || d < c->hi) {
            int64_t p_hi = ub[d] < r ? ub[d] : r;
            int64_t p_lo = 1;
            if (!c->hi_inf) {
                int64_t room = c->hi - d;
                int64_t need = r / room + (r % room != 0);
                if (need > p_lo) p_lo = need;
            }
            if (c->lo > d + 1) {
                int64_t cut = r - (c->lo - d - 1);    /* rest must be >= lo-d-1 ones */
                if (cut < p_hi) p_hi = cut;
            }
            if (p_hi >= p_lo) {
                if ((size_t)d + 1 >= cap) {
                    size_t nc = cap * 2;
                    int64_t* g1 = realloc(cur, nc * sizeof(int64_t));
                    if (g1) cur = g1;
                    int64_t* g2 = g1 ? realloc(ub, nc * sizeof(int64_t)) : NULL;
                    if (g2) ub = g2;
                    if (!g1 || !g2) { c->error = true; break; }
                    cap = nc;
                }
                cur[d] = p_hi;
                ub[d] = p_hi - 1;
                r -= p_hi;
                d++;
                ub[d] = p_hi;                       /* non-increasing */
                descend = true;
            }
        }
        if (descend) continue;
        if (d == 0) break;                          /* backtrack */
        d--;
        r += cur[d];
    }

    for (size_t v = 0; v <= ncache; v++) if (node[v]) expr_free(node[v]);
    free(node); free(cur); free(ub);
    expr_free(lhead);
}

/* The same walk over GMP integers, for a bignum n (IntegerPartitions[10^20,
 * {2}]). Only reached with a length bound or a first-m limit -- otherwise the
 * output could not exist -- so speed is not the point; agreement is. */
static void ip_range_mpz(PartCtx* c, const mpz_t n, const mpz_t smax_in) {
    mpz_t smax, r, t, lim;
    mpz_inits(smax, r, t, lim, NULL);
    mpz_set(smax, smax_in);
    if (mpz_cmp(smax, n) > 0) mpz_set(smax, n);
    if (mpz_sgn(n) <= 0 || mpz_sgn(smax) < 1) { mpz_clears(smax, r, t, lim, NULL); return; }

    size_t cap = 64;
    mpz_t* cur = malloc(cap * sizeof(mpz_t));
    mpz_t* ub  = malloc(cap * sizeof(mpz_t));
    if (!cur || !ub) { free(cur); free(ub); mpz_clears(smax, r, t, lim, NULL); c->error = true; return; }
    for (size_t i = 0; i < cap; i++) { mpz_init(cur[i]); mpz_init(ub[i]); }

    int64_t d = 0;
    mpz_set(r, n);
    mpz_set(ub[0], smax);
    for (;;) {
        bool descend = false;
        if (mpz_sgn(r) == 0) {
            if (ip_len_ok(c, d)) {
                Expr* part = expr_new_function(expr_new_symbol(SYM_List), NULL, (size_t)d);
                if (part)
                    for (int64_t k = 0; k < d; k++)
                        part->data.function.args[k] =
                            expr_bigint_normalize(expr_new_bigint_from_mpz(cur[k]));
                ip_push(c, part);
                if (c->stop || c->error) break;
            }
        } else if (c->hi_inf || d < c->hi) {
            /* p_hi = min(ub[d], r, r - (lo - d - 1)); p_lo = max(1, ceil(r / (hi - d))) */
            mpz_set(t, mpz_cmp(ub[d], r) < 0 ? ub[d] : r);
            if (c->lo > d + 1) {
                ip_i64_mpz(lim, c->lo - d - 1);
                mpz_sub(lim, r, lim);
                if (mpz_cmp(lim, t) < 0) mpz_set(t, lim);
            }
            mpz_set_ui(lim, 1);
            if (!c->hi_inf) {
                mpz_t room;
                mpz_init(room);
                ip_i64_mpz(room, c->hi - d);
                mpz_cdiv_q(room, r, room);
                if (mpz_cmp(room, lim) > 0) mpz_set(lim, room);
                mpz_clear(room);
            }
            if (mpz_cmp(t, lim) >= 0) {
                if ((size_t)d + 1 >= cap) {
                    size_t nc = cap * 2;
                    mpz_t* g1 = realloc(cur, nc * sizeof(mpz_t));
                    if (g1) cur = g1;
                    mpz_t* g2 = g1 ? realloc(ub, nc * sizeof(mpz_t)) : NULL;
                    if (g2) ub = g2;
                    if (!g1 || !g2) { c->error = true; break; }
                    for (size_t i = cap; i < nc; i++) { mpz_init(cur[i]); mpz_init(ub[i]); }
                    cap = nc;
                }
                mpz_set(cur[d], t);
                mpz_sub_ui(ub[d], t, 1);
                mpz_sub(r, r, t);
                d++;
                mpz_set(ub[d], t);
                descend = true;
            }
        }
        if (descend) continue;
        if (d == 0) break;
        d--;
        mpz_add(r, r, cur[d]);
    }

    for (size_t i = 0; i < cap; i++) { mpz_clear(cur[i]); mpz_clear(ub[i]); }
    free(cur); free(ub);
    mpz_clears(smax, r, t, lim, NULL);
}

/* ----- diagnostics ----------------------------------------------------- */

/* "IntegerPartitions[a0, a1, ...]" for the first `nargs` arguments of res. */
static char* ip_call_str(const Expr* res, size_t nargs) {
    size_t cap = 64, len = 0;
    char* buf = malloc(cap);
    if (!buf) return NULL;
    len += (size_t)snprintf(buf, cap, "IntegerPartitions[");
    for (size_t i = 0; i < nargs; i++) {
        char* a = expr_to_string(res->data.function.args[i]);
        const char* piece = a ? a : "?";
        /* ", " + piece, plus room for the closing ']' and NUL */
        size_t need = len + 2 + strlen(piece) + 2;
        if (need > cap) {
            while (need > cap) cap *= 2;
            char* g = realloc(buf, cap);
            if (!g) { free(buf); free(a); return NULL; }
            buf = g;
        }
        if (i) { buf[len++] = ','; buf[len++] = ' '; }
        memcpy(buf + len, piece, strlen(piece));
        len += strlen(piece);
        free(a);
    }
    buf[len++] = ']';
    buf[len] = '\0';
    return buf;
}

/* Emit IntegerPartitions::<tag> for the whole call: `fmt` has one %s for it. */
static void ip_msg_call(const Expr* res, size_t nargs, const char* tag, const char* fmt) {
    char* call = ip_call_str(res, nargs);
    mth_message("IntegerPartitions", tag, fmt, call ? call : "IntegerPartitions[...]");
    free(call);
}

/* ----- argument parsing ------------------------------------------------ */

/* An Integer or BigInt, clamped to the int64 range. */
static bool ip_read_int(const Expr* e, int64_t* val) {
    if (e->type == EXPR_INTEGER) { *val = e->data.integer; return true; }
    if (e->type == EXPR_BIGINT) {
        *val = mpz_sgn(e->data.bigint) < 0 ? -INT64_MAX : INT64_MAX;
        return true;
    }
    return false;
}

/* Position 2: All | Infinity | k >= 0 | {kmin} | {kmin, kmax} | {kmin, kmax, dk}
 * with kmin >= 0, kmax >= 0 or Infinity, dk a nonzero integer. Fills the
 * length set; *empty when it admits no length at all. false = ::nninfseq. */
static bool ip_parse_kspec(const Expr* ks, PartCtx* c, bool* empty) {
    c->lo = 0; c->hi = 0; c->anchor = 0; c->step = 1; c->hi_inf = true;
    *empty = false;
    int64_t v;
    if (!ks || ip_is_sym(ks, SYM_All) || ip_is_sym(ks, SYM_Infinity)) return true;
    if (ip_read_int(ks, &v)) {
        if (v < 0) return false;
        c->hi = v; c->hi_inf = false;
        return true;
    }
    if (!ip_is_list(ks)) return false;
    size_t kc = ks->data.function.arg_count;
    Expr** ka = ks->data.function.args;
    if (kc < 1 || kc > 3) return false;
    int64_t kmin, kmax = 0, dk = 1;
    bool kmax_inf = false;
    if (!ip_read_int(ka[0], &kmin) || kmin < 0) return false;
    if (kc == 1) kmax = kmin;
    else if (ip_is_sym(ka[1], SYM_Infinity)) kmax_inf = true;
    else if (!ip_read_int(ka[1], &kmax) || kmax < 0) return false;
    if (kc == 3 && (!ip_read_int(ka[2], &dk) || dk == 0)) return false;

    if (dk > 0) {
        c->lo = kmin; c->hi = kmax; c->hi_inf = kmax_inf;
        c->anchor = kmin; c->step = dk;
        if (!kmax_inf && kmin > kmax) *empty = true;
    } else {                                   /* counting down from kmin */
        if (kmax_inf || kmin < kmax) { *empty = true; return true; }
        c->lo = kmax; c->hi = kmin; c->hi_inf = false;
        c->anchor = kmin; c->step = dk == -INT64_MAX - 1 ? INT64_MAX : -dk;
    }
    return true;
}

/* Position 4 as a Take sequence {a, b, s}; *none for None (nothing). false =
 * ::seq. A bare m > 0 is {1, m, 1}, m < 0 is {m, -1, 1}, {a} is {a, a, 1}. */
typedef struct { int64_t a, b, s; bool all, none; } IpTake;

static bool ip_parse_take(const Expr* me, IpTake* t) {
    t->a = 1; t->b = -1; t->s = 1; t->all = false; t->none = false;
    int64_t v;
    if (!me || ip_is_sym(me, SYM_All)) { t->all = true; return true; }
    if (ip_is_sym(me, SYM_None)) { t->none = true; return true; }
    if (ip_read_int(me, &v)) {
        if (v > 0)      { t->a = 1; t->b = v; }
        else if (v < 0) { t->a = v; t->b = -1; }
        else            { t->none = true; }        /* m = 0: nothing */
        return true;
    }
    if (!ip_is_list(me)) return false;
    size_t n = me->data.function.arg_count;
    Expr** a = me->data.function.args;
    if (n < 1 || n > 3) return false;
    if (!ip_read_int(a[0], &t->a)) return false;
    t->b = t->a;
    if (n >= 2 && !ip_read_int(a[1], &t->b)) return false;
    if (n == 3 && (!ip_read_int(a[2], &t->s) || t->s == 0)) return false;
    return true;
}

/* Apply a Take sequence to the collected rows (which it consumes), clamping
 * an index that runs past either end with an IntegerPartitions::take warning
 * and answering what is there -- {20} of 11 rows is {}, {2, 30} is rows 2..11.
 * An index 0 gives {} silently. Returns the selected rows' count in *nsel
 * and the rows in a fresh array. */
static Expr** ip_apply_take(const Expr* res, PartCtx* c, const IpTake* t, size_t* nsel) {
    int64_t T = (int64_t)c->out_n;
    *nsel = 0;
    Expr** sel = malloc((c->out_n ? c->out_n : 1) * sizeof(Expr*));
    bool* used = calloc(c->out_n ? c->out_n : 1, sizeof(bool));
    if (!sel || !used) { free(sel); free(used); c->error = true; return NULL; }

    if (t->a != 0 && t->b != 0) {
        int64_t a = t->a > 0 ? t->a : T + 1 + t->a;
        int64_t b = t->b > 0 ? t->b : T + 1 + t->b;
        bool warn = false, any = true;
        if (t->s > 0) {
            if (a < 1) { warn = true; a = 1; }
            if (b > T) { warn = true; b = T; }
            if (a > T || b < 1) { warn = true; any = false; }
        } else {
            if (a > T) { warn = true; a = T; }
            if (b < 1) { warn = true; b = 1; }
            if (a < 1 || b > T) { warn = true; any = false; }
        }
        if (warn) {
            char* call = ip_call_str(res, 3);
            mth_message("IntegerPartitions", "take",
                        "Warning: not all elements were found when attempting to take the "
                        "sequence {%" PRId64 ", %" PRId64 ", %" PRId64 "} from %s, which has "
                        "length %" PRId64 ".", t->a, t->b, t->s,
                        call ? call : "IntegerPartitions[...]", T);
            free(call);
        }
        if (any) {
            for (int64_t i = a; t->s > 0 ? i <= b : i >= b; i += t->s) {
                sel[(*nsel)++] = c->out[i - 1];
                used[i - 1] = true;
                if ((t->s > 0 && i > INT64_MAX - t->s) || (t->s < 0 && i < INT64_MIN - t->s)) break;
            }
        }
    }
    for (size_t i = 0; i < c->out_n; i++) if (!used[i]) expr_free(c->out[i]);
    free(used);
    free(c->out);
    c->out = NULL;
    c->out_n = c->out_cap = 0;
    return sel;
}

void partitions_init(void);

Expr* builtin_integerpartitions(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc < 1 || argc > 4) {
        mth_message("IntegerPartitions", "argb", "IntegerPartitions called with %zu argument%s; between 1 and 4 arguments are expected.", argc, argc == 1 ? "" : "s");
        return NULL;
    }

    /* --- validate every position first, reporting the first bad one --- */
    mpq_t N;
    if (!ip_to_mpq(args[0], N)) {
        ip_msg_call(res, argc, "rat", "Position 1 of %s must be a rational number.");
        return NULL;
    }

    PartCtx ctx;
    memset(&ctx, 0, sizeof ctx);
    bool k_empty = false;
    if (!ip_parse_kspec(argc >= 2 ? args[1] : NULL, &ctx, &k_empty)) {
        ip_msg_call(res, argc, "nninfseq",
                    "Position 2 of %s must be All, Infinity, a non-negative integer, or a List "
                    "whose first element (required) is a non-negative integer, second element "
                    "(optional) is a non-negative integer or Infinity, and third element "
                    "(optional) is a nonzero integer.");
        mpq_clear(N);
        return NULL;
    }

    /* sspec: All (Range[n]), a rational s (Range[s]) or a list of rationals. */
    bool range = true;
    mpq_t smax;
    mpq_init(smax);
    mpq_set(smax, N);
    if (argc >= 3 && !ip_is_sym(args[2], SYM_All)) {
        mpq_t s;
        if (ip_to_mpq(args[2], s)) {
            mpq_set(smax, s);
            mpq_clear(s);
        } else if (ip_is_list(args[2])) {
            range = false;
            size_t r = args[2]->data.function.arg_count;
            for (size_t i = 0; i < r && !range; i++) {
                mpq_t q;
                if (!ip_to_mpq(args[2]->data.function.args[i], q)) { range = true; break; }
                mpq_clear(q);
            }
            if (range) {                       /* a non-rational entry */
                ip_msg_call(res, argc, "rseq", "Position 3 of %s must be All, a rational number, or a list of rational numbers.");
                mpq_clear(smax); mpq_clear(N);
                return NULL;
            }
        } else {
            ip_msg_call(res, argc, "rseq", "Position 3 of %s must be All, a rational number, or a list of rational numbers.");
            mpq_clear(smax); mpq_clear(N);
            return NULL;
        }
    }

    IpTake take;
    if (!ip_parse_take(argc >= 4 ? args[3] : NULL, &take)) {
        ip_msg_call(res, argc, "seq", "Position 4 of %s must be All, None, an integer, or a List of 1, 2, or 3 integers, with the third (if present) nonzero.");
        mpq_clear(smax); mpq_clear(N);
        return NULL;
    }

    /* --- explicit sspec: the reversed part list (display priority order) --- */
    mpq_t* parts = NULL;
    size_t nparts = 0;
    bool has_pos = false, has_neg = false, has_zero = false;
    if (!range) {
        size_t r = args[2]->data.function.arg_count;
        Expr** sa = args[2]->data.function.args;
        parts = malloc((r ? r : 1) * sizeof(mpq_t));
        if (!parts) { mpq_clear(smax); mpq_clear(N); return NULL; }
        for (size_t i = 0; i < r; i++) {
            ip_to_mpq(sa[r - 1 - i], parts[i]);     /* validated above */
            int sg = mpq_sgn(parts[i]);
            if (sg > 0) has_pos = true; else if (sg < 0) has_neg = true; else has_zero = true;
        }
        nparts = r;
    }

    /* --- infinite-result guard (::undef) --- */
    if (!range && !k_empty && ctx.hi_inf && (has_zero || (has_pos && has_neg))) {
        ip_msg_call(res, argc, "undef", "%s contains partitions that are undefined because they are infinitely large.");
        for (size_t i = 0; i < nparts; i++) mpq_clear(parts[i]);
        free(parts);
        mpq_clear(smax); mpq_clear(N);
        return NULL;
    }

    /* A positive Take window bounds the enumeration: rows past b are never
     * needed (and when there are fewer, the ::take check sees the true total). */
    if (!take.all && !take.none && take.s > 0 && take.a > 0 && take.b >= take.a)
        ctx.limit = take.b;

    /* --- enumerate --- */
    bool* suffix_nonneg = NULL;
    bool* suffix_nonpos = NULL;
    int64_t* counts = NULL;
    if (!k_empty && !take.none) {
        if (range) {
            /* The parts are the integers 1..floor(smax); n must be an integer. */
            if (mpz_cmp_ui(mpq_denref(N), 1) == 0) {
                mpz_t fl;
                mpz_init(fl);
                mpz_fdiv_q(fl, mpq_numref(smax), mpq_denref(smax));
                int64_t n64, s64;
                if (ip_mpz_i64(mpq_numref(N), &n64)) {
                    if (!ip_mpz_i64(fl, &s64)) s64 = mpz_sgn(fl) > 0 ? INT64_MAX : 0;
                    ip_range_i64(&ctx, n64, s64);
                } else if (!ctx.hi_inf || ctx.limit > 0) {
                    ip_range_mpz(&ctx, mpq_numref(N), fl);
                } else {
                    ctx.error = true;               /* no bound: no answer */
                }
                mpz_clear(fl);
            }
        } else {
            suffix_nonneg = nparts ? malloc(nparts * sizeof(bool)) : NULL;
            suffix_nonpos = nparts ? malloc(nparts * sizeof(bool)) : NULL;
            counts        = nparts ? malloc(nparts * sizeof(int64_t)) : NULL;
            if (nparts && (!suffix_nonneg || !suffix_nonpos || !counts)) {
                ctx.error = true;
            } else {
                for (size_t k = nparts; k > 0; k--) {
                    size_t i = k - 1;
                    int sg = mpq_sgn(parts[i]);
                    bool nn = (sg >= 0), np = (sg <= 0);
                    if (i + 1 < nparts) {
                        nn = nn && suffix_nonneg[i + 1];
                        np = np && suffix_nonpos[i + 1];
                    }
                    suffix_nonneg[i] = nn;
                    suffix_nonpos[i] = np;
                }
                ctx.parts = parts;
                ctx.suffix_nonneg = suffix_nonneg;
                ctx.suffix_nonpos = suffix_nonpos;
                ctx.counts = counts;
                ctx.nparts = nparts;
                if (!ip_fast(&ctx, N) && !ip_counts_i64(&ctx, N)) ip_recurse(&ctx, 0, N, 0);
            }
        }
    }

    free(suffix_nonneg);
    free(suffix_nonpos);
    free(counts);
    for (size_t i = 0; i < nparts; i++) mpq_clear(parts[i]);
    free(parts);
    mpq_clear(smax);
    mpq_clear(N);

    if (ctx.error) {
        for (size_t i = 0; i < ctx.out_n; i++) expr_free(ctx.out[i]);
        free(ctx.out);
        return NULL;
    }

    Expr** sel = ctx.out;
    size_t sel_n = ctx.out_n;
    if (take.none) {
        sel_n = 0;
    } else if (!take.all) {
        sel = ip_apply_take(res, &ctx, &take, &sel_n);
        if (!sel) return NULL;
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), sel, sel_n);
    free(sel);
    /* 2*10^5 partitions are ~3.6*10^6 nodes of literal data: mark them as a
     * fixed point rather than have the evaluator re-walk every one. */
    eval_stamp_list_literal(result);
    return result;
}

/* ===================================================================== */
/* PartitionsP[n] — the partition-counting function p(n).                 */
/*                                                                        */
/* A table, then two engines (see partitions_count):                      */
/*                                                                        */
/*  0. n < PARTITIONS_TABLE_N: a copy from a table filled at startup.     */
/*                                                                        */
/*  1. partitionsp_recurrence — Euler's pentagonal-number-theorem         */
/*     recurrence, exact GMP integers. O(n^1.5) time and O(n*sqrt n) bits */
/*     of memory. Fills the table; the fallback without MPFR.             */
/*                                                                        */
/*  2. partitionsp_hrr — the non-recursive Hardy-Ramanujan-Rademacher     */
/*     exact formula evaluated in MPFR. O(n^{3/4}) cosine evaluations and */
/*     O(sqrt n)-bit working precision, so memory stays tiny. Used for    */
/*     large n, where the recurrence's table would be prohibitive.        */
/*     Reference: F. Johansson, "Efficient implementation of the          */
/*     Hardy-Ramanujan-Rademacher formula" (arXiv:1205.5991).             */
/* ===================================================================== */

/* p(n) and q(n) below this come from the startup tables; from here on the
 * HRR series is the faster engine for both (with the per-term precision of
 * hrr_q_eval, q's series overtakes the recurrence near n = 700). The
 * documented big values (p(1024), p(2048), p(4096)) exercise the HRR path. */
#define PARTITIONS_TABLE_N 1000

/* --- Engine 1: pentagonal recurrence ----------------------------------- */

/* Compute p(n) for n >= 0 into the freshly-initialised result `out`.
 * Returns 0 on success, -1 on allocation failure (out left at 0/uninit). */
/* Fill p[from..n] (initialised, zero) given p[0..from). */
static void partitionsp_fill(mpz_t* p, unsigned long from, unsigned long n) {
    if (from == 0) { mpz_set_ui(p[0], 1); from = 1; }
    for (unsigned long m = from; m <= n; m++) {
        /* p[m] accumulates the alternating pentagonal sum. */
        for (unsigned long k = 1;; k++) {
            /* g_k = k(3k-1)/2 ; once g_k > m no further k contributes. */
            unsigned long g1 = k * (3 * k - 1) / 2;
            if (g1 > m) break;
            int sign = (k & 1UL) ? 1 : -1; /* +1 for odd k, -1 for even k */
            if (sign > 0) mpz_add(p[m], p[m], p[m - g1]);
            else          mpz_sub(p[m], p[m], p[m - g1]);

            unsigned long g2 = k * (3 * k + 1) / 2;
            if (g2 <= m) {
                if (sign > 0) mpz_add(p[m], p[m], p[m - g2]);
                else          mpz_sub(p[m], p[m], p[m - g2]);
            }
        }
    }
}

int partitionsp_recurrence(unsigned long n, mpz_t out) {
    mpz_t* p = malloc((n + 1) * sizeof(mpz_t));
    if (!p) return -1;
    for (unsigned long m = 0; m <= n; m++) mpz_init(p[m]);
    partitionsp_fill(p, 0, n);
    mpz_set(out, p[n]);
    for (unsigned long m = 0; m <= n; m++) mpz_clear(p[m]);
    free(p);
    return 0;
}

/* --- Engine 2: Hardy-Ramanujan-Rademacher (MPFR) ----------------------- */
#ifdef USE_MPFR

/* W_k(n) = sqrt(3/k) * A_k(n), where A_k(n) is the Rademacher exponential
 * sum (Johansson eq. 1.6). We never form A_k(n) on its own: pulling the
 * sqrt(3/k) into Selberg's identity (eq. 2.1) cancels its sqrt(k/3) prefactor,
 * leaving, for k >= 3, the bare cosine sum
 *     sum over l in [0, 2k) with (3l^2+l)/2 == -n (mod k) of
 *         (-1)^l cos( pi (6l+1) / (6k) ).
 * Computed with Johansson's Algorithm 1: the residue m = ((3l^2+l)/2 + n)
 * mod k and its first difference r are advanced with integer adds only, so a
 * cosine is evaluated solely at the O(sqrt k) solutions. k = 1 and k = 2 are
 * the closed forms A_1 = 1, A_2 = (-1)^n scaled by sqrt(3/k). */
static void hrr_Wk(mpfr_t out, unsigned long k, unsigned long n,
                   const mpfr_t pi, mpfr_prec_t prec) {
    if (k == 1) {                       /* sqrt(3) * 1 */
        mpfr_sqrt_ui(out, 3, MPFR_RNDN);
        return;
    }
    if (k == 2) {                       /* sqrt(3/2) * (-1)^n */
        mpfr_set_ui(out, 3, MPFR_RNDN);
        mpfr_div_ui(out, out, 2, MPFR_RNDN);
        mpfr_sqrt(out, out, MPFR_RNDN);
        if (n & 1UL) mpfr_neg(out, out, MPFR_RNDN);
        return;
    }

    mpfr_t term, ang;
    mpfr_init2(term, prec);
    mpfr_init2(ang, prec);
    mpfr_set_zero(out, 1);

    unsigned long m = n % k;            /* l = 0 gives (3*0+0)/2 + n == n */
    unsigned long r = 2;               /* first difference of (3l^2+l)/2 */
    for (unsigned long l = 0; l < 2 * k; l++) {
        if (m == 0) {
            /* ang = pi * (6l + 1) / (6k) */
            mpfr_mul_ui(ang, pi, 6 * l + 1, MPFR_RNDN);
            mpfr_div_ui(ang, ang, 6 * k, MPFR_RNDN);
            mpfr_cos(term, ang, MPFR_RNDN);
            if (l & 1UL) mpfr_sub(out, out, term, MPFR_RNDN);
            else         mpfr_add(out, out, term, MPFR_RNDN);
        }
        m += r;                         /* next residue (m, r < k => one wrap) */
        if (m >= k) m -= k;
        r += 3;
        if (r >= k) r -= k;
    }

    mpfr_clear(term);
    mpfr_clear(ang);
}

/* Rademacher truncation bound M(n,N) (Johansson eq. 1.8); |p(n)-partial|<M.
 *     M = 44 pi^2/(225 sqrt 3) * N^{-1/2}
 *       + pi sqrt2/75 * (N/(n-1))^{1/2} * sinh( (pi/N) sqrt(2n/3) ).
 * Computed at the supplied precision; requires n >= 2. */
static void hrr_remainder_bound(mpfr_t out, unsigned long n, unsigned long N,
                                const mpfr_t pi, mpfr_prec_t prec) {
    mpfr_t a, b, t;
    mpfr_init2(a, prec);
    mpfr_init2(b, prec);
    mpfr_init2(t, prec);

    /* term 1: 44 pi^2/(225 sqrt 3) * N^{-1/2} */
    mpfr_mul(a, pi, pi, MPFR_RNDN);
    mpfr_mul_ui(a, a, 44, MPFR_RNDN);
    mpfr_sqrt_ui(t, 3, MPFR_RNDN);
    mpfr_mul_ui(t, t, 225, MPFR_RNDN);
    mpfr_div(a, a, t, MPFR_RNDN);
    mpfr_sqrt_ui(t, N, MPFR_RNDN);
    mpfr_div(a, a, t, MPFR_RNDN);

    /* term 2: pi sqrt2/75 * sqrt(N/(n-1)) * sinh( (pi/N) sqrt(2n/3) ) */
    mpfr_sqrt_ui(b, 2, MPFR_RNDN);
    mpfr_mul(b, b, pi, MPFR_RNDN);
    mpfr_div_ui(b, b, 75, MPFR_RNDN);
    mpfr_set_ui(t, N, MPFR_RNDN);
    mpfr_div_ui(t, t, n - 1, MPFR_RNDN);
    mpfr_sqrt(t, t, MPFR_RNDN);
    mpfr_mul(b, b, t, MPFR_RNDN);
    /* sinh argument */
    mpfr_set_ui(t, 2 * n, MPFR_RNDN);
    mpfr_div_ui(t, t, 3, MPFR_RNDN);
    mpfr_sqrt(t, t, MPFR_RNDN);
    mpfr_mul(t, t, pi, MPFR_RNDN);
    mpfr_div_ui(t, t, N, MPFR_RNDN);
    mpfr_sinh(t, t, MPFR_RNDN);
    mpfr_mul(b, b, t, MPFR_RNDN);

    mpfr_add(out, a, b, MPFR_RNDN);
    mpfr_clear(a);
    mpfr_clear(b);
    mpfr_clear(t);
}

/* Evaluate the HRR formula (Johansson eq. 1.4) at a fixed precision `prec`
 * with `N` terms, rounding the real sum to the nearest integer in `out`.
 * On success stores in `dist` the distance from the real sum to that integer
 * (a sanity margin; the caller retries at higher precision if it is large). */
static void hrr_eval(unsigned long n, unsigned long N, mpfr_prec_t prec,
                     mpz_t out, mpfr_t dist) {
    mpfr_t pi, C, x, U, sh, W, term, sum, rounded;
    mpfr_init2(pi, prec);
    mpfr_init2(C, prec);
    mpfr_init2(x, prec);
    mpfr_init2(U, prec);
    mpfr_init2(sh, prec);
    mpfr_init2(W, prec);
    mpfr_init2(term, prec);
    mpfr_init2(sum, prec);
    mpfr_init2(rounded, prec);

    mpfr_const_pi(pi, MPFR_RNDN);

    /* C = (pi/6) sqrt(24n - 1) */
    mpfr_set_ui(C, 24, MPFR_RNDN);
    mpfr_mul_ui(C, C, n, MPFR_RNDN);
    mpfr_sub_ui(C, C, 1, MPFR_RNDN);
    mpfr_sqrt(C, C, MPFR_RNDN);
    mpfr_mul(C, C, pi, MPFR_RNDN);
    mpfr_div_ui(C, C, 6, MPFR_RNDN);

    mpfr_set_zero(sum, 1);
    for (unsigned long k = 1; k <= N; k++) {
        hrr_Wk(W, k, n, pi, prec);
        if (mpfr_zero_p(W)) continue;       /* A_k(n) = 0 ~ half the time */

        /* x = C / k ; U(x) = cosh(x) - sinh(x)/x */
        mpfr_div_ui(x, C, k, MPFR_RNDN);
        mpfr_cosh(U, x, MPFR_RNDN);
        mpfr_sinh(sh, x, MPFR_RNDN);
        mpfr_div(sh, sh, x, MPFR_RNDN);
        mpfr_sub(U, U, sh, MPFR_RNDN);

        mpfr_mul(term, W, U, MPFR_RNDN);
        mpfr_add(sum, sum, term, MPFR_RNDN);
    }

    /* multiply by the common factor 4/(24n - 1) */
    mpfr_mul_ui(sum, sum, 4, MPFR_RNDN);
    mpfr_set_ui(x, 24, MPFR_RNDN);
    mpfr_mul_ui(x, x, n, MPFR_RNDN);
    mpfr_sub_ui(x, x, 1, MPFR_RNDN);
    mpfr_div(sum, sum, x, MPFR_RNDN);

    mpfr_round(rounded, sum);
    mpfr_sub(dist, sum, rounded, MPFR_RNDN);
    mpfr_abs(dist, dist, MPFR_RNDN);
    mpfr_get_z(out, rounded, MPFR_RNDN);

    mpfr_clear(pi);
    mpfr_clear(C);
    mpfr_clear(x);
    mpfr_clear(U);
    mpfr_clear(sh);
    mpfr_clear(W);
    mpfr_clear(term);
    mpfr_clear(sum);
    mpfr_clear(rounded);
}

/* Compute p(n) for n >= 2 via the HRR formula into `out`. Returns 0 on
 * success, -1 if it could not converge to a confident integer (the caller
 * then falls back to the recurrence). */
int partitionsp_hrr(unsigned long n, mpz_t out) {
    if (n < 2) return -1;               /* M(n,N) divides by n-1 */

    /* Number of terms: grow N from ceil(sqrt n) until the Rademacher bound
     * M(n,N) < 1/4, guaranteeing the rounded sum is exact (eq. 1.8). */
    unsigned long N = (unsigned long)ceil(sqrt((double)n));
    if (N < 1) N = 1;

    /* Working precision: the result has ~ C/ln2 bits where
     * C = (pi/6) sqrt(24n-1); add guard bits for accumulation/rounding. */
    double Cd = (M_PI / 6.0) * sqrt((double)24.0 * (double)n - 1.0);
    mpfr_prec_t prec = (mpfr_prec_t)(Cd / 0.6931471805599453) + 64;
    if (prec < 64) prec = 64;

    /* Pick N with a low-precision evaluation of the bound. */
    {
        mpfr_t pi, bound, quarter;
        mpfr_init2(pi, 64);
        mpfr_init2(bound, 64);
        mpfr_init2(quarter, 64);
        mpfr_const_pi(pi, MPFR_RNDN);
        mpfr_set_d(quarter, 0.25, MPFR_RNDN);
        for (unsigned long guard = 0; guard < 1000000UL; guard++) {
            hrr_remainder_bound(bound, n, N, pi, 64);
            if (mpfr_less_p(bound, quarter)) break;
            N++;
        }
        mpfr_clear(pi);
        mpfr_clear(bound);
        mpfr_clear(quarter);
    }

    /* Evaluate, retrying at higher precision if the rounding margin is thin
     * (it should not be: truncation < 1/4 and rounding error is tiny). */
    mpfr_t dist, margin;
    mpfr_init2(dist, 64);
    mpfr_init2(margin, 64);
    mpfr_set_d(margin, 0.25, MPFR_RNDN);
    int rc = -1;
    for (int attempt = 0; attempt < 4; attempt++) {
        hrr_eval(n, N, prec, out, dist);
        if (mpfr_less_p(dist, margin)) { rc = 0; break; }
        prec *= 2;
    }
    mpfr_clear(dist);
    mpfr_clear(margin);
    return rc;
}

#else  /* !USE_MPFR */

/* Without MPFR the HRR engine is unavailable; signal fallback. */
int partitionsp_hrr(unsigned long n, mpz_t out) {
    (void)n; (void)out;
    return -1;
}

#endif /* USE_MPFR */

/* --- builtin ----------------------------------------------------------- */

/* Emit `PartitionsP::argx: PartitionsP called with N arguments; 1 argument
 * is expected.` to stderr and return NULL. */
static Expr* pp_emit_argx(size_t argc) {
    mth_message("PartitionsP", "argx", "PartitionsP called with %zu argument%s; 1 argument is expected.", argc, argc == 1 ? "" : "s");
    return NULL;
}

/* PartitionsP[n] — number of unrestricted partitions of the integer n.
 * Reads `res` only; returns NULL (unevaluated) for symbolic, non-integer or
 * out-of-range arguments. Negative n gives 0. */
Expr* builtin_partitionsp(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return pp_emit_argx(argc);

    return partitions_scalar(false, res);
}

/* ===================================================================== */
/* PartitionsQ[n] — q(n), the number of partitions of n into DISTINCT     */
/* parts (equivalently, into odd parts; OEIS A000009).                    */
/*                                                                        */
/* The same table, then two engines (see partitions_count):               */
/*                                                                        */
/*  1. partitionsq_recurrence — an exact GMP recurrence derived from the  */
/*     Euler identity prod(1-x^k) prod(1+x^k) = prod(1-x^{2k}). It has    */
/*     the same generalized-pentagonal shape as partitionsp_recurrence    */
/*     plus an inhomogeneous term from prod(1-x^{2k}). O(n^1.5) time and  */
/*     O(n sqrt n) bits of memory. Fills the table; the MPFR-less fallback. */
/*                                                                        */
/*  2. partitionsq_hrr — the Hardy-Ramanujan-Rademacher / Hagis exact     */
/*     convergent series for q(n), evaluated in MPFR. Tiny working memory */
/*     (O(sqrt n)-bit precision); used for large n where the recurrence's */
/*     table would be prohibitive. References: P. Hagis, "Partitions into */
/*     odd summands", Amer. J. Math. 85 (1963) 213-222; Rademacher,       */
/*     "Topics in Analytic Number Theory"; eta-quotient Rademacher        */
/*     invariants (Sussman arXiv:1710.03415, Guadalupe arXiv:2404.02367). */
/* ===================================================================== */

/* --- Engine 1: pentagonal-convolution recurrence ----------------------- */

/* Compute q(n) for n >= 0 into the freshly-initialised result `out`.
 * Returns 0 on success, -1 on allocation failure.
 *
 * From prod(1-x^k) * Q(x) = prod(1-x^{2k}) with Q(x) = sum q(n) x^n:
 *     q(m) = r(m) + sum_{k>=1} (-1)^{k+1} [ q(m-g1_k) + q(m-g2_k) ]
 * where g1_k = k(3k-1)/2, g2_k = k(3k+1)/2 are the generalized pentagonal
 * numbers (identical to the p(n) recurrence) and the inhomogeneous term
 *     r(m) = (-1)^k  when m = 2*g_k (i.e. m even and m/2 is generalized
 *                    pentagonal with index k),  else 0
 * is the coefficient of x^m in prod(1-x^{2k}). */
/* Fill q[from..n] (initialised, zero) given q[0..from). */
static void partitionsq_fill(mpz_t* q, unsigned long from, unsigned long n) {
    if (from == 0) { mpz_set_ui(q[0], 1); from = 1; }
    for (unsigned long m = from; m <= n; m++) {
        /* Homogeneous generalized-pentagonal sum (as in partitionsp). */
        for (unsigned long k = 1;; k++) {
            unsigned long g1 = k * (3 * k - 1) / 2;
            if (g1 > m) break;
            int sign = (k & 1UL) ? 1 : -1; /* (-1)^{k+1}: +1 for odd k */
            if (sign > 0) mpz_add(q[m], q[m], q[m - g1]);
            else          mpz_sub(q[m], q[m], q[m - g1]);

            unsigned long g2 = k * (3 * k + 1) / 2;
            if (g2 <= m) {
                if (sign > 0) mpz_add(q[m], q[m], q[m - g2]);
                else          mpz_sub(q[m], q[m], q[m - g2]);
            }
        }

        /* Inhomogeneous term r(m) from prod(1-x^{2k}): nonzero only when
         * m is even and m/2 is a generalized pentagonal number. */
        if ((m & 1UL) == 0) {
            unsigned long pm = m / 2;
            for (unsigned long k = 1;; k++) {
                unsigned long a1 = k * (3 * k - 1) / 2;
                if (a1 > pm) break;
                unsigned long a2 = k * (3 * k + 1) / 2;
                if (a1 == pm || a2 == pm) {
                    if (k & 1UL) mpz_sub_ui(q[m], q[m], 1); /* (-1)^k = -1 */
                    else         mpz_add_ui(q[m], q[m], 1); /* (-1)^k = +1 */
                    break;
                }
            }
        }
    }
}

int partitionsq_recurrence(unsigned long n, mpz_t out) {
    mpz_t* q = malloc((n + 1) * sizeof(mpz_t));
    if (!q) return -1;
    for (unsigned long m = 0; m <= n; m++) mpz_init(q[m]);
    partitionsq_fill(q, 0, n);
    mpz_set(out, q[n]);
    for (unsigned long m = 0; m <= n; m++) mpz_clear(q[m]);
    free(q);
    return 0;
}

/* --- Small-n tables ---------------------------------------------------- */

/* p(n) and q(n) below PARTITIONS_TABLE_N are filled once at startup (well under a millisecond: ~65k additions of
 * numbers of at most 32 digits) and never change, so lookups need no lock.
 * A loop over small n -- Do[PartitionsP[i], {i, 300}] -- used to rerun the
 * O(n^1.5) recurrence per call; it is now a copy. The prefix that fits in
 * int64 (p up to n = 405, q up to n = 769) backs the machine helpers below. */
static mpz_t   g_ptab[PARTITIONS_TABLE_N];
static mpz_t   g_qtab[PARTITIONS_TABLE_N];
static int64_t g_p64[PARTITIONS_TABLE_N];
static int64_t g_q64[PARTITIONS_TABLE_N];
static int64_t g_p64_n, g_q64_n;          /* entries [0, n) fit in int64 */
static bool    g_tables_ready;

static void partitions_tables_build(void) {
    if (g_tables_ready) return;
    for (unsigned long m = 0; m < PARTITIONS_TABLE_N; m++) mpz_init(g_ptab[m]);
    for (unsigned long m = 0; m < PARTITIONS_TABLE_N; m++) mpz_init(g_qtab[m]);
    partitionsp_fill(g_ptab, 0, PARTITIONS_TABLE_N - 1);
    partitionsq_fill(g_qtab, 0, PARTITIONS_TABLE_N - 1);
    g_p64_n = 0;
    while (g_p64_n < PARTITIONS_TABLE_N && ip_mpz_i64(g_ptab[g_p64_n], &g_p64[g_p64_n]))
        g_p64_n++;
    g_q64_n = 0;
    while (g_q64_n < PARTITIONS_TABLE_N && ip_mpz_i64(g_qtab[g_q64_n], &g_q64[g_q64_n]))
        g_q64_n++;
    g_tables_ready = true;
}

bool partitions_p_i64(int64_t n, int64_t* out) {
    if (n < 0) { *out = 0; return true; }
    if (!g_tables_ready || n >= g_p64_n) return false;
    *out = g_p64[n];
    return true;
}

bool partitions_q_i64(int64_t n, int64_t* out) {
    if (n < 0) { *out = 0; return true; }
    if (!g_tables_ready || n >= g_q64_n) return false;
    *out = g_q64[n];
    return true;
}

/* p(n) or q(n) for an n >= 0 into `out` (initialised): the table, then the
 * HRR series, then (no MPFR, or no convergence) the recurrence. */
static int partitions_count(bool distinct, int64_t n, mpz_t out) {
    if (g_tables_ready && n < PARTITIONS_TABLE_N) {
        mpz_set(out, distinct ? g_qtab[n] : g_ptab[n]);
        return 0;
    }
    int rc = distinct ? partitionsq_hrr((unsigned long)n, out)
                      : partitionsp_hrr((unsigned long)n, out);
    if (rc != 0)
        rc = distinct ? partitionsq_recurrence((unsigned long)n, out)
                      : partitionsp_recurrence((unsigned long)n, out);
    return rc;
}

/* The shared scalar body of PartitionsP / PartitionsQ: a machine answer from
 * the int64 helper (the one the packed kernel and Compile use), else GMP. */
static Expr* partitions_scalar(bool distinct, Expr* res) {
    Expr* arg = res->data.function.args[0];
    /* A negative big integer has no partitions (PartitionsP[-10^20] is 0,
     * as for every negative n); a positive one has an astronomically large
     * answer. */
    if (arg->type == EXPR_BIGINT && mpz_sgn(arg->data.bigint) < 0) return expr_new_integer(0);
    if (arg->type != EXPR_INTEGER) return NULL;
    int64_t n = arg->data.integer, v;
    if (distinct ? partitions_q_i64(n, &v) : partitions_p_i64(n, &v))
        return expr_new_integer(v);
    /* The engines index by unsigned long, 32 bits on LLP64 Windows. */
    if ((uint64_t)n > (uint64_t)ULONG_MAX) return NULL;
    mpz_t out;
    mpz_init(out);
    if (partitions_count(distinct, n, out) != 0) { mpz_clear(out); return NULL; }
    Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(out));
    mpz_clear(out);
    return result;
}

/* --- Engine 2: Hardy-Ramanujan-Rademacher / Hagis (MPFR) --------------- */
#ifdef USE_MPFR

/* gcd for the coprimality test in the character sum. */
static unsigned long pq_gcd(unsigned long a, unsigned long b) {
    while (b) { unsigned long t = a % b; a = b; b = t; }
    return a;
}

/* Modified Bessel function I_1(x) for real x >= 0, into `out`.
 * MPFR provides only the unmodified Bessel J/Y, so we sum the everywhere-
 * positive power series I_1(x) = sum_{m>=0} (x/2)^{2m+1} / (m! (m+1)!),
 * advanced by the term ratio (x/2)^2 / (m(m+1)). All terms are positive,
 * so there is no cancellation even for large x; we stop once a term is
 * `prec`+ guard bits below the running sum. */
static void i1_mpfr(mpfr_t out, const mpfr_t x, mpfr_prec_t prec) {
    if (mpfr_zero_p(x)) { mpfr_set_zero(out, 1); return; }

    mpfr_t half, hx2, term, sum;
    mpfr_init2(half, prec);
    mpfr_init2(hx2, prec);
    mpfr_init2(term, prec);
    mpfr_init2(sum, prec);

    mpfr_div_2ui(half, x, 1, MPFR_RNDN);   /* x/2 */
    mpfr_mul(hx2, half, half, MPFR_RNDN);  /* (x/2)^2 */
    mpfr_set(term, half, MPFR_RNDN);       /* m = 0 term: x/2 */
    mpfr_set(sum, term, MPFR_RNDN);

    for (unsigned long m = 1; m < 1000000UL; m++) {
        mpfr_mul(term, term, hx2, MPFR_RNDN);
        mpfr_div_ui(term, term, m, MPFR_RNDN);
        mpfr_div_ui(term, term, m + 1, MPFR_RNDN);
        mpfr_add(sum, sum, term, MPFR_RNDN);
        if (mpfr_zero_p(term)) break;
        if (mpfr_get_exp(sum) - mpfr_get_exp(term) > (mpfr_exp_t)prec + 16) break;
    }

    mpfr_set(out, sum, MPFR_RNDN);
    mpfr_clear(half);
    mpfr_clear(hx2);
    mpfr_clear(term);
    mpfr_clear(sum);
}

/* Dedekind sum s(h,k) as an exact rational, for 0 <= h < k, gcd(h,k)=1.
 * Computed by reciprocity in O(log k) steps:
 *     s(h,k) = -1/4 + (1/12)(h/k + k/h + 1/(hk)) - s(k mod h, h),
 * with base case s(0,k) = 0. */
static void pq_dedekind_sum(mpq_t out, unsigned long h, unsigned long k) {
    if (h == 0) { mpq_set_ui(out, 0, 1); return; }

    mpq_t acc, tmp, rec;
    mpq_init(acc);
    mpq_init(tmp);
    mpq_init(rec);

    /* acc = h/k + k/h + 1/(h k) */
    mpq_set_ui(acc, h, k);          mpq_canonicalize(acc);
    mpq_set_ui(tmp, k, h);          mpq_canonicalize(tmp);
    mpq_add(acc, acc, tmp);
    mpq_set_ui(tmp, 1, h * k);      mpq_canonicalize(tmp);
    mpq_add(acc, acc, tmp);

    /* acc *= 1/12 ; acc += -1/4 */
    mpq_set_ui(tmp, 1, 12);         mpq_canonicalize(tmp);
    mpq_mul(acc, acc, tmp);
    mpq_set_si(tmp, -1, 4);         mpq_canonicalize(tmp);
    mpq_add(acc, acc, tmp);

    pq_dedekind_sum(rec, k % h, h);
    mpq_sub(out, acc, rec);

    mpq_clear(acc);
    mpq_clear(tmp);
    mpq_clear(rec);
}

/* 6k s(h, k) for gcd(h, k) = 1. It is an integer (Rademacher) of size at
 * most k^2/2, so the Dedekind sum is run by reciprocity in double and the
 * product rounded: for the k this series reaches the rounding error is
 * below 10^-6. A result not within 1/8 of an integer takes the exact mpq
 * route instead, so the shortcut can never return a wrong value. */
static int64_t pq_dedekind6k(unsigned long h, unsigned long k) {
    double s = 0.0, sign = 1.0;
    unsigned long a = h % k, b = k;
    while (a != 0) {                 /* s(a,b) = T(a,b) - s(b mod a, a) */
        double ad = (double)a, bd = (double)b;
        s += sign * ((ad / bd + bd / ad + 1.0 / (ad * bd)) / 12.0 - 0.25);
        sign = -sign;
        unsigned long t = b % a;
        b = a;
        a = t;
    }
    double v = 6.0 * (double)k * s;
    double r = floor(v + 0.5);
    if (fabs(v - r) < 0.125) return (int64_t)r;
    mpq_t q;
    mpq_init(q);
    pq_dedekind_sum(q, h, k);
    mpz_mul_ui(mpq_numref(q), mpq_numref(q), 6 * k);
    mpq_canonicalize(q);                       /* now an integer */
    int64_t out = (int64_t)mpz_get_d(mpq_numref(q));
    mpq_clear(q);
    return out;
}

/* The Hagis character sum (the analogue of Rademacher's A_k(n) for q(n)):
 *     A_k(n) = sum_{h mod k, gcd(h,k)=1} cos( pi(s(h,k) - s(2h,k)) - 2 pi n h/k ).
 * k is odd, so gcd(2h, k) = 1 whenever gcd(h, k) = 1. Every angle is the
 * exact rational pi R / (6k) with
 *     R = 6k s(h,k) - 6k s(2h,k) - 12 (n h mod k)   reduced mod 12k,
 * so no precision is spent on argument reduction. With `prec` == 0 the sum is
 * taken in double (the caller has checked the term is too small for its
 * error to matter); otherwise in MPFR at `prec` bits. */
static void partitionsq_Ak(mpfr_t out, unsigned long k, unsigned long n,
                           mpfr_prec_t prec) {
    if (k == 1) { mpfr_set_ui(out, 1, MPFR_RNDN); return; }
    unsigned long nmod = n % k;
    int64_t m12 = 12 * (int64_t)k;
    if (prec == 0) {
        double acc = 0.0;
        for (unsigned long h = 1; h < k; h++) {
            if (pq_gcd(h, k) != 1) continue;
            int64_t R = pq_dedekind6k(h, k) - pq_dedekind6k((2 * h) % k, k)
                      - 12 * (int64_t)((nmod * h) % k);
            R %= m12;
            if (R < 0) R += m12;
            acc += cos(M_PI * (double)R / (6.0 * (double)k));
        }
        mpfr_set_d(out, acc, MPFR_RNDN);
        return;
    }
    mpfr_t pi, ang, c;
    mpfr_inits2(prec, pi, ang, c, (mpfr_ptr)0);
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_set_zero(out, 1);
    for (unsigned long h = 1; h < k; h++) {
        if (pq_gcd(h, k) != 1) continue;
        int64_t R = pq_dedekind6k(h, k) - pq_dedekind6k((2 * h) % k, k)
                  - 12 * (int64_t)((nmod * h) % k);
        R %= m12;
        if (R < 0) R += m12;
        mpfr_mul_ui(ang, pi, (unsigned long)R, MPFR_RNDN);
        mpfr_div_ui(ang, ang, 6 * k, MPFR_RNDN);
        mpfr_cos(c, ang, MPFR_RNDN);
        mpfr_add(out, out, c, MPFR_RNDN);
    }
    mpfr_clears(pi, ang, c, (mpfr_ptr)0);
}

/* Evaluate the Hagis series for q(n) over odd k = 1,3,...,<= Kmax, rounding
 * the real sum to the nearest integer in `out`:
 *     q(n) = (pi/sqrt(24n+1)) sum_{k odd} (1/k) A_k(n) I_1( pi sqrt(48n+2)/(12k) ).
 * Each term gets only the precision it needs: |term_k| <= pref e^(B/k), so
 * term k carries about log2(pref) + B/(k ln 2) bits above the units place.
 * Terms whose double-precision error (|A_k| <= k cosines of 2^-50 each)
 * stays below 2^-12 / Kmax are summed in double; the rest in MPFR at their
 * own magnitude plus guard bits. The total is accumulated at `prec`.
 * Stores in `dist` the distance from the sum to that integer and in
 * `last_mag` the magnitude of the last term added (a tail estimate). */
static void hrr_q_eval(unsigned long n, unsigned long Kmax, mpfr_prec_t prec,
                       mpz_t out, mpfr_t dist, mpfr_t last_mag) {
    mpfr_t pi, pref, B, arg, Ak, I1, term, sum, rounded, tmp;
    mpfr_inits2(prec, pi, pref, B, sum, rounded, tmp, (mpfr_ptr)0);
    mpfr_inits2(64, arg, Ak, I1, term, (mpfr_ptr)0);

    mpfr_const_pi(pi, MPFR_RNDN);

    /* pref = pi / sqrt(24n + 1) */
    mpfr_set_ui(tmp, 24, MPFR_RNDN);
    mpfr_mul_ui(tmp, tmp, n, MPFR_RNDN);
    mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN);
    mpfr_sqrt(tmp, tmp, MPFR_RNDN);
    mpfr_div(pref, pi, tmp, MPFR_RNDN);

    /* B = pi sqrt(48n + 2) / 12 ; the Bessel argument is B/k */
    mpfr_set_ui(B, 48, MPFR_RNDN);
    mpfr_mul_ui(B, B, n, MPFR_RNDN);
    mpfr_add_ui(B, B, 2, MPFR_RNDN);
    mpfr_sqrt(B, B, MPFR_RNDN);
    mpfr_mul(B, B, pi, MPFR_RNDN);
    mpfr_div_ui(B, B, 12, MPFR_RNDN);

    double Bd = mpfr_get_d(B, MPFR_RNDN);
    double lpref = log2(mpfr_get_d(pref, MPFR_RNDN));
    double lK = log2((double)Kmax + 1.0);

    mpfr_set_zero(sum, 1);
    mpfr_set_zero(last_mag, 1);
    for (unsigned long k = 1; k <= Kmax; k += 2) {
        double mag = lpref + Bd / (double)k / 0.6931471805599453;  /* log2 bound */
        bool dbl = mag < 38.0 - lK;
        mpfr_prec_t pk = dbl ? 64 : (mpfr_prec_t)(mag + log2((double)k) + lK + 40.0);
        if (pk < 64) pk = 64;
        if (pk > prec) pk = prec;
        mpfr_set_prec(Ak, pk);
        mpfr_set_prec(arg, pk);
        mpfr_set_prec(I1, pk);
        mpfr_set_prec(term, pk);

        partitionsq_Ak(Ak, k, n, dbl ? 0 : pk);
        if (mpfr_zero_p(Ak)) continue;

        mpfr_div_ui(arg, B, k, MPFR_RNDN);
        i1_mpfr(I1, arg, pk);

        /* term = (1/k) A_k(n) I_1(B/k) */
        mpfr_mul(term, Ak, I1, MPFR_RNDN);
        mpfr_div_ui(term, term, k, MPFR_RNDN);
        mpfr_add(sum, sum, term, MPFR_RNDN);

        mpfr_abs(last_mag, term, MPFR_RNDN);
    }

    mpfr_mul(sum, sum, pref, MPFR_RNDN);
    mpfr_mul(last_mag, last_mag, pref, MPFR_RNDN);

    mpfr_round(rounded, sum);
    mpfr_sub(dist, sum, rounded, MPFR_RNDN);
    mpfr_abs(dist, dist, MPFR_RNDN);
    mpfr_get_z(out, rounded, MPFR_RNDN);

    mpfr_clears(pi, pref, B, sum, rounded, tmp, arg, Ak, I1, term, (mpfr_ptr)0);
}


/* Compute q(n) for n >= 2 via the Hagis series into `out`. Returns 0 on
 * success, -1 if it could not converge to a confident integer (the caller
 * then falls back to the recurrence).
 *
 * Unlike p(n), no clean closed-form Rademacher remainder bound with
 * explicit constants is available for this series, so termination is
 * convergence-driven: we sum ~sqrt(n) odd-k terms and accept the rounding
 * only when the sum lies within 1/4 of an integer and the last term is
 * below 1/8 (so the unsummed tail cannot flip the rounding), growing the
 * term count and precision otherwise. Correctness across the supported
 * range is pinned by the exhaustive HRR-equals-recurrence cross-check in
 * the unit tests. */
int partitionsq_hrr(unsigned long n, mpz_t out) {
    if (n < 2) return -1;

    /* Working precision: q(n) ~ exp(pi sqrt(n/3)) has ~ pi sqrt(n/3)/ln2
     * bits; add generous guard bits for accumulation and rounding. */
    double e = M_PI * sqrt((double)n / 3.0);
    mpfr_prec_t prec = (mpfr_prec_t)(e / 0.6931471805599453) + 96;
    if (prec < 96) prec = 96;

    unsigned long base = (unsigned long)ceil(sqrt((double)n));
    unsigned long Kmax = 2 * base + 9;     /* ~ sqrt(n) odd terms to start */

    mpfr_t dist, last_mag, q14, q18;
    mpfr_init2(dist, 64);
    mpfr_init2(last_mag, 64);
    mpfr_init2(q14, 64);
    mpfr_init2(q18, 64);
    mpfr_set_d(q14, 0.25, MPFR_RNDN);
    mpfr_set_d(q18, 0.125, MPFR_RNDN);

    int rc = -1;
    for (int attempt = 0; attempt < 16; attempt++) {
        hrr_q_eval(n, Kmax | 1UL, prec, out, dist, last_mag);
        if (mpfr_less_p(dist, q14) && mpfr_less_p(last_mag, q18)) {
            rc = 0;
            break;
        }
        Kmax = Kmax + Kmax / 2 + 4;          /* add more terms */
        if (!mpfr_less_p(dist, q14)) prec *= 2; /* and precision if margin thin */
    }

    mpfr_clear(dist);
    mpfr_clear(last_mag);
    mpfr_clear(q14);
    mpfr_clear(q18);
    return rc;
}

#else  /* !USE_MPFR */

/* Without MPFR the HRR engine is unavailable; signal fallback. */
int partitionsq_hrr(unsigned long n, mpz_t out) {
    (void)n; (void)out;
    return -1;
}

#endif /* USE_MPFR */

/* Emit `PartitionsQ::argx: PartitionsQ called with N arguments; 1 argument
 * is expected.` to stderr and return NULL. */
static Expr* qq_emit_argx(size_t argc) {
    mth_message("PartitionsQ", "argx", "PartitionsQ called with %zu argument%s; 1 argument is expected.", argc, argc == 1 ? "" : "s");
    return NULL;
}

/* PartitionsQ[n] — number of partitions of n into distinct parts.
 * Reads `res` only; returns NULL (unevaluated) for symbolic, non-integer or
 * big-integer arguments. Negative n gives 0. */
Expr* builtin_partitionsq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return qq_emit_argx(argc);
    return partitions_scalar(true, res);
}

void partitions_init(void) {
    partitions_tables_build();
    symtab_add_builtin("IntegerPartitions", builtin_integerpartitions);
    symtab_get_def("IntegerPartitions")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("PartitionsP", builtin_partitionsp);
    symtab_get_def("PartitionsP")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);

    symtab_add_builtin("PartitionsQ", builtin_partitionsq);
    symtab_get_def("PartitionsQ")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);
}

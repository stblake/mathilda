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
 *   IntegerPartitions[n, {kmin, kmax, dk}]   kmin, kmin+dk, ... parts
 *   IntegerPartitions[n, kspec, sspec]       parts drawn only from sspec
 *   IntegerPartitions[n, kspec, sspec, m]    first m (m>0) or last |m| (m<0)
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* M_PI is POSIX, not C99; glibc hides it under -std=c99. */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

/* floor(a / b) clamped to [0, LONG_MAX]; b must be non-zero. */
static long ip_floor_div(const mpq_t a, const mpq_t b) {
    mpq_t q;
    mpq_init(q);
    mpq_div(q, a, b);
    mpz_t num, den, fl;
    mpz_inits(num, den, fl, NULL);
    mpq_get_num(num, q);
    mpq_get_den(den, q);
    mpz_fdiv_q(fl, num, den);          /* floor */
    long r;
    if (mpz_sgn(fl) < 0)            r = 0;
    else if (!mpz_fits_slong_p(fl)) r = LONG_MAX;
    else                            r = mpz_get_si(fl);
    mpz_clears(num, den, fl, NULL);
    mpq_clear(q);
    return r;
}

static bool ip_is_sym(const Expr* e, const char* s) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == s;
}

/* ----- enumeration context --------------------------------------------- */

typedef struct {
    mpq_t*  parts;          /* P: reversed sspec (display priority order)   */
    bool*   suffix_nonneg;  /* all parts[j>=i] have sign >= 0               */
    bool*   suffix_nonpos;  /* all parts[j>=i] have sign <= 0               */
    long*   counts;         /* current count chosen per part                */
    size_t  nparts;

    long    kmin, kmax, dk; /* length constraints                           */
    bool    kmax_inf;

    Expr**  out;            /* collected partitions, in order               */
    size_t  out_n, out_cap;
    long    limit;          /* stop after this many (first-m); 0 = unlimited*/
    bool    stop;
    bool    error;          /* allocation failure                           */
} PartCtx;

static bool ip_len_ok(const PartCtx* c, long len) {
    if (len < c->kmin) return false;
    if (!c->kmax_inf && len > c->kmax) return false;
    long step = c->dk > 0 ? c->dk : 1;
    return ((len - c->kmin) % step) == 0;
}

/* Materialise the current count vector as List[...] and append it. */
static void ip_emit(PartCtx* c, long len) {
    Expr** elems = NULL;
    if (len > 0) {
        elems = malloc((size_t)len * sizeof(Expr*));
        if (!elems) { c->error = true; return; }
        size_t idx = 0;
        for (size_t i = 0; i < c->nparts; i++) {
            for (long k = 0; k < c->counts[i]; k++) {
                elems[idx++] = ip_mpq_to_expr(c->parts[i]);
            }
        }
    }
    Expr* part = expr_new_function(expr_new_symbol(SYM_List), elems, (size_t)len);
    free(elems);

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

/* Assign a count to parts[i], recursing in descending-count order. */
static void ip_recurse(PartCtx* c, size_t i, const mpq_t remaining, long len) {
    if (c->stop || c->error) return;

    if (i == c->nparts) {
        if (mpq_sgn(remaining) == 0 && ip_len_ok(c, len)) ip_emit(c, len);
        return;
    }

    const mpq_t* p = (const mpq_t*)&c->parts[i];
    int psgn = mpq_sgn(c->parts[i]);
    int rsgn = mpq_sgn(remaining);

    /* Length budget from a finite kmax. */
    long budget = -1;                       /* -1 == unlimited */
    if (!c->kmax_inf) {
        budget = c->kmax - len;
        if (budget < 0) budget = 0;
    }

    /* Tight numeric bound: only valid when the remaining parts are single
     * signed in p's direction, so contributions cannot overshoot and recover. */
    long nb = -1;                           /* -1 == unbounded */
    if (psgn > 0 && c->suffix_nonneg[i]) {
        nb = (rsgn <= 0) ? 0 : ip_floor_div(remaining, *p);
    } else if (psgn < 0 && c->suffix_nonpos[i]) {
        nb = (rsgn >= 0) ? 0 : ip_floor_div(remaining, *p);
    }

    long cmax;
    if (nb >= 0 && budget >= 0)      cmax = (nb < budget) ? nb : budget;
    else if (nb >= 0)                cmax = nb;
    else if (budget >= 0)            cmax = budget;
    else { c->error = true; return; }       /* unreachable: undef gate guards */

    mpq_t step, rem2;
    mpq_init(step);
    mpq_init(rem2);
    for (long cc = cmax; cc >= 0; cc--) {
        mpq_set_si(step, cc, 1);
        mpq_mul(step, step, *p);            /* cc * p */
        mpq_sub(rem2, remaining, step);     /* remaining - cc*p */
        c->counts[i] = cc;
        ip_recurse(c, i + 1, rem2, len + cc);
        if (c->stop || c->error) break;
    }
    mpq_clear(step);
    mpq_clear(rem2);
}

/* ----- diagnostics ----------------------------------------------------- */

/* "IntegerPartitions[a0,a1,...]" for the first `nargs` arguments of res. */
static char* ip_call_str(const Expr* res, size_t nargs) {
    size_t cap = 64, len = 0;
    char* buf = malloc(cap);
    if (!buf) return NULL;
    len += (size_t)snprintf(buf, cap, "IntegerPartitions[");
    for (size_t i = 0; i < nargs; i++) {
        char* a = expr_to_string(res->data.function.args[i]);
        const char* piece = a ? a : "?";
        size_t need = len + strlen(piece) + 2; /* piece + ',' or ']' + NUL */
        if (need > cap) {
            while (need > cap) cap *= 2;
            char* g = realloc(buf, cap);
            if (!g) { free(buf); free(a); return NULL; }
            buf = g;
        }
        if (i) buf[len++] = ',';
        memcpy(buf + len, piece, strlen(piece));
        len += strlen(piece);
        free(a);
    }
    buf[len++] = ']';
    buf[len] = '\0';
    return buf;
}

/* ----- argument parsing ------------------------------------------------ */

/* Read a length value (Integer/BigInt) or Infinity. Returns false for
 * anything else. Huge bigints clamp to LONG_MAX. */
static bool ip_read_len(const Expr* e, long* val, bool* inf) {
    if (ip_is_sym(e, SYM_Infinity)) { *inf = true; *val = 0; return true; }
    if (e->type == EXPR_INTEGER) { *inf = false; *val = (long)e->data.integer; return true; }
    if (e->type == EXPR_BIGINT) {
        *inf = false;
        *val = mpz_fits_slong_p(e->data.bigint) ? mpz_get_si(e->data.bigint)
             : (mpz_sgn(e->data.bigint) < 0 ? LONG_MIN : LONG_MAX);
        return true;
    }
    return false;
}

void partitions_init(void);

Expr* builtin_integerpartitions(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc < 1 || argc > 4) {
        fprintf(stderr,
                "IntegerPartitions::argb: IntegerPartitions called with %zu "
                "argument%s; between 1 and 4 arguments are expected.\n",
                argc, argc == 1 ? "" : "s");
        return NULL;
    }

    /* --- n --- */
    mpq_t N;
    if (!ip_to_mpq(args[0], N)) return NULL;   /* symbolic / real: unevaluated */

    /* --- kspec --- */
    long kmin = 0, kmax = 0, dk = 1;
    bool kmax_inf = true;
    if (argc >= 2) {
        Expr* ks = args[1];
        if (ip_is_sym(ks, SYM_All) || ip_is_sym(ks, SYM_Infinity)) {
            kmin = 0; kmax_inf = true; dk = 1;
        } else if (ks->type == EXPR_INTEGER || ks->type == EXPR_BIGINT) {
            long v; bool inf;
            ip_read_len(ks, &v, &inf);
            kmin = 0; kmax = v; kmax_inf = false; dk = 1;
        } else if (ks->type == EXPR_FUNCTION &&
                   ip_is_sym(ks->data.function.head, SYM_List)) {
            size_t kc = ks->data.function.arg_count;
            Expr** ka = ks->data.function.args;
            long v0 = 0, v1 = 0, v2 = 1; bool i0 = false, i1 = false, i2 = false;
            if (kc == 1) {
                if (!ip_read_len(ka[0], &v0, &i0)) { mpq_clear(N); return NULL; }
                if (i0) { kmin = 0; kmax_inf = true; }
                else    { kmin = kmax = v0; kmax_inf = false; }
                dk = 1;
            } else if (kc == 2 || kc == 3) {
                if (!ip_read_len(ka[0], &v0, &i0) ||
                    !ip_read_len(ka[1], &v1, &i1)) { mpq_clear(N); return NULL; }
                if (kc == 3 && !ip_read_len(ka[2], &v2, &i2)) { mpq_clear(N); return NULL; }
                if (i0) { mpq_clear(N); return NULL; }  /* kmin must be finite */
                kmin = v0;
                kmax_inf = i1;
                kmax = v1;
                dk = (kc == 3) ? v2 : 1;
                if (dk <= 0) dk = 1;
            } else {
                mpq_clear(N); return NULL;
            }
        } else {
            mpq_clear(N); return NULL;          /* invalid kspec: unevaluated */
        }
    }

    /* --- sspec: build P (reversed) directly --- */
    mpq_t* parts = NULL;
    size_t nparts = 0;
    bool has_pos = false, has_neg = false, has_zero = false;

    bool sspec_all = (argc < 3) || ip_is_sym(args[2], SYM_All);
    if (sspec_all) {
        /* Range[n] reversed = floor(n), floor(n)-1, ..., 1 (empty if n<1). */
        mpz_t fl;
        mpz_init(fl);
        mpz_fdiv_q(fl, mpq_numref(N), mpq_denref(N));   /* floor(n) */
        if (mpz_sgn(fl) >= 1 && mpz_fits_slong_p(fl)) {
            long f = mpz_get_si(fl);
            nparts = (size_t)f;
            parts = malloc(nparts * sizeof(mpq_t));
            if (!parts) { mpz_clear(fl); mpq_clear(N); return NULL; }
            for (size_t i = 0; i < nparts; i++) {
                mpq_init(parts[i]);
                mpq_set_si(parts[i], f - (long)i, 1);   /* f, f-1, ..., 1 */
            }
            has_pos = (nparts > 0);
        }
        mpz_clear(fl);
    } else if (args[2]->type == EXPR_FUNCTION &&
               ip_is_sym(args[2]->data.function.head, SYM_List)) {
        size_t r = args[2]->data.function.arg_count;
        Expr** sa = args[2]->data.function.args;
        if (r > 0) {
            parts = malloc(r * sizeof(mpq_t));
            if (!parts) { mpq_clear(N); return NULL; }
        }
        for (size_t i = 0; i < r; i++) {
            /* reversed: P[i] = sspec[r-1-i] */
            if (!ip_to_mpq(sa[r - 1 - i], parts[i])) {
                for (size_t j = 0; j < i; j++) mpq_clear(parts[j]);
                free(parts);
                mpq_clear(N);
                return NULL;                    /* non-numeric s_i: unevaluated */
            }
            int sg = mpq_sgn(parts[i]);
            if (sg > 0) has_pos = true;
            else if (sg < 0) has_neg = true;
            else has_zero = true;
        }
        nparts = r;
    } else {
        mpq_clear(N); return NULL;               /* invalid sspec: unevaluated */
    }

    /* --- m --- */
    bool m_all = (argc < 4);
    long m = 0;
    if (argc >= 4) {
        Expr* me = args[3];
        if (ip_is_sym(me, SYM_All) || ip_is_sym(me, SYM_Infinity)) {
            m_all = true;
        } else if (me->type == EXPR_INTEGER || me->type == EXPR_BIGINT) {
            long v; bool inf;
            ip_read_len(me, &v, &inf);
            m = v;
        } else {
            for (size_t i = 0; i < nparts; i++) mpq_clear(parts[i]);
            free(parts);
            mpq_clear(N);
            return NULL;                          /* invalid m: unevaluated */
        }
    }

    /* --- infinite-result guard (::undef) --- */
    if (kmax_inf && (has_zero || (has_pos && has_neg))) {
        char* call = ip_call_str(res, argc);
        fprintf(stderr,
                "IntegerPartitions::undef: %s contains partitions that are "
                "undefined because they are infinitely large.\n",
                call ? call : "IntegerPartitions[...]");
        free(call);
        for (size_t i = 0; i < nparts; i++) mpq_clear(parts[i]);
        free(parts);
        mpq_clear(N);
        return NULL;
    }

    /* --- precompute suffix sign flags --- */
    bool* suffix_nonneg = nparts ? malloc(nparts * sizeof(bool)) : NULL;
    bool* suffix_nonpos = nparts ? malloc(nparts * sizeof(bool)) : NULL;
    long* counts        = nparts ? malloc(nparts * sizeof(long)) : NULL;
    if (nparts && (!suffix_nonneg || !suffix_nonpos || !counts)) {
        free(suffix_nonneg); free(suffix_nonpos); free(counts);
        for (size_t i = 0; i < nparts; i++) mpq_clear(parts[i]);
        free(parts);
        mpq_clear(N);
        return NULL;
    }
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

    /* --- enumerate --- */
    PartCtx ctx;
    ctx.parts = parts;
    ctx.suffix_nonneg = suffix_nonneg;
    ctx.suffix_nonpos = suffix_nonpos;
    ctx.counts = counts;
    ctx.nparts = nparts;
    ctx.kmin = kmin;
    ctx.kmax = kmax;
    ctx.dk = dk;
    ctx.kmax_inf = kmax_inf;
    ctx.out = NULL;
    ctx.out_n = 0;
    ctx.out_cap = 0;
    ctx.limit = (!m_all && m > 0) ? m : 0;     /* first-m early stop */
    ctx.stop = false;
    ctx.error = false;

    if (!(!m_all && m == 0)) {                  /* m == 0 -> empty, skip work */
        ip_recurse(&ctx, 0, N, 0);
    }

    /* free working state */
    free(suffix_nonneg);
    free(suffix_nonpos);
    free(counts);
    for (size_t i = 0; i < nparts; i++) mpq_clear(parts[i]);
    free(parts);
    mpq_clear(N);

    if (ctx.error) {
        for (size_t i = 0; i < ctx.out_n; i++) expr_free(ctx.out[i]);
        free(ctx.out);
        return NULL;
    }

    /* --- apply m (slice / ::take) --- */
    size_t total = ctx.out_n;
    Expr** sel = ctx.out;
    size_t sel_n = total;

    if (!m_all && m == 0) {
        sel = NULL; sel_n = 0;                  /* {} */
    } else if (!m_all && m > 0) {
        if (total < (size_t)m) {
            char* call = ip_call_str(res, 3);
            fprintf(stderr,
                    "IntegerPartitions::take: Warning: not all elements were "
                    "found when attempting to take the sequence {1,%ld,1} from "
                    "%s, which has length %zu.\n",
                    m, call ? call : "IntegerPartitions[...]", total);
            free(call);
        }
        /* out already holds the first min(m,total) in order */
    } else if (!m_all && m < 0) {
        size_t want = (size_t)(-(long long)m);
        if (want >= total) {
            if (want > total) {
                char* call = ip_call_str(res, 3);
                fprintf(stderr,
                        "IntegerPartitions::take: Warning: not all elements were "
                        "found when attempting to take the sequence {%ld,-1,1} "
                        "from %s, which has length %zu.\n",
                        m, call ? call : "IntegerPartitions[...]", total);
                free(call);
            }
        } else {
            size_t drop = total - want;
            for (size_t i = 0; i < drop; i++) expr_free(ctx.out[i]);
            memmove(ctx.out, ctx.out + drop, want * sizeof(Expr*));
            sel_n = want;
        }
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), sel, sel_n);
    free(ctx.out);
    return result;
}

/* ===================================================================== */
/* PartitionsP[n] — the partition-counting function p(n).                 */
/*                                                                        */
/* Two engines, dispatched by size (see builtin_partitionsp):             */
/*                                                                        */
/*  1. partitionsp_recurrence — Euler's pentagonal-number-theorem         */
/*     recurrence, exact GMP integers. O(n^1.5) time and O(n*sqrt n) bits */
/*     of memory. Simple and robust; used for small n.                    */
/*                                                                        */
/*  2. partitionsp_hrr — the non-recursive Hardy-Ramanujan-Rademacher     */
/*     exact formula evaluated in MPFR. O(n^{3/4}) cosine evaluations and */
/*     O(sqrt n)-bit working precision, so memory stays tiny. Used for    */
/*     large n, where the recurrence's table would be prohibitive.        */
/*     Reference: F. Johansson, "Efficient implementation of the          */
/*     Hardy-Ramanujan-Rademacher formula" (arXiv:1205.5991).             */
/* ===================================================================== */

/* Threshold (in n) above which the HRR engine is preferred. Below it the
 * recurrence is faster and needs no MPFR. The documented big values
 * (p(1024), p(2048), p(4096)) consequently exercise the HRR path. */
#define PARTITIONSP_HRR_THRESHOLD 1000

/* --- Engine 1: pentagonal recurrence ----------------------------------- */

/* Compute p(n) for n >= 0 into the freshly-initialised result `out`.
 * Returns 0 on success, -1 on allocation failure (out left at 0/uninit). */
int partitionsp_recurrence(unsigned long n, mpz_t out) {
    mpz_t* p = malloc((n + 1) * sizeof(mpz_t));
    if (!p) return -1;

    for (unsigned long m = 0; m <= n; m++) mpz_init(p[m]);
    mpz_set_ui(p[0], 1);

    for (unsigned long m = 1; m <= n; m++) {
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
    fprintf(stderr,
            "PartitionsP::argx: PartitionsP called with %zu argument%s; "
            "1 argument is expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/* PartitionsP[n] — number of unrestricted partitions of the integer n.
 * Reads `res` only; returns NULL (unevaluated) for symbolic, non-integer or
 * out-of-range arguments. Negative n gives 0. */
Expr* builtin_partitionsp(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return pp_emit_argx(argc);

    Expr* arg = res->data.function.args[0];

    /* A big-integer argument is astronomically large; we cannot build a table
     * of that size. Leave the call unevaluated. */
    if (!expr_is_integer_like(arg) || arg->type == EXPR_BIGINT) return NULL;

    int64_t n = arg->data.integer;
    if (n < 0) return expr_new_integer(0); /* p(n) = 0 for n < 0 */

    mpz_t out;
    mpz_init(out);

    int rc = -1;
    if (n >= PARTITIONSP_HRR_THRESHOLD) {
        rc = partitionsp_hrr((unsigned long)n, out);
        if (rc != 0)                       /* fall back if HRR unavailable */
            rc = partitionsp_recurrence((unsigned long)n, out);
    } else {
        rc = partitionsp_recurrence((unsigned long)n, out);
    }
    if (rc != 0) {
        mpz_clear(out);
        return NULL;
    }

    Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(out));
    mpz_clear(out);
    return result;
}

/* ===================================================================== */
/* PartitionsQ[n] — q(n), the number of partitions of n into DISTINCT     */
/* parts (equivalently, into odd parts; OEIS A000009).                    */
/*                                                                        */
/* Two engines, dispatched by size (see builtin_partitionsq):             */
/*                                                                        */
/*  1. partitionsq_recurrence — an exact GMP recurrence derived from the  */
/*     Euler identity prod(1-x^k) prod(1+x^k) = prod(1-x^{2k}). It has    */
/*     the same generalized-pentagonal shape as partitionsp_recurrence    */
/*     plus an inhomogeneous term from prod(1-x^{2k}). O(n^1.5) time and  */
/*     O(n sqrt n) bits of memory. Simple and robust; used for small n.   */
/*                                                                        */
/*  2. partitionsq_hrr — the Hardy-Ramanujan-Rademacher / Hagis exact     */
/*     convergent series for q(n), evaluated in MPFR. Tiny working memory */
/*     (O(sqrt n)-bit precision); used for large n where the recurrence's */
/*     table would be prohibitive. References: P. Hagis, "Partitions into */
/*     odd summands", Amer. J. Math. 85 (1963) 213-222; Rademacher,       */
/*     "Topics in Analytic Number Theory"; eta-quotient Rademacher        */
/*     invariants (Sussman arXiv:1710.03415, Guadalupe arXiv:2404.02367). */
/* ===================================================================== */

/* Threshold (in n) above which the HRR engine is preferred. Below it the
 * recurrence is faster and needs no MPFR. */
#define PARTITIONSQ_HRR_THRESHOLD 1000

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
int partitionsq_recurrence(unsigned long n, mpz_t out) {
    mpz_t* q = malloc((n + 1) * sizeof(mpz_t));
    if (!q) return -1;

    for (unsigned long m = 0; m <= n; m++) mpz_init(q[m]);
    mpz_set_ui(q[0], 1);

    for (unsigned long m = 1; m <= n; m++) {
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

    mpz_set(out, q[n]);
    for (unsigned long m = 0; m <= n; m++) mpz_clear(q[m]);
    free(q);
    return 0;
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

/* The Hagis character sum (the analog of the Rademacher A_k(n) for q(n)):
 *     A_k(n) = sum_{h mod k, gcd(h,k)=1} cos( pi(s(h,k) - s(2h,k)) - 2pi n h/k )
 * which is real (the imaginary parts cancel). k is odd, so gcd(2h,k)=1
 * whenever gcd(h,k)=1. A_1(n) = 1 (the single residue h = 0). */
static void partitionsq_Ak(mpfr_t out, unsigned long k, unsigned long n,
                           mpfr_prec_t prec) {
    mpfr_t pi, ang, c;
    mpfr_init2(pi, prec);
    mpfr_init2(ang, prec);
    mpfr_init2(c, prec);
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_set_zero(out, 1);

    mpq_t s1, s2, sd;
    mpq_init(s1);
    mpq_init(s2);
    mpq_init(sd);

    unsigned long nmod = n % k;
    for (unsigned long h = 0; h < k; h++) {
        if (pq_gcd(h, k) != 1) continue;   /* h = 0 only contributes for k = 1 */

        pq_dedekind_sum(s1, h, k);
        pq_dedekind_sum(s2, (2 * h) % k, k);
        mpq_sub(sd, s1, s2);

        /* ang = pi * (s(h,k) - s(2h,k)) */
        mpfr_set_q(ang, sd, MPFR_RNDN);
        mpfr_mul(ang, ang, pi, MPFR_RNDN);

        /* ang -= 2 pi ((n h) mod k) / k */
        unsigned long nh = (nmod * h) % k;
        mpfr_mul_ui(c, pi, 2 * nh, MPFR_RNDN);
        mpfr_div_ui(c, c, k, MPFR_RNDN);
        mpfr_sub(ang, ang, c, MPFR_RNDN);

        mpfr_cos(c, ang, MPFR_RNDN);
        mpfr_add(out, out, c, MPFR_RNDN);
    }

    mpq_clear(s1);
    mpq_clear(s2);
    mpq_clear(sd);
    mpfr_clear(pi);
    mpfr_clear(ang);
    mpfr_clear(c);
}

/* Evaluate the Hagis series for q(n) over odd k = 1,3,...,<= Kmax at the
 * given precision, rounding the real sum to the nearest integer in `out`:
 *     q(n) = (pi/sqrt(24n+1)) sum_{k odd} (1/k) A_k(n) I_1( pi sqrt(48n+2)/(12k) ).
 * Stores in `dist` the distance from the sum to that integer and in
 * `last_mag` the magnitude of the last term added (a tail estimate); the
 * caller uses both to decide whether to trust the rounding. */
static void hrr_q_eval(unsigned long n, unsigned long Kmax, mpfr_prec_t prec,
                       mpz_t out, mpfr_t dist, mpfr_t last_mag) {
    mpfr_t pi, pref, B, arg, Ak, I1, term, sum, rounded, tmp;
    mpfr_init2(pi, prec);
    mpfr_init2(pref, prec);
    mpfr_init2(B, prec);
    mpfr_init2(arg, prec);
    mpfr_init2(Ak, prec);
    mpfr_init2(I1, prec);
    mpfr_init2(term, prec);
    mpfr_init2(sum, prec);
    mpfr_init2(rounded, prec);
    mpfr_init2(tmp, prec);

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

    mpfr_set_zero(sum, 1);
    mpfr_set_zero(last_mag, 1);
    for (unsigned long k = 1; k <= Kmax; k += 2) {
        partitionsq_Ak(Ak, k, n, prec);
        if (mpfr_zero_p(Ak)) continue;

        mpfr_div_ui(arg, B, k, MPFR_RNDN);
        i1_mpfr(I1, arg, prec);

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

    mpfr_clear(pi);
    mpfr_clear(pref);
    mpfr_clear(B);
    mpfr_clear(arg);
    mpfr_clear(Ak);
    mpfr_clear(I1);
    mpfr_clear(term);
    mpfr_clear(sum);
    mpfr_clear(rounded);
    mpfr_clear(tmp);
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
    fprintf(stderr,
            "PartitionsQ::argx: PartitionsQ called with %zu argument%s; "
            "1 argument is expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/* PartitionsQ[n] — number of partitions of n into distinct parts.
 * Reads `res` only; returns NULL (unevaluated) for symbolic, non-integer or
 * big-integer arguments. Negative n gives 0. */
Expr* builtin_partitionsq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return qq_emit_argx(argc);

    Expr* arg = res->data.function.args[0];

    /* A big-integer argument is astronomically large; leave unevaluated. */
    if (!expr_is_integer_like(arg) || arg->type == EXPR_BIGINT) return NULL;

    int64_t n = arg->data.integer;
    if (n < 0) return expr_new_integer(0); /* q(n) = 0 for n < 0 */

    mpz_t out;
    mpz_init(out);

    int rc = -1;
    if (n >= PARTITIONSQ_HRR_THRESHOLD) {
        rc = partitionsq_hrr((unsigned long)n, out);
        if (rc != 0)                       /* fall back if HRR unavailable */
            rc = partitionsq_recurrence((unsigned long)n, out);
    } else {
        rc = partitionsq_recurrence((unsigned long)n, out);
    }
    if (rc != 0) {
        mpz_clear(out);
        return NULL;
    }

    Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(out));
    mpz_clear(out);
    return result;
}

void partitions_init(void) {
    symtab_add_builtin("IntegerPartitions", builtin_integerpartitions);
    symtab_get_def("IntegerPartitions")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("PartitionsP", builtin_partitionsp);
    symtab_get_def("PartitionsP")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);

    symtab_add_builtin("PartitionsQ", builtin_partitionsq);
    symtab_get_def("PartitionsQ")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);
}

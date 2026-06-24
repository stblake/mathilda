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
#include <gmp.h>

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
        e->data.function.head->data.symbol == SYM_Rational &&
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
    return e && e->type == EXPR_SYMBOL && e->data.symbol == s;
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

void partitions_init(void) {
    symtab_add_builtin("IntegerPartitions", builtin_integerpartitions);
    symtab_get_def("IntegerPartitions")->attributes |= ATTR_PROTECTED;
}

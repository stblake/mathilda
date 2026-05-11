#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "expr.h"
#include "arithmetic.h"
#include "sym_intern.h"
#include "sym_names.h"
#include <stdbool.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

// Create/allocate a new integer expression.
Expr* expr_new_integer(int64_t value) {
    Expr* e = malloc(sizeof(Expr));
    if (!e) return NULL;
    e->type = EXPR_INTEGER;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    e->data.integer = value;
    return e;
}

// Create/allocate a new real (double) expression.
Expr* expr_new_real(double value) {
    Expr* e = malloc(sizeof(Expr));
    if (!e) return NULL;
    e->type = EXPR_REAL;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    e->data.real = value;
    return e;
}

// Create/allocate a new symbol expression.
//
// The symbol name is funneled through the global interner so that two
// symbols with the same name share the same `const char*`. This makes:
//   - expr_eq() on symbols a pointer compare,
//   - expr_copy() of a symbol a pointer copy,
//   - expr_free() of a symbol a no-op for the name field.
//
// The cast to `char*` is intentional: `data.symbol` retains its existing
// type for ABI stability, but the memory is owned by the interner and
// must never be freed by Expr code.
Expr* expr_new_symbol(const char* name) {
    Expr* e = malloc(sizeof(Expr));
    if (!e) return NULL;
    e->type = EXPR_SYMBOL;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    e->data.symbol = (char*)intern_symbol(name);
    if (!e->data.symbol) {
        free(e);
        return NULL;
    }
    return e;
}

// Create/allocate a new string expression.
Expr* expr_new_string(const char* str) {
    Expr* e = malloc(sizeof(Expr));
    if (!e) return NULL;

    e->type = EXPR_STRING;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    e->data.string = strdup(str);
    if (!e->data.string) {
        free(e);
        return NULL;
    }
    return e;
}

// Create/allocate an expression: h[arg1, arg2, ...]
Expr* expr_new_function(Expr* head, Expr** args, size_t arg_count) {
    Expr* e = malloc(sizeof(Expr));
    if (!e) return NULL;

    e->type = EXPR_FUNCTION;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    e->data.function.head = head;
    if (arg_count > 0) {
        e->data.function.args = calloc(arg_count, sizeof(Expr*));
        if (!e->data.function.args) {
            free(e);
            return NULL;
        }
        if (args) memcpy(e->data.function.args, args, sizeof(Expr*) * arg_count);
    } else {
        e->data.function.args = NULL;
    }
    e->data.function.arg_count = arg_count;
    return e;
}

// BigInt constructors
Expr* expr_new_bigint_from_mpz(const mpz_t val) {
    Expr* e = malloc(sizeof(Expr));
    if (!e) return NULL;
    e->type = EXPR_BIGINT;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    mpz_init_set(e->data.bigint, val);
    return e;
}

Expr* expr_new_bigint_from_int64(int64_t val) {
    Expr* e = malloc(sizeof(Expr));
    if (!e) return NULL;
    e->type = EXPR_BIGINT;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    mpz_init_set_si(e->data.bigint, val);
    return e;
}

Expr* expr_new_bigint_from_str(const char* str) {
    Expr* e = malloc(sizeof(Expr));
    if (!e) return NULL;
    e->type = EXPR_BIGINT;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    if (mpz_init_set_str(e->data.bigint, str, 10) == -1) {
        mpz_clear(e->data.bigint);
        free(e);
        return NULL;
    }
    return e;
}

#ifdef USE_MPFR
/* MPFR constructors. All allocate an Expr and initialize the payload
 * `mpfr_t` at the requested precision; the caller owns the result and
 * should free it with `expr_free`, which calls `mpfr_clear`. */
Expr* expr_new_mpfr_bits(mpfr_prec_t bits) {
    Expr* e = malloc(sizeof(Expr));
    if (!e) return NULL;
    e->type = EXPR_MPFR;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    mpfr_init2(e->data.mpfr, bits);
    mpfr_set_zero(e->data.mpfr, +1);
    return e;
}
Expr* expr_new_mpfr_from_d(double v, mpfr_prec_t bits) {
    Expr* e = expr_new_mpfr_bits(bits);
    if (e) mpfr_set_d(e->data.mpfr, v, MPFR_RNDN);
    return e;
}
Expr* expr_new_mpfr_from_si(long v, mpfr_prec_t bits) {
    Expr* e = expr_new_mpfr_bits(bits);
    if (e) mpfr_set_si(e->data.mpfr, v, MPFR_RNDN);
    return e;
}
Expr* expr_new_mpfr_from_mpz(const mpz_t z, mpfr_prec_t bits) {
    Expr* e = expr_new_mpfr_bits(bits);
    if (e) mpfr_set_z(e->data.mpfr, z, MPFR_RNDN);
    return e;
}
Expr* expr_new_mpfr_from_str(const char* str, mpfr_prec_t bits) {
    Expr* e = expr_new_mpfr_bits(bits);
    if (!e) return NULL;
    if (mpfr_set_str(e->data.mpfr, str, 10, MPFR_RNDN) != 0) {
        mpfr_clear(e->data.mpfr);
        free(e);
        return NULL;
    }
    return e;
}
Expr* expr_new_mpfr_move(mpfr_t src) {
    Expr* e = malloc(sizeof(Expr));
    if (!e) { mpfr_clear(src); return NULL; }
    e->type = EXPR_MPFR;
    e->refcount = 1;
    e->last_evaluated_at = 0;
    /* mpfr_t is an array-type alias for __mpfr_struct[1]; `memcpy` moves
     * the header — MPFR's documentation (mpfr_swap, GMP's mpz semantics)
     * sanctions this as long as the source is then treated as uninit. */
    memcpy(e->data.mpfr, src, sizeof(e->data.mpfr));
    return e;
}
Expr* expr_new_mpfr_copy(const mpfr_t src) {
    Expr* e = expr_new_mpfr_bits(mpfr_get_prec(src));
    if (e) mpfr_set(e->data.mpfr, src, MPFR_RNDN);
    return e;
}
#endif

void expr_to_mpz(const Expr* e, mpz_t out) {
    if (e->type == EXPR_INTEGER) {
        mpz_init_set_si(out, e->data.integer);
    } else { // EXPR_BIGINT
        mpz_init_set(out, e->data.bigint);
    }
}

bool expr_is_integer_like(const Expr* e) {
    return e && (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT);
}

bool expr_is_numeric_like(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:
        case EXPR_REAL:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        case EXPR_FUNCTION:
            /* Rational[n,d] and Complex[re,im] with numeric components.
             * Use the bigint-aware rational check so rationals whose
             * components have overflowed int64 still fold via the GMP
             * fallback in multiply_numbers / add_numbers (otherwise
             * Times[BigInt, Rational[1, BigInt]] survives unsimplified
             * and breaks exact polynomial division). */
            if (is_rational_like((Expr*)e)) return true;
            {
                Expr *re, *im;
                if (is_complex((Expr*)e, &re, &im)) {
                    return expr_is_numeric_like(re) && expr_is_numeric_like(im);
                }
            }
            return false;
        default:
            return false;
    }
}

Expr* expr_bigint_normalize(Expr* e) {
    if (e->type == EXPR_BIGINT && mpz_fits_slong_p(e->data.bigint)) {
        long val = mpz_get_si(e->data.bigint);
        Expr* result = expr_new_integer((int64_t)val);
        /* Refcount-safe: dec-ref the BigInt; if shared, the original
         * stays alive for other holders (each will normalize on their
         * own); if unique (refcount==1, the common case), expr_free
         * physically clears the mpz_t and frees the node. */
        expr_free(e);
        return result;
    }
    return e;
}

/* Inc-ref. NULL passes through. */
Expr* expr_ref(Expr* e) {
    if (e) e->refcount++;
    return e;
}

/* M3 phase-2 copy-on-write helper. See header doc. */
Expr* expr_unshare(Expr* e) {
    if (!e) return NULL;
    if (e->refcount == 1) return e;

    /* Build a one-level private copy. Children stay shared (inc-ref). */
    Expr* fresh = malloc(sizeof(Expr));
    if (!fresh) return e;  /* OOM: fall back to (still-shared) original */
    fresh->type = e->type;
    fresh->refcount = 1;
    /* Reset the eval timestamp on unshare. The caller is preparing to
     * mutate this node; even if they do not, treating the node as
     * "never evaluated" only costs one extra evaluation, whereas
     * inheriting the timestamp risks a false cache hit on a node that
     * the caller subsequently rewrites. */
    fresh->last_evaluated_at = 0;

    switch (e->type) {
        case EXPR_INTEGER:
            fresh->data.integer = e->data.integer;
            break;
        case EXPR_REAL:
            fresh->data.real = e->data.real;
            break;
        case EXPR_SYMBOL:
            fresh->data.symbol = e->data.symbol;  /* interned */
            break;
        case EXPR_STRING:
            fresh->data.string = e->data.string ? strdup(e->data.string) : NULL;
            break;
        case EXPR_BIGINT:
            mpz_init(fresh->data.bigint);
            mpz_set(fresh->data.bigint, e->data.bigint);
            break;
#ifdef USE_MPFR
        case EXPR_MPFR:
            mpfr_init2(fresh->data.mpfr, mpfr_get_prec(e->data.mpfr));
            mpfr_set(fresh->data.mpfr, e->data.mpfr, MPFR_RNDN);
            break;
#endif
        case EXPR_FUNCTION:
            fresh->data.function.head = expr_copy(e->data.function.head);
            fresh->data.function.arg_count = e->data.function.arg_count;
            if (e->data.function.arg_count > 0) {
                fresh->data.function.args = malloc(sizeof(Expr*) * e->data.function.arg_count);
                if (!fresh->data.function.args) {
                    if (fresh->data.function.head) expr_free(fresh->data.function.head);
                    free(fresh);
                    return e;
                }
                for (size_t i = 0; i < e->data.function.arg_count; i++) {
                    fresh->data.function.args[i] = expr_copy(e->data.function.args[i]);
                }
            } else {
                fresh->data.function.args = NULL;
            }
            break;
        default:
            break;
    }

    /* Drop the caller's ref on the original. */
    expr_free(e);
    return fresh;
}

// Deallocate an expression. Decrements the refcount; only physically
// destroys the node (and its children, for FUNCTION) when refcount
// drops to 0. Atoms are eligible to have refcount > 1 starting in M3
// (see expr_copy). FUNCTION nodes always have refcount == 1 today
// because expr_copy still deep-copies them.
void expr_free(Expr* e) {
    if (!e) return;
    if (e->refcount > 1) {
        e->refcount--;
        return;
    }
    /* refcount == 1 (or, defensively, 0 — should not happen, but treat
     * as last reference). Physically free. */

    switch (e->type) {
        case EXPR_SYMBOL:
            /* data.symbol is interned (owned by sym_intern); never free. */
            break;
        case EXPR_STRING:
            if (e->data.string) free(e->data.string);
            break;
        case EXPR_FUNCTION:
            if (e->data.function.head) expr_free(e->data.function.head);
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (e->data.function.args && e->data.function.args[i]) {
                    expr_free(e->data.function.args[i]);
                }
            }
            if (e->data.function.args) free(e->data.function.args);
            break;
        case EXPR_BIGINT:
            mpz_clear(e->data.bigint);
            break;
#ifdef USE_MPFR
        case EXPR_MPFR:
            mpfr_clear(e->data.mpfr);
            break;
#endif
        default:
            break;
    }
    free(e);
}



// Logically copy an expression. M3 phase-2: every node type, INCLUDING
// FUNCTION, is shared via inc-ref. Mutating helpers must call
// expr_unshare() to obtain a private (refcount==1) version before
// rewriting fields in place. The audit (commit notes) verified that
// the only sites that mutated a possibly-shared FUNCTION node lived
// in print.c; those now wrap their expr_copy() in expr_unshare().
// Other mutators (eval_flatten_args, flatten_sequences, builtin_plus
// numeric-contagion args[i] writes, parse.c arg appending, match.c
// sequence-binding writes, core.c QuotientRemainder zeroing, etc.)
// all act on freshly-constructed FUNCTION nodes whose refcount is
// guaranteed to be 1 at the point of mutation.
Expr* expr_copy(Expr* e) {
    if (!e) return NULL;
    e->refcount++;
    return e;
}

bool expr_eq(const Expr* a, const Expr* b) {
    if (a == b) return true;
    if (!a || !b) return false;
    if (a->type != b->type) {
        // Cross-type equality: EXPR_INTEGER vs EXPR_BIGINT
        if ((a->type == EXPR_INTEGER && b->type == EXPR_BIGINT) ||
            (a->type == EXPR_BIGINT && b->type == EXPR_INTEGER)) {
            mpz_t va, vb;
            expr_to_mpz(a, va);
            expr_to_mpz(b, vb);
            bool result = (mpz_cmp(va, vb) == 0);
            mpz_clear(va);
            mpz_clear(vb);
            return result;
        }
        return false;
    }

    switch (a->type) {
        case EXPR_INTEGER:
            return a->data.integer == b->data.integer;
        case EXPR_REAL:
            if (isnan(a->data.real) && isnan(b->data.real)) return true;
            return a->data.real == b->data.real;
        case EXPR_SYMBOL:
            /* Interned: same name guarantees same pointer. */
            return a->data.symbol == b->data.symbol;
        case EXPR_STRING:
            return strcmp(a->data.string, b->data.string) == 0;
        case EXPR_FUNCTION: {
            if (!expr_eq(a->data.function.head, b->data.function.head)) return false;
            if (a->data.function.arg_count != b->data.function.arg_count) return false;
            for (size_t i = 0; i < a->data.function.arg_count; i++) {
                if (!expr_eq(a->data.function.args[i], b->data.function.args[i])) return false;
            }
            return true;
        }
        case EXPR_BIGINT:
            return mpz_cmp(a->data.bigint, b->data.bigint) == 0;
#ifdef USE_MPFR
        case EXPR_MPFR:
            /* Equal iff same precision AND same value (matches SameQ
             * semantics: 1.`20 and 1.`30 are not SameQ even though
             * their values agree). */
            if (mpfr_get_prec(a->data.mpfr) != mpfr_get_prec(b->data.mpfr)) return false;
            return mpfr_equal_p(a->data.mpfr, b->data.mpfr) != 0;
#endif
    }
    return false;
}

static double get_numeric_value(const Expr* e) {
    if (e->type == EXPR_INTEGER) return (double)e->data.integer;
    if (e->type == EXPR_BIGINT) return mpz_get_d(e->data.bigint);
    if (e->type == EXPR_REAL) return e->data.real;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_get_d(e->data.mpfr, MPFR_RNDN);
#endif
    int64_t n, d;
    if (is_rational((Expr*)e, &n, &d)) return (double)n / d;
    return 0;
}

static int string_compare_canonical(const char* s1, const char* s2) {
    const char* p1 = s1;
    const char* p2 = s2;
    while (*p1 && *p2) {
        char l1 = tolower((unsigned char)*p1);
        char l2 = tolower((unsigned char)*p2);
        if (l1 != l2) return (l1 < l2) ? -1 : 1;
        p1++; p2++;
    }
    if (*p1) return 1;
    if (*p2) return -1;
    
    p1 = s1; p2 = s2;
    while (*p1 && *p2) {
        if (*p1 != *p2) {
            if (islower((unsigned char)*p1) && !islower((unsigned char)*p2)) return -1;
            if (!islower((unsigned char)*p1) && islower((unsigned char)*p2)) return 1;
            return (*p1 < *p2) ? -1 : 1;
        }
        p1++; p2++;
    }
    return 0;
}


/* True for the atomic numeric kinds that have an extractable double value:
 * Integer, BigInt, Real, Rational, MPFR. False for symbolic numeric
 * constants like Sqrt[2-Sqrt[2]] which canonical-sort as numbers but have
 * no meaningful get_numeric_value(). */
static bool is_atomic_numeric(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (is_rational((Expr*)e, NULL, NULL)) return true;
    return false;
}

/* Return true iff `e` contains a Symbol with the interned name `sym`
 * anywhere in its tree (heads and args included). */
static bool expr_contains_symbol(const Expr* e, const char* sym) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol == sym;
    if (e->type != EXPR_FUNCTION) return false;
    if (expr_contains_symbol(e->data.function.head, sym)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (expr_contains_symbol(e->data.function.args[i], sym)) return true;
    }
    return false;
}

/* Polynomial degree of `e` in the variable `v`, as a double:
 *   0.0          -- e does not contain v.
 *   k > 0        -- e is v (k=1), Power[v, k] for positive integer or
 *                   rational k, or a Times product whose direct factors
 *                   sum to k.
 *   +inf         -- e contains v in a non-polynomial position: inside a
 *                   compound that isn't bare-v / Power[v, positive rational]
 *                   (e.g. Sin[v], Plus[..., v, ...], Power[v, negative],
 *                   Power[non-v-containing-v, _]).
 *
 * Only DIRECT factors of a Times contribute to the sum; a non-poly
 * direct factor that itself contains v (e.g. Sin[v]) saturates the
 * result at +inf.  This matches Mathematica's canonical order:
 *   - `Sin[x] + x` displays as `x + Sin[x]` (poly-in-x before
 *      non-poly-in-x).
 *   - `Sqrt[x] + x^2` displays as `Sqrt[x] + x^2` (deg 1/2 < deg 2).
 *   - `1/x + 1` displays as `1 + 1/x` (negative-exp x treated as +inf
 *      so it sorts last).
 *   - `(3 + 6 a + 3 a^2) x` sorts by deg_x = 1 regardless of the Plus
 *      coefficient's internal a-degree. */
static double expr_poly_degree(const Expr* e, const char* v) {
    if (!e) return 0.0;
    if (e->type == EXPR_SYMBOL) return (e->data.symbol == v) ? 1.0 : 0.0;
    if (e->type != EXPR_FUNCTION) return 0.0;
    if (e->data.function.head->type != EXPR_SYMBOL) {
        return expr_contains_symbol(e, v) ? INFINITY : 0.0;
    }
    const char* h = e->data.function.head->data.symbol;
    if (h == SYM_Power && e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp  = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL && base->data.symbol == v) {
            if (exp->type == EXPR_INTEGER) {
                int64_t k = exp->data.integer;
                if (k < 0) return INFINITY;
                return (double)k;
            }
            int64_t num, den;
            if (is_rational((Expr*)exp, &num, &den) && den > 0) {
                if (num < 0) return INFINITY;
                return (double)num / (double)den;
            }
            return INFINITY;
        }
        return (expr_contains_symbol(base, v) || expr_contains_symbol(exp, v))
               ? INFINITY : 0.0;
    }
    if (h == SYM_Times) {
        double sum = 0.0;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            double d = expr_poly_degree(e->data.function.args[i], v);
            if (isinf(d)) return INFINITY;
            sum += d;
        }
        return sum;
    }
    return expr_contains_symbol(e, v) ? INFINITY : 0.0;
}

/* Append every distinct Symbol-name pointer appearing in `e` to `out`
 * (up to capacity `cap`), via a depth-first walk that visits heads and
 * args.  Returns the new logical length of `out`.  Truncates silently
 * if the buffer fills -- pathological inputs with hundreds of distinct
 * symbols lose precision in the canonical order but won't crash. */
static size_t collect_symbols_in(const Expr* e, const char** out,
                                 size_t count, size_t cap) {
    if (!e || count >= cap) return count;
    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < count; i++) {
            if (out[i] == e->data.symbol) return count;
        }
        out[count++] = e->data.symbol;
        return count;
    }
    if (e->type != EXPR_FUNCTION) return count;
    count = collect_symbols_in(e->data.function.head, out, count, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        count = collect_symbols_in(e->data.function.args[i], out, count, cap);
    }
    return count;
}

/* qsort comparator that sorts char* pointers in REVERSE lexicographic
 * order, so the lex-LAST variable name ends up first in the sorted
 * array.  Used to drive the degree-vector comparison: in Mathematica's
 * canonical order, the lex-last variable (e.g. `y` in {a, b, x, y}) is
 * the most significant dimension. */
static int symbol_reverse_cmp(const void* pa, const void* pb) {
    const char* a = *(const char* const*)pa;
    const char* b = *(const char* const*)pb;
    return strcmp(b, a);
}

/* Canonical comparison of two expressions, modeled on Mathematica's
 * `Order`.  Public entry point used by Sort, OrderedQ, and the
 * Orderless-attribute argument-sorting pipeline in eval.c.
 *
 * Order of precedence:
 *   1. Numeric atoms (Integer, Real, Rational, BigInt, MPFR) come first
 *      and compare by value.  A numeric vs non-numeric pair places the
 *      numeric first.
 *   2. Strings come next, compared lexicographically.
 *   3. Everything else (Symbols and Functions) is compared by polynomial
 *      degree vector:
 *        - Collect every Symbol name appearing in either expression.
 *        - Sort those names in REVERSE alphabetical order.
 *        - For each name v, compute the polynomial degree of each side
 *          in v (see expr_poly_degree).
 *        - Compare the degree vectors lexicographically ASCENDING.
 *      This is what makes `Plus[a x^2, b x]` sort as `b x + a x^2`
 *      (deg_x = 1 < deg_x = 2) and what gives the multivariate
 *      grevlex-like display order on `Expand[(1 + x y^2 + y x^3)^3]`.
 *   4. When the degree vectors are equal, fall back to a structural
 *      tiebreak: bare Symbol before compound, then by head, arity, and
 *      args (the args themselves are compared with this same canonical
 *      comparator, so the tiebreak recursively re-enters polynomial
 *      ordering for nested expressions). */
int expr_compare(const Expr* a, const Expr* b) {
    if (a == b) return 0;
    if (!a) return -1;
    if (!b) return 1;

    /* 1. Numeric atoms. */
    bool a_atomic = is_atomic_numeric(a);
    bool b_atomic = is_atomic_numeric(b);
    if (a_atomic && b_atomic) {
        if (a->type == EXPR_INTEGER && b->type == EXPR_INTEGER) {
            if (a->data.integer < b->data.integer) return -1;
            if (a->data.integer > b->data.integer) return 1;
            return 0;
        }
        if (expr_is_integer_like(a) && expr_is_integer_like(b)) {
            mpz_t ma2, mb2;
            expr_to_mpz(a, ma2);
            expr_to_mpz(b, mb2);
            int cmp = mpz_cmp(ma2, mb2);
            mpz_clear(ma2);
            mpz_clear(mb2);
            return cmp;
        }
        double va = get_numeric_value(a);
        double vb = get_numeric_value(b);
        if (va < vb) return -1;
        if (va > vb) return 1;
        if (a->type != b->type) return (int)a->type - (int)b->type;
        return 0;
    }
    if (a_atomic) return -1;
    if (b_atomic) return 1;

    /* 2. Strings. */
    if (a->type == EXPR_STRING && b->type == EXPR_STRING) {
        return string_compare_canonical(a->data.string, b->data.string);
    }
    if (a->type == EXPR_STRING) return -1;
    if (b->type == EXPR_STRING) return 1;

    /* 3. Polynomial degree vector comparison.  VAR_CAP=64 covers any
     * realistic input; pathological inputs degrade gracefully (the
     * comparator stays total because both sides see the same truncated
     * variable set). */
    enum { VAR_CAP = 64 };
    const char* vars[VAR_CAP];
    size_t count = 0;
    count = collect_symbols_in(a, vars, count, VAR_CAP);
    count = collect_symbols_in(b, vars, count, VAR_CAP);

    qsort(vars, count, sizeof(const char*), symbol_reverse_cmp);

    for (size_t i = 0; i < count; i++) {
        double da = expr_poly_degree(a, vars[i]);
        double db = expr_poly_degree(b, vars[i]);
        if (da < db) return -1;
        if (da > db) return 1;
    }

    /* 4. Structural tiebreak.  At this point both expressions share an
     * identical polynomial signature; differences live in coefficients
     * (atomic numerics) or non-poly factors (Plus/Sin/etc. that contain
     * the same variables in the same non-poly way). */
    if (a->type == EXPR_SYMBOL && b->type == EXPR_SYMBOL) {
        if (a->data.symbol == b->data.symbol) return 0;
        return strcmp(a->data.symbol, b->data.symbol);
    }
    if (a->type == EXPR_SYMBOL) return -1;
    if (b->type == EXPR_SYMBOL) return 1;

    if (a->type == EXPR_FUNCTION && b->type == EXPR_FUNCTION) {
        int head_cmp = expr_compare(a->data.function.head, b->data.function.head);
        if (head_cmp != 0) return head_cmp;
        size_t ac = a->data.function.arg_count;
        size_t bc = b->data.function.arg_count;
        if (ac != bc) return (ac < bc) ? -1 : 1;
        for (size_t i = 0; i < ac; i++) {
            int c = expr_compare(a->data.function.args[i], b->data.function.args[i]);
            if (c != 0) return c;
        }
        return 0;
    }

    return 0;
}

uint64_t expr_hash(const Expr* e) {
    if (!e) return 0;
    uint64_t h = 14695981039346656037ULL;
    const uint64_t prime = 1099511628211ULL;

    h ^= (uint64_t)e->type;
    h *= prime;

    switch (e->type) {
        case EXPR_INTEGER: {
            uint64_t v;
            memcpy(&v, &e->data.integer, 8);
            h ^= v; h *= prime;
            break;
        }
        case EXPR_REAL: {
            uint64_t v;
            memcpy(&v, &e->data.real, 8);
            h ^= v; h *= prime;
            break;
        }
        case EXPR_SYMBOL:
        case EXPR_STRING: {
            const char* s = e->data.symbol;
            if (s) {
                while (*s) {
                    h ^= (uint8_t)(*s++);
                    h *= prime;
                }
            }
            break;
        }
        case EXPR_FUNCTION: {
            h ^= expr_hash(e->data.function.head);
            h *= prime;
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                h ^= expr_hash(e->data.function.args[i]);
                h *= prime;
            }
            break;
        }
        case EXPR_BIGINT: {
            h ^= (uint64_t)(mpz_sgn(e->data.bigint) + 2);
            h *= prime;
            size_t nlimbs = mpz_size(e->data.bigint);
            for (size_t i = 0; i < nlimbs; i++) {
                h ^= mpz_getlimbn(e->data.bigint, i);
                h *= prime;
            }
            break;
        }
#ifdef USE_MPFR
        case EXPR_MPFR: {
            /* Hash precision + IEEE double approximation. Not perfectly
             * collision-free across precisions with equal double value,
             * but it's a hash, not an equality. */
            h ^= (uint64_t)mpfr_get_prec(e->data.mpfr);
            h *= prime;
            double d = mpfr_get_d(e->data.mpfr, MPFR_RNDN);
            uint64_t v; memcpy(&v, &d, 8);
            h ^= v; h *= prime;
            break;
        }
#endif
    }
    return h;
}


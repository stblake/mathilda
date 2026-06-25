
// Mathilda

#ifndef EXPR_H
#define EXPR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdarg.h>
#include <math.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

typedef enum {
    EXPR_INTEGER,
    EXPR_REAL,
    EXPR_SYMBOL,
    EXPR_STRING,
    EXPR_FUNCTION,
    EXPR_BIGINT
#ifdef USE_MPFR
    , EXPR_MPFR                /* arbitrary-precision real (MPFR) */
#endif
} ExprType;

typedef struct Expr {
    ExprType type;
    /* Reference count for shared-ownership / copy-on-write semantics
     * (M3 milestone). Every successful expr_new_* / expr_copy returns a
     * node with refcount == 1; expr_ref bumps it; expr_free decrements
     * and only physically frees on transition to 0. Atoms are eligible
     * to be shared today (see expr_copy); compound (FUNCTION) nodes
     * remain deep-copied for now and therefore always have refcount==1
     * in current code paths. */
    unsigned refcount;
    /* M3 phase-3 evaluation timestamp. The value of `eval_clock_get()` at
     * the moment this node was last evaluated to a fixed point. Fresh
     * (parser- or builtin-constructed) nodes carry 0, which never matches
     * the live clock (clock starts at 1) so they always evaluate on the
     * first call. Once `evaluate(e)` reaches fixed point it stamps the
     * result with the current clock; a subsequent `evaluate(e)` with the
     * same clock returns an inc-ref'd view immediately, skipping the
     * outer fixed-point loop and all of `evaluate_step`. The clock is
     * bumped by every symbol-table mutation (Set, SetDelayed, Clear,
     * SetAttributes, ClearAttributes), so any user definition change
     * cleanly invalidates every cached evaluation in one shot. The
     * field is benign metadata and is intentionally NOT considered by
     * `expr_eq` / `expr_hash`; updating it on a shared node is safe. */
    uint64_t last_evaluated_at;
    union {
        int64_t integer;
        double real;
        char* symbol;
        char* string;
        struct {
            struct Expr* head;
            struct Expr** args;
            size_t arg_count;
        } function;
        mpz_t bigint;
#ifdef USE_MPFR
        mpfr_t mpfr;          /* carries its own precision in bits */
#endif
    } data;
} Expr;

Expr* expr_new_integer(int64_t value);
Expr* expr_new_real(double value);
Expr* expr_new_symbol(const char* name);
Expr* expr_new_string(const char* str);
Expr* expr_new_function(Expr* head, Expr** args, size_t arg_count);
void expr_free(Expr* e);
Expr* expr_copy(Expr* e);
/* Inc-ref `e` (no copy) and return the same pointer. Pairs with
 * expr_free, which dec-refs. NULL is passed through. Use this in places
 * that previously called expr_copy on an immutable atom. */
Expr* expr_ref(Expr* e);

/* M3 phase-2 copy-on-write helper. Returns a refcount==1 Expr that is
 * logically equal to `e`, consuming the input ref. If `e->refcount` is
 * already 1, returns `e` unchanged (zero work). Otherwise builds a
 * one-level private copy: a fresh top-level node with a private args[]
 * array (FUNCTION) or owned mpz_t/mpfr_t/string payload (others), with
 * children inc-ref'd via expr_copy. The caller may then mutate the
 * returned node's *direct* fields (args[i], arg_count, head, etc.).
 * To mutate deeper than one level, unshare each interior node along
 * the path. NULL passes through. */
Expr* expr_unshare(Expr* e);

#include <stdbool.h>
bool expr_eq(const Expr* a, const Expr* b);
int expr_compare(const Expr* a, const Expr* b);
uint64_t expr_hash(const Expr* e);

/* BigInt constructors */
Expr* expr_new_bigint_from_mpz(const mpz_t val);
Expr* expr_new_bigint_from_int64(int64_t val);
Expr* expr_new_bigint_from_str(const char* str);

#ifdef USE_MPFR
/* MPFR constructors. Each allocates an Expr whose payload `mpfr_t` is
 * initialized to the requested precision. Caller must `expr_free` when
 * done, which calls `mpfr_clear`. */
Expr* expr_new_mpfr_bits(mpfr_prec_t bits);                       /* zero */
Expr* expr_new_mpfr_from_d(double v, mpfr_prec_t bits);
Expr* expr_new_mpfr_from_si(long v, mpfr_prec_t bits);
Expr* expr_new_mpfr_from_mpz(const mpz_t z, mpfr_prec_t bits);
Expr* expr_new_mpfr_from_str(const char* str, mpfr_prec_t bits);
/* Build an Expr taking ownership of `src`. The mpfr_t is moved, not copied;
 * afterwards the caller must not touch `src`. Precision is inherited. */
Expr* expr_new_mpfr_move(mpfr_t src);
/* Copy constructor: new Expr with an independent mpfr_t at the same
 * precision as `src`. */
Expr* expr_new_mpfr_copy(const mpfr_t src);
#endif

/* Portable strdup: `strdup` is POSIX, not C99, and glibc hides it under
 * -std=c99.  Returns a malloc'd copy of `s` (NULL if `s` is NULL). Caller
 * frees with free(). Use this instead of strdup throughout the codebase. */
char* mathilda_strdup(const char* s);

/* Portable "this static function/variable may be unused" marker. Wraps the
 * GNU/Clang attribute behind a guard so strict C99 elsewhere still builds
 * (SPEC §10 forbids an unguarded __attribute__). */
#if defined(__GNUC__) || defined(__clang__)
#define MATHILDA_MAYBE_UNUSED __attribute__((unused))
#else
#define MATHILDA_MAYBE_UNUSED
#endif

/* Helpers used by arithmetic modules */
void  expr_to_mpz(const Expr* e, mpz_t out);
bool  expr_is_integer_like(const Expr* e);
/* True if `e` represents any concrete number: Integer, BigInt, Real,
 * Rational[n,d], Complex with numeric parts, or (with USE_MPFR) MPFR. */
bool  expr_is_numeric_like(const Expr* e);
Expr* expr_bigint_normalize(Expr* e);

#endif // EXPR_H

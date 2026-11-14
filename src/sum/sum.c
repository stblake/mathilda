/*
 * sum.c -- Sum dispatcher for Mathilda.
 *
 * Sum is HoldAll: the summation variable and bounds must be held so that the
 * iterator is not prematurely evaluated against an outer binding (exactly as
 * Table/Do hold their iterator specs).
 *
 * Responsibilities of this file (Stage 0):
 *   - strip trailing options (Method -> "...", etc.);
 *   - rewrite multiple iterators Sum[f, s1, ..., sk] into nested single-spec
 *     sums (outer-depends-on-inner bounds come for free);
 *   - finite explicit expansion: when a range resolves to a finite span of
 *     integers, or the spec iterates an explicit list, bind the variable and
 *     fold the evaluated terms with Plus;
 *   - otherwise (symbolic bounds, Infinity, or the indefinite form Sum[f,i])
 *     run a Method cascade over the context-qualified sub-algorithms
 *     Sum`Polynomial, Sum`Geometric, Sum`Gosper.  Each sub-builtin returns the
 *     closed form (definite: Sum`M[f,i,imin,imax]; indefinite: Sum`M[f,i]) or
 *     comes back unevaluated to signal "fall through".  When all stages fall
 *     through the Sum[...] is returned unevaluated (held).
 *
 * Adding a later stage is purely additive: a new src/sum/sum_*.c file, one
 * try_* line in the cascade, and one *_init() call in sum_init().
 *
 * Memory contract: builtin_sum takes ownership of res but must not free it
 * (the evaluator owns it).  Every Expr* allocated here is freed on all paths.
 */

#include "sum.h"
#include "term/term_progress.h"   /* the REPL spinner's loop progress */
#include "checked_int.h"   /* ci_add_i64: integer index step */
#include "sum_internal.h"   /* sum_body_has_opaque_index — shared with sum_gosper.c */
#include "ndarray.h"   /* is_ndarray — see sum_body_is_array */
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "iter.h"
#include "attr.h"
#include "arithmetic.h"   /* is_rational, make_complex */
#include "sym_names.h"
#include "message.h"   /* mth_message: Sum::div through the Quiet/Check funnel */
#include "compile/autocompile.h"   /* machine fast path for finite numeric sums */
#include "rationalize.h"   /* internal_contains_inexact */
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>   /* creal/cimag on the real-iterator accumulator */

/* Upper bound on the number of terms a finite explicit expansion will add
 * before giving up (a runaway guard, not a correctness limit). */
#define SUM_MAX_FINITE_TERMS 100000000LL
/* Finite ranges up to this many terms are enumerated before any closed form is
 * tried -- the reference's own switch point (it enumerates Sum[x^i, {i, 0, 10^6}]
 * term by term and tries symbolic summation only beyond). */
#define SUM_ENUMERATE_MAX 1000000.0

/* Ceiling for the compiled machine fast path (sum_try_compiled).  Far above
 * SUM_MAX_FINITE_TERMS because the compiled accumulator allocates nothing per
 * term — the enumeration cap is a memory bound (an Expr* array), which the
 * compiled loop does not pay.  This is purely a runaway/typo guard so a stray
 * {i, 1, 10^18} cannot spin; it is NOT a correctness limit and may be tuned. */
#define SUM_MAX_COMPILED_TERMS 1000000000000LL

/* ------------------------------------------------------------------ */
/*  Small helpers                                                      */
/* ------------------------------------------------------------------ */

/* True if e is a function with the given symbol head. */
static bool head_is_name(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol.name, name) == 0;
}

/* A result "came back unevaluated" if its head is still the stage we called
 * (no rule fired). */
static bool result_is_unresolved(const Expr* result, const char* head_name) {
    return head_is_name(result, head_name);
}

/* True if the step is the integer 1.  The closed-form stages
 * (Sum`Polynomial / Sum`Geometric / Sum`Gosper) take no step argument and
 * assume a unit step, so they may only be consulted when di == 1; a non-unit
 * step (e.g. {i, 1, n, 2}) must not be dispatched to them or it would silently
 * produce the wrong, step-1 closed form. */
static bool is_unit_step(const Expr* di) {
    return di && di->type == EXPR_INTEGER && di->data.integer == 1;
}

/* ------------------------------------------------------------------ */
/*  Method option                                                      */
/* ------------------------------------------------------------------ */

typedef enum {
    SUM_METHOD_AUTOMATIC = 0,
    SUM_METHOD_POLYNOMIAL,
    SUM_METHOD_GEOMETRIC,
    SUM_METHOD_GOSPER,
    SUM_METHOD_RATIONAL,
    SUM_METHOD_EULER,
    SUM_METHOD_ALTERNATING,
    SUM_METHOD_TRIGONOMETRIC,
    SUM_METHOD_INVALID
} SumMethod;

/* Is arg an option of the form sym -> rhs (Rule or RuleDelayed)? */
static bool is_option(const Expr* arg) {
    if (arg->type != EXPR_FUNCTION) return false;
    if (arg->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = arg->data.function.head->data.symbol.name;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (arg->data.function.arg_count != 2) return false;
    return arg->data.function.args[0]->type == EXPR_SYMBOL;
}

/* Map a Method -> "Name" option to a SumMethod.  Returns SUM_METHOD_AUTOMATIC
 * for any non-Method option (those are ignored in this pass). */
static SumMethod parse_method_option(const Expr* opt) {
    Expr* lhs = opt->data.function.args[0];
    Expr* rhs = opt->data.function.args[1];
    if (lhs->data.symbol.name != SYM_Method) return SUM_METHOD_AUTOMATIC;
    if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic)
        return SUM_METHOD_AUTOMATIC;
    if (rhs->type != EXPR_STRING) return SUM_METHOD_INVALID;
    const char* m = rhs->data.string;
    if (strcmp(m, "Automatic") == 0)                 return SUM_METHOD_AUTOMATIC;
    if (strcmp(m, "Polynomial") == 0)                return SUM_METHOD_POLYNOMIAL;
    if (strcmp(m, "PolynomialExponential") == 0)     return SUM_METHOD_GEOMETRIC;
    if (strcmp(m, "Geometric") == 0)                 return SUM_METHOD_GEOMETRIC;
    if (strcmp(m, "Gosper") == 0)                    return SUM_METHOD_GOSPER;
    if (strcmp(m, "HypergeometricTermGosper") == 0)  return SUM_METHOD_GOSPER;
    if (strcmp(m, "Rational") == 0)                  return SUM_METHOD_RATIONAL;
    if (strcmp(m, "Euler") == 0)                     return SUM_METHOD_EULER;
    if (strcmp(m, "Alternating") == 0)               return SUM_METHOD_ALTERNATING;
    if (strcmp(m, "Trigonometric") == 0)             return SUM_METHOD_TRIGONOMETRIC;
    return SUM_METHOD_INVALID;
}

/* ------------------------------------------------------------------ */
/*  Stage calls                                                        */
/* ------------------------------------------------------------------ */

/* Build head[f, i, imin, imax] (definite) and evaluate it. */
static Expr* call_stage_def(const char* head, Expr* f, Expr* var,
                            Expr* imin, Expr* imax) {
    Expr* args[4] = { expr_copy(f), expr_copy(var),
                      expr_copy(imin), expr_copy(imax) };
    Expr* call = expr_new_function(expr_new_symbol(head), args, 4);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Build head[f, i] (indefinite) and evaluate it. */
static Expr* call_stage_indef(const char* head, Expr* f, Expr* var) {
    Expr* args[2] = { expr_copy(f), expr_copy(var) };
    Expr* call = expr_new_function(expr_new_symbol(head), args, 2);
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* Try one definite stage; NULL means "fall through". */
static Expr* try_def(const char* head, Expr* f, Expr* var,
                     Expr* imin, Expr* imax) {
    Expr* r = call_stage_def(head, f, var, imin, imax);
    if (!r) return NULL;
    if (result_is_unresolved(r, head)) { expr_free(r); return NULL; }
    return r;
}

/* Try one indefinite stage; NULL means "fall through". */
static Expr* try_indef(const char* head, Expr* f, Expr* var) {
    Expr* r = call_stage_indef(head, f, var);
    if (!r) return NULL;
    if (result_is_unresolved(r, head)) { expr_free(r); return NULL; }
    return r;
}

/* Run the Automatic cascade (or a single strict stage) for a definite sum. */
static Expr* dispatch_def(SumMethod method, Expr* f, Expr* var,
                          Expr* imin, Expr* imax) {
    switch (method) {
        case SUM_METHOD_AUTOMATIC: {
            Expr* r = try_def("Sum`Polynomial", f, var, imin, imax);
            if (!r) r = try_def("Sum`Geometric", f, var, imin, imax);
            if (!r) r = try_def("Sum`Gosper", f, var, imin, imax);
            if (!r) r = try_def("Sum`Alternating", f, var, imin, imax);
            if (!r) r = try_def("Sum`Trigonometric", f, var, imin, imax);
            if (!r) r = try_def("Sum`Dirichlet", f, var, imin, imax);
            if (!r) r = try_def("Sum`ZetaSeries", f, var, imin, imax);
            if (!r) r = try_def("Sum`Rational", f, var, imin, imax);
            if (!r) r = try_def("Sum`LogZeta", f, var, imin, imax);
            if (!r) r = try_def("Sum`LogRational", f, var, imin, imax);
            if (!r) r = try_def("Sum`Euler", f, var, imin, imax);
            if (!r) r = try_def("Sum`EulerNonlinear", f, var, imin, imax);
            if (!r) r = try_def("Sum`Hypergeometric", f, var, imin, imax);
            if (!r) r = try_def("Sum`Lambert", f, var, imin, imax);
            return r;
        }
        case SUM_METHOD_POLYNOMIAL:   return try_def("Sum`Polynomial", f, var, imin, imax);
        case SUM_METHOD_GEOMETRIC:    return try_def("Sum`Geometric",  f, var, imin, imax);
        case SUM_METHOD_GOSPER:       return try_def("Sum`Gosper",     f, var, imin, imax);
        case SUM_METHOD_RATIONAL:     return try_def("Sum`Rational",   f, var, imin, imax);
        case SUM_METHOD_EULER:        return try_def("Sum`Euler",      f, var, imin, imax);
        case SUM_METHOD_ALTERNATING:  return try_def("Sum`Alternating", f, var, imin, imax);
        case SUM_METHOD_TRIGONOMETRIC: return try_def("Sum`Trigonometric", f, var, imin, imax);
        default:                    return NULL;
    }
}

/* Run the Automatic cascade (or a single strict stage) for an indefinite sum. */
static Expr* dispatch_indef(SumMethod method, Expr* f, Expr* var) {
    switch (method) {
        case SUM_METHOD_AUTOMATIC: {
            Expr* r = try_indef("Sum`Polynomial", f, var);
            if (!r) r = try_indef("Sum`Geometric", f, var);
            if (!r) r = try_indef("Sum`Gosper", f, var);
            return r;
        }
        case SUM_METHOD_POLYNOMIAL: return try_indef("Sum`Polynomial", f, var);
        case SUM_METHOD_GEOMETRIC:  return try_indef("Sum`Geometric",  f, var);
        case SUM_METHOD_GOSPER:     return try_indef("Sum`Gosper",     f, var);
        case SUM_METHOD_RATIONAL:   return NULL;   /* infinite-only stage */
        default:                    return NULL;
    }
}

/* ------------------------------------------------------------------ */
/*  Finite explicit expansion                                          */
/* ------------------------------------------------------------------ */

/* Build Plus[terms...] from a count of owned term Expr*, evaluate, and free
 * the array.  An empty sum is 0. */
static Expr* fold_plus(Expr** terms, size_t count) {
    if (count == 0) { free(terms); return expr_new_integer(0); }
    Expr* plus = expr_new_function(expr_new_symbol(SYM_Plus), terms, count);
    free(terms);
    Expr* r = evaluate(plus);
    expr_free(plus);
    return r;
}

/* Iterate an explicit list spec, binding var to each element and collecting
 * the evaluated body.  Returns the folded Plus.  Never fails. */
static Expr* expand_list(Expr* f, Expr* var, Expr* list) {
    size_t n = list->data.function.arg_count;
    Expr** terms = malloc(sizeof(Expr*) * (n ? n : 1));
    size_t count = 0;

    IterShadow saved = iter_spec_shadow(var);
    int ptok = mth_progress_push(SYM_Sum, (long)n);
    for (size_t i = 0; i < n; i++) {
        mth_progress_step(ptok, (long)i);
        iter_bind(var, list->data.function.args[i]);
        terms[count++] = evaluate(f);
    }
    mth_progress_pop(ptok);
    iter_spec_restore(var, saved);

    return fold_plus(terms, count);
}

/*
 * Iterate a numeric arithmetic-progression range [imin : imax : di], binding
 * var to the exact value each step and collecting the evaluated body.  Returns
 * the folded Plus, or NULL if the span exceeds SUM_MAX_FINITE_TERMS.
 *
 * Exact value stepping mirrors Do: curr_e is advanced via Plus[curr_e, di]
 * while a synchronised double drives the termination test.
 */
static Expr* expand_range_first(Expr* f, Expr* var, Expr* imin, Expr* imax, Expr* di,
                                double min_val, double max_val, double di_val,
                                bool is_real, Expr* first);
static Expr* expand_range(Expr* f, Expr* var, Expr* imin, Expr* imax, Expr* di,
                          double min_val, double max_val, double di_val,
                          bool is_real) {
    return expand_range_first(f, var, imin, imax, di, min_val, max_val, di_val,
                              is_real, NULL);
}

/* expand_range with the body's value at the first index already computed (the
 * compiled path's inexactness probe, owned here; NULL when none): each term is
 * evaluated exactly once, as in the reference, so a summand with side effects
 * runs once per index and a costly one is not paid for twice. */
static Expr* expand_range_first(Expr* f, Expr* var, Expr* imin, Expr* imax, Expr* di,
                                double min_val, double max_val, double di_val,
                                bool is_real, Expr* first) {
    size_t cap = 16, count = 0;
    Expr** terms = malloc(sizeof(Expr*) * cap);

    IterShadow saved = iter_spec_shadow(var);
    double val = min_val;
    Expr* curr_e = expr_copy(imin);
    bool overflow = false;

    int ptok = mth_progress_push(SYM_Sum, mth_progress_range_count(min_val, max_val, di_val));
    while (iter_range_continue(is_real, /*is_inf=*/false, curr_e, imax,
                               val, max_val, di_val)) {
        mth_progress_step(ptok, (long)count);
        if ((int64_t)count >= SUM_MAX_FINITE_TERMS) { overflow = true; break; }

        if (count == cap) { cap *= 2; terms = realloc(terms, sizeof(Expr*) * cap); }
        if (first) {
            terms[count++] = first;
            first = NULL;
        } else {
            Expr* i_val = is_real ? expr_new_real(val) : expr_copy(curr_e);
            iter_bind(var, i_val);
            terms[count++] = evaluate(f);
            expr_free(i_val);
        }

        /* advance curr_e and val: machine integers add directly (the common
         * case, once per term), anything else through Plus */
        Expr* next_e = NULL;
        int64_t nx;
        if (curr_e->type == EXPR_INTEGER && di->type == EXPR_INTEGER &&
            !ci_add_i64(curr_e->data.integer, di->data.integer, &nx)) {
            next_e = expr_new_integer(nx);
        } else {
            Expr* nargs[2] = { expr_copy(curr_e), expr_copy(di) };
            Expr* nexpr = expr_new_function(expr_new_symbol(SYM_Plus), nargs, 2);
            next_e = evaluate(nexpr);
            expr_free(nexpr);
        }
        expr_free(curr_e);
        curr_e = next_e;
        /* Refresh `val` from the exact value -- a symbolic one (Pi + 1) through
         * its machine value -- or, failing that, by the double recurrence, so
         * the loop always advances. */
        if (is_real || !iter_real_value(curr_e, &val)) val += di_val;
    }
    mth_progress_pop(ptok);
    expr_free(curr_e);
    expr_free(first);   /* unused only on an empty range */
    iter_spec_restore(var, saved);

    if (overflow) {
        for (size_t i = 0; i < count; i++) expr_free(terms[i]);
        free(terms);
        return NULL;
    }
    return fold_plus(terms, count);
}

/* ------------------------------------------------------------------ */
/*  Single-spec handler                                                */
/* ------------------------------------------------------------------ */

/* Handle Sum[f, spec] for one already-isolated spec, with the parsed method.
 * Returns the summed value, or NULL to leave the Sum[...] unevaluated. */
/* Above this many terms, keep trying the closed form even for an array-valued
 * body: the expansion is O(terms x elements) and the closed form, when it works,
 * is O(1) in the number of terms. Below it the expansion is certainly cheap, and
 * the attempt is certainly not. */
#define SUM_ARRAY_EXPAND_MAX 256

/* WHY BODY SIZE DECIDES THE PATH. The closed-form stages (polynomial /
 * geometric / hypergeometric antidifference, and Rational -> Cancel /
 * PolynomialGCD / Expand / FLINT) cost time proportional to the SIZE OF THE
 * BODY, not to the number of terms, because they take the body apart
 * symbolically. On a scalar summand that is nothing. On a large one it is
 * ruinous, and it FAILS on a body that has no closed form, then falls through to
 * the enumeration that should have run first:
 *   - A 256x256 packed grid measured **214 ms per Sum**, versus 0.256 ms to
 *     expand: that fixed cost was the whole of the vectorised Game of Life
 *     benchmark, Sum[RotateLeft[m,{i,j}], {i,-1,1}, {j,-1,1}], 463 ms vs 0.81 ms
 *     for the same nine terms as an explicit Plus (and unpacked as well).
 *   - A 20-D rotated Rastrigin objective, Sum[z[[i]]^2 - 10 Cos[2 Pi z[[i]]],
 *     {i, 1, 20}] with z = q.vars, took MINUTES to build: each z[[i]] carries
 *     the whole dotted vector, so every stage churned over it and failed.
 *     Enumerating the 20 terms is instant.
 * Only SHORT ranges skip the closed form, where the enumeration is bounded and
 * cheap no matter how big the body is. A long range keeps trying the closed
 * form -- Sum[m, {i, 1, 10^5}] for a constant matrix m answers in 0.53 ms that
 * way and would cost 10^5 buffer adds by expansion.
 *
 * SUM_BODY_CLOSED_FORM_MAX is the node budget below which a body is cheap to
 * decompose. Ordinary scalar summands (i^2, r^i, i x[i], small polynomials --
 * tens of nodes) stay under it and on the closed-form-first path unchanged. */
#define SUM_BODY_CLOSED_FORM_MAX 48

/* Node count of the body AS THE CLOSED-FORM STAGES WILL SEE IT, capped at
 * `limit`. Two twists over a plain node walk, both essential because Sum is
 * HoldAll and the body arrives unevaluated:
 *   - a SYMBOL is charged the size of its OwnValue (the tree the cascade gets
 *     once it evaluates the body). A summand like Part[z, i] is three held
 *     nodes, but z is a length-n list, so the cascade sees Part[<the whole
 *     list>, i]; charging z's value is what exposes that. (+1 per symbol before
 *     recursing keeps a cyclic OwnValue like a := a + 1 from looping -- acc rises
 *     every hop and the cap stops it.)
 *   - an NDArray counts as the whole budget: cheap as a single node, but
 *     expensive to take apart, and spliced in directly by helpers that bind
 *     their grid argument by substitution (the Game of Life case).
 * Capped, so this is O(min(size, limit)): a handful of nodes for an ordinary
 * summand -- and, unlike the old array probe, it never EVALUATES the body, so a
 * purely arithmetic Sum pays only a short pointer walk, not a spare evaluation. */
static int64_t sum_body_expanded_size(const Expr* e, int64_t limit, int64_t acc) {
    if (!e || acc > limit) return acc;
    if (is_ndarray(e)) return limit + 1;
    if (e->type == EXPR_SYMBOL) {
        SymbolDef* d = symtab_lookup(e->data.symbol.name);
        if (d && d->own_values && d->own_values->replacement)
            return sum_body_expanded_size(d->own_values->replacement, limit, acc + 1);
        return acc + 1;
    }
    if (e->type != EXPR_FUNCTION) return acc + 1;
    acc += 1;
    acc = sum_body_expanded_size(e->data.function.head, limit, acc);
    for (size_t i = 0; i < e->data.function.arg_count && acc <= limit; i++)
        acc = sum_body_expanded_size(e->data.function.args[i], limit, acc);
    return acc;
}

static bool sum_body_is_expensive(const Expr* f) {
    return sum_body_expanded_size(f, SUM_BODY_CLOSED_FORM_MAX, 0)
               > SUM_BODY_CLOSED_FORM_MAX;
}

/* Integer-membership predicates: True/False for an INTEGER argument, but a
 * definite (index-oblivious) value for a SYMBOLIC one -- EvenQ[k]/PrimeQ[k]/
 * IntegerQ[k] fold to False, SquareFreeQ[k] to True. */
static bool sum_head_is_index_predicate(const char* h) {
    return h == SYM_EvenQ || h == SYM_OddQ
        || h == SYM_PrimeQ || strcmp(h, "CompositeQ") == 0
        || h == SYM_PrimePowerQ || h == SYM_IntegerQ
        || h == SYM_SquareFreeQ || h == SYM_CoprimeQ
        || h == SYM_Divisible;
}

static bool expr_mentions_symbol(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return strcmp(e->data.symbol.name, name) == 0;
    if (e->type != EXPR_FUNCTION) return false;
    if (expr_mentions_symbol(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (expr_mentions_symbol(e->data.function.args[i], name)) return true;
    return false;
}

/* True if the HELD body applies an integer-membership predicate to a
 * subexpression that mentions the iterator -- `If[EvenQ[k], k, -k]`,
 * `Boole[PrimeQ[k]]`, and the like.
 *
 * The predicate folds to a definite Boolean for the SYMBOLIC index, so the
 * conditional around it collapses to ONE branch, and the closed-form cascade --
 * which evaluates the body symbolically before telescoping it -- would sum the
 * wrong summand exactly (`Sum[If[EvenQ[k], k, -k], {k, 1, 10}]` telescoped
 * `Sum[-k]` = -55 instead of enumerating to 5).  A finite range must enumerate
 * such a body; a predicate applied only to CONSTANTS (`PrimeQ[7]`) is fine,
 * because it folds identically per term, so the iterator-mention test gates it. */
static bool sum_body_has_index_predicate(const Expr* e, const char* ivar) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && sum_head_is_index_predicate(e->data.function.head->data.symbol.name)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (expr_mentions_symbol(e->data.function.args[i], ivar)) return true;
    }
    /* The index as a Part position ({a, b}[[j]]): at a symbolic j the Part
     * stays unevaluated and the rest of the body is evaluated against that
     * opaque object -- D[a^2, {a, b}[[j]]] folds to 0 -- so the closed form
     * is of a different summand.  Enumerate, as Mathematica does. */
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Part) {
        for (size_t i = 1; i < e->data.function.arg_count; i++)
            if (expr_mentions_symbol(e->data.function.args[i], ivar)) return true;
    }
    if (sum_body_has_index_predicate(e->data.function.head, ivar)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (sum_body_has_index_predicate(e->data.function.args[i], ivar)) return true;
    return false;
}

/* True when `e` is a machine (inexact) number: a Real, or a Complex with at
 * least one Real component.  This is the correctness gate for the compiled fast
 * path: only when the body evaluates to a machine number at a representative
 * index does inexact-contagion through Plus guarantee the interpreter's own
 * total is a machine number too, so the compiled machine sum agrees to machine
 * precision.  An exact Integer / Rational / BigInt / symbolic result must keep
 * the interpreter path to stay exact. */
static bool sum_body_probe_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol.name, SYM_Complex) == 0) {
        const Expr* re = e->data.function.args[0];
        const Expr* im = e->data.function.args[1];
        return (re && re->type == EXPR_REAL) || (im && im->type == EXPR_REAL);
    }
    return false;
}

/*
 * Machine fast path for a finite numeric sum whose body is inexact.  Returns the
 * scalar total as a fresh Expr, or NULL to fall through to the interpreter's
 * enumeration (expand_range).  Two strategies, one correctness gate (the probe):
 *
 *   A. Integer iterator (!is_real): compile the WHOLE Sum[body, {i, lo, hi[,di]}]
 *      as a zero-variable program (the compiler lowers it to one native
 *      accumulation loop; compile_emit_ctrl.c) and read back the scalar total in
 *      a single VM run — no per-index evaluate(), no per-term allocation.
 *
 *   B. Real iterator (is_real, e.g. {i, 0., 10., 0.1}): the whole-Sum lowering
 *      needs integer bounds, so instead compile the BODY as f(i) and accumulate
 *      in C, driving iteration with the interpreter's own iter_range_continue so
 *      the same index points are visited.  On any non-finite point, abandon.
 *
 * The probe (evaluate the body once at the first index) is what keeps exact sums
 * exact: an exact result declines the fast path, leaving expand_range to produce
 * the exact Integer/Rational/BigInt total, or the empty-range identity 0.
 */
static Expr* sum_try_compiled(Expr* f, const IterSpec* s,
                              double min_val, double max_val, double di_val,
                              bool is_real, Expr** probe_out) {
    *probe_out = NULL;
    if (!autocompile_enabled() || s->var->type != EXPR_SYMBOL) return NULL;

    /* Empty range: let expand_range fold the exact identity 0 (a compiled loop
     * would answer 0. — a spurious inexactness). */
    bool nonempty = (di_val > 0) ? (min_val <= max_val + 1e-14)
                                 : (min_val >= max_val - 1e-14);
    if (!nonempty) return NULL;

    /* Runaway guard only — well above SUM_MAX_FINITE_TERMS (see the macro). */
    double nterms = (max_val - min_val) / di_val + 1.0;
    if (!(nterms >= 1.0) || nterms > (double)SUM_MAX_COMPILED_TERMS) return NULL;

    /* Probe the body at the first index for inexactness (see sum_body_probe_inexact).
     * Shadow/restore the iterator because Sum is HoldAll — mirrors expand_range. */
    IterShadow saved = iter_spec_shadow(s->var);
    Expr* first_val = is_real ? expr_new_real(min_val) : expr_copy(s->imin);
    iter_bind(s->var, first_val);
    Expr* probe = evaluate(f);
    bool inexact = sum_body_probe_inexact(probe);
    expr_free(first_val);
    iter_spec_restore(s->var, saved);
    /* The probe IS the first term: hand it to the enumeration rather than
     * evaluate the body at the first index a second time. */
    *probe_out = probe;
    if (!inexact) return NULL;

    if (!is_real) {
        /* Strategy A: only for exact integer bounds/step — the compiler's Sum
         * lowering requires CT_INT bounds.  A non-integer (BigInt/Rational) bound
         * would just bail inside the compiler; skip it to avoid a wasted attempt. */
        if (s->imin->type != EXPR_INTEGER || s->imax->type != EXPR_INTEGER
            || s->di->type != EXPR_INTEGER)
            return NULL;

        int nspec = (s->di->data.integer == 1) ? 3 : 4;
        Expr* spec_args[4];
        spec_args[0] = expr_copy(s->var);
        spec_args[1] = expr_copy(s->imin);
        spec_args[2] = expr_copy(s->imax);
        if (nspec == 4) spec_args[3] = expr_copy(s->di);
        Expr* spec = expr_new_function(expr_new_symbol(SYM_List), spec_args, nspec);
        Expr* whole_args[2] = { expr_copy(f), spec };
        Expr* whole = expr_new_function(expr_new_symbol(SYM_Sum), whole_args, 2);
        Expr* total = autocompile_eval_closed(whole);   /* NULL unless CT_REAL/CT_COMPLEX */
        expr_free(whole);
        return total;
    }

    /* Strategy B: compile the body as f(i), accumulate the machine terms in C. */
    AutoCompiled* ac = autocompile_new(f, (const Expr* const*)&s->var, 1);
    if (!ac) return NULL;   /* body outside the compilable subset → enumerate */
    double _Complex acc = 0.0;
    double val = min_val;
    int64_t count = 0;
    bool ok = true;
    while (iter_range_continue(/*is_real=*/true, /*is_inf=*/false, NULL, NULL,
                               val, max_val, di_val)) {
        if (count >= SUM_MAX_COMPILED_TERMS) { ok = false; break; }
        double _Complex y;
        if (!autocompiled_eval_complex(ac, &val, &y)) { ok = false; break; }
        acc += y;
        val += di_val;
        count++;
    }
    autocompiled_free(ac);
    if (!ok) return NULL;   /* non-finite point → fall back to the interpreter */
    double re = creal(acc), im = cimag(acc);
    if (!isfinite(re) || !isfinite(im)) return NULL;
    return im == 0.0 ? expr_new_real(re)
                     : make_complex(expr_new_real(re), expr_new_real(im));
}

/*
 * A closed form of a finite numeric range must be as exact as its terms.  The
 * stages see the body at a SYMBOLIC index, where an explicit numericalisation
 * can vanish -- N[1/k] is 1/k -- so Sum[N[1/k], {k, 1, 10^7}] would come back
 * as the exact EulerGamma + PolyGamma[0, 10000001] although every term is a
 * machine number.  When the closed form is exact, evaluate the body at the
 * first index: an inexact term rejects it (consumed; NULL), and the caller's
 * compiled machine loop sums the terms, as the reference does.  An inexact
 * closed form (Sum[1.0 k^2, ...]) or an exact term keeps it.
 */
static Expr* sum_keep_inexact(Expr* cf, Expr* f, const IterSpec* s) {
    if (!cf || internal_contains_inexact(cf)) return cf;
    IterShadow saved = iter_spec_shadow(s->var);
    iter_bind(s->var, s->imin);
    Expr* probe = evaluate(f);
    iter_spec_restore(s->var, saved);
    bool inexact = internal_contains_inexact(probe);
    expr_free(probe);
    if (inexact) { expr_free(cf); return NULL; }
    return cf;
}

/*
 * Closed form of a finite sum with an exact, non-unit step,
 *
 *     Sum[f(i), {i, imin, imax, di}]  =  Sum[f(imin + di j), {j, 0, M}],
 *     M = Floor[(imax - imin)/di],
 *
 * the reindexing that turns any arithmetic progression into a unit-step range
 * the closed-form stages accept (they assume di == 1).  It is exact for every
 * body, so a long stepped range costs what its unit-step twin costs instead of
 * one evaluation (and one big-number addition) per term.  The body is
 * evaluated once at the symbolic index -- what the stages do themselves -- with
 * the iterator shadowed, then reindexed.  NULL when M is not a non-negative
 * integer (an empty range must fold to 0 by enumeration), no stage applies, or
 * the answer still mentions the dummy index.
 */
static Expr* sum_closed_form_stepped(SumMethod method, Expr* f, const IterSpec* s) {
    if (s->var->type != EXPR_SYMBOL) return NULL;
    Expr* span = eval_and_free(expr_new_function(expr_new_symbol(SYM_Floor),
        (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Plus),
                           (Expr*[]){ expr_copy(s->imax),
                                      expr_new_function(expr_new_symbol(SYM_Times),
                                          (Expr*[]){ expr_new_integer(-1),
                                                     expr_copy(s->imin) }, 2) }, 2),
                       expr_new_function(expr_new_symbol(SYM_Power),
                           (Expr*[]){ expr_copy(s->di), expr_new_integer(-1) }, 2) }, 2) }, 1));
    bool span_ok = expr_is_integer_like(span);
    if (span_ok) {
        span_ok = span->type == EXPR_INTEGER ? span->data.integer >= 0
                                             : mpz_sgn(span->data.bigint) >= 0;
    }
    if (!span_ok) { expr_free(span); return NULL; }
    const char* jname = "Sum`Private`step$";
    Expr* j = expr_new_symbol(jname);

    IterShadow saved = iter_spec_shadow(s->var);
    Expr* body = evaluate(f);
    /* i -> imin + di j */
    Expr* at = expr_new_function(expr_new_symbol(SYM_Plus),
                   (Expr*[]){ expr_copy(s->imin),
                              expr_new_function(expr_new_symbol(SYM_Times),
                                  (Expr*[]){ expr_copy(s->di), expr_copy(j) }, 2) }, 2);
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
                     (Expr*[]){ expr_copy(s->var), at }, 2);
    Expr* reindexed = eval_and_free(expr_new_function(expr_new_symbol(SYM_ReplaceAll),
                          (Expr*[]){ body, rule }, 2));
    Expr* zero = expr_new_integer(0);
    Expr* cf = dispatch_def(method, reindexed, j, zero, span);
    iter_spec_restore(s->var, saved);

    if (cf && iter_expr_mentions_symbol(cf, jname)) { expr_free(cf); cf = NULL; }
    expr_free(zero);
    expr_free(span);
    expr_free(reindexed);
    expr_free(j);
    return cf;
}

/* ------------------------------------------------------------------ */
/*  Poles inside the range                                             */
/* ------------------------------------------------------------------ */

#define SUM_POLE_SCAN_MAX 64     /* distinct integer poles inspected */

typedef enum {
    SUM_POLE_NONE,      /* no integer pole on the range: closed forms apply  */
    SUM_POLE_VALUE,     /* a term is non-finite: the sum is that value       */
    SUM_POLE_HOLD       /* a pole may lie on the range (symbolic bound): hold */
} SumPoleVerdict;

static bool sum_nonfinite(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL)
        return e->data.symbol.name == SYM_ComplexInfinity
            || e->data.symbol.name == SYM_Indeterminate;
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_DirectedInfinity;
}

/*
 * Does a term of Sum[f, {var, imin, imax}] (unit step) sit on a pole?
 *
 * Every closed-form stage sums an analytic continuation -- F(imax+1) - F(imin),
 * Zeta[s, imin - rho], a pFq read off t(imin) -- that is blind to a pole of the
 * summand INSIDE the range.  Sum[1/k^2, {k, 0, Infinity}] came out Pi^2/6 and
 * Sum[1/(k - 3)^2, {k, 1, Infinity}] 5/4 + Pi^2/6, where the k = 0 / k = 3
 * term is 1/0; Sum[1/(k (k + 1)), {k, -1, n}] gave -1 - 1/(1 + n).
 *
 * The integer poles are read off the denominator's linear factors
 * (sum_integer_poles).  For a concrete integer imin:
 *   - a pole p >= imin with imax = Infinity or a concrete imax >= p: the term
 *     there is evaluated (with its Power::infy) and, when it is non-finite,
 *     becomes the value of the sum -- what enumerating the range gives;
 *   - a pole p >= imin with a symbolic imax: whether the range reaches p
 *     depends on the bound, so the Sum is held, as the reference does.
 * A pole the term itself removes (a finite value at p) is no pole.  A symbolic
 * or infinite imin is left to the stages (Sum[1/k^2, {k, m, Infinity}] is
 * PolyGamma[1, m] there, as in the reference).  The body is expensive-checked
 * by the caller; Together + Factor of a small body's denominator is cheap.
 */
static SumPoleVerdict sum_pole_verdict(Expr* f, Expr* var, Expr* imin, Expr* imax,
                                       Expr** value) {
    *value = NULL;
    if (imin->type != EXPR_INTEGER || var->type != EXPR_SYMBOL) return SUM_POLE_NONE;
    /* Upper bound class: unbounded (Infinity, or a big integer past every
     * int64 pole), a machine integer, or anything else (symbolic: hold). */
    bool imax_inf = (imax->type == EXPR_SYMBOL && imax->data.symbol.name == SYM_Infinity)
                 || (imax->type == EXPR_BIGINT && mpz_sgn(imax->data.bigint) > 0);
    bool imax_int = imax->type == EXPR_INTEGER;

    IterShadow saved = iter_spec_shadow(var);
    Expr* fe = evaluate(f);
    int64_t poles[SUM_POLE_SCAN_MAX];
    size_t np = sum_integer_poles(fe, var, poles, SUM_POLE_SCAN_MAX);
    SumPoleVerdict v = SUM_POLE_NONE;
    for (size_t i = 0; i < np && v == SUM_POLE_NONE; i++) {
        int64_t p = poles[i];
        if (p < imin->data.integer) continue;
        if (imax_int && p > imax->data.integer) continue;
        if (!imax_inf && !imax_int) { v = SUM_POLE_HOLD; break; }
        Expr* pe = expr_new_integer(p);
        Expr* t = sum_subst(fe, var, pe);
        expr_free(pe);
        if (sum_nonfinite(t)) { *value = t; v = SUM_POLE_VALUE; }
        else if (t) expr_free(t);
    }
    expr_free(fe);
    iter_spec_restore(var, saved);
    return v;
}

/* VerifyConvergence (default True) for the Sum being evaluated; the builtin
 * sets it from the option and restores it, so a nested Sum sees its own. */
static bool g_sum_verify_convergence = true;

/* A divergent infinite series is held with Sum::div, as the reference does,
 * instead of being handed to a stage that would sum its analytic
 * continuation (Sum[2^k, ...] = -1).  The term is examined at a symbolic,
 * shadowed index; UNKNOWN (a symbolic ratio x^k, k^-s) lets the stages run. */
static bool sum_series_diverges(Expr* f, Expr* var) {
    if (!g_sum_verify_convergence || var->type != EXPR_SYMBOL) return false;
    IterShadow saved = iter_spec_shadow(var);
    Expr* fe = evaluate(f);
    SumConvergence c = sum_series_convergence(fe, var);
    expr_free(fe);
    iter_spec_restore(var, saved);
    if (c != SUM_CONV_NO) return false;
    mth_message("Sum", "div", "Sum does not converge.");
    return true;
}

static Expr* sum_one_spec(Expr* f, Expr* spec, SumMethod method) {
    /* Indefinite form Sum[f, i]: spec is a bare symbol. */
    if (spec->type == EXPR_SYMBOL) {
        Expr* var = spec;
        return dispatch_indef(method, f, var);
    }

    IterSpec s;
    if (!iter_spec_parse_lattice(spec, &s)) return NULL;

    /* A single-element iterator {c}: c copies of the body (c real), or the
     * indefinite form (c a free symbol); see iter_count_spec_for_sum.  A
     * bare count without braces stays unevaluated, as in the reference. */
    if (s.kind == ITER_KIND_COUNT) {
        Expr* conv = NULL;
        IterCountKind k = iter_count_spec_for_sum(spec, &s, "Sum", "Sum`Private`count$", &conv);
        iter_spec_free(&s);
        Expr* r = NULL;
        if (k == ITER_COUNT_INDEFINITE) {
            r = dispatch_indef(method, f, conv);
        } else if (k == ITER_COUNT_RANGE) {
            r = sum_one_spec(f, conv, method);
            /* the dummy never occurs in f, so an answer that mentions it is an
             * unevaluated range form: keep the original call instead */
            if (r && iter_expr_mentions_symbol(r, conv->data.function.args[0]->data.symbol.name)) {
                expr_free(r);
                r = NULL;
            }
        }
        if (conv) expr_free(conv);
        return r;
    }

    if (s.kind == ITER_KIND_LIST) {
        Expr* r = expand_list(f, s.var, s.list);
        iter_spec_free(&s);
        return r;
    }

    /* RANGE: finite numeric span -> explicit expansion; otherwise closed form. */
    double min_val, max_val, di_val;
    bool is_real, is_inf;
    bool numeric = iter_spec_resolve_numeric(&s, /*allow_inf=*/true,
                                             &min_val, &max_val, &di_val,
                                             &is_real, &is_inf);
    /* An indexed iterator variable (v[1]) is only enumerated: the closed-form
     * and compiled paths work on a symbol. A symbolic or infinite range over
     * one is left unevaluated. */
    if (iter_var_is_indexed(s.var)) {
        Expr* r = (numeric && !is_inf)
                      ? expand_range(f, s.var, s.imin, s.imax, s.di, min_val, max_val, di_val, is_real)
                      : NULL;
        iter_spec_free(&s);
        return r;
    }
    if (numeric && !is_inf) {
        /* Closed form first: for a unit-step, non-empty, integer-bounded
         * range a polynomial / geometric / hypergeometric body telescopes to
         * F(imax+1) - F(imin), whose cost is independent of the span width --
         * far cheaper than evaluating the body once per term.  This is exact
         * because the antidifference identity holds for any imin <= imax+1.
         * Guards (all required):
         *   - !is_real:  integer bounds/step, no float edge cases;
         *   - di_val==1: the unit-step antidifference is invalid for di != 1;
         *   - min<=max:  empty ranges must fold to 0 via expansion, not the
         *                telescoping form (which would give a wrong value).
         * The iterator is shadowed because Sum is HoldAll: an outer binding of
         * var would otherwise leak into the held body and the stage args. */
        double nterms = (max_val - min_val) / di_val + 1.0;
        /* Skip the closed-form cascade for a SHORT range whose body is
         * expensive to take apart symbolically -- a large symbolic tree or an
         * array (see sum_body_is_expensive and SUM_ARRAY_EXPAND_MAX). The
         * cascade's cost scales with body size while the enumeration below is
         * bounded and cheap, so on such a body the cascade only churns and
         * fails before falling through to it. */
        bool short_expensive_body =
            nterms <= (double)SUM_ARRAY_EXPAND_MAX && sum_body_is_expensive(f);
        /* A body whose symbolic evaluation collapses an index predicate (EvenQ,
         * PrimeQ, ...) has NO valid closed form -- the cascade would telescope
         * the wrong summand -- so it must enumerate regardless of range length. */
        bool index_pred_body =
            sum_body_has_index_predicate(f, s.var->data.symbol.name);

        /* A short range whose body contains an opaque indexed accessor (v[i] for
         * an undefined indexed variable -- the shape an optimiser objective over
         * Table[v[i], ...] produces) has no closed form, and the cascade only
         * churns Simplify on it before falling through here anyway.  Enumerate it
         * directly.  Mirrors short_expensive_body; the Gosper stage has the same
         * guard for the long-range path. */
        bool short_opaque_index_body =
            nterms <= (double)SUM_ARRAY_EXPAND_MAX
            && sum_body_has_opaque_index(f, s.var);

        /* A range whose end points are symbolic lattice points (Pi .. 3 + Pi,
         * from iter_normalize_bounds) is enumerated when short, as Mathematica
         * does: the closed form is equally exact but hands back a factored
         * 2 (3 + 2 Pi) where the term-by-term sum is 6 + 4 Pi. */
        bool short_symbolic_lattice =
            nterms <= (double)SUM_ARRAY_EXPAND_MAX
            && (!expr_is_integer_like(s.imin) || !expr_is_integer_like(s.imax));
        /* A finite range of at most SUM_ENUMERATE_MAX terms is enumerated, as
         * the reference does ("i is assigned a sequence of values, with f
         * evaluated for each"): the closed form evaluates the body at a
         * SYMBOLIC index, which is wrong for a body that takes the index apart
         * -- Sum[Length[Delete[{2, 3, 4}, m]], {m, 1, 3}] was 9 (it is 6),
         * Sum[dt @@ Delete[list, m], ...] collapsed dt @@ Delete[list, m] to
         * dt[list, m] -- and needlessly changes the form of an ordinary one:
         * Sum[x^i, {i, 0, 5}] is 1 + x + ... + x^5, not (x^6 - 1)/(x - 1).
         * Longer ranges keep the closed form first. */
        bool enumerate_first = nterms <= (double)SUM_ENUMERATE_MAX;
        if (!is_real && di_val == 1.0 && min_val <= max_val && !enumerate_first
            && !short_expensive_body && !index_pred_body && !short_symbolic_lattice
            && !short_opaque_index_body) {
            Expr* pv = NULL;
            SumPoleVerdict pole = sum_pole_verdict(f, s.var, s.imin, s.imax, &pv);
            if (pole == SUM_POLE_VALUE) { iter_spec_free(&s); return pv; }
            IterShadow saved = iter_spec_shadow(s.var);
            Expr* cf = dispatch_def(method, f, s.var, s.imin, s.imax);
            iter_spec_restore(s.var, saved);
            cf = sum_keep_inexact(cf, f, &s);
            if (cf) { iter_spec_free(&s); return cf; }
        }
        /* The same for a long range with an exact non-unit step: the closed
         * form of the reindexed unit-step sum (sum_closed_form_stepped). */
        if (!is_real && di_val != 1.0 && !enumerate_first
            && !short_expensive_body && !index_pred_body && !short_symbolic_lattice
            && !short_opaque_index_body) {
            Expr* cf = sum_keep_inexact(sum_closed_form_stepped(method, f, &s), f, &s);
            if (cf) { iter_spec_free(&s); return cf; }
        }
        /* No closed form: before the slow evaluate()-per-term enumeration, try
         * the compiled machine loop for an inexact body (exact bodies decline). */
        Expr* first = NULL;
        Expr* comp = sum_try_compiled(f, &s, min_val, max_val, di_val, is_real, &first);
        if (comp) { expr_free(first); iter_spec_free(&s); return comp; }
        Expr* r = expand_range_first(f, s.var, s.imin, s.imax, s.di,
                                     min_val, max_val, di_val, is_real, first);
        if (r) { iter_spec_free(&s); return r; }
        /* Span too large to enumerate.  For an ordinary body the closed form is
         * exact and worth trying; for an index-predicate body it would be WRONG,
         * so leave the Sum unevaluated rather than telescope a collapsed summand. */
        if (index_pred_body) { iter_spec_free(&s); return NULL; }
        /* span too large: fall through to closed form */
    }

    /* Symbolic / Infinity / over-wide ranges: only the closed-form stages
     * remain, and they assume a unit step.  A non-unit step has no step-aware
     * closed form here, so leave the Sum held rather than emit the wrong
     * step-1 result. */
    if (!is_unit_step(s.di)) { iter_spec_free(&s); return NULL; }

    if (!sum_body_is_expensive(f)) {
        Expr* pv = NULL;
        SumPoleVerdict pole = sum_pole_verdict(f, s.var, s.imin, s.imax, &pv);
        if (pole != SUM_POLE_NONE) { iter_spec_free(&s); return pv; }   /* NULL holds */
    }
    if (s.imax->type == EXPR_SYMBOL && s.imax->data.symbol.name == SYM_Infinity
        && !sum_body_is_expensive(f)
        && sum_series_diverges(f, s.var)) { iter_spec_free(&s); return NULL; }

    Expr* r = dispatch_def(method, f, s.var, s.imin, s.imax);
    iter_spec_free(&s);
    return r;
}

/* ------------------------------------------------------------------ */
/*  Entry point                                                        */
/* ------------------------------------------------------------------ */

Expr* builtin_sum(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2)
        return NULL;

    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;

    /* Split trailing options from the leading body+specs. */
    SumMethod method = SUM_METHOD_AUTOMATIC;
    size_t nspecs_end = argc;            /* one past the last non-option arg */
    bool verify = true;
    for (size_t i = argc; i > 1; i--) {
        if (is_option(a[i - 1])) {
            const Expr* o = a[i - 1];
            if (o->data.function.args[0]->type == EXPR_SYMBOL
                && o->data.function.args[0]->data.symbol.name == SYM_VerifyConvergence) {
                const Expr* v = o->data.function.args[1];
                verify = !(v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_False);
                nspecs_end = i - 1;
                continue;
            }
            SumMethod m = parse_method_option(a[i - 1]);
            if (m == SUM_METHOD_INVALID) return NULL;
            if (m != SUM_METHOD_AUTOMATIC) method = m;
            nspecs_end = i - 1;
        } else {
            break;
        }
    }
    /* a[0] = body; a[1 .. nspecs_end-1] = iterator specs. */
    if (nspecs_end < 2) return NULL;
    size_t nspecs = nspecs_end - 1;

    Expr* f = a[0];

    /* Multiple iterators: rewrite Sum[f, s1, ..., sk, opts] as
     *   Sum[Sum[f, sk, opts], s1, ..., s_{k-1}, opts]
     * and re-evaluate.  The inner sum is fully evaluated under each binding of
     * the outer variables, so outer bounds may depend on inner variables. */
    if (nspecs > 1) {
        size_t nopts = argc - nspecs_end;

        /* inner = Sum[f, sk, opts...] */
        Expr** in = malloc(sizeof(Expr*) * (2 + nopts));
        in[0] = expr_copy(f);
        in[1] = expr_copy(a[nspecs_end - 1]);
        for (size_t k = 0; k < nopts; k++) in[2 + k] = expr_copy(a[nspecs_end + k]);
        Expr* inner = expr_new_function(expr_new_symbol(SYM_Sum), in, 2 + nopts);
        free(in);

        /* outer = Sum[inner, s1, ..., s_{k-1}, opts...] */
        size_t outc = 1 + (nspecs - 1) + nopts;
        Expr** out = malloc(sizeof(Expr*) * outc);
        size_t w = 0;
        out[w++] = inner;
        for (size_t i = 1; i < nspecs; i++) out[w++] = expr_copy(a[i]);
        for (size_t k = 0; k < nopts; k++) out[w++] = expr_copy(a[nspecs_end + k]);
        Expr* outer = expr_new_function(expr_new_symbol(SYM_Sum), out, outc);
        free(out);

        Expr* r = evaluate(outer);
        expr_free(outer);
        return r;
    }

    bool saved_verify = g_sum_verify_convergence;
    g_sum_verify_convergence = verify;
    Expr* r = sum_one_spec(f, a[1], method);
    g_sum_verify_convergence = saved_verify;
    return r;
}

void sum_init(void) {
    symtab_add_builtin("Sum", builtin_sum);
    symtab_get_def("Sum")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("Sum",
        "Sum[f, {i, imax}] gives the sum of f for i from 1 to imax. "
        "Sum[f, {i, imin, imax}], Sum[f, {i, imin, imax, di}] and "
        "Sum[f, {i, {i1, i2, ...}}] use the standard iterator forms; multiple "
        "iterators give nested sums. Sum[f, i] gives the indefinite sum "
        "(antidifference). Symbolic and infinite sums are evaluated in closed "
        "form via Method -> \"Polynomial\" | \"Geometric\" | \"Gosper\".");

    /* Sub-algorithm stages (registered as Sum`Polynomial, Sum`Geometric,
     * Sum`Gosper). Added incrementally; the cascade tolerates absent stages. */
    void sum_polynomial_init(void);
    sum_polynomial_init();
    void sum_geometric_init(void);
    sum_geometric_init();
    void sum_gosper_init(void);
    sum_gosper_init();
    void sum_hypergeometric_init(void);
    sum_hypergeometric_init();
    void sum_lambert_init(void);
    sum_lambert_init();
    void sum_rational_init(void);
    sum_rational_init();
    void sum_euler_init(void);
    sum_euler_init();
    void sum_euler_nonlinear_init(void);
    sum_euler_nonlinear_init();
    void sum_alternating_init(void);
    sum_alternating_init();
    void sum_trigonometric_init(void);
    sum_trigonometric_init();
    void sum_dirichlet_init(void);
    sum_dirichlet_init();
    void sum_zetaseries_init(void);
    sum_zetaseries_init();
    void sum_logzeta_init(void);
    sum_logzeta_init();
    void sum_logrational_init(void);
    sum_logrational_init();
}

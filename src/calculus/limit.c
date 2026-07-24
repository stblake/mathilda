/* ============================================================================
 * limit.c -- Symbolic limits for Mathilda.
 * ============================================================================
 *
 * Implements the Mathematica-style Limit built-in per the pipeline sketched
 * in plans/limit_candidate_spec.md. The architecture is a layered dispatcher;
 * each layer either resolves the limit and short-circuits or passes the
 * problem down to the next layer:
 *
 *     Layer 0 -- Interface normalization (three calling forms, Direction).
 *     Layer 1 -- Cheap structural fast paths.
 *     Layer 2 -- Series-based evaluation (leverages Series[] natively).
 *     Layer 3 -- Rational-function dispatch (P(x)/Q(x) short-cuts).
 *     Layer 5 -- L'Hospital + logarithmic reduction heuristics.
 *     Layer 6 -- Bound analysis (Interval[] returns).
 *
 * The series layer in Mathilda is powerful enough that it subsumes most
 * classical DELIMITER cases. L'Hospital is reserved for those shapes
 * where Series cannot compute a useful expansion (unknown heads, non-
 * analytic inputs, etc.).
 *
 * Memory conventions follow Mathilda standards: every helper that returns
 * an Expr* returns a freshly-allocated tree owned by the caller. The
 * top-level built-in returns a newly-allocated result on success (the
 * evaluator frees the original `res` for us on a non-NULL return) or
 * NULL to leave `res` unevaluated. In particular we never free `res`
 * ourselves -- that would be a double-free against src/eval.c.
 *
 * The module is intentionally layered with small single-purpose helpers
 * so new test failures can be addressed by extending or swapping a single
 * layer rather than re-plumbing the whole pipeline.
 * ========================================================================= */

#include "limit.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "common.h"
#include "print.h"       /* expr_to_string -- Method::method diagnostic */
#include "rationalize.h"
#include "arithmetic.h"  /* arith_warnings_mute_push/pop -- silences the
                          * Power::infy / Infinity::indet messages that our
                          * internal probes would otherwise emit while poking
                          * at candidate sub-expressions. */
#include "sym_names.h"
#include "gruntz.h"       /* Gruntz mrv-algorithm limit engine (layer_gruntz) */
/* Note: Series and D are invoked symbolically (through the evaluator),
 * not via direct C calls, so series.h / deriv.h are intentionally not
 * included here. Adding the Series and Derivative symbols to the symbol
 * table is the responsibility of series_init / deriv_init, which run
 * before limit_init in core_init(). */

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------- */
/* Configuration                                                           */
/* ---------------------------------------------------------------------- */

/* Starting order for series-based evaluation; doubled up to the cap. */
#define LIMIT_SERIES_START_ORDER   4
#define LIMIT_SERIES_MAX_ORDER    32

/* Maximum recursion depth across nested Limit invocations (iterated, log
 * reduction, 1/t substitution, L'Hospital, etc.). Pathological inputs
 * abort into `unevaluated` rather than blowing the C stack. */
#define LIMIT_MAX_DEPTH           24

/* L'Hospital guardrails. */
#define LIMIT_LH_MAX_ITERATIONS   12
#define LIMIT_LH_MAX_GROWTH       3

/* ---------------------------------------------------------------------- */
/* Internal direction encoding                                             */
/*                                                                         */
/* Two-sided = 0; approach from above (x -> x*+) = +1; from below = -1.    */
/* Note: the public `Direction -> ...` option uses Mathematica's sign      */
/* convention (which is opposite of MACSYMA's); the flip is applied at     */
/* parse time, so the rest of the pipeline only ever sees these values.    */
/* A distinct tag is used for Complexes, where one-sided splitting is      */
/* meaningless.                                                            */
/* ---------------------------------------------------------------------- */
#define LIMIT_DIR_TWOSIDED   0
#define LIMIT_DIR_FROMABOVE  1
#define LIMIT_DIR_FROMBELOW (-1)
#define LIMIT_DIR_COMPLEX    2
/* LIMIT_DIR_REALS is the *explicit* real-line two-sided mode. It differs
 * from the implicit TWOSIDED default at sign-disagreement poles: where
 * TWOSIDED returns ComplexInfinity (matching the complex-plane fall-back
 * Mathematica used for unflagged two-sided limits on rational functions),
 * Direction -> Reals asks for the on-reals answer, which is Indeterminate
 * because the one-sided limits disagree in sign. */
#define LIMIT_DIR_REALS      3
/* Approach along +I direction (upper-half-plane): Direction -> I tells
 * us to pick the "other" branch at a branch cut on the negative real
 * axis. For Sqrt and Log the sign of the imaginary part flips relative
 * to the principal branch. */
#define LIMIT_DIR_IMAGINARY  4

/* ---------------------------------------------------------------------- */
/* Method selection (the Method -> m option)                               */
/* ---------------------------------------------------------------------- */
/* The compute_limit cascade is a sequence of strategy layers. By default  */
/* (Automatic) all of them run in order. A named Method restricts the      */
/* *outermost* compute_limit call to a single strategy group; if that      */
/* group produces no answer the whole Limit is left unevaluated. The       */
/* restriction is deliberately confined to the top-level call (see the     */
/* depth==1 gate in compute_limit) so recursive sub-limits -- one-sided    */
/* probes, L'Hospital, Abs splitting -- keep the full cascade.             */
#define LIMIT_M_AUTOMATIC     0
#define LIMIT_M_SUBSTITUTION  1
#define LIMIT_M_RATIONAL      2
#define LIMIT_M_ASYMPTOTIC    3
#define LIMIT_M_BOUNDED       4
#define LIMIT_M_SERIES        5
#define LIMIT_M_LHOSPITAL     6
#define LIMIT_M_GRUNTZ        7

/* ---------------------------------------------------------------------- */
/* LimitCtx -- threaded through the pipeline                               */
/* ---------------------------------------------------------------------- */
typedef struct {
    Expr* x;       /* limit variable, borrowed */
    Expr* point;   /* limit point, borrowed */
    int   dir;     /* one of LIMIT_DIR_* */
    int   depth;   /* recursion depth guard */
    int   method;  /* one of LIMIT_M_*; enforced only at depth==1 */
} LimitCtx;

/* ---------------------------------------------------------------------- */
/* Tiny Expr builders                                                      */
/* ---------------------------------------------------------------------- */

static Expr* mk_int(int64_t v)           { return expr_new_integer(v); }
static Expr* mk_sym(const char* s)       { return expr_new_symbol(s); }

static Expr* mk_fn1(const char* n, Expr* a) {
    Expr* args[1] = { a };
    return expr_new_function(mk_sym(n), args, 1);
}

static Expr* mk_fn2(const char* n, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(mk_sym(n), args, 2);
}

static Expr* mk_times(Expr* a, Expr* b)  { return mk_fn2("Times", a, b); }
static Expr* mk_neg(Expr* a)             { return mk_times(mk_int(-1), a); }

/* Evaluate + free the source. We rely on the fact that evaluate() copies
 * its input internally. */
static Expr* simp(Expr* e) {
    if (!e) return NULL;
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

/* ---------------------------------------------------------------------- */
/* Predicates                                                              */
/* ---------------------------------------------------------------------- */

static bool is_sym(Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

static int literal_sign(Expr* e);  /* forward */

static bool is_lit_zero(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_REAL)    return e->data.real == 0.0;
    if (e->type == EXPR_BIGINT)  return mpz_sgn(e->data.bigint) == 0;
    return false;
}

/* is_infinity_sym is shared via arithmetic.h (identical semantics). */

static bool is_neg_infinity(Expr* e) {
    /* Canonical form of -Infinity is Times[-1, Infinity]. */
    if (!head_is(e, SYM_Times) || e->data.function.arg_count != 2) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    bool has_inf = is_infinity_sym(a) || is_infinity_sym(b);
    Expr* other = is_infinity_sym(a) ? b : a;
    if (!has_inf) return false;
    if (other->type == EXPR_INTEGER && other->data.integer < 0) return true;
    if (other->type == EXPR_REAL && other->data.real < 0.0)     return true;
    return false;
}

static bool is_complex_infinity(Expr* e) {
    return is_sym(e, "ComplexInfinity");
}

static bool is_indeterminate(Expr* e) {
    return is_sym(e, "Indeterminate");
}

static bool is_directed_infinity(Expr* e) {
    return head_is(e, SYM_DirectedInfinity);
}

/* Returns true for any flavour of infinity or undefined value anywhere
 * inside `e`. The fast paths use this to refuse an answer that "looks"
 * finite but still has an un-simplified Infinity / Indeterminate buried
 * in a sub-expression (e.g. 3 + 6/(Infinity^2 - 2), which Mathilda does
 * not fold). */
static bool is_divergent(Expr* e) {
    if (!e) return true;
    if (is_infinity_sym(e) || is_complex_infinity(e) ||
        is_indeterminate(e) || is_neg_infinity(e) ||
        is_directed_infinity(e)) return true;
    if (e->type == EXPR_FUNCTION) {
        if (is_divergent(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (is_divergent(e->data.function.args[i])) return true;
        }
    }
    return false;
}

/* Structural "does the expression contain the target subtree anywhere?". */
static bool expr_contains(Expr* e, Expr* target) {
    if (!e) return false;
    if (expr_eq(e, target)) return true;
    if (e->type == EXPR_FUNCTION) {
        if (expr_contains(e->data.function.head, target)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (expr_contains(e->data.function.args[i], target)) return true;
        }
    }
    return false;
}

static bool free_of(Expr* e, Expr* target) { return !expr_contains(e, target); }

/* A function head is "known" if it has a C builtin, OwnValues, or DownValues
 * in the symbol table. Anything else (e.g. FractionalPart when that module
 * is not implemented; user symbol `f` with no definition) is treated as
 * opaque: we cannot assume continuity, so Limit refuses to substitute a
 * limit point into it. We additionally treat a curated list of
 * discontinuous heads (Floor, Ceiling, Sign, ...) as opaque even when a
 * builtin exists, because plain substitution would silently pick one
 * side's value at a jump and produce a wrong-looking "clean" answer. */
static bool is_discontinuous_head(const char* name) {
    return strcmp(name, "Floor") == 0
        || strcmp(name, "Ceiling") == 0
        || strcmp(name, "Round") == 0
        || strcmp(name, "FractionalPart") == 0
        || strcmp(name, "IntegerPart") == 0
        || strcmp(name, "Sign") == 0
        || strcmp(name, "UnitStep") == 0
        || strcmp(name, "HeavisideTheta") == 0
        || strcmp(name, "KroneckerDelta") == 0
        || strcmp(name, "DiscreteDelta") == 0
        || strcmp(name, "Piecewise") == 0
        || strcmp(name, "Boole") == 0
        || strcmp(name, "Mod") == 0
        || strcmp(name, "Quotient") == 0;
}

static bool is_known_head_symbol(const char* name) {
    if (is_discontinuous_head(name)) return false;
    SymbolDef* def = symtab_get_def(name);
    return def->builtin_func || def->down_values || def->own_values;
}

/* True iff `e` applies an opaque head to any sub-expression containing `x`.
 * Used by the top-level dispatcher to bail out for shapes like f[x] or
 * FractionalPart[x^2] Sin[x] where no layer has a hope of making progress
 * and returning a symbolic value would be unsafe. */
static bool contains_opaque_head_over(Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* name = e->data.function.head->data.symbol.name;
        if (!is_known_head_symbol(name)) {
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (expr_contains(e->data.function.args[i], x)) return true;
            }
        }
    }
    if (contains_opaque_head_over(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_opaque_head_over(e->data.function.args[i], x)) return true;
    }
    return false;
}

/* Leaf count (used for L'Hospital complexity-growth guardrail). */
static int64_t leaf_count(Expr* e) {
    if (!e) return 0;
    if (e->type != EXPR_FUNCTION) return 1;
    int64_t c = leaf_count(e->data.function.head);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        c += leaf_count(e->data.function.args[i]);
    }
    return c;
}

/* ---------------------------------------------------------------------- */
/* Substitution                                                            */
/* ---------------------------------------------------------------------- */

/* ReplaceAll helper: returns evaluated expr with `from` -> `to` everywhere.
 * Caller owns the returned Expr*; inputs are untouched. */
static Expr* subst_eval(Expr* e, Expr* from, Expr* to) {
    Expr* rule = mk_fn2("Rule", expr_copy(from), expr_copy(to));
    Expr* ra   = mk_fn2("ReplaceAll", expr_copy(e), rule);
    return simp(ra);
}

/* ---------------------------------------------------------------------- */
/* Forward declarations                                                    */
/* ---------------------------------------------------------------------- */
static Expr* compute_limit(Expr* f, LimitCtx* ctx);
static Expr* layer1_fast_paths(Expr* f, LimitCtx* ctx);
static Expr* layer2_series(Expr* f, LimitCtx* ctx);
static Expr* layer3_rational(Expr* f, LimitCtx* ctx);
static Expr* layer5_lhospital(Expr* f, LimitCtx* ctx);
static Expr* layer5_log_reduction(Expr* f, LimitCtx* ctx);
static Expr* layer6_bounded(Expr* f, LimitCtx* ctx);
static Expr* layer_compose_at_infinity(Expr* f, LimitCtx* ctx);
static Expr* layer_bounded_envelope(Expr* f, LimitCtx* ctx);
static Expr* layer_log_merge(Expr* f, LimitCtx* ctx);
static Expr* layer_log_of_finite(Expr* f, LimitCtx* ctx);
static Expr* layer_plus_termwise(Expr* f, LimitCtx* ctx);
static Expr* layer_plus_split_convergent(Expr* f, LimitCtx* ctx);
static Expr* layer_abs_rewrite(Expr* f, LimitCtx* ctx);
static Expr* layer_gruntz(Expr* f, LimitCtx* ctx);
static Expr* rewrite_reciprocal_trig(Expr* e);
static Expr* magnitude_upper_bound(Expr* e, Expr* x, bool var_abs);
static bool  contains_bounded_head(Expr* e);
#define LIMIT_UNKNOWN_GROWTH INT64_MAX
static int64_t growth_exponent_upper(Expr* e, Expr* x);

/* ---------------------------------------------------------------------- */
/* Reciprocal trig normalization                                           */
/*                                                                         */
/* Rewrite Csc/Sec/Cot (and hyperbolic twins) in terms of Sin/Cos/Sinh/    */
/* Cosh. This must run before any layer that substitutes the limit point  */
/* into f, because Mathilda's evaluator aggressively folds                  */
/*     0 * Csc[0] = 0 * ComplexInfinity -> 0                               */
/* whereas the equivalent `0 / Sin[0] = 0/0` survives as an indeterminate */
/* 0/0 form that Series / L'Hospital can actually handle.                  */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Hyperbolic -> exponential rewrite for limits at Infinity                */
/*                                                                         */
/* Sinh[z]   = (E^z - E^-z) / 2                                            */
/* Cosh[z]   = (E^z + E^-z) / 2                                            */
/* Tanh[z]   = (E^z - E^-z) / (E^z + E^-z)                                 */
/*                                                                         */
/* Rewriting these before the Series layer matters at Infinity, where the */
/* Exp form exposes the dominant E^z factor while the Sinh/Cosh form is   */
/* asymptotically opaque. Concretely it resolves (1 + Sinh[x])/Exp[x] at  */
/* Infinity to 1/2: after the rewrite the numerator becomes               */
/* 1 + E^x/2 - E^-x/2 and division by E^x cancels to E^-x + 1/2 - ...,    */
/* which the termwise-Plus layer then folds to 1/2.                       */
/* ---------------------------------------------------------------------- */
/* True iff e is a single-argument application of a hyperbolic head that has a
 * well-defined value at Infinity (so the compose-at-infinity layer can fold it
 * directly instead of going through the exp rewrite). */
static bool is_bare_hyperbolic_call(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 1)
        return false;
    return head_is(e, SYM_Sinh) || head_is(e, SYM_Cosh) || head_is(e, SYM_Tanh) ||
           head_is(e, SYM_Coth) || head_is(e, SYM_Sech) || head_is(e, SYM_Csch);
}

static Expr* rewrite_hyperbolic_to_exp(Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    size_t n = e->data.function.arg_count;
    if (e->data.function.head->type == EXPR_SYMBOL && n == 1) {
        const char* hn = e->data.function.head->data.symbol.name;
        Expr* z = rewrite_hyperbolic_to_exp(e->data.function.args[0]);
        if (hn == SYM_Sinh) {
            /* (E^z - E^-z)/2 */
            Expr* ez  = mk_fn1("Exp", expr_copy(z));
            Expr* enz = mk_fn1("Exp", mk_neg(expr_copy(z)));
            expr_free(z);
            Expr* diff = mk_fn2("Plus", ez, mk_neg(enz));
            return mk_fn2("Times", mk_fn2("Power", mk_int(2), mk_int(-1)), diff);
        }
        if (hn == SYM_Cosh) {
            Expr* ez  = mk_fn1("Exp", expr_copy(z));
            Expr* enz = mk_fn1("Exp", mk_neg(expr_copy(z)));
            expr_free(z);
            Expr* sum = mk_fn2("Plus", ez, enz);
            return mk_fn2("Times", mk_fn2("Power", mk_int(2), mk_int(-1)), sum);
        }
        if (hn == SYM_Tanh) {
            Expr* ez1 = mk_fn1("Exp", expr_copy(z));
            Expr* enz1 = mk_fn1("Exp", mk_neg(expr_copy(z)));
            Expr* num = mk_fn2("Plus", ez1, mk_neg(enz1));
            Expr* ez2 = mk_fn1("Exp", expr_copy(z));
            Expr* enz2 = mk_fn1("Exp", mk_neg(expr_copy(z)));
            expr_free(z);
            Expr* den = mk_fn2("Plus", ez2, enz2);
            return mk_fn2("Times", num, mk_fn2("Power", den, mk_int(-1)));
        }
        expr_free(z);
    }

    Expr* head = rewrite_hyperbolic_to_exp(e->data.function.head);
    Expr** args = (Expr**)malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        args[i] = rewrite_hyperbolic_to_exp(e->data.function.args[i]);
    }
    Expr* out = expr_new_function(head, args, n);
    free(args);
    return out;
}

static Expr* rewrite_reciprocal_trig(Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    size_t n = e->data.function.arg_count;

    if (e->data.function.head->type == EXPR_SYMBOL && n == 1) {
        const char* hn = e->data.function.head->data.symbol.name;
        /* Recurse into the argument first so nested reciprocal-trig
         * expressions are all rewritten. */
        Expr* z = rewrite_reciprocal_trig(e->data.function.args[0]);
        if      (hn == SYM_Csc) {
            Expr* r = mk_fn2("Power", mk_fn1("Sin", z), mk_int(-1));
            return r;
        }
        else if (hn == SYM_Sec) {
            Expr* r = mk_fn2("Power", mk_fn1("Cos", z), mk_int(-1));
            return r;
        }
        else if (hn == SYM_Cot) {
            Expr* r = mk_times(mk_fn1("Cos", expr_copy(z)),
                               mk_fn2("Power", mk_fn1("Sin", z), mk_int(-1)));
            return r;
        }
        else if (hn == SYM_Csch) {
            Expr* r = mk_fn2("Power", mk_fn1("Sinh", z), mk_int(-1));
            return r;
        }
        else if (hn == SYM_Sech) {
            Expr* r = mk_fn2("Power", mk_fn1("Cosh", z), mk_int(-1));
            return r;
        }
        else if (hn == SYM_Coth) {
            Expr* r = mk_times(mk_fn1("Cosh", expr_copy(z)),
                               mk_fn2("Power", mk_fn1("Sinh", z), mk_int(-1)));
            return r;
        }
        else if (hn == SYM_Tan) {
            /* Rewriting Tan -> Sin/Cos lets the Series layer expand
             * around Cos-zeros (e.g. Tan[3x] at x = Pi/2) where a direct
             * Series on Tan would hit the Tan[Pi/2] = ComplexInfinity
             * pole and fold the product to 0 via 0 * ComplexInfinity. */
            Expr* r = mk_times(mk_fn1("Sin", expr_copy(z)),
                               mk_fn2("Power", mk_fn1("Cos", z), mk_int(-1)));
            return r;
        }
        else if (hn == SYM_Tanh) {
            Expr* r = mk_times(mk_fn1("Sinh", expr_copy(z)),
                               mk_fn2("Power", mk_fn1("Cosh", z), mk_int(-1)));
            return r;
        }
        expr_free(z);
    }

    /* Generic recursion: rebuild the function node with rewritten children. */
    Expr* head = rewrite_reciprocal_trig(e->data.function.head);
    Expr** args = (Expr**)malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        args[i] = rewrite_reciprocal_trig(e->data.function.args[i]);
    }
    Expr* out = expr_new_function(head, args, n);
    free(args);
    return out;
}

/* ---------------------------------------------------------------------- */
/* Bounded envelope                                                        */
/*                                                                         */
/* `magnitude_upper_bound(e, x)` returns an expression that dominates      */
/* |e(x)| pointwise for every real x. The bound uses:                     */
/*     |Sin[g]|, |Cos[g]|, |Tanh[g]|   <= 1                                */
/*     |ArcTan[g]|, |ArcCot[g]|        <= Pi/2                             */
/*     |a + b + ...|                   <= |a| + |b| + ...                  */
/*     |a b ...|                       =  |a| |b| ...                      */
/*     |a^n| for integer n >= 0        =  |a|^n                            */
/* Anything else is wrapped in Abs[...]. The bound is useful only when    */
/* the caller can show that its limit is zero: then by squeeze f -> 0.    */
/* ---------------------------------------------------------------------- */
static bool contains_bounded_head(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Sin) || head_is(e, SYM_Cos) ||
        head_is(e, SYM_Tanh) || head_is(e, SYM_ArcTan) ||
        head_is(e, SYM_ArcCot)) return true;
    if (contains_bounded_head(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_bounded_head(e->data.function.args[i])) return true;
    }
    return false;
}

static Expr* magnitude_upper_bound(Expr* e, Expr* x, bool var_abs) {
    if (!e) return mk_int(0);
    if (free_of(e, x)) return mk_fn1("Abs", expr_copy(e));
    if (e->type != EXPR_FUNCTION) {
        /* Bare variable. With var_abs the true magnitude |x| = Abs[x] is
         * returned, so the bound is sound at a finite point too (used by
         * envelope_bounded_power for shapes like x Sin[1/x] at x -> 0).
         * Without it we return `x` (not Abs[x]): the +Infinity-only general
         * squeeze relies on `x` being eventually positive, and keeping Abs
         * around causes Limit[1/Abs[x], x->Infinity] to bottom out through
         * L'Hospital on a non-polynomial denominator with no useful
         * progress. */
        return var_abs ? mk_fn1("Abs", expr_copy(e)) : expr_copy(e);
    }

    if (head_is(e, SYM_Sin) || head_is(e, SYM_Cos) || head_is(e, SYM_Tanh)) {
        return mk_int(1);
    }
    if (head_is(e, SYM_ArcTan) || head_is(e, SYM_ArcCot)) {
        return mk_fn2("Times", mk_fn2("Power", mk_int(2), mk_int(-1)),
                               mk_sym("Pi"));
    }

    if (head_is(e, SYM_Plus)) {
        Expr* sum = mk_int(0);
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            sum = mk_fn2("Plus", sum, magnitude_upper_bound(
                                            e->data.function.args[i], x, var_abs));
        }
        return simp(sum);
    }

    if (head_is(e, SYM_Times)) {
        Expr* prod = mk_int(1);
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            prod = mk_times(prod, magnitude_upper_bound(
                                        e->data.function.args[i], x, var_abs));
        }
        return simp(prod);
    }

    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp  = e->data.function.args[1];
        if (exp->type == EXPR_INTEGER) {
            int64_t n = exp->data.integer;
            if (n >= 0) {
                Expr* bb = magnitude_upper_bound(base, x, var_abs);
                return simp(mk_fn2("Power", bb, mk_int(n)));
            }
            /* negative integer: |1/base^|n|| = 1/|base|^|n| -- keep Abs wrap. */
            Expr* bb = magnitude_upper_bound(base, x, var_abs);
            return simp(mk_fn2("Power", bb, mk_int(n)));
        }
        /* Exp[bounded]: bound by E. */
        if (is_sym(base, "E") && contains_bounded_head(exp) &&
            head_is(exp, SYM_Sin)) {
            return mk_sym("E");
        }
        return mk_fn1("Abs", expr_copy(e));
    }

    /* For Log, Exp, Sqrt etc. when free_of check missed them, default to
     * Abs wrapping (the limit of Abs[Log[...]] etc. is the true
     * magnitude; if it's not zero we fall through anyway). */
    return mk_fn1("Abs", expr_copy(e));
}

/* ---------------------------------------------------------------------- */
/* Layer 0 -- Interface normalization                                      */
/*                                                                         */
/* The calling forms are unpacked in builtin_limit() (below). This section */
/* provides the Direction-option parser.                                   */
/* ---------------------------------------------------------------------- */

/* True iff `e` is a literal imaginary unit or a purely imaginary constant
 * (e.g. I, 2 I, -I, Complex[0, k]). Such values as a Direction argument
 * select approach along the imaginary axis, which we currently route
 * through LIMIT_DIR_COMPLEX. */
static bool is_imaginary_direction(Expr* e) {
    if (is_sym(e, "I")) return true;
    if (head_is(e, SYM_Complex) && e->data.function.arg_count == 2) {
        Expr* re = e->data.function.args[0];
        return (re->type == EXPR_INTEGER && re->data.integer == 0) ||
               (re->type == EXPR_REAL && re->data.real == 0.0);
    }
    if (head_is(e, SYM_Times)) {
        /* k * I, -I, etc. */
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (is_imaginary_direction(e->data.function.args[i])) return true;
        }
    }
    return false;
}

/* Translate a user-facing Direction value to the internal direction tag.
 * Returns true on success; false for an unrecognised / symbolic value.    */
static bool parse_direction(Expr* dir, int* out) {
    if (!dir) { *out = LIMIT_DIR_TWOSIDED; return true; }
    if (is_sym(dir, "Reals"))     { *out = LIMIT_DIR_REALS;     return true; }
    if (is_sym(dir, "TwoSided"))  { *out = LIMIT_DIR_TWOSIDED;  return true; }
    if (is_sym(dir, "Complexes")) { *out = LIMIT_DIR_COMPLEX;   return true; }
    if (dir->type == EXPR_STRING) {
        if (strcmp(dir->data.string, "TwoSided")   == 0) { *out = LIMIT_DIR_TWOSIDED;  return true; }
        /* Mathematica sign convention: "FromAbove" == Direction -> -1,
         * "FromBelow" == Direction -> +1. Internally +1 = FromAbove,
         * -1 = FromBelow (so the math layer never has to remember the flip). */
        if (strcmp(dir->data.string, "FromAbove")  == 0) { *out = LIMIT_DIR_FROMABOVE; return true; }
        if (strcmp(dir->data.string, "FromBelow")  == 0) { *out = LIMIT_DIR_FROMBELOW; return true; }
    }
    if (dir->type == EXPR_INTEGER) {
        /* Apply the sign flip once and only here. */
        if (dir->data.integer == -1) { *out = LIMIT_DIR_FROMABOVE; return true; }
        if (dir->data.integer == +1) { *out = LIMIT_DIR_FROMBELOW; return true; }
    }
    /* Imaginary direction (I, k I for positive k, Complex[0, positive]):
     * approach the branch point from the upper half plane. Tagged with
     * LIMIT_DIR_IMAGINARY so the branch-cut post-pass in builtin_limit
     * can flip the imaginary part for Sqrt/Log etc. */
    if (is_imaginary_direction(dir)) { *out = LIMIT_DIR_IMAGINARY; return true; }
    return false;
}

/* Translate a user-facing Method value to the internal method tag. Accepts
 * the symbol Automatic and the six strategy-group names (as strings).
 * Returns true on success; false for an unrecognised value.               */
static bool parse_method(Expr* m, int* out) {
    if (!m) { *out = LIMIT_M_AUTOMATIC; return true; }
    if (is_sym(m, "Automatic")) { *out = LIMIT_M_AUTOMATIC; return true; }
    if (is_sym(m, "Gruntz"))    { *out = LIMIT_M_GRUNTZ;    return true; }
    if (m->type == EXPR_STRING) {
        const char* s = m->data.string;
        if (strcmp(s, "Automatic")        == 0) { *out = LIMIT_M_AUTOMATIC;    return true; }
        if (strcmp(s, "Substitution")     == 0) { *out = LIMIT_M_SUBSTITUTION; return true; }
        if (strcmp(s, "RationalFunction") == 0) { *out = LIMIT_M_RATIONAL;     return true; }
        if (strcmp(s, "Asymptotic")       == 0) { *out = LIMIT_M_ASYMPTOTIC;   return true; }
        if (strcmp(s, "Bounded")          == 0) { *out = LIMIT_M_BOUNDED;      return true; }
        if (strcmp(s, "Series")           == 0) { *out = LIMIT_M_SERIES;       return true; }
        if (strcmp(s, "LHospital")        == 0) { *out = LIMIT_M_LHOSPITAL;    return true; }
        if (strcmp(s, "Gruntz")           == 0) { *out = LIMIT_M_GRUNTZ;       return true; }
    }
    return false;
}

/* ---------------------------------------------------------------------- */
/* Direction helpers                                                       */
/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
/* Infinity/sign helpers                                                   */
/* ---------------------------------------------------------------------- */

/* Sign of a numeric literal: +1 / 0 / -1, or 0 if not decidable. */
static int literal_sign(Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER) {
        if (e->data.integer > 0) return +1;
        if (e->data.integer < 0) return -1;
        return 0;
    }
    if (e->type == EXPR_REAL) {
        if (e->data.real > 0.0) return +1;
        if (e->data.real < 0.0) return -1;
        return 0;
    }
    if (e->type == EXPR_BIGINT) {
        int s = mpz_sgn(e->data.bigint);
        if (s > 0) return +1;
        if (s < 0) return -1;
        return 0;
    }
    /* Rational[n, d]: sign follows numerator's. */
    if (head_is(e, SYM_Rational) && e->data.function.arg_count == 2) {
        return literal_sign(e->data.function.args[0]);
    }
    /* Positive real constants. */
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        if (nm == SYM_Pi || nm == SYM_E || nm == SYM_EulerGamma ||
            nm == SYM_GoldenRatio || nm == SYM_Catalan || nm == SYM_Degree)
            return +1;
        return 0;
    }
    /* Power[b, e]: a positive base raised to a real exponent is positive
     * (Sqrt[2] = 2^(1/2), 1/Sqrt[2] = 2^(-1/2), ...). A negative base with
     * an integer exponent follows the parity of the exponent; any other
     * shape (negative base, fractional exponent -> complex) is undecidable. */
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* b = e->data.function.args[0];
        Expr* ex = e->data.function.args[1];
        int sb = literal_sign(b);
        if (sb > 0) return +1;
        if (sb < 0 && ex->type == EXPR_INTEGER)
            return (ex->data.integer % 2 == 0) ? +1 : -1;
        return 0;
    }
    /* Times[...]: product of the signs of every factor. Undecidable (0) as
     * soon as any factor's sign is unknown -- covers 1/(2 Sqrt[2]) ArcTan...
     * coefficients whose sign the old leading-factor-only test missed. */
    if (head_is(e, SYM_Times) && e->data.function.arg_count > 0) {
        int prod = 1;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            int s = literal_sign(e->data.function.args[i]);
            if (s == 0) return 0;
            prod *= s;
        }
        return prod;
    }
    /* Plus[...]: decidable only when every summand shares one nonzero sign. */
    if (head_is(e, SYM_Plus) && e->data.function.arg_count > 0) {
        int common = 0;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            int s = literal_sign(e->data.function.args[i]);
            if (s == 0) return 0;
            if (common == 0) common = s;
            else if (common != s) return 0;
        }
        return common;
    }
    /* Log[c] for a positive real literal c: sign follows c vs 1.
     * This is the narrow shape the log-reduction layer leaves behind
     * after folding b^(g(x)) = Exp[g(x) Log[b]]: whether the Limit is
     * +Infinity or -Infinity depends on Sign[Log[b]], which the series
     * sign test would otherwise miss (Log[3] is symbolic to literal_sign
     * without this arm). */
    if (head_is(e, SYM_Log) && e->data.function.arg_count == 1) {
        Expr* a = e->data.function.args[0];
        if (a->type == EXPR_INTEGER) {
            if (a->data.integer >  1) return +1;
            if (a->data.integer == 1) return  0;
            if (a->data.integer >  0) return -1;  /* 0 < c < 1: impossible for Integer, no-op */
        } else if (a->type == EXPR_REAL) {
            if (a->data.real > 1.0) return +1;
            if (a->data.real == 1.0) return 0;
            if (a->data.real > 0.0) return -1;
        } else if (a->type == EXPR_BIGINT) {
            if (mpz_cmp_ui(a->data.bigint, 1) > 0) return +1;
            if (mpz_cmp_ui(a->data.bigint, 1) == 0) return 0;
        } else if (head_is(a, SYM_Rational) && a->data.function.arg_count == 2) {
            /* Positive rational: compare numerator and denominator. */
            Expr* rn = a->data.function.args[0];
            Expr* rd = a->data.function.args[1];
            if (rn->type == EXPR_INTEGER && rd->type == EXPR_INTEGER &&
                rd->data.integer > 0 && rn->data.integer > 0) {
                if (rn->data.integer > rd->data.integer) return +1;
                if (rn->data.integer < rd->data.integer) return -1;
                return 0;
            }
        }
    }
    return 0;
}

/* Produce +Infinity or -Infinity from a sign. */
static Expr* signed_infinity(int sign) {
    if (sign >= 0) return mk_sym("Infinity");
    return mk_neg(mk_sym("Infinity"));
}

/* ---------------------------------------------------------------------- */
/* Continuous substitution (Layer 1c)                                      */
/*                                                                         */
/* Substitute the limit point into f and check whether the result is a     */
/* clean, finite expression free of the limit variable. We accept any      */
/* answer that is not divergent and does not still mention x.              */
/* ---------------------------------------------------------------------- */
/* True iff `e` contains an exponential form Power[base, exp] with the
 * limit variable appearing in `exp`. These shapes (x^x, (1+Ax)^(1/x),
 * etc.) are classic indeterminate 1^inf / 0^0 / inf^0 forms that the
 * continuous-substitution fast path cannot handle correctly because
 * Mathilda's arithmetic happily folds 1^ComplexInfinity to 1. We refuse
 * the fast path and let the log-reduction layer deal with them. */
static bool has_var_in_exponent(Expr* e, Expr* x) {
    if (!e) return false;
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* exp = e->data.function.args[1];
        if (expr_contains(exp, x)) return true;
    }
    if (e->type == EXPR_FUNCTION) {
        if (has_var_in_exponent(e->data.function.head, x)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (has_var_in_exponent(e->data.function.args[i], x)) return true;
        }
    }
    return false;
}

/* True iff `e` is a numeric literal (integer, bigint, real, or Rational).
 * Such a limit point is "plain old number": if substituting it produces a
 * finite value with no residual variable, the expression is analytic
 * there and we can return the substituted value directly, skipping the
 * expensive Together / Numerator / Denominator pipeline used by the
 * generic continuous-substitution path. */
static bool is_numeric_literal_point(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL ||
        e->type == EXPR_BIGINT) return true;
    if (head_is(e, SYM_Rational) && e->data.function.arg_count == 2) return true;
    /* A symbolic `Times[-1, <numeric>]` (e.g. user wrote `-5`) is also a
     * numeric literal once evaluated; the caller has already evaluated
     * the point, so we only need to recognise the canonical forms. */
    if (head_is(e, SYM_Times) && e->data.function.arg_count == 2 &&
        e->data.function.args[0]->type == EXPR_INTEGER &&
        e->data.function.args[0]->data.integer == -1) {
        return is_numeric_literal_point(e->data.function.args[1]);
    }
    return false;
}

/* Cheap fast path: if the limit point is a plain numeric literal and a
 * single substitution + evaluate produces a clean result (not divergent,
 * no residual limit variable), return it. This avoids ever running
 * Together on expressions like `Log[1 - (Log[Exp[z]/z - 1] + Log[z])/z]/z`
 * at z = 100, where Together's sub-expression normalisation can spin in
 * the evaluator even though the input is analytic at the point. */
/* True iff `e` contains a `Power[base, exp]` subterm whose exponent
 * diverges when the limit point is substituted. Divergent exponents are
 * the classic 1^inf / 0^0 / inf^0 indeterminate seeds, and Mathilda
 * arithmetic folds them to "clean" (but wrong) values -- we must refuse
 * the direct-substitution fast path in that case. */
static bool has_divergent_exponent_at(Expr* e, Expr* x, Expr* point) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* exp = e->data.function.args[1];
        if (expr_contains(exp, x)) {
            Expr* exp_at = subst_eval(exp, x, point);
            bool bad = is_divergent(exp_at);
            expr_free(exp_at);
            if (bad) return true;
        }
    }
    if (has_divergent_exponent_at(e->data.function.head, x, point)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_divergent_exponent_at(e->data.function.args[i], x, point)) return true;
    }
    return false;
}

static Expr* try_numeric_point_substitution(Expr* f, LimitCtx* ctx) {
    if (!is_numeric_literal_point(ctx->point)) return NULL;

    /* Refuse when any Power[_, exp] has x in the exponent AND exp diverges
     * at the point. That's the 1^inf / 0^0 / inf^0 indeterminate family
     * that Mathilda's arithmetic silently folds to a plausible-looking but
     * wrong answer (e.g. 1^ComplexInfinity -> 1). */
    if (has_divergent_exponent_at(f, ctx->x, ctx->point)) return NULL;

    /* Denominator-at-point zero check. If Together's denominator vanishes
     * at the point we have a potential 0/0 or non-zero/0 shape; either
     * way, direct substitution is unsafe -- defer. */
    Expr* tog = simp(mk_fn1("Together", expr_copy(f)));
    Expr* den = simp(mk_fn1("Denominator", expr_copy(tog)));
    Expr* den_at = subst_eval(den, ctx->x, ctx->point);
    expr_free(den);
    bool den_bad = is_lit_zero(den_at) || is_divergent(den_at);
    expr_free(den_at);
    if (den_bad) { expr_free(tog); return NULL; }

    /* Substitute into the *Together-normalised* form. This matters when
     * the cancel-first form has a non-zero denominator at the point while
     * the original has a 0/0 shape that Mathilda's arithmetic folds to 0.
     * Example: (r Sin[t])/(r (Cos[t]+Sin[t])) at r = 0 -- Together
     * cancels the r and we get Sin[t]/(Cos[t]+Sin[t]) without the 0/0
     * trap. */
    Expr* sub = subst_eval(tog, ctx->x, ctx->point);
    expr_free(tog);
    if (!sub) return NULL;
    if (is_divergent(sub) || expr_contains(sub, ctx->x)) {
        expr_free(sub);
        return NULL;
    }
    return sub;
}

static Expr* try_continuous_substitution(Expr* f, LimitCtx* ctx) {
    /* Skip for non-finite limit points; those need a different path. */
    if (is_infinity_sym(ctx->point) || is_neg_infinity(ctx->point) ||
        is_complex_infinity(ctx->point)) {
        return NULL;
    }
    /* Refuse exponential-indeterminate shapes; see has_var_in_exponent. */
    if (has_var_in_exponent(f, ctx->x)) return NULL;

    /* Guard rails. Mathilda's evaluator aggressively folds 0 * ComplexInfinity
     * to 0 and 1/0 to ComplexInfinity, which would produce misleading
     * "clean" answers for shapes like Sin[x]/x at x=0. Rule out anything
     * that could be a division by zero, log of zero, or similar:
     *
     *   - Together-normalize f to a single quotient N/D and evaluate the
     *     denominator at the point. If D(point) = 0, the expression is
     *     not continuous there and we bail to the next layer.
     *   - Rejecting "looks clean" is not enough: we also need to refuse
     *     inputs whose closed-form substitution hit Power::infy along the
     *     way, so we scan intermediate sub-expressions too.
     */
    Expr* tog = simp(mk_fn1("Together", expr_copy(f)));
    Expr* num = simp(mk_fn1("Numerator",   expr_copy(tog)));
    Expr* den = simp(mk_fn1("Denominator", expr_copy(tog)));
    expr_free(tog);
    Expr* den_at = subst_eval(den, ctx->x, ctx->point);
    bool den_bad = is_lit_zero(den_at) || is_divergent(den_at);
    expr_free(den_at);
    if (den_bad) {
        expr_free(num); expr_free(den);
        return NULL;
    }
    /* Denominator is safe: compose N(point)/D(point). We use the original
     * expression for the substitution so that trivially-evaluable heads
     * (Sin, Cos, ...) at the point still simplify exactly. */
    Expr* num_at = subst_eval(num, ctx->x, ctx->point);
    Expr* den_at2 = subst_eval(den, ctx->x, ctx->point);
    expr_free(num); expr_free(den);
    if (is_divergent(num_at) || is_divergent(den_at2) ||
        expr_contains(num_at, ctx->x) || expr_contains(den_at2, ctx->x)) {
        expr_free(num_at); expr_free(den_at2);
        return NULL;
    }
    Expr* s = simp(mk_fn2("Divide", num_at, den_at2));
    if (is_divergent(s) || expr_contains(s, ctx->x)) {
        expr_free(s);
        return NULL;
    }
    return s;
}

/* ---------------------------------------------------------------------- */
/* Layer 1 -- Fast paths                                                   */
/* ---------------------------------------------------------------------- */
static Expr* layer1_fast_paths(Expr* f, LimitCtx* ctx) {
    /* free of x -> constant expression, pass through. */
    if (free_of(f, ctx->x)) return expr_copy(f);

    /* f is exactly the variable: answer is the target point. */
    if (expr_eq(f, ctx->x)) return expr_copy(ctx->point);

    /* Layer 1a: cheap numeric-point substitution. Skips Together entirely,
     * so expressions that happen to be analytic at a numeric point return
     * in microseconds even if Together would have triggered a pathological
     * evaluator loop. */
    Expr* r = try_numeric_point_substitution(f, ctx);
    if (r) return r;

    /* Layer 1b: generic continuous substitution via Together +
     * Numerator/Denominator zero-check. Handles non-literal points and
     * shapes where the naive substitution would emit Power::infy. */
    return try_continuous_substitution(f, ctx);
}

/* ---------------------------------------------------------------------- */
/* Layer 2 -- Series-based evaluation                                      */
/*                                                                         */
/* Call Series[f, {x, point, k}] with increasing order k until we can      */
/* read off a leading term. Handles Infinity / MInfinity by substitution.  */
/* ---------------------------------------------------------------------- */

/* Try to interpret `s` as a SeriesData and read off the leading term.
 * Returns a freshly-allocated Expr* for the limit on success, or NULL if
 * the input is not a SeriesData / the direction disagreement is
 * ambiguous. Fills *unbounded_sign with -1/0/+1 when the answer is a
 * directional infinity (caller may re-sign). */
static Expr* read_leading_term_limit(Expr* s, LimitCtx* ctx) {
    /* s may be free of x (Series returned a constant) -- just return it. */
    if (!s) return NULL;
    if (free_of(s, ctx->x)) return expr_copy(s);

    /* Unwrap the Series may be wrapped in Times[prefactor, SeriesData[...]]. */
    Expr* prefactor = NULL;
    Expr* sd = s;
    if (head_is(s, SYM_Times)) {
        /* find the SeriesData factor (at most one) */
        Expr* found = NULL;
        Expr* pf = mk_int(1);
        for (size_t i = 0; i < s->data.function.arg_count; i++) {
            Expr* a = s->data.function.args[i];
            if (head_is(a, SYM_SeriesData) && !found) {
                found = a;
            } else {
                pf = simp(mk_times(pf, expr_copy(a)));
            }
        }
        if (found) {
            sd = found;
            prefactor = pf;
        } else {
            expr_free(pf);
        }
    }

    if (!head_is(sd, SYM_SeriesData) || sd->data.function.arg_count != 6) {
        /* Not a recognised series shape; if the expression is still in a
         * closed form free of x we already handled that above. */
        if (prefactor) expr_free(prefactor);
        return NULL;
    }

    /* SeriesData[var, x0, {a0, a1, ...}, nmin, nmax, den] */
    Expr* var   = sd->data.function.args[0];
    Expr* x0    = sd->data.function.args[1];
    Expr* coefs = sd->data.function.args[2];
    Expr* nmin_e= sd->data.function.args[3];
    Expr* den_e = sd->data.function.args[5];

    if (!head_is(coefs, SYM_List) || nmin_e->type != EXPR_INTEGER ||
        den_e->type != EXPR_INTEGER) {
        if (prefactor) expr_free(prefactor);
        return NULL;
    }
    int64_t nmin = nmin_e->data.integer;
    int64_t den  = den_e->data.integer;
    (void)var; (void)x0;

    /* Find first nonzero coefficient index. */
    size_t k = 0;
    size_t kn = coefs->data.function.arg_count;
    while (k < kn && is_lit_zero(coefs->data.function.args[k])) k++;
    if (k == kn) {
        /* All coefficients are literally zero within the computed range.
         * The leading behaviour is below the O-term; we can't read off a
         * value without expanding further. The caller will retry at a
         * higher order. */
        if (prefactor) expr_free(prefactor);
        return NULL;
    }

    int64_t leading_num = nmin + (int64_t)k;
    Expr* leading_coef = expr_copy(coefs->data.function.args[k]);

    /* A genuine power-series coefficient cannot depend on the expansion
     * variable. Mathilda's Series does not always enforce this -- for
     * asymptotic shapes like Log[E^x - E^a] at x -> a it emits a Log[x-a]
     * term as the "constant coefficient" (the logarithmic part of an
     * asymptotic expansion). Treating that as the limit value gives a
     * result that still depends on x, which is structurally wrong.
     * Escalating the Series order won't remove the residual -- the
     * Log-term is part of the expansion's form, not a truncation
     * artefact -- so we bail out of the whole series layer by setting the
     * variable to NULL and returning. The caller's escalation loop will
     * move on to subsequent attempts. We can't easily tell layer2_series
     * to skip further k values without a side-channel flag; instead, to
     * avoid an expensive retry at k=32 (Series on these shapes is O(n^3)
     * in the expansion order), we detect the same residual early. */
    if (expr_contains(leading_coef, ctx->x)) {
        /* For a leading exponent of exactly 0 the residual x is the
         * logarithmic-constant part of an asymptotic expansion (e.g.
         * ExpIntegralEi[x] at x -> 0, whose x^0 coefficient is
         * EulerGamma + Log[x]). Every higher-order term carries a strictly
         * positive power of x and vanishes at the point, so the limit of f
         * equals the limit of this leading coefficient. Recurse to resolve
         * it; only accept an answer that no longer mentions x. If the
         * recursion can't make progress we fall through to the bail below,
         * preserving the old behaviour for genuinely-stuck shapes. */
        if (leading_num == 0) {
            LimitCtx sub = *ctx; sub.depth += 1;
            Expr* lim_c = compute_limit(leading_coef, &sub);
            if (lim_c && !expr_contains(lim_c, ctx->x)) {
                expr_free(leading_coef);
                if (prefactor) lim_c = simp(mk_times(prefactor, lim_c));
                return lim_c;
            }
            if (lim_c) expr_free(lim_c);
        }
        expr_free(leading_coef);
        if (prefactor) expr_free(prefactor);
        return NULL;
    }

    /* Decide the limit from the leading exponent leading_num/den.
     *    > 0  -> 0
     *    = 0  -> leading coefficient
     *    < 0  -> directional infinity
     * All three cases apply whether the expansion is around a finite point
     * or around infinity (we normalize to t -> 0+ before calling Series).
     */
    Expr* result = NULL;
    if (leading_num > 0) {
        result = mk_int(0);
        expr_free(leading_coef);
    } else if (leading_num == 0) {
        result = leading_coef;
    } else {
        /* Infinity with a sign. For one-sided at 0+, sign = sign(coef).
         * For 0-, multiply by (-1)^leading_num (when integer); if the
         * exponent has a nontrivial denominator the side is complex and
         * we only answer if the direction is FromAbove. */
        int coef_sign = literal_sign(leading_coef);
        if (coef_sign == 0) {
            /* Unknown sign -- only meaningful when the coefficient is a
             * pure constant (no residual limit variable). If x survived
             * into the coefficient (e.g. Log[x] in Log[x]/x at x=0), a
             * `DirectedInfinity[Log[x]]` answer is misleading; bail so the
             * remaining layers get a shot. */
            if (expr_contains(leading_coef, ctx->x)) {
                expr_free(leading_coef);
                if (prefactor) expr_free(prefactor);
                return NULL;
            }
            result = mk_fn1("DirectedInfinity", leading_coef);
        } else {
            int side_factor = +1;
            if (ctx->dir == LIMIT_DIR_FROMBELOW) {
                /* (-1)^leading_num if integer; odd -> flip, even -> keep */
                if (den == 1) {
                    side_factor = (leading_num % 2 == 0) ? +1 : -1;
                } else {
                    /* Fractional exponent on the negative side: branch cut. */
                    expr_free(leading_coef);
                    if (prefactor) expr_free(prefactor);
                    return NULL;
                }
            } else if (ctx->dir == LIMIT_DIR_TWOSIDED) {
                /* Default (un-qualified) two-sided: for an odd-order integer
                 * pole the two one-sided limits disagree in sign. Pick the
                 * complex-plane fall-back (ComplexInfinity). For non-integer
                 * exponents (sqrt-type), the negative side is off the real
                 * line anyway and we bail to let the caller retry with a
                 * sharper layer. */
                if (den == 1 && (leading_num % 2 != 0)) {
                    expr_free(leading_coef);
                    result = mk_sym("ComplexInfinity");
                    goto apply_prefactor;
                }
            } else if (ctx->dir == LIMIT_DIR_REALS) {
                /* Explicit Direction -> Reals. Unlike the default TWOSIDED
                 * mode, this specifically asks for the real-line answer: a
                 * pole with disagreeing signs on the two sides has no real
                 * limit, so return Indeterminate. Non-integer exponents hit
                 * a branch cut on the negative side; also Indeterminate. */
                if (den != 1 || (leading_num % 2 != 0)) {
                    expr_free(leading_coef);
                    result = mk_sym("Indeterminate");
                    goto apply_prefactor;
                }
            } else if (ctx->dir == LIMIT_DIR_COMPLEX) {
                /* Direction -> Complexes asks for the radial / all-
                 * direction answer. Any isolated pole (integer or
                 * fractional order) is ComplexInfinity. */
                expr_free(leading_coef);
                result = mk_sym("ComplexInfinity");
                goto apply_prefactor;
            }
            expr_free(leading_coef);
            result = signed_infinity(coef_sign * side_factor);
        }
    }

apply_prefactor:
    if (prefactor) {
        /* Fold the prefactor into the result. For Infinity we'd like to
         * keep a sensible sign; let the evaluator's arithmetic handle it. */
        result = simp(mk_times(prefactor, result));
    }
    return result;
}

static Expr* layer2_series(Expr* f, LimitCtx* ctx) {
    /* We rely on Series's native at-Infinity handling: Series[f, {x,
     * Infinity, k}] produces SeriesData[Power[x, -1], 0, ...] -- i.e.
     * the expansion variable is 1/x. The leading-term reader only looks
     * at the exponent and its sign, so no additional substitution is
     * required for +Infinity.
     *
     * For -Infinity we can't feed Series directly (it only recognises
     * the bare Infinity symbol as a special expansion point), so we
     * substitute x -> -y first, reduce to a +Infinity problem in y, and
     * recurse. The direction flips accordingly. */
    Expr* f_use = f;
    Expr* x_use = ctx->x;
    Expr* x0_use = ctx->point;
    Expr* x0_owned = NULL;
    Expr* f_owned = NULL;
    Expr* y_sym = NULL;
    int effective_dir = ctx->dir;

    if (is_infinity_sym(ctx->point)) {
        effective_dir = LIMIT_DIR_FROMABOVE;
    } else if (is_neg_infinity(ctx->point)) {
        y_sym = mk_sym("$LimitNegInfVar$");
        Expr* neg_y = mk_neg(expr_copy(y_sym));
        f_owned = subst_eval(f, ctx->x, neg_y);
        expr_free(neg_y);
        f_use = f_owned;
        x_use = y_sym;
        x0_owned = mk_sym("Infinity");
        x0_use = x0_owned;
        effective_dir = LIMIT_DIR_FROMABOVE;
    }

    Expr* result = NULL;
    for (int64_t k = LIMIT_SERIES_START_ORDER;
         k <= LIMIT_SERIES_MAX_ORDER;
         k *= 2) {
        Expr* spec = expr_new_function(
            mk_sym("List"),
            (Expr*[]){ expr_copy(x_use), expr_copy(x0_use), mk_int(k) }, 3);
        Expr* call = mk_fn2("Series", expr_copy(f_use), spec);
        Expr* s    = simp(call);
        LimitCtx leaf = { x_use, x0_use, effective_dir, ctx->depth, ctx->method };
        result = read_leading_term_limit(s, &leaf);
        /* If the series expansion still mentions the limit variable in
         * its coefficients (asymptotic shapes like Log[E^x-E^a] at x->a
         * whose leading "coefficient" is Log[x-a]), escalating to higher
         * order won't help and is expensive (O(k^3) for common heads).
         * Short-circuit and hand off to the later layers. */
        bool stuck = (result == NULL) && s && expr_contains(s, x_use);
        expr_free(s);
        if (result) break;
        if (stuck) break;
    }

    if (f_owned)  expr_free(f_owned);
    if (x0_owned) expr_free(x0_owned);
    if (y_sym)    expr_free(y_sym);
    return result;
}

/* ---------------------------------------------------------------------- */
/* Layer 3 -- Rational function short-cut                                  */
/*                                                                         */
/* Handles the classical P(x)/Q(x) cases where the series layer would be   */
/* overkill. We lean on Together / Numerator / Denominator rather than     */
/* re-implementing polynomial arithmetic here.                             */
/* ---------------------------------------------------------------------- */
/* Degree of a univariate polynomial in `x`: use CoefficientList which
 * is a native builtin in Mathilda. Returns -1 on failure. */
static int poly_degree_in(Expr* p, Expr* x) {
    Expr* call = mk_fn2("CoefficientList", expr_copy(p), expr_copy(x));
    Expr* cl = simp(call);
    int deg = -1;
    if (head_is(cl, SYM_List) && cl->data.function.arg_count > 0) {
        /* Walk from the tail back to the first non-zero coefficient.
         * CoefficientList returns {c0, c1, ..., cn}, i.e. ordered by
         * ascending power. */
        for (size_t i = cl->data.function.arg_count; i > 0; i--) {
            if (!is_lit_zero(cl->data.function.args[i - 1])) { deg = (int)(i - 1); break; }
        }
    }
    expr_free(cl);
    return deg;
}

static Expr* poly_leading_coeff(Expr* p, Expr* x) {
    int deg = poly_degree_in(p, x);
    if (deg < 0) return NULL;
    Expr* call = mk_fn2("CoefficientList", expr_copy(p), expr_copy(x));
    Expr* cl = simp(call);
    Expr* out = NULL;
    if (head_is(cl, SYM_List) && (int)cl->data.function.arg_count > deg) {
        out = expr_copy(cl->data.function.args[deg]);
    }
    expr_free(cl);
    return out;
}

static bool is_polynomial_in(Expr* p, Expr* x) {
    Expr* call = mk_fn2("PolynomialQ", expr_copy(p), expr_copy(x));
    Expr* r = simp(call);
    bool ok = is_sym(r, "True");
    expr_free(r);
    return ok;
}

static Expr* layer3_rational(Expr* f, LimitCtx* ctx) {
    /* Together-normalize and extract numerator/denominator. */
    Expr* together = simp(mk_fn1("Together", expr_copy(f)));
    Expr* num = simp(mk_fn1("Numerator",   expr_copy(together)));
    Expr* den = simp(mk_fn1("Denominator", expr_copy(together)));
    expr_free(together);

    if (!is_polynomial_in(num, ctx->x) || !is_polynomial_in(den, ctx->x)) {
        expr_free(num); expr_free(den);
        return NULL;
    }

    Expr* result = NULL;

    if (is_infinity_sym(ctx->point) || is_neg_infinity(ctx->point)) {
        int dn = poly_degree_in(num, ctx->x);
        int dd = poly_degree_in(den, ctx->x);
        if (dn < 0 || dd < 0) goto done;

        if (dn < dd) {
            result = mk_int(0);
        } else if (dn > dd) {
            Expr* ln = poly_leading_coeff(num, ctx->x);
            Expr* ld = poly_leading_coeff(den, ctx->x);
            if (!ln || !ld) { if (ln) expr_free(ln); if (ld) expr_free(ld); goto done; }
            int sn = literal_sign(ln);
            int sd = literal_sign(ld);
            int parity = 1;
            if (is_neg_infinity(ctx->point)) {
                /* As x -> -infty, sign of x^(dn-dd) depends on parity of dn-dd. */
                parity = ((dn - dd) % 2 == 0) ? 1 : -1;
            }
            if (sn != 0 && sd != 0) {
                result = signed_infinity(sn * sd * parity);
            } else {
                /* Give up on symbolic leading coefficients here. */
                expr_free(ln); expr_free(ld);
                goto done;
            }
            expr_free(ln); expr_free(ld);
        } else {
            /* Equal degree: ratio of leading coefficients. */
            Expr* ln = poly_leading_coeff(num, ctx->x);
            Expr* ld = poly_leading_coeff(den, ctx->x);
            if (!ln || !ld) { if (ln) expr_free(ln); if (ld) expr_free(ld); goto done; }
            result = simp(mk_fn2("Divide", ln, ld));
        }
        goto done;
    }

    /* Finite limit point: substitute and, if we hit 0/0, cancel and
     * recurse. We cap iterations to avoid surprises on pathological
     * inputs; the outer dispatcher falls through to the series layer if
     * we return NULL here. */
    Expr* num_at = subst_eval(num, ctx->x, ctx->point);
    Expr* den_at = subst_eval(den, ctx->x, ctx->point);
    bool num_zero = is_lit_zero(num_at);
    bool den_zero = is_lit_zero(den_at);

    if (!den_zero) {
        result = simp(mk_fn2("Divide", expr_copy(num_at), expr_copy(den_at)));
    } else if (num_zero && den_zero) {
        /* 0/0 -- let the series / L'Hospital layers handle it. */
        result = NULL;
    } else {
        /* 0 in denominator, nonzero numerator: the pole's order parity
         * determines the answer. Defer to the series layer unconditionally
         * -- it inspects the leading exponent and produces Infinity for an
         * even-order pole (e.g. 1/(x-2)^2 -> +Infinity) and ComplexInfinity
         * for an odd-order pole with a two-sided direction. */
        (void)num_at;
        result = NULL;
    }
    expr_free(num_at); expr_free(den_at);

done:
    expr_free(num);
    expr_free(den);
    return result;
}

/* ---------------------------------------------------------------------- */
/* Layer 5 -- L'Hospital + logarithmic reduction                           */
/* ---------------------------------------------------------------------- */

/* Apply L'Hospital's rule repeatedly. We only engage this layer when
 * the pointwise evaluation is a recognised indeterminate form (0/0 or
 * Indeterminate). Strict growth of the leaf count in three consecutive
 * iterations is treated as "not making progress". */
static Expr* layer5_lhospital(Expr* f, LimitCtx* ctx) {
    /* Require f = N/D structurally via Together + Numerator/Denominator. */
    Expr* together = simp(mk_fn1("Together", expr_copy(f)));
    Expr* num = simp(mk_fn1("Numerator",   expr_copy(together)));
    Expr* den = simp(mk_fn1("Denominator", expr_copy(together)));
    expr_free(together);

    /* Denominator must not be constant 1 -- that means `f` isn't written
     * as a quotient and L'Hospital doesn't directly apply. */
    if (free_of(den, ctx->x)) {
        expr_free(num); expr_free(den);
        return NULL;
    }

    int64_t last_leaves = leaf_count(num) + leaf_count(den);
    int grown = 0;

    Expr* result = NULL;
    for (int iter = 0; iter < LIMIT_LH_MAX_ITERATIONS; iter++) {
        Expr* n_at = subst_eval(num, ctx->x, ctx->point);
        Expr* d_at = subst_eval(den, ctx->x, ctx->point);
        bool n_zero = is_lit_zero(n_at);
        bool d_zero = is_lit_zero(d_at);
        /* Infinity^2, 3*Infinity etc. all count as "diverging to infinity"
         * for the purposes of classifying this as an inf/inf indeterminate
         * form -- sign / scaling doesn't matter here since we're about to
         * differentiate. */
        bool n_inf  = is_divergent(n_at) && !n_zero;
        bool d_inf  = is_divergent(d_at) && !d_zero;

        if (!n_zero && !d_zero && !n_inf && !d_inf) {
            /* Quotient determinate (and clean); compute and return. */
            result = simp(mk_fn2("Divide", expr_copy(n_at), expr_copy(d_at)));
            expr_free(n_at); expr_free(d_at);
            break;
        }
        expr_free(n_at); expr_free(d_at);

        if (!((n_zero && d_zero) || (n_inf && d_inf))) {
            /* 0/finite or finite/0 -- not L'Hospital territory. */
            break;
        }

        /* Differentiate both. */
        Expr* dn = simp(mk_fn2("D", num, expr_copy(ctx->x)));
        Expr* dd = simp(mk_fn2("D", den, expr_copy(ctx->x)));
        num = dn;
        den = dd;

        int64_t leaves = leaf_count(num) + leaf_count(den);
        if (leaves > last_leaves) {
            grown++;
            if (grown >= LIMIT_LH_MAX_GROWTH) break;
        } else {
            grown = 0;
        }
        last_leaves = leaves;
    }

    expr_free(num);
    expr_free(den);
    return result;
}

/* Logarithmic reduction for exponential indeterminate forms
 *     f = f1^g1 * f2^g2 * ...    ->    Exp[Limit[Sum[g_i Log[f_i]]]]
 *
 * Fires whenever any Power factor in f has the limit variable in its
 * exponent (not just a bare Power at the top level -- that covered the
 * (1+a/x)^(bx) family but missed shapes like ((2+3x)/x)^x / 2^x where
 * the indeterminate behaviour comes from a *product* of Powers).
 *
 * We construct `log_expr = Sum[g_i Log[f_i]]` by walking the Times args,
 * reducing each Power[f_i, g_i] (or bare factor) to its log contribution.
 * Other factors are handled via a catch-all `Log[factor]` term -- the
 * evaluator usually simplifies these to closed form when they are free
 * of x.                                                                  */
static Expr* log_contribution(Expr* factor) {
    if (head_is(factor, SYM_Power) && factor->data.function.arg_count == 2) {
        /* Log[a^b] = b Log[a]. */
        return mk_times(expr_copy(factor->data.function.args[1]),
                        mk_fn1("Log", expr_copy(factor->data.function.args[0])));
    }
    return mk_fn1("Log", expr_copy(factor));
}

/* Post-process a log-limit into a finite/zero/Infinity answer.
 * Returns NULL when the value is ambiguous (ComplexInfinity, Indeterminate). */
static Expr* exp_of_limit(Expr* lim_log) {
    if (is_infinity_sym(lim_log))   return mk_sym("Infinity");
    if (is_neg_infinity(lim_log))   return mk_int(0);
    if (is_complex_infinity(lim_log) || is_indeterminate(lim_log)) return NULL;
    return simp(mk_fn1("Exp", expr_copy(lim_log)));
}

static Expr* layer5_log_reduction(Expr* f, LimitCtx* ctx) {
    if (!has_var_in_exponent(f, ctx->x)) return NULL;

    /* Only fire when f is a "pure power product": a top-level Power, or
     * a top-level Times whose factors are each a constant or a Power.
     * This keeps us away from expressions like (Cos[x]-1)/(Exp[x^2]-1)
     * where `has_var_in_exponent` is true only because of an Exp buried
     * inside a Plus. Log-reducing those shapes would split a determinate
     * 0/0 into two divergent Log terms and trigger infinite Series
     * escalation. */
    bool top_ok = false;
    if (head_is(f, SYM_Power) && f->data.function.arg_count == 2) {
        top_ok = true;
    } else if (head_is(f, SYM_Times)) {
        top_ok = true;
        for (size_t i = 0; i < f->data.function.arg_count; i++) {
            Expr* a = f->data.function.args[i];
            bool is_factor_power = head_is(a, SYM_Power) &&
                                   a->data.function.arg_count == 2;
            if (!is_factor_power && !free_of(a, ctx->x)) { top_ok = false; break; }
        }
    }
    if (!top_ok) return NULL;

    /* Build log_expr = Sum[log_contribution(factor)]. If f itself is a
     * Power head we treat it as a single-factor product. */
    Expr* log_expr = NULL;
    if (head_is(f, SYM_Times)) {
        log_expr = mk_int(0);
        for (size_t i = 0; i < f->data.function.arg_count; i++) {
            log_expr = mk_fn2("Plus", log_expr,
                              log_contribution(f->data.function.args[i]));
        }
    } else {
        log_expr = log_contribution(f);
    }
    Expr* log_s = simp(log_expr);

    LimitCtx sub = *ctx; sub.depth += 1;
    Expr* lim_log = compute_limit(log_s, &sub);
    expr_free(log_s);

    if (!lim_log) return NULL;
    if (expr_contains(lim_log, ctx->x)) { expr_free(lim_log); return NULL; }

    Expr* out = exp_of_limit(lim_log);
    expr_free(lim_log);
    return out;
}

/* ---------------------------------------------------------------------- */
/* Layer 6 -- Bound analysis (Interval returns)                            */
/*                                                                         */
/* We cover one concrete case for now: a bounded head (Sin, Cos) of an    */
/* expression that diverges in the limit -- e.g. Sin[1/x] at x=0. This     */
/* is the main shape exercised by the test suite.                          */
/* ---------------------------------------------------------------------- */
static Expr* layer6_bounded(Expr* f, LimitCtx* ctx) {
    if (!head_is(f, SYM_Sin) && !head_is(f, SYM_Cos)) return NULL;
    if (f->data.function.arg_count != 1) return NULL;

    /* Check that the inner argument has no limit at this point (diverges).
     * The simplest positive check: inner argument has ComplexInfinity or
     * Infinity-like behaviour, or we failed to compute a limit on it. */
    Expr* inner = f->data.function.args[0];
    if (free_of(inner, ctx->x)) return NULL; /* handled elsewhere */

    /* An unbounded-oscillation at the limit point does not converge. The
     * earlier Interval[{-1,1}] return was useful as a bound but it wasn't
     * a *limit* -- Mathematica returns Indeterminate in these shapes
     * (Sin[1/x] at 0, Sin[x] at Infinity). Match that behaviour: the caller
     * only reaches this layer after the squeeze / substitution / Series
     * paths have already failed, so we know the oscillation is real. */
    return mk_sym("Indeterminate");
}

/* Recursively true iff any node of `e` is a function application whose head is
 * the interned symbol `head_sym`. Used to decide whether `head[Infinity]`
 * actually reduced (head gone) or stayed unevaluated. */
static bool contains_head_symbol(Expr* e, const char* head_sym) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == head_sym) return true;
    if (contains_head_symbol(e->data.function.head, head_sym)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_head_symbol(e->data.function.args[i], head_sym)) return true;
    return false;
}

/* ---------------------------------------------------------------------- */
/* Layer -- f[g(x)] at +/-Infinity via the head's own value at Infinity     */
/*                                                                         */
/* Generalises the old ArcTan/ArcCot special case. For a single-argument   */
/* f[inner] whose inner argument diverges to +/-Infinity, apply f to that   */
/* infinite limit and let the builtin fold it:                             */
/*     Erf[Infinity] = 1,   Tanh[Infinity] = 1,   ArcTan[Infinity] = Pi/2,  */
/*     Exp[Infinity] = Infinity,   Gamma[Infinity] = Infinity,   ...        */
/* This is exact for every f that self-evaluates at Infinity -- precisely   */
/* the functions that have a well-defined value there. Series cannot see    */
/* through such heads at an infinite inner argument, so we intercept early. */
/* Oscillatory heads (Sin, Cos) leave f[Infinity] unevaluated, so the       */
/* residual-head guard rejects them and we fall through to the bounded/     */
/* Interval layers (which return Indeterminate). Restricting to +/-Infinity  */
/* (not ComplexInfinity) avoids direction-ambiguous folds.                  */
/* ---------------------------------------------------------------------- */
static Expr* layer_compose_at_infinity(Expr* f, LimitCtx* ctx) {
    if (f->type != EXPR_FUNCTION || f->data.function.arg_count != 1) return NULL;
    Expr* head = f->data.function.head;
    if (!head || head->type != EXPR_SYMBOL) return NULL;
    Expr* inner = f->data.function.args[0];
    if (free_of(inner, ctx->x)) return NULL;

    LimitCtx sub = *ctx; sub.depth += 1;
    Expr* lim_inner = compute_limit(inner, &sub);
    if (!lim_inner) return NULL;
    if (expr_contains(lim_inner, ctx->x)) { expr_free(lim_inner); return NULL; }

    bool inner_real_inf = is_infinity_sym(lim_inner) || is_neg_infinity(lim_inner);
    /* A directionless (ComplexInfinity) or directed (DirectedInfinity[d])
     * inner limit is also foldable when the head has a definite value there,
     * e.g. ArcTan[DirectedInfinity[I Sqrt[2]]] = Pi/2. */
    bool inner_dir_inf = is_complex_infinity(lim_inner) || is_directed_infinity(lim_inner);
    if (!inner_real_inf && !inner_dir_inf) {
        expr_free(lim_inner);
        return NULL;
    }

    /* Apply the head to the infinite inner limit and let the builtin evaluate
     * it. mk_fn1 adopts lim_inner (via expr_new_function). */
    Expr* val = simp(mk_fn1(head->data.symbol.name, lim_inner));

    /* Accept only when the head actually resolved (no residual head[...] and no
     * pending x). Rejects Sin[Infinity], Cos[Infinity], ... that stay symbolic. */
    bool resolved = val && !contains_head_symbol(val, head->data.symbol.name) &&
                    !expr_contains(val, ctx->x);
    /* For a directionless / complex inner infinity, additionally require an
     * unambiguous value: reject Indeterminate and any residual infinity
     * (ArcTan[ComplexInfinity] = Indeterminate must NOT be accepted as a
     * limit). A clean finite value or a real signed Infinity is genuine. */
    if (resolved && inner_dir_inf &&
        (is_indeterminate(val) || is_complex_infinity(val) ||
         is_directed_infinity(val))) {
        resolved = false;
    }
    if (resolved) return val;
    expr_free(val);
    return NULL;
}

/* ---------------------------------------------------------------------- */
/* Layer -- Bounded envelope (squeeze theorem)                             */
/*                                                                         */
/* Try to show |f| -> 0 by bounding |f| pointwise with an expression       */
/* whose limit is easier to compute. Fires when f actually contains a      */
/* bounded oscillatory head; otherwise no value is added over the other    */
/* layers. Returns 0 on success, NULL otherwise. This gives us the         */
/* Sin[h(x)]/p(x), (1 +/- Cos[x])/x, x Sin[x]/(5 + x^2) family at infinity */
/* as well as x^2 Sin[1/x] at 0.                                           */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Bounded base raised to a divergent positive power                       */
/*                                                                         */
/* Handles f = base^exp where `base` carries a bounded oscillatory head    */
/* (Sin/Cos/Tanh/ArcTan/ArcCot) and exp -> +Infinity at the limit point.   */
/* For real exp, |base^exp| = |base|^exp, and because exp is eventually     */
/* positive the magnitude is bounded above by B^exp, with B the *constant*  */
/* magnitude bound of the base. When B^exp -> 0 the squeeze forces f -> 0   */
/* (the value is complex where base < 0, but its magnitude still vanishes,  */
/* so the two-sided complex limit is 0). Unlike the general envelope this   */
/* uses no bare-variable magnitude reasoning, so it is valid at a finite    */
/* point as well as at Infinity.                                            */
/* Example: (Sin[1/x]/2)^(1/x^2) at x -> 0  =>  <= (1/2)^(1/x^2) -> 0.      */
/* ---------------------------------------------------------------------- */
/* True iff `e` is a concrete real number in [0, 1). Symbolic values fail
 * the comparison (stay unevaluated, not True), so only comparable literals
 * -- Integer, Real, Rational -- can satisfy it. */
static bool lit_in_unit_interval(Expr* e) {
    if (!e) return false;
    Expr* lo = simp(mk_fn2("GreaterEqual", expr_copy(e), mk_int(0)));
    Expr* hi = simp(mk_fn2("Less",         expr_copy(e), mk_int(1)));
    bool ok = is_sym(lo, "True") && is_sym(hi, "True");
    expr_free(lo); expr_free(hi);
    return ok;
}

/* True iff `e` is a concrete real number strictly greater than 1. */
static bool lit_gt_one(Expr* e) {
    if (!e) return false;
    Expr* hi = simp(mk_fn2("Greater", expr_copy(e), mk_int(1)));
    bool ok = is_sym(hi, "True");
    expr_free(hi);
    return ok;
}

static Expr* envelope_bounded_power(Expr* f, LimitCtx* ctx) {
    if (!head_is(f, SYM_Power) || f->data.function.arg_count != 2) return NULL;
    Expr* base = f->data.function.args[0];
    Expr* exp  = f->data.function.args[1];

    /* The base must be the oscillatory factor; the exponent must vary. */
    if (!contains_bounded_head(base)) return NULL;
    if (free_of(exp, ctx->x)) return NULL;

    /* Require exp -> +Infinity. This makes exp eventually positive (so the
     * squeeze direction is fixed): B^exp -> 0 forces base^exp -> 0, and a
     * base bounded below by a constant > 1 forces base^exp -> Infinity. A
     * -Infinity or finite exponent limit is rejected: the inequalities flip
     * and neither conclusion holds. */
    LimitCtx sub = *ctx; sub.depth += 1;
    Expr* lim_exp = compute_limit(exp, &sub);
    if (!lim_exp) return NULL;
    bool pos_inf = is_infinity_sym(lim_exp);
    expr_free(lim_exp);
    if (!pos_inf) return NULL;

    /* --- Upper-bound squeeze to 0 ---------------------------------------
     * B(x) >= |base(x)| pointwise (var_abs makes the bare-variable bound
     * sound at a finite point, e.g. |x Sin[1/x]/2| <= Abs[x]/2). When
     * Limit[B] = L with 0 <= L < 1 then eventually B < r < 1, so
     * |base^exp| = |base|^exp <= B^exp <= r^exp -> 0 (exp -> +Infinity).
     * The two-sided complex limit is then 0 (magnitude vanishes even where
     * base < 0). Firing on Limit[B] rather than Limit[B^exp] lets a
     * *shrinking* x-dependent bound (Limit[B] = 0) succeed, not just a
     * constant one. */
    Expr* B = magnitude_upper_bound(base, ctx->x, /*var_abs=*/true);
    if (B) {
        Expr* lim_B = compute_limit(B, &sub);
        expr_free(B);
        if (lim_B) {
            bool decays = lit_in_unit_interval(lim_B);
            expr_free(lim_B);
            if (decays) return mk_int(0);
        }
    }

    /* --- Lower-bound blow-up to Infinity --------------------------------
     * For base = c + g(x) with c an x-free real constant and g the bounded
     * oscillatory remainder, |base| >= |c| - U where U >= |g|. If
     * |c| - U -> L > 1 and c > 0, then base is eventually positive and
     * base >= L > 1, so base^exp -> +Infinity (exp -> +Infinity). */
    if (head_is(base, SYM_Plus)) {
        Expr* c    = mk_int(0);   /* sum of x-free terms   */
        Expr* rest = mk_int(0);   /* sum of x-bearing terms */
        for (size_t i = 0; i < base->data.function.arg_count; i++) {
            Expr* t = base->data.function.args[i];
            if (free_of(t, ctx->x))
                c    = mk_fn2("Plus", c,    expr_copy(t));
            else
                rest = mk_fn2("Plus", rest, expr_copy(t));
        }
        c = simp(c); rest = simp(rest);

        if (literal_sign(c) > 0) {
            Expr* U     = magnitude_upper_bound(rest, ctx->x, /*var_abs=*/true);
            Expr* lower = simp(mk_fn2("Plus", mk_fn1("Abs", expr_copy(c)),
                                              mk_neg(U)));
            Expr* lim_L = compute_limit(lower, &sub);
            expr_free(lower);
            bool blows_up = lim_L && lit_gt_one(lim_L);
            if (lim_L) expr_free(lim_L);
            if (blows_up) { expr_free(c); expr_free(rest); return mk_sym("Infinity"); }
        }
        expr_free(c); expr_free(rest);
    }

    return NULL;
}

/* ---------------------------------------------------------------------- */
static Expr* layer_bounded_envelope(Expr* f, LimitCtx* ctx) {
    if (!contains_bounded_head(f)) return NULL;

    /* Bounded base raised to a divergent positive power. Valid at finite
     * points too, so it runs before the +Infinity-only general squeeze. */
    Expr* bp = envelope_bounded_power(f, ctx);
    if (bp) return bp;

    /* General squeeze is +Infinity-only: magnitude_upper_bound returns `x`
     * (not Abs[x]) for the bare variable, which is only an actual upper
     * bound on |x| when x is positive. Negative-infinity limits must first
     * be transformed via x -> -y. */
    if (!is_infinity_sym(ctx->point)) return NULL;

    Expr* bound = magnitude_upper_bound(f, ctx->x, /*var_abs=*/false);
    if (!bound) return NULL;

    /* Guard: if the bound itself still mentions a bounded head, the
     * magnitude reduction failed to replace the oscillator by its constant
     * envelope -- it survived buried inside an Abs[Log[..]] / Abs[E^(..)]
     * default wrap. Recursing compute_limit on such a bound re-enters this
     * layer and re-wraps another Abs each level down to the depth cap (an
     * expensive simp per level), so bail: the squeeze cannot make progress.
     * This is what makes Limit[x Log[Sin[x]/2], x->Infinity] terminate,
     * and thereby unblocks the (Sin[x]/2)^x / (1+Sin[x]/x)^(x^2) families. */
    if (contains_bounded_head(bound)) { expr_free(bound); return NULL; }

    LimitCtx sub = *ctx; sub.depth += 1;
    Expr* lim_bound = compute_limit(bound, &sub);
    expr_free(bound);
    if (!lim_bound) return NULL;

    bool zero = is_lit_zero(lim_bound);
    expr_free(lim_bound);
    if (zero) return mk_int(0);
    return NULL;
}

/* ---------------------------------------------------------------------- */
/* Layer -- Log / polynomial merge at infinity                             */
/*                                                                         */
/* Rewrites `Sum(Log[g_i]) + Sum(h_j)` into a single `Log[∏ g_i * ∏ Exp[h_j]]`  */
/* and re-computes the limit. Meant for shapes like `-x + Log[2 + E^x]` at  */
/* infinity, where the individual summands diverge but after Expand the    */
/* combined argument has a finite positive limit (here 1).                 */
/* ---------------------------------------------------------------------- */
static Expr* layer_log_merge(Expr* f, LimitCtx* ctx) {
    if (!is_infinity_sym(ctx->point) && !is_neg_infinity(ctx->point)) return NULL;
    if (!head_is(f, SYM_Plus)) return NULL;
    size_t n = f->data.function.arg_count;

    bool has_x_log = false;
    for (size_t i = 0; i < n; i++) {
        Expr* t = f->data.function.args[i];
        if (head_is(t, SYM_Log) && t->data.function.arg_count == 1 &&
            !free_of(t, ctx->x)) { has_x_log = true; break; }
    }
    if (!has_x_log) return NULL;

    Expr** factors = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* t = f->data.function.args[i];
        if (head_is(t, SYM_Log) && t->data.function.arg_count == 1) {
            factors[i] = expr_copy(t->data.function.args[0]);
        } else {
            factors[i] = mk_fn1("Exp", expr_copy(t));
        }
    }
    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times), factors, n);
    free(factors);
    Expr* prod_s = simp(prod);
    Expr* prod_expanded = simp(mk_fn1("Expand", prod_s));
    Expr* new_f = simp(mk_fn1("Log", prod_expanded));

    /* Require a structural change to avoid unbounded recursion. */
    if (expr_eq(new_f, f)) { expr_free(new_f); return NULL; }

    LimitCtx sub = *ctx; sub.depth += 1;
    Expr* r = compute_limit(new_f, &sub);
    expr_free(new_f);
    return r;
}

/* ---------------------------------------------------------------------- */
/* Gruntz-lite: Log of a sum where one summand dominates                   */
/*                                                                         */
/* For `Log[a + b + ...]` at +Infinity where one summand dominates the    */
/* rest (strictly larger growth_exponent), rewrite as                      */
/*   Log[dom * (1 + rest/dom)] = Log[dom] + Log[1 + rest/dom]              */
/* and recurse. The `rest/dom` terms go to 0 so `Log[1 + small]` folds to  */
/* its Taylor series `small + ...`. Handles the iterated-log shapes that  */
/* the generic Series / L'Hospital layers can't peel.                     */
/* ---------------------------------------------------------------------- */
static Expr* layer_log_sum_gruntz(Expr* f, LimitCtx* ctx) {
    if (!is_infinity_sym(ctx->point)) return NULL;
    if (!head_is(f, SYM_Log) || f->data.function.arg_count != 1) return NULL;
    Expr* inner = f->data.function.args[0];
    if (!head_is(inner, SYM_Plus)) return NULL;
    size_t n = inner->data.function.arg_count;
    if (n < 2) return NULL;

    /* Find the unique dominant summand by growth_exponent. */
    int64_t max_g = 0;
    int max_idx = -1;
    int max_count = 0;
    bool unknown_seen = false;
    for (size_t i = 0; i < n; i++) {
        int64_t g = growth_exponent_upper(inner->data.function.args[i], ctx->x);
        if (g == LIMIT_UNKNOWN_GROWTH) { unknown_seen = true; break; }
        if (max_idx < 0 || g > max_g) { max_g = g; max_idx = (int)i; max_count = 1; }
        else if (g == max_g) max_count++;
    }
    if (unknown_seen || max_count != 1 || max_g <= 0) return NULL;

    /* Build rest = Plus of the other summands; then the rewrite
     * Log[dom + rest] = Log[dom] + Log[1 + rest/dom]. Recurse via the
     * outer dispatcher so Log-of-finite, Series, etc. can fold the
     * Log[1 + rest/dom] term (which goes to 0). */
    Expr* dom = expr_copy(inner->data.function.args[max_idx]);
    size_t rc = n - 1;
    Expr** rest_args = calloc(rc, sizeof(Expr*));
    size_t k = 0;
    for (size_t i = 0; i < n; i++) {
        if ((int)i == max_idx) continue;
        rest_args[k++] = expr_copy(inner->data.function.args[i]);
    }
    Expr* rest;
    if (rc == 1) { rest = rest_args[0]; free(rest_args); }
    else { rest = expr_new_function(mk_sym("Plus"), rest_args, rc); free(rest_args); }

    Expr* ratio = simp(mk_fn2("Times", rest,
                              mk_fn2("Power", expr_copy(dom), mk_int(-1))));
    Expr* log1p = mk_fn1("Log", mk_fn2("Plus", mk_int(1), ratio));
    Expr* log_dom = mk_fn1("Log", dom);
    Expr* rewritten = simp(mk_fn2("Plus", log_dom, log1p));

    LimitCtx sub = *ctx; sub.depth += 1;
    Expr* r = compute_limit(rewritten, &sub);
    expr_free(rewritten);
    return r;
}

/* Higher-level Gruntz shortcut: when the OUTER expression is a quotient
 * and both numerator and denominator have top-level Log[sum] shapes at
 * +Infinity, running `layer_log_sum_gruntz` on each side and letting
 * the termwise layer pick up the asymptotic cancellation produces the
 * iterated-log limits. This is called from the top dispatcher before
 * the generic Series layer. */
static Expr* layer_gruntz_iterated_log(Expr* f, LimitCtx* ctx) {
    if (!is_infinity_sym(ctx->point)) return NULL;
    if (!head_is(f, SYM_Times) && !head_is(f, SYM_Log)) return NULL;

    /* Walk f and apply the Log[sum] rewrite to every eligible Log[sum]
     * subexpression. ReplaceAll-style traversal driven by
     * layer_log_sum_gruntz can't be expressed directly, so we emulate
     * a one-level rewrite by computing the limit of the transformed
     * expression: if any rewrite fires and the result is determinate,
     * we win; otherwise we bail. */
    /* Simpler: only attempt when f has the specific shape Log[sum]
     * multiplied by something free-of-x, or Log[sum] directly. */
    if (head_is(f, SYM_Log)) {
        Expr* r = layer_log_sum_gruntz(f, ctx);
        if (r) return r;
    }
    return NULL;
}

/* ---------------------------------------------------------------------- */
/* Layer -- Log of an expression with a finite non-zero limit              */
/*                                                                         */
/* If f = Log[g] and Limit[g] is a finite value c, return Log[c]. The      */
/* evaluator will produce the correct value for any real c (including      */
/* -Infinity for c = 0). Refuses divergent inner limits, which are picked  */
/* up elsewhere.                                                           */
/* ---------------------------------------------------------------------- */
static Expr* layer_log_of_finite(Expr* f, LimitCtx* ctx) {
    if (!head_is(f, SYM_Log) || f->data.function.arg_count != 1) return NULL;
    Expr* g = f->data.function.args[0];
    if (free_of(g, ctx->x)) return NULL;

    LimitCtx sub = *ctx; sub.depth += 1;
    Expr* lim_g = compute_limit(g, &sub);
    if (!lim_g) return NULL;
    if (is_divergent(lim_g) || expr_contains(lim_g, ctx->x)) {
        expr_free(lim_g); return NULL;
    }
    return simp(mk_fn1("Log", lim_g));
}

/* Upper bound on the growth exponent of `e` viewed as a function of x
 * as x -> Infinity. Returns 0 for constants and bounded heads
 * (Sin, Cos, Tanh, ArcTan, ArcCot of anything), 1 for x, and adds /
 * multiplies through Plus/Times/Power. Returns INT64_MAX to signal "I
 * don't know" -- the caller must then refuse. Used by the dominant-term
 * classifier for Plus at infinity: a term with strictly larger growth
 * than everything else drives the sum to its own limit. */
static int64_t growth_exponent_upper(Expr* e, Expr* x) {
    if (!e) return 0;
    if (free_of(e, x)) return 0;
    if (expr_eq(e, x)) return 1;
    if (e->type != EXPR_FUNCTION) return LIMIT_UNKNOWN_GROWTH;

    /* Bounded heads. */
    if (head_is(e, SYM_Sin) || head_is(e, SYM_Cos) || head_is(e, SYM_Tanh) ||
        head_is(e, SYM_ArcTan) || head_is(e, SYM_ArcCot)) {
        return 0;
    }

    if (head_is(e, SYM_Plus)) {
        int64_t m = 0;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            int64_t g = growth_exponent_upper(e->data.function.args[i], x);
            if (g == LIMIT_UNKNOWN_GROWTH) return LIMIT_UNKNOWN_GROWTH;
            if (g > m) m = g;
        }
        return m;
    }
    if (head_is(e, SYM_Times)) {
        int64_t total = 0;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            int64_t g = growth_exponent_upper(e->data.function.args[i], x);
            if (g == LIMIT_UNKNOWN_GROWTH) return LIMIT_UNKNOWN_GROWTH;
            total += g;
        }
        return total;
    }
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* b = e->data.function.args[0];
        Expr* p = e->data.function.args[1];
        if (p->type == EXPR_INTEGER && p->data.integer >= 0) {
            int64_t g = growth_exponent_upper(b, x);
            if (g == LIMIT_UNKNOWN_GROWTH) return LIMIT_UNKNOWN_GROWTH;
            return g * p->data.integer;
        }
        if (p->type == EXPR_INTEGER && p->data.integer < 0) {
            /* Negative-power term: treat as bounded (<= const) when the
             * base diverges. Conservative: if base is polynomial in x
             * with positive growth, then 1/base -> 0, safely < any
             * diverging term. */
            int64_t g = growth_exponent_upper(b, x);
            if (g > 0 && g != LIMIT_UNKNOWN_GROWTH) return 0;
        }
    }

    /* Log at infinity is sub-polynomial -- treat as 0 for ordering but
     * we don't know if it ever wins against a positive polynomial
     * degree. Conservative: return 0 so Log[x] alongside x is clearly
     * dominated by x. */
    if (head_is(e, SYM_Log) && e->data.function.arg_count == 1) {
        int64_t g = growth_exponent_upper(e->data.function.args[0], x);
        if (g > 0 && g != LIMIT_UNKNOWN_GROWTH) return 0;
    }

    return LIMIT_UNKNOWN_GROWTH;
}

/* True iff `e` is structurally bounded over the reals: every reachable
 * sub-expression's magnitude is dominated by a constant, via the family
 * |Sin|, |Cos|, |Tanh|, |ArcTan|, |ArcCot| <= const, products of
 * bounded things, sums of bounded things, constants. A non-exhaustive
 * but conservative check -- used by the dominant-term Plus layer to
 * allow a clean Infinity answer when one term diverges and the others
 * are merely bounded (not convergent). */
static bool is_structurally_bounded(Expr* e, Expr* x) {
    if (!e) return true;
    if (free_of(e, x)) return true;   /* constant in x */
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL ||
        e->type == EXPR_BIGINT || e->type == EXPR_STRING) return true;
    if (e->type == EXPR_SYMBOL) return !expr_eq(e, x);
    if (e->type != EXPR_FUNCTION) return false;

    /* Heads that are bounded regardless of argument. */
    if (head_is(e, SYM_Sin) || head_is(e, SYM_Cos) || head_is(e, SYM_Tanh) ||
        head_is(e, SYM_ArcTan) || head_is(e, SYM_ArcCot)) {
        return true;
    }

    /* Plus / Times of bounded things is bounded. */
    if (head_is(e, SYM_Plus) || head_is(e, SYM_Times)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (!is_structurally_bounded(e->data.function.args[i], x)) return false;
        }
        return true;
    }
    /* Power[bounded, constant_non_negative_int]. */
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* b = e->data.function.args[0];
        Expr* p = e->data.function.args[1];
        if (p->type == EXPR_INTEGER && p->data.integer >= 0) {
            return is_structurally_bounded(b, x);
        }
    }
    return false;
}

/* ---------------------------------------------------------------------- */
/* Layer -- Term-wise sum at infinity                                      */
/*                                                                         */
/* For a Plus at ±Infinity, compute each term's limit.                     */
/*   - If every term has a finite limit, return the sum.                   */
/*   - If exactly one term has a +Infinity / -Infinity limit and the       */
/*     remaining terms are structurally bounded (Sin/Cos/Tanh/ArcTan       */
/*     applied to anything) or have finite limits, the dominant term       */
/*     wins: return its limit. This is the `x^2 + x Sin[x^2]` -> Infinity */
/*     and `x + Sin[x]` -> Infinity family.                                */
/*   - Otherwise refuse (NULL).                                            */
/* ---------------------------------------------------------------------- */
static Expr* layer_plus_termwise(Expr* f, LimitCtx* ctx) {
    if (!head_is(f, SYM_Plus)) return NULL;
    bool at_infinity = is_infinity_sym(ctx->point) || is_neg_infinity(ctx->point);
    size_t n = f->data.function.arg_count;
    if (n == 0) return NULL;

    /* First, try the growth-exponent dominant-term shortcut. If one
     * term has strictly larger growth than every other term and its
     * individual limit is +/- Infinity, return that -- every other
     * term is at most O(x^(lower)) which the dominant absorbs. This
     * reasoning is asymptotic, so it only applies at infinity. */
    if (at_infinity) {
    int64_t max_g = 0;
    int max_idx = -1;
    int max_count = 0;
    bool unknown_seen = false;
    for (size_t i = 0; i < n; i++) {
        int64_t g = growth_exponent_upper(f->data.function.args[i], ctx->x);
        if (g == LIMIT_UNKNOWN_GROWTH) { unknown_seen = true; break; }
        if (max_idx < 0 || g > max_g) { max_g = g; max_idx = (int)i; max_count = 1; }
        else if (g == max_g) max_count++;
    }
    if (!unknown_seen && max_idx >= 0 && max_count == 1 && max_g > 0) {
        LimitCtx sub = *ctx; sub.depth += 1;
        Expr* lim_dom = compute_limit(f->data.function.args[max_idx], &sub);
        if (lim_dom && (is_infinity_sym(lim_dom) || is_neg_infinity(lim_dom))) {
            return lim_dom;
        }
        if (lim_dom) expr_free(lim_dom);
    }
    } /* end at_infinity growth-exponent shortcut */

    Expr** terms = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) terms[i] = NULL;

    int dominant_count = 0;
    int dominant_idx = -1;
    Expr* dominant_lim = NULL;

    for (size_t i = 0; i < n; i++) {
        Expr* t = f->data.function.args[i];
        if (free_of(t, ctx->x)) {
            terms[i] = expr_copy(t);
            continue;
        }
        LimitCtx sub = *ctx; sub.depth += 1;
        Expr* lim_t = compute_limit(t, &sub);
        /* Three outcomes for the term: finite (keep), ±Infinity
         * (potential dominant), or bounded-but-unresolved (admissible
         * alongside a single dominant term). */
        bool is_pos_inf = lim_t && is_infinity_sym(lim_t);
        bool is_neg_inf = lim_t && is_neg_infinity(lim_t);
        bool is_signed_inf = is_pos_inf || is_neg_inf;
        bool ok_finite = lim_t && !is_divergent(lim_t) &&
                         !expr_contains(lim_t, ctx->x);
        if (is_signed_inf) {
            dominant_count++;
            dominant_idx = (int)i;
            if (dominant_lim) expr_free(dominant_lim);
            dominant_lim = lim_t;
            continue;
        }
        if (ok_finite) { terms[i] = lim_t; continue; }
        /* Unresolved or non-signed-infinity divergent. Accept as
         * structurally-bounded only if we have a dominant term elsewhere
         * -- we don't know the actual limit value, but its magnitude is
         * capped, so it doesn't change the dominant term's infinity. */
        if (lim_t) expr_free(lim_t);
        if (is_structurally_bounded(t, ctx->x)) {
            /* Placeholder only usable if a dominant term eventually
             * appears. Tag via a negative dominant_count signal. */
            terms[i] = mk_int(0);
            dominant_idx = -1000;  /* flag: saw a placeholder */
            continue;
        }
        /* Can't classify -> bail. */
        for (size_t j = 0; j < n; j++) if (terms[j]) expr_free(terms[j]);
        if (dominant_lim) expr_free(dominant_lim);
        free(terms);
        return NULL;
    }

    if (dominant_count == 1) {
        /* Dominant-term win: ignore the finite / bounded terms. */
        for (size_t j = 0; j < n; j++) if (terms[j]) expr_free(terms[j]);
        free(terms);
        (void)dominant_idx;
        return dominant_lim;
    }
    if (dominant_count > 1) {
        /* Two or more ±Infinity; we'd have to rank them. Not handled
         * here; let later layers try. */
        for (size_t j = 0; j < n; j++) if (terms[j]) expr_free(terms[j]);
        if (dominant_lim) expr_free(dominant_lim);
        free(terms);
        return NULL;
    }
    /* dominant_count == 0: normal termwise sum. But if we had a
     * bounded-unresolved placeholder and no dominant term actually
     * fired, the sum is not valid -- bail. */
    if (dominant_idx == -1000) {
        for (size_t j = 0; j < n; j++) if (terms[j]) expr_free(terms[j]);
        free(terms);
        return NULL;
    }
    /* dominant_count == 0 and no bounded placeholder: every term has a
     * finite limit. The limit of a sum equals the sum of the limits whenever
     * each summand converges -- valid at a finite point as well as at
     * infinity. Any summand that diverges (ComplexInfinity, an unresolved
     * oscillation, ...) has already caused a bail above, so a cancellation
     * shape like 1/x - Cot[x] never reaches here and is left to Series; only
     * genuinely convergent sums are folded. This closes branch-point
     * boundaries such as  -2 Sqrt[r]/(1+r) + 2 ArcTan[Sqrt[r]]  at r -> oo,
     * where each summand converges but Series cannot expand. */
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, n);
    free(terms);
    return simp(sum);
}

/* ---------------------------------------------------------------------- */
/* Layer -- Plus split: convergent terms + jointly-convergent remainder    */
/*                                                                         */
/* layer_plus_termwise fails when two-or-more summands each diverge, even   */
/* if their divergences cancel.  The archetype is the real log-part of a    */
/* rational-function antiderivative,                                        */
/*                                                                         */
/*     ... + a ArcTan[p] + a ArcTan[q]                                      */
/*         - b Log[1 - c x + x^2] + b Log[1 + c x + x^2]                    */
/*                                                                         */
/* at x -> +/-Infinity: the two ArcTan terms converge, but each Log term    */
/* individually -> +/-Infinity so termwise sees two dominant infinities and */
/* bails -- even though the Log pair together -> 0 (their ratio -> 1).      */
/*                                                                         */
/* The rescue rests on the identity                                        */
/*                                                                         */
/*     lim (A + B) = (lim A) + (lim B)   whenever  lim A  exists (finite),  */
/*                                                                         */
/* so it is always valid to (1) sum the summands whose INDIVIDUAL limit is  */
/* a finite value, then (2) take the limit of the leftover group as a       */
/* SINGLE sub-expression and add it.  The leftover group -- handed back to  */
/* compute_limit intact -- reaches layer_log_merge / Series, which cancel   */
/* the mutual divergences that no single term's limit exposes.  Requires at */
/* least one finite ("easy") summand, so the recursion runs on a strictly   */
/* smaller expression (termination), and at least one leftover summand.     */
/* ---------------------------------------------------------------------- */

/* A limit value is directly usable as a finite contribution iff it is a
 * concrete value free of the limit variable, not a divergence sentinel, and
 * not an Interval (an oscillation envelope has no single limit). */
static bool limit_value_is_finite(Expr* v, Expr* x) {
    return v && !is_divergent(v) && !expr_contains(v, x) &&
           !head_is(v, SYM_Interval);
}

static Expr* layer_plus_split_convergent(Expr* f, LimitCtx* ctx) {
    if (!head_is(f, SYM_Plus)) return NULL;
    size_t n = f->data.function.arg_count;
    if (n < 2) return NULL;

    Expr** easy = malloc(sizeof(Expr*) * n);   /* owned finite limits */
    Expr** hard = malloc(sizeof(Expr*) * n);   /* borrowed leftover terms */
    size_t ne = 0, nh = 0;

    for (size_t i = 0; i < n; i++) {
        Expr* t = f->data.function.args[i];
        if (free_of(t, ctx->x)) { easy[ne++] = expr_copy(t); continue; }
        LimitCtx sub = *ctx; sub.depth += 1;
        Expr* lim_t = compute_limit(t, &sub);
        if (limit_value_is_finite(lim_t, ctx->x)) {
            easy[ne++] = lim_t;                 /* finite -> peel off */
        } else {
            if (lim_t) expr_free(lim_t);
            hard[nh++] = t;                     /* borrowed leftover */
        }
    }

    /* Need >=1 finite term (recursion strictly shrinks) and a leftover. */
    if (ne == 0 || nh == 0) {
        for (size_t i = 0; i < ne; i++) expr_free(easy[i]);
        free(easy); free(hard);
        return NULL;
    }

    /* Limit of the leftover group, taken together so cancellations show. */
    Expr* sub_f;
    if (nh == 1) {
        sub_f = expr_copy(hard[0]);
    } else {
        Expr** ha = malloc(sizeof(Expr*) * nh);
        for (size_t i = 0; i < nh; i++) ha[i] = expr_copy(hard[i]);
        sub_f = expr_new_function(expr_new_symbol(SYM_Plus), ha, nh);
        free(ha);
    }
    free(hard);

    LimitCtx sub = *ctx; sub.depth += 1;
    Expr* H = compute_limit(sub_f, &sub);
    expr_free(sub_f);

    if (!H || expr_contains(H, ctx->x)) {
        for (size_t i = 0; i < ne; i++) expr_free(easy[i]);
        free(easy);
        if (H) expr_free(H);
        return NULL;
    }

    /* (sum of finite limits) + H.  simp folds finite+finite, and
     * finite + (+/-Infinity) -> that infinity, so a genuinely divergent
     * leftover still surfaces correctly rather than being masked. */
    Expr** all = malloc(sizeof(Expr*) * (ne + 1));
    for (size_t i = 0; i < ne; i++) all[i] = easy[i];
    all[ne] = H;
    free(easy);
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), all, ne + 1);
    free(all);
    return simp(sum);
}

/* ---------------------------------------------------------------------- */
/* Layer -- Atom substitution for Power-in-x-exponent shapes               */
/*                                                                         */
/* Replaces a uniform Power[b, e(x)] subterm with a fresh symbol u, then   */
/* recurses on Limit[f_sub, u -> Limit[atom, x -> point]]. Works for       */
/* shapes like                                                             */
/*                                                                         */
/*     (-1 + 3^(2/x)) / (1 + 3^(2/x))     at x -> 0 one-sided              */
/*                                                                         */
/* where Series cannot expand around the essential singularity of the      */
/* inner Power, but the expression IS polynomial in the atom itself and   */
/* the atom has a clean one-sided limit (+Infinity for x -> 0+, 0 for     */
/* x -> 0-). The rational-function layer then closes out u -> Infinity.   */
/*                                                                         */
/* The atom must appear in every x-dependent position of f; if x survives */
/* outside the atom we can't reason about f by u alone. This is exactly    */
/* the condition that subst_eval(f, atom, u) leaves no x behind. We only  */
/* try when the original pipeline (Series + L'Hospital) has already       */
/* failed, so the extra cost is amortised over genuinely hard inputs.     */
/* ---------------------------------------------------------------------- */
static Expr* find_mrv_power(Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* exp = e->data.function.args[1];
        /* Require the exponent to depend on x AND the base to NOT depend
         * on x; otherwise the "atom" is actually a moving target and the
         * substitution trick gives the wrong answer. */
        Expr* base = e->data.function.args[0];
        if (expr_contains(exp, x) && free_of(base, x)) return e;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        Expr* r = find_mrv_power(e->data.function.args[i], x);
        if (r) return r;
    }
    return NULL;
}

static Expr* layer_atom_substitute(Expr* f, LimitCtx* ctx) {
    /* Narrow gating: the trigger shape is a Power[constant, e(x)] where
     * e diverges at the point. If no such subtree already exists in f we
     * save an expensive Together call and bail immediately. */
    if (!find_mrv_power(f, ctx->x)) return NULL;

    /* Numeric literal limit points are the working scope. Symbolic points
     * (Limit[..., x -> a]) don't produce a clean "atom -> infinity" shape;
     * attempting Together here can also cost seconds on exponential
     * subexpressions like E^x - E^a. */
    if (!is_numeric_literal_point(ctx->point)) return NULL;

    /* Depth guard: Together on a recursively-substituted f can explode in
     * size. Keep this layer to the top few recursion levels. */
    if (ctx->depth > 3) return NULL;

    /* Need a Together-normalised view: the trigger shape (e.g.
     * (3^(1/x) - 3^(-1/x))/(3^(1/x) + 3^(-1/x))) only exposes a single
     * atom after Together factors out the shared 3^(-1/x) and collapses
     * to (-1 + 3^(2/x))/(1 + 3^(2/x)). */
    Expr* g = simp(mk_fn1("Together", expr_copy(f)));

    Expr* atom = find_mrv_power(g, ctx->x);
    if (!atom) { expr_free(g); return NULL; }

    /* Compute the atom's limit with the current direction. If the atom's
     * limit is divergent *or* the atom itself doesn't resolve, we can't
     * make progress by substituting (a dangling symbol with no value would
     * defeat the whole point). */
    LimitCtx sub = *ctx; sub.depth += 1;
    Expr* atom_copy = expr_copy(atom);   /* atom is borrowed from g */
    Expr* atom_lim = compute_limit(atom_copy, &sub);
    expr_free(atom_copy);
    if (!atom_lim) { expr_free(g); return NULL; }
    if (expr_contains(atom_lim, ctx->x)) { expr_free(atom_lim); expr_free(g); return NULL; }

    /* Substitute atom -> u. We want a symbol that cannot clash with any
     * user variable in f; the `$` prefix convention in Mathilda marks
     * system-local symbols. */
    Expr* u_sym = mk_sym("$LimitAtomU$");
    Expr* f_sub = subst_eval(g, atom, u_sym);
    expr_free(g);

    /* If x survives outside the atom, the atom trick is invalid. */
    if (expr_contains(f_sub, ctx->x)) {
        expr_free(f_sub); expr_free(u_sym); expr_free(atom_lim);
        return NULL;
    }

    /* Run the u-limit. Use TWOSIDED inside -- u takes a one-sided value
     * (e.g. Infinity or 0) driven by the direction on x, but once we're
     * asking for Limit in u the approach is along the real line of u
     * itself. */
    LimitCtx sub2 = { u_sym, atom_lim, LIMIT_DIR_TWOSIDED, ctx->depth, ctx->method };
    Expr* result = compute_limit(f_sub, &sub2);
    expr_free(f_sub); expr_free(u_sym); expr_free(atom_lim);
    return result;
}

/* ---------------------------------------------------------------------- */
/* Layer -- Two-sided disagreement probe                                   */
/*                                                                         */
/* When a two-sided limit at a finite numeric point has reached the end of */
/* the pipeline without resolving, compute the one-sided limits from both  */
/* sides and compare. This catches shapes like                             */
/*                                                                         */
/*     (3^(1/x) - 3^(-1/x)) / (3^(1/x) + 3^(-1/x))    at x -> 0            */
/*                                                                         */
/* where Series / L'Hospital can't expand around the essential singularity */
/* but the one-sided limits (+1 from above, -1 from below) are both        */
/* reachable through the standard layers (the divergent Power collapses to */
/* 0 or Infinity on one side and the quotient simplifies).                 */
/*                                                                         */
/* Gating is deliberately narrow: only at a finite numeric point, only     */
/* with dir == TWOSIDED, only when f has the variable in some exponent     */
/* (the failure mode we're after), and only near the top of the recursion */
/* (depth cap prevents a combinatorial explosion in deeply nested Limits). */
/* ---------------------------------------------------------------------- */

/* ------------------------------------------------------------------ */
/* Abs rewrite layer                                                  */
/*                                                                    */
/* For one-sided limits, rewrite Abs[g(x)] -> g(x) or -g(x) according */
/* to the sign of g approaching the limit point in the given          */
/* direction. For two-sided limits where the rewrite would give       */
/* different answers on the two sides, return Indeterminate. This     */
/* layer catches the canonical Limit[Abs[g(x)]/h(x), x -> a, ...]     */
/* shapes where the kink at the zero of g would otherwise reach       */
/* L'Hospital / Series and produce Derivative[1][Abs][a] (a non-      */
/* numeric symbolic value that those layers cannot classify and       */
/* mistakenly accept as "determinate").                               */
/*                                                                    */
/* Sign analysis for a g that vanishes at a finite point uses the     */
/* iterated derivative test: if g^(k)(a) is the first nonzero         */
/* derivative, then g(a + t) = g^(k)(a)/k! * t^k + O(t^(k+1)), so     */
/*   FromAbove (t > 0+):  sign(g) = sign(g^(k)(a))                    */
/*   FromBelow (t < 0):   sign(g) = sign(g^(k)(a)) * (-1)^k           */
/* ------------------------------------------------------------------ */
static int sign_at_finite_zero(Expr* g, LimitCtx* ctx) {
    Expr* current = expr_copy(g);
    int sign = 0;
    for (int order = 1; order <= 5; order++) {
        Expr* d = simp(mk_fn2("D", current, expr_copy(ctx->x)));
        current = d;
        Expr* at_pt = subst_eval(current, ctx->x, ctx->point);
        int s = literal_sign(at_pt);
        bool zero = is_lit_zero(at_pt);
        bool divergent = is_divergent(at_pt);
        bool x_residual = expr_contains(at_pt, ctx->x);
        expr_free(at_pt);
        if (s != 0) {
            int side = +1;
            if (ctx->dir == LIMIT_DIR_FROMBELOW) {
                side = (order % 2 == 0) ? +1 : -1;
            }
            sign = s * side;
            break;
        }
        if (zero) continue;          /* derivative still vanishes -- try next order */
        if (divergent || x_residual) break; /* unresolved symbolic -- give up */
        break;                       /* unsigned literal we can't classify -- give up */
    }
    expr_free(current);
    return sign;
}

/* Use a recursive one-sided Limit on g to classify its sign. Falls back
 * to 0 (undetermined) when the inner limit cannot be resolved or the
 * value still depends on x. */
static int sign_via_inner_limit(Expr* g, LimitCtx* ctx) {
    LimitCtx sub = { ctx->x, ctx->point, ctx->dir, ctx->depth, ctx->method };
    Expr* g_lim = compute_limit(g, &sub);
    if (!g_lim) return 0;
    int sgn = 0;
    if (is_infinity_sym(g_lim)) sgn = +1;
    else if (is_neg_infinity(g_lim)) sgn = -1;
    else if (!expr_contains(g_lim, ctx->x)) sgn = literal_sign(g_lim);
    expr_free(g_lim);
    return sgn;
}

/* Returns +1 / -1 / 0 (undetermined). Only meaningful for one-sided dir. */
static int sign_near_point(Expr* g, LimitCtx* ctx) {
    if (ctx->dir != LIMIT_DIR_FROMABOVE && ctx->dir != LIMIT_DIR_FROMBELOW) {
        return 0;
    }

    Expr* g_at = subst_eval(g, ctx->x, ctx->point);
    if (is_infinity_sym(g_at)) { expr_free(g_at); return +1; }
    if (is_neg_infinity(g_at)) { expr_free(g_at); return -1; }
    int s = literal_sign(g_at);
    if (s != 0) { expr_free(g_at); return s; }
    bool zero = is_lit_zero(g_at);
    bool divergent = is_divergent(g_at);
    bool x_residual = expr_contains(g_at, ctx->x);
    expr_free(g_at);

    /* g(point) is a literal zero: leading-order derivative test
     * (only applies at finite points; at infinity we'd need a 1/x
     * substitution which is delegated to the inner-limit fallback below). */
    if (zero && !is_infinity_sym(ctx->point) && !is_neg_infinity(ctx->point) &&
        !is_complex_infinity(ctx->point)) {
        int s2 = sign_at_finite_zero(g, ctx);
        if (s2 != 0) return s2;
    }

    /* Either g(point) is divergent (ComplexInfinity, DirectedInfinity,
     * Indeterminate, ...), or naive substitution left an x-residual, or
     * the leading-order test was inconclusive. Recurse: compute the
     * proper one-sided Limit of g. This handles Abs[1/x] (1/x diverges
     * with a definite one-sided sign), Abs[Tan[x]] near Pi/2, Abs[Log[x]]
     * near 0+, etc. */
    if (divergent || x_residual || zero) {
        return sign_via_inner_limit(g, ctx);
    }
    return 0;
}

/* Walk e, rewriting each x-dependent Abs[g] to g or -g per direction
 * sign. Sets *changed if any rewrite fires. */
static Expr* rewrite_abs_in_expr(Expr* e, LimitCtx* ctx, bool* changed) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    size_t n = e->data.function.arg_count;
    Expr* head = rewrite_abs_in_expr(e->data.function.head, ctx, changed);
    Expr** args = (Expr**)malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        args[i] = rewrite_abs_in_expr(e->data.function.args[i], ctx, changed);
    }
    Expr* current = expr_new_function(head, args, n);
    free(args);

    if (current->data.function.head->type == EXPR_SYMBOL && n == 1 &&
        current->data.function.head->data.symbol.name == SYM_Abs &&
        expr_contains(current->data.function.args[0], ctx->x)) {
        Expr* g = current->data.function.args[0];
        int s = sign_near_point(g, ctx);
        if (s == +1) {
            Expr* gcopy = expr_copy(g);
            expr_free(current);
            *changed = true;
            return gcopy;
        }
        if (s == -1) {
            Expr* neg = simp(mk_neg(expr_copy(g)));
            expr_free(current);
            *changed = true;
            return neg;
        }
    }
    return current;
}

static bool contains_abs_over(Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Abs &&
        e->data.function.arg_count == 1 &&
        expr_contains(e->data.function.args[0], x)) {
        return true;
    }
    if (contains_abs_over(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_abs_over(e->data.function.args[i], x)) return true;
    }
    return false;
}

static Expr* layer_abs_rewrite(Expr* f, LimitCtx* ctx) {
    if (!contains_abs_over(f, ctx->x)) return NULL;

    /* Fast path: f is exactly Abs[g(x)]. Abs is continuous everywhere
     * (including the kink at zero), so Limit[Abs[g], x->a] = Abs[Limit[g,
     * x->a]] for every direction. Resolves shapes like Abs[Tan[x]] near
     * Pi/2, Abs[1/x] near 0, Abs[Log[x]] near 0+, etc., where the inner
     * limit diverges and the structural rewrite path would produce an
     * intermediate -g that the rest of the pipeline cannot simplify. */
    if (head_is(f, SYM_Abs) && f->data.function.arg_count == 1 &&
        expr_contains(f->data.function.args[0], ctx->x)) {
        LimitCtx sub = { ctx->x, ctx->point, ctx->dir, ctx->depth, ctx->method };
        Expr* g_lim = compute_limit(f->data.function.args[0], &sub);
        if (g_lim) {
            if (is_infinity_sym(g_lim) || is_neg_infinity(g_lim) ||
                is_complex_infinity(g_lim) || is_directed_infinity(g_lim)) {
                expr_free(g_lim);
                return mk_sym("Infinity");
            }
            if (is_indeterminate(g_lim)) {
                return g_lim;
            }
            if (!expr_contains(g_lim, ctx->x)) {
                return simp(mk_fn1("Abs", g_lim));
            }
            expr_free(g_lim);
        }
        /* Inner limit unresolved; fall through to the structural rewrite
         * which may still resolve via direction-aware sign analysis. */
    }

    /* One-sided: rewrite x-dependent Abs[g] using the direction's sign,
     * then recurse on the (Abs-free, or strictly simpler) rewritten form. */
    if (ctx->dir == LIMIT_DIR_FROMABOVE || ctx->dir == LIMIT_DIR_FROMBELOW) {
        bool changed = false;
        Expr* rewritten = rewrite_abs_in_expr(f, ctx, &changed);
        if (!changed) {
            expr_free(rewritten);
            return NULL;
        }
        Expr* result = compute_limit(rewritten, ctx);
        expr_free(rewritten);
        return result;
    }

    /* Two-sided / Reals: probe the two sides. If they agree we have the
     * common limit; if they disagree the two-sided limit is genuinely
     * Indeterminate. Only fires when an x-dependent Abs is present so we
     * don't reroute every two-sided limit through the disagree pair. */
    if (ctx->dir == LIMIT_DIR_TWOSIDED || ctx->dir == LIMIT_DIR_REALS) {
        if (ctx->depth > LIMIT_MAX_DEPTH - 4) return NULL;
        LimitCtx left  = { ctx->x, ctx->point, LIMIT_DIR_FROMBELOW, ctx->depth, ctx->method };
        LimitCtx right = { ctx->x, ctx->point, LIMIT_DIR_FROMABOVE, ctx->depth, ctx->method };
        Expr* L = compute_limit(f, &left);
        if (!L) return NULL;
        if (expr_contains(L, ctx->x) || is_indeterminate(L)) {
            expr_free(L);
            return NULL;
        }
        Expr* R = compute_limit(f, &right);
        if (!R) { expr_free(L); return NULL; }
        if (expr_contains(R, ctx->x) || is_indeterminate(R)) {
            expr_free(L); expr_free(R);
            return NULL;
        }
        if (expr_eq(L, R)) {
            expr_free(R);
            return L;
        }
        expr_free(L); expr_free(R);
        return mk_sym("Indeterminate");
    }

    return NULL;
}

static Expr* layer_onesided_disagree(Expr* f, LimitCtx* ctx) {
    if (ctx->dir != LIMIT_DIR_TWOSIDED) return NULL;
    if (!is_numeric_literal_point(ctx->point)) return NULL;
    if (is_divergent(ctx->point)) return NULL;
    /* Only fire for shapes that defeat the closed-form layers. Having x in
     * an exponent is the canonical trigger (3^(1/x), (1+1/x)^x, and the
     * family); it also keeps the cost in check -- ordinary rationals and
     * analytic points never reach this layer because the earlier pipeline
     * resolves them. */
    if (!has_var_in_exponent(f, ctx->x)) return NULL;
    /* Depth guard so recursive inner Limits don't expand the pair on every
     * step. The two recursive compute_limit calls below run with a FROM*
     * direction, so the early dir!=TWOSIDED check above prevents re-entry. */
    if (ctx->depth > 3) return NULL;

    LimitCtx left  = { ctx->x, ctx->point, LIMIT_DIR_FROMBELOW, ctx->depth, ctx->method };
    LimitCtx right = { ctx->x, ctx->point, LIMIT_DIR_FROMABOVE, ctx->depth, ctx->method };

    Expr* L = compute_limit(f, &left);
    if (!L) return NULL;
    if (expr_contains(L, ctx->x)) { expr_free(L); return NULL; }

    Expr* R = compute_limit(f, &right);
    if (!R) { expr_free(L); return NULL; }
    if (expr_contains(R, ctx->x)) { expr_free(L); expr_free(R); return NULL; }

    /* Both sides agree -> that's the two-sided value. */
    if (expr_eq(L, R)) { expr_free(R); return L; }

    /* Both sides well-defined but distinct -> Indeterminate. The one-sided
     * outputs are kept long enough to confirm they are both not divergent
     * (an Indeterminate on either side would be a wash and we should stay
     * unevaluated; Infinity on exactly one side is a genuine disagreement). */
    bool L_indet = is_indeterminate(L);
    bool R_indet = is_indeterminate(R);
    expr_free(L); expr_free(R);
    if (L_indet || R_indet) return NULL;
    return mk_sym("Indeterminate");
}

/* ---------------------------------------------------------------------- */
/* Layer -- constant-factor linearity                                      */
/*                                                                         */
/* Limit[c f(x), x -> a] = c Limit[f(x), x -> a] for any factor c free of  */
/* the limit variable. Runs late in the cascade -- only after the analytic */
/* layers have failed -- so it never pre-empts a direct evaluation. It      */
/* rescues composite shapes like 2 ArcTan[g(x)] whose constant multiple    */
/* keeps the single-argument layer_compose_at_infinity from firing on the  */
/* enclosing Times (also fixes Limit[3 ArcTan[t], t -> Infinity]).         */
/* ---------------------------------------------------------------------- */
static Expr* layer_constant_factor(Expr* f, LimitCtx* ctx) {
    if (!head_is(f, SYM_Times)) return NULL;
    size_t n = f->data.function.arg_count;

    Expr* c = mk_int(1);
    Expr** var_factors = (Expr**)malloc(n * sizeof(Expr*));
    size_t nv = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* a = f->data.function.args[i];
        if (free_of(a, ctx->x)) c = simp(mk_times(c, expr_copy(a)));
        else var_factors[nv++] = expr_copy(a);
    }

    /* Need a genuine split (at least one constant and one x-dependent
     * factor) and a finite, nonzero constant so that c * (a possibly
     * divergent inner limit) stays well-defined -- 0 * Infinity would be
     * an indeterminate form we must not fabricate. */
    if (nv == 0 || nv == n || is_divergent(c) || is_lit_zero(c)) {
        expr_free(c);
        for (size_t i = 0; i < nv; i++) expr_free(var_factors[i]);
        free(var_factors);
        return NULL;
    }

    Expr* rest = (nv == 1)
        ? var_factors[0]
        : expr_new_function(expr_new_symbol(SYM_Times), var_factors, nv);
    free(var_factors);

    LimitCtx sub = *ctx; sub.depth += 1;
    Expr* lim_rest = compute_limit(rest, &sub);
    expr_free(rest);
    if (!lim_rest || expr_contains(lim_rest, ctx->x) || is_indeterminate(lim_rest)) {
        expr_free(c);
        if (lim_rest) expr_free(lim_rest);
        return NULL;
    }
    return simp(mk_times(c, lim_rest));
}

/* ---------------------------------------------------------------------- */
/* Gruntz mrv-algorithm layer. Thin adapter around gruntz_limit(): the      */
/* engine works at a fixed limit point/direction, so we hand it the ctx's   */
/* variable, point, and direction. Only real directions are supported;      */
/* complex/branch-cut requests fall through. Returns a fresh Expr* or NULL. */
/* ---------------------------------------------------------------------- */
static Expr* layer_gruntz(Expr* f, LimitCtx* ctx) {
    int gdir;
    switch (ctx->dir) {
        case LIMIT_DIR_FROMABOVE: gdir =  1; break;
        case LIMIT_DIR_FROMBELOW: gdir = -1; break;
        case LIMIT_DIR_TWOSIDED:
        case LIMIT_DIR_REALS:     gdir =  0; break;
        default:                  return NULL; /* complex / imaginary: skip */
    }
    return gruntz_limit(f, ctx->x, ctx->point, gdir, ctx->depth);
}

/* ---------------------------------------------------------------------- */
/* Top-level dispatch                                                      */
/* ---------------------------------------------------------------------- */
static Expr* compute_limit(Expr* f_in, LimitCtx* ctx) {
    if (ctx->depth >= LIMIT_MAX_DEPTH) return NULL;
    ctx->depth++;

    /* Refuse early if f applies an undefined/opaque head to anything
     * involving the limit variable. Without a continuity assumption we
     * cannot simplify `f[x]` to `f[a]`, and none of the analytic layers
     * (Series, L'Hospital, log reduction) can make progress either. */
    if (contains_opaque_head_over(f_in, ctx->x)) {
        ctx->depth--;
        return NULL;
    }

    /* At an infinite point, a *bare* single hyperbolic call (Tanh[g], Coth[g],
     * Sech[g], Csch[g], Sinh[g], Cosh[g]) is left completely intact -- both the
     * reciprocal-trig normalization and the exp rewrite are skipped -- so
     * layer_compose_at_infinity can fold it via the head's own value at
     * Infinity (Tanh[Infinity]=1, Sech[Infinity]=0, Sinh[Infinity]=Infinity...).
     * rewrite_reciprocal_trig would otherwise turn Tanh/Coth/Sech/Csch into
     * Sinh/Cosh ratios, which the exp rewrite then expands into an
     * asymptotically opaque quotient that the Series/termwise layers miss. */
    bool at_inf = is_infinity_sym(ctx->point) || is_neg_infinity(ctx->point);
    bool bare_hyp_at_inf = at_inf && is_bare_hyperbolic_call(f_in);

    /* Normalize reciprocal trig up-front. This is cheap (one tree walk + one
     * evaluate) and rescues the continuous-substitution path on shapes like
     * x Csc[x], x Cot[a x], Sec[2x] (1 - Tan[x]), etc. */
    Expr* rewritten = bare_hyp_at_inf ? expr_copy(f_in)
                                      : rewrite_reciprocal_trig(f_in);

    /* At Infinity the Sinh/Cosh/Tanh series kernels are misleading because
     * their leading behaviour is dominated by a single Exp term. Rewrite them
     * in exponential form so Series can see the dominant growth and the
     * term-wise Plus layer can cancel decaying tails. We also Expand so shapes
     * like `E^(-x) (1 + (E^x - E^-x)/2)` fold into `1/2 - E^(-2x)/2 + E^(-x)`,
     * which the term-wise Plus layer can directly sum. */
    if (at_inf && !bare_hyp_at_inf) {
        Expr* r2 = rewrite_hyperbolic_to_exp(rewritten);
        expr_free(rewritten);
        rewritten = simp(mk_fn1("Expand", r2));
    }

    Expr* f = simp(rewritten);

    /* Re-check for opaque / discontinuous heads AFTER evaluation. User-
     * defined wrappers (`h[n_] := Ceiling[n]`) look "known" to the pre-simp
     * gate — h has DownValues — but the expression they unfold to may
     * involve a discontinuous head (Ceiling, Floor, Sign, ...) applied to
     * something that still depends on x. If we let the continuous-
     * substitution fast path run we would silently pick one side's value
     * at a jump. Bail out here so the caller sees an unevaluated Limit.
     * The check is also useful against `Ceiling[g[x]]` written directly
     * with a finite numeric point: the naive substitution would produce
     * `Ceiling[finite]` and then return it as if analytic. */
    if (contains_opaque_head_over(f, ctx->x)) {
        expr_free(f);
        ctx->depth--;
        return NULL;
    }

    Expr* r = NULL;

    /* Method gating. `only` is 0 (unrestricted) unless a specific Method
     * was requested AND this is the outermost call (depth==1). Restricting
     * only at the top level keeps recursive sub-limits on the full cascade.
     * TRY runs a layer if it belongs to the selected group, and short-
     * circuits the whole function on the first success. */
    int only = (ctx->method != LIMIT_M_AUTOMATIC && ctx->depth == 1)
               ? ctx->method : 0;
    #define TRY(grp, call) do {                                  \
        if (!only || only == (grp)) {                            \
            r = (call);                                          \
            if (r) { expr_free(f); ctx->depth--; return r; }     \
        } } while (0)

    /* Layer 1: structural fast paths, including continuous substitution. */
    TRY(LIMIT_M_SUBSTITUTION, layer1_fast_paths(f, ctx));

    /* Abs[g(x)] direction-aware rewrite. Runs early so the kink at the
     * zero of g is resolved before L'Hospital / Series can latch onto
     * Derivative[1][Abs][...] as a "clean" answer. */
    TRY(LIMIT_M_SUBSTITUTION, layer_abs_rewrite(f, ctx));

    /* f[g(x)] at +/-Infinity via the head's own value there (Erf, Tanh,
     * ArcTan, Exp, Gamma, ...). Generalises the old ArcTan/ArcCot rule. */
    TRY(LIMIT_M_ASYMPTOTIC, layer_compose_at_infinity(f, ctx));

    /* Layer 3: rational function short-cut (P(x)/Q(x) classical forms). */
    TRY(LIMIT_M_RATIONAL, layer3_rational(f, ctx));

    /* Log[g(x)] at infinity when g has a finite inner limit. Runs before
     * Series because Series can miss Log[1 + decay] at infinity shapes. */
    TRY(LIMIT_M_ASYMPTOTIC, layer_log_of_finite(f, ctx));

    /* Gruntz-lite: Log[sum] at Infinity with a unique dominant summand. */
    TRY(LIMIT_M_ASYMPTOTIC, layer_gruntz_iterated_log(f, ctx));

    /* Log + linear merge at infinity for shapes like -x + Log[2 + E^x],
     * which individually diverge but combine to a finite Log. */
    TRY(LIMIT_M_ASYMPTOTIC, layer_log_merge(f, ctx));

    /* Plus at infinity where every summand has a finite limit; sum them. */
    TRY(LIMIT_M_ASYMPTOTIC, layer_plus_termwise(f, ctx));

    /* Layer 5.3: f^g log reduction. Must run before Series because
     * Series does not have a kernel for Power[base, exp] where exp is
     * itself an expression in x (e.g. (1 + a/x)^(b x)). */
    TRY(LIMIT_M_ASYMPTOTIC, layer5_log_reduction(f, ctx));

    /* Bounded envelope before Series. Series often blows up on bounded
     * heads at infinity (Sin[t^2] has no Taylor expansion at infinity),
     * so we try to squeeze to 0 first. Gating is enforced inside the
     * layer (currently only fires at +Infinity). */
    TRY(LIMIT_M_BOUNDED, layer_bounded_envelope(f, ctx));

    /* Layer 2: series-based evaluation -- the workhorse. */
    TRY(LIMIT_M_SERIES, layer2_series(f, ctx));

    /* Gruntz's mrv algorithm. Reached exclusively under Method->"Gruntz",
     * and as a last-resort fallback under Automatic (fires only when the
     * cheaper series/log layers above returned NULL). Handles the hard
     * cancellation-heavy nested-exponential exp-log limits. */
    TRY(LIMIT_M_GRUNTZ, layer_gruntz(f, ctx));

    /* Layer 5.1: L'Hospital's rule with guardrails. */
    TRY(LIMIT_M_LHOSPITAL, layer5_lhospital(f, ctx));

    /* Layer 6: bounded-oscillation Interval. */
    TRY(LIMIT_M_BOUNDED, layer6_bounded(f, ctx));

    /* Last-resort Plus split: peel off the summands with a finite individual
     * limit and recurse on the leftover group so its mutual divergences can
     * cancel (e.g. the  a ArcTan[p] + a ArcTan[q] - b Log[u] + b Log[v]  real
     * log-part of a rational antiderivative at +/-Infinity). Placed after
     * Series / L'Hospital so it only fires on Plus shapes those cannot fold,
     * never pre-empting a cleaner whole-expression evaluation. */
    TRY(LIMIT_M_ASYMPTOTIC, layer_plus_split_convergent(f, ctx));

    /* Atom substitution: rewrite a Power[b, e(x)] subterm as a fresh
     * symbol u and recurse on the u-limit. Catches shapes like
     * (-1 + 3^(2/x))/(1 + 3^(2/x)) at x -> 0 one-sided, where Series
     * hits the essential singularity but the ratio is polynomial-in-u. */
    TRY(LIMIT_M_SUBSTITUTION, layer_atom_substitute(f, ctx));

    /* Constant-factor linearity: Limit[c f, x->a] = c Limit[f, x->a]. Runs
     * late so it only rescues shapes the analytic layers left unresolved
     * (e.g. 2 ArcTan[g] where compose-at-infinity can't see through the
     * enclosing Times). */
    TRY(LIMIT_M_SUBSTITUTION, layer_constant_factor(f, ctx));

    /* Last-resort one-sided disagreement probe for two-sided limits at a
     * finite numeric point with an x-bearing exponent. Returns
     * Indeterminate for the 3^(1/x)-family shapes. */
    TRY(LIMIT_M_SUBSTITUTION, layer_onesided_disagree(f, ctx));

    #undef TRY

    expr_free(f);
    ctx->depth--;
    return NULL; /* leave unevaluated */
}

/* ---------------------------------------------------------------------- */
/* Interface normalization                                                 */
/* ---------------------------------------------------------------------- */

/* Given a single Rule[x, a], fill *out_var / *out_point and return true.
 * The caller retains ownership of rule; the out-pointers are borrowed
 * into rule's children (callers must not free them separately). */
static bool split_rule(Expr* rule, Expr** out_var, Expr** out_point) {
    if (!head_is(rule, SYM_Rule) || rule->data.function.arg_count != 2) return false;
    *out_var   = rule->data.function.args[0];
    *out_point = rule->data.function.args[1];
    return true;
}

/* Extract the `Assumptions -> ...` option (if any). Returns the raw
 * option value (borrowed from opts) or NULL when not supplied. */
static Expr* find_assumptions_opt(Expr** opts, size_t nopts) {
    for (size_t i = 0; i < nopts; i++) {
        Expr* o = opts[i];
        if (head_is(o, SYM_Rule) && o->data.function.arg_count == 2 &&
            is_sym(o->data.function.args[0], "Assumptions")) {
            return o->data.function.args[1];
        }
        if (head_is(o, SYM_RuleDelayed) && o->data.function.arg_count == 2 &&
            is_sym(o->data.function.args[0], "Assumptions")) {
            return o->data.function.args[1];
        }
    }
    return NULL;
}

/* Parse an Abs-comparison assumption clause:
 *   Less[Abs[expr], c],     LessEqual[Abs[expr], c]
 *   Greater[Abs[expr], c],  GreaterEqual[Abs[expr], c]
 *   Equal[Abs[expr], c]
 * (and the swapped variants where the constant appears on the left)
 * On success: *out_expr -> NON-owned pointer to the Abs argument,
 *             *out_op   -> -1 (<), -2 (<=), +1 (>), +2 (>=), 0 (==)
 *             *out_const -> NON-owned pointer to the constant.
 * Returns false otherwise. The constant is required to be numeric so a
 * caller can compare it against 1 (the threshold separating contraction
 * from blow-up for x^n at infinity). */
static bool assumption_abs_compare(Expr* a, Expr** out_expr, int* out_op,
                                   Expr** out_const) {
    if (!a || a->type != EXPR_FUNCTION || a->data.function.arg_count != 2) return false;
    const char* head_name = NULL;
    if (a->data.function.head->type == EXPR_SYMBOL) {
        head_name = a->data.function.head->data.symbol.name;
    }
    if (!head_name) return false;

    int op;
    bool can_swap = false;
    if (!strcmp(head_name, "Less"))           { op = -1; can_swap = true; }
    else if (!strcmp(head_name, "LessEqual")) { op = -2; can_swap = true; }
    else if (!strcmp(head_name, "Greater"))   { op = +1; can_swap = true; }
    else if (!strcmp(head_name, "GreaterEqual")) { op = +2; can_swap = true; }
    else if (!strcmp(head_name, "Equal"))     { op = 0;  can_swap = true; }
    else return false;

    Expr* lhs = a->data.function.args[0];
    Expr* rhs = a->data.function.args[1];

    /* Try lhs = Abs[expr], rhs = constant. */
    if (head_is(lhs, SYM_Abs) && lhs->data.function.arg_count == 1 &&
        expr_is_numeric_like(rhs)) {
        *out_expr  = lhs->data.function.args[0];
        *out_op    = op;
        *out_const = rhs;
        return true;
    }
    /* Try the swapped form: lhs = const, rhs = Abs[expr]. The op flips
     * sense (Less <-> Greater, LessEqual <-> GreaterEqual). */
    if (can_swap && head_is(rhs, SYM_Abs) && rhs->data.function.arg_count == 1 &&
        expr_is_numeric_like(lhs)) {
        *out_expr  = rhs->data.function.args[0];
        if (op == -1) *out_op = +1;
        else if (op == -2) *out_op = +2;
        else if (op == +1) *out_op = -1;
        else if (op == +2) *out_op = -2;
        else *out_op = 0;            /* Equal stays Equal */
        *out_const = lhs;
        return true;
    }
    return false;
}

/* Compare a numeric `c` against integer 1. Returns -1 (c < 1), 0 (c == 1),
 * +1 (c > 1). Used to dispatch Limit[B^var, var->Infinity] under an
 * assumption Abs[B] R c. */
static int compare_numeric_to_one(Expr* c) {
    if (c->type == EXPR_INTEGER) {
        if (c->data.integer < 1) return -1;
        if (c->data.integer > 1) return +1;
        return 0;
    }
    if (c->type == EXPR_BIGINT) {
        int s = mpz_cmp_ui(c->data.bigint, 1);
        return (s < 0) ? -1 : (s > 0) ? +1 : 0;
    }
    if (c->type == EXPR_REAL) {
        if (c->data.real < 1.0) return -1;
        if (c->data.real > 1.0) return +1;
        return 0;
    }
    int64_t n, d;
    if (is_rational(c, &n, &d)) {
        /* sign(n - d) when d > 0 (canonical Rational guarantees d > 0). */
        if (n < d) return -1;
        if (n > d) return +1;
        return 0;
    }
    return 0;            /* unknown */
}

/* If f matches Power[base, lim_var] with base free of lim_var, set
 * *out_base and return true. Otherwise return false. */
static bool match_power_in_var(Expr* f, Expr* lim_var, Expr** out_base) {
    if (f && f->type == EXPR_FUNCTION &&
        is_sym(f->data.function.head, "Power") &&
        f->data.function.arg_count == 2 &&
        expr_eq(f->data.function.args[1], lim_var) &&
        free_of(f->data.function.args[0], lim_var)) {
        *out_base = f->data.function.args[0];
        return true;
    }
    return false;
}

/* Dispatch Limit[Power[B, var], var -> ±Infinity] under an assumption
 * that pins Abs[B] against 1. Returns a fresh result Expr* on a clean
 * dispatch, or NULL when the assumption does not apply (caller falls
 * back to the standard Limit machinery).
 *
 * Decision table (var -> +Infinity):
 *   Abs[B] < 1    -> 0
 *   Abs[B] <= 1   -> 0  if strict-comparison is sound; here we conservatively
 *                       use Indeterminate at the |B| == 1 boundary, so we
 *                       only reduce strict <.
 *   Abs[B] > 1    -> ComplexInfinity
 *   Abs[B] >= 1   -> ditto: only reduce strict >.
 *   Abs[B] == 1   -> Indeterminate
 *
 * For var -> -Infinity, swap (< and > inversion).
 */
static Expr* limit_power_under_abs_assumption(Expr* f, Expr* lim_var,
                                              Expr* point, Expr* assumption) {
    Expr* base = NULL;
    if (!match_power_in_var(f, lim_var, &base)) return NULL;

    Expr* a_expr = NULL;
    int op = 0;
    Expr* a_const = NULL;
    if (!assumption_abs_compare(assumption, &a_expr, &op, &a_const)) return NULL;

    /* The assumption must constrain the limit's base. */
    if (!expr_eq(a_expr, base)) return NULL;

    int c_vs_one = compare_numeric_to_one(a_const);

    bool to_pos_inf = is_infinity_sym(point);
    bool to_neg_inf = is_neg_infinity(point);
    if (!to_pos_inf && !to_neg_inf) return NULL;

    /* Translate (op, c vs 1) into a verdict on |B| relative to 1 under
     * the assumption:
     *   verdict = -1 (|B| < 1),  +1 (|B| > 1),  0 (|B| == 1),
     *             else 99 (no determinate verdict from this assumption). */
    int verdict = 99;
    if (op == 0 && c_vs_one == 0) verdict = 0;     /* Abs[B] == 1 */
    else if (op == -1 && c_vs_one <= 0) verdict = -1;   /* Abs[B] < 1 */
    else if (op == -2 && c_vs_one < 0)  verdict = -1;   /* Abs[B] <= c < 1 */
    else if (op == +1 && c_vs_one >= 0) verdict = +1;   /* Abs[B] > c >= 1 */
    else if (op == +2 && c_vs_one > 0)  verdict = +1;   /* Abs[B] >= c > 1 */
    if (verdict == 99) return NULL;

    if (to_neg_inf) verdict = -verdict;            /* invert for x^(-Infinity) */

    if (verdict == -1) return mk_int(0);
    if (verdict == +1) return mk_sym("ComplexInfinity");
    return mk_sym("Indeterminate");                /* boundary |B|==1 */
}

/* Extract the `Direction -> ...` option (if any) from a list of option
 * arguments. Returns the raw option value (borrowed from opts) or NULL
 * when no Direction option was supplied. */
static Expr* find_direction_opt(Expr** opts, size_t nopts) {
    for (size_t i = 0; i < nopts; i++) {
        Expr* o = opts[i];
        if (head_is(o, SYM_Rule) && o->data.function.arg_count == 2 &&
            is_sym(o->data.function.args[0], "Direction")) {
            return o->data.function.args[1];
        }
        if (head_is(o, SYM_RuleDelayed) && o->data.function.arg_count == 2 &&
            is_sym(o->data.function.args[0], "Direction")) {
            return o->data.function.args[1];
        }
    }
    return NULL;
}

/* Extract the `Method -> ...` option (if any) from a list of option
 * arguments. Returns the raw option value (borrowed from opts) or NULL
 * when no Method option was supplied. */
static Expr* find_method_opt(Expr** opts, size_t nopts) {
    for (size_t i = 0; i < nopts; i++) {
        Expr* o = opts[i];
        if ((head_is(o, SYM_Rule) || head_is(o, SYM_RuleDelayed)) &&
            o->data.function.arg_count == 2 &&
            is_sym(o->data.function.args[0], "Method")) {
            return o->data.function.args[1];
        }
    }
    return NULL;
}

/* Handle Limit[f, {x1 -> a1, ..., xn -> an}] by iterated right-to-left
 * folding: the innermost (rightmost) rule's Limit is computed first. */
static Expr* run_iterated(Expr* f, Expr* rule_list, int dir, int method, int depth) {
    Expr* current = expr_copy(f);
    size_t n = rule_list->data.function.arg_count;
    for (size_t i = n; i-- > 0; ) {
        Expr* rule = rule_list->data.function.args[i];
        Expr *var, *point;
        if (!split_rule(rule, &var, &point)) { expr_free(current); return NULL; }
        LimitCtx ctx = { var, point, dir, depth, method };
        Expr* next = compute_limit(current, &ctx);
        expr_free(current);
        if (!next) return NULL;
        current = next;
    }
    return current;
}

/* ---------------------------------------------------------------------- */
/* Multivariate limits: polar/spherical path-dependence analysis           */
/* ---------------------------------------------------------------------- */
/* All points are {0,0,...} or all are {Infinity,Infinity,...}. The
 * n-dimensional joint limit gets an (n-1)-parameter angular substitution
 * and we compute the resulting single-variable limit in r. If the answer
 * is free of the angles, that IS the joint limit. If it depends on
 * angles, we sample a handful of axis and diagonal directions; if they
 * all agree we return that value, otherwise we return Indeterminate.  */

/* True iff every entry of `points` matches `val` in canonical form
 * (integer 0 for origin, or the Infinity symbol for +Infinity). */
static bool all_points_are(Expr* points, int kind) {
    for (size_t i = 0; i < points->data.function.arg_count; i++) {
        Expr* p = points->data.function.args[i];
        if (kind == 0) {
            if (!is_lit_zero(p)) return false;
        } else {
            if (!is_infinity_sym(p)) return false;
        }
    }
    return true;
}

/* Compute Limit[f_polar, r -> r0, Direction -> "FromAbove"]. r is a
 * fresh symbol; r0 is 0 (origin) or Infinity. Returns a new Expr* or
 * NULL if the inner limit cannot be resolved. */
static Expr* limit_r_fromabove(Expr* f_polar, Expr* r_sym, int kind) {
    Expr* point = (kind == 0) ? mk_int(0) : mk_sym("Infinity");
    LimitCtx sub = { r_sym, point, LIMIT_DIR_FROMABOVE, 0, LIMIT_M_AUTOMATIC };
    Expr* res = compute_limit(f_polar, &sub);
    expr_free(point);
    return res;
}

/* For every multi-argument angle substitution we make, check that the
 * substituted `f` is free of the original variables. Mathilda's
 * ReplaceAll doesn't recurse through evaluation so we always simp().
 *
 * Order of substitution matters: we must replace all variables
 * simultaneously rather than one at a time, otherwise e.g. replacing
 * `x -> r Cos[t]` then `y -> r Sin[t]` would leave the expression
 * referencing the fresh `r` and `t` while `y` has only just vanished.
 * We ReplaceAll with a list of rules. */
static Expr* subst_all(Expr* f, Expr** vars, Expr** vals, size_t n) {
    Expr** rules = calloc(n, sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        rules[i] = mk_fn2("Rule", expr_copy(vars[i]), expr_copy(vals[i]));
    }
    Expr* rule_list = expr_new_function(mk_sym("List"), rules, n);
    free(rules);
    Expr* ra = mk_fn2("ReplaceAll", expr_copy(f), rule_list);
    return simp(ra);
}

/* Polar substitution at origin or infinity for 2D. Returns the r-limit
 * result as an expression in `angle_sym`, or NULL. */
static Expr* polar_2d_limit(Expr* f, Expr** vars, int kind,
                            Expr* r_sym, Expr* t_sym) {
    /* x = r Cos[t], y = r Sin[t] */
    Expr* rcos = simp(mk_times(expr_copy(r_sym), mk_fn1("Cos", expr_copy(t_sym))));
    Expr* rsin = simp(mk_times(expr_copy(r_sym), mk_fn1("Sin", expr_copy(t_sym))));
    Expr* vals[2] = { rcos, rsin };
    Expr* f_polar = subst_all(f, vars, vals, 2);
    expr_free(rcos); expr_free(rsin);

    /* Simplify before the r-limit: cancelling the common r-powers exposes
     * shapes like ArcTan[(r Sin[t])^2/((r Cos[t])^2 + (r Cos[t])^3)] as
     * ArcTan[Tan[t]^2/(1 + r Cos[t])], removing the buried 0/0 that
     * Mathilda's arithmetic would otherwise fold to a spurious 0. */
    Expr* f_simp = simp(mk_fn1("Simplify", f_polar));
    Expr* rlim = limit_r_fromabove(f_simp, r_sym, kind);
    expr_free(f_simp);
    return rlim;
}

/* Spherical substitution for 3D joint limits at the origin. */
static Expr* polar_3d_limit(Expr* f, Expr** vars, int kind,
                            Expr* r_sym, Expr* t_sym, Expr* p_sym) {
    /* x = r Sin[p] Cos[t], y = r Sin[p] Sin[t], z = r Cos[p] */
    Expr* sp = mk_fn1("Sin", expr_copy(p_sym));
    Expr* cp = mk_fn1("Cos", expr_copy(p_sym));
    Expr* ct = mk_fn1("Cos", expr_copy(t_sym));
    Expr* st = mk_fn1("Sin", expr_copy(t_sym));
    Expr* vx = simp(mk_times(expr_copy(r_sym), mk_times(expr_copy(sp), expr_copy(ct))));
    Expr* vy = simp(mk_times(expr_copy(r_sym), mk_times(expr_copy(sp), expr_copy(st))));
    Expr* vz = simp(mk_times(expr_copy(r_sym), expr_copy(cp)));
    expr_free(sp); expr_free(cp); expr_free(ct); expr_free(st);
    Expr* vals[3] = { vx, vy, vz };
    Expr* f_sph = subst_all(f, vars, vals, 3);
    expr_free(vx); expr_free(vy); expr_free(vz);

    /* Simplify first (see polar_2d_limit) so buried 0/0 shapes cancel. */
    Expr* f_simp = simp(mk_fn1("Simplify", f_sph));
    Expr* rlim = limit_r_fromabove(f_simp, r_sym, kind);
    expr_free(f_simp);
    return rlim;
}

/* Sample a handful of concrete directions and return Indeterminate if
 * they disagree, the common value if they all match, or NULL if we
 * cannot resolve. Used as a backstop when the polar r-limit depends on
 * the angle parameters. */
static Expr* sample_joint_limit(Expr* f, Expr** vars, size_t n, int kind) {
    /* Direction vectors. For kind=0 (origin) each coordinate can move
     * alone (axis directions) or together (diagonals) because "fixing at
     * 0" is meaningful -- a coordinate that's exactly zero while the
     * others approach zero is the canonical axis limit. For kind=1
     * (infinity) every coordinate must actually approach +Infinity, so
     * we only emit directions with all-positive entries: the all-ones
     * direction plus a couple of skewed diagonals (2,1,1...) and
     * (1,2,1...). */
    int max_dirs = 10;
    int (*dirs)[3] = calloc(max_dirs, sizeof(*dirs));
    int nd = 0;
    if (kind == 0) {
        for (size_t i = 0; i < n && nd < max_dirs; i++) {
            for (size_t j = 0; j < n; j++) dirs[nd][j] = (j == i) ? 1 : 0;
            nd++;
        }
        if (n == 2 && nd + 2 <= max_dirs) {
            dirs[nd][0] = 1; dirs[nd][1] = 1; nd++;
            dirs[nd][0] = 1; dirs[nd][1] = -1; nd++;
        }
        if (n == 3 && nd + 2 <= max_dirs) {
            dirs[nd][0] = 1; dirs[nd][1] = 1; dirs[nd][2] = 1; nd++;
            dirs[nd][0] = 1; dirs[nd][1] = 1; dirs[nd][2] = -1; nd++;
        }
    } else {
        /* kind == 1: all-positive directions. */
        for (size_t j = 0; j < n; j++) dirs[nd][j] = 1;
        nd++;
        for (size_t i = 0; i < n && nd < max_dirs; i++) {
            for (size_t j = 0; j < n; j++) dirs[nd][j] = (j == i) ? 2 : 1;
            nd++;
        }
    }

    /* Pick a fresh scalar parameter. */
    Expr* t_sym = mk_sym("$LimitPathT$");
    Expr* values[3] = { NULL, NULL, NULL };

    Expr* common = NULL;
    bool ok = true;
    for (int d = 0; d < nd && ok; d++) {
        for (size_t j = 0; j < n; j++) {
            int c = dirs[d][j];
            Expr* v;
            if (c == 0 && kind == 0) {
                /* Pin this coordinate at 0; the path approaches the
                 * origin only along the other axis. */
                v = mk_int(0);
            } else {
                /* x_j = c * t. Origin: t -> 0+. Infinity: t -> Infinity
                 * (all c values are positive by construction above). */
                v = simp(mk_times(mk_int(c), expr_copy(t_sym)));
            }
            values[j] = v;
        }
        bool skip = false;
        for (size_t j = 0; j < n; j++) if (!values[j]) { skip = true; break; }
        if (!skip) {
            Expr* f_path = subst_all(f, vars, values, n);
            Expr* point = (kind == 0) ? mk_int(0) : mk_sym("Infinity");
            LimitCtx sub = { t_sym, point, LIMIT_DIR_FROMABOVE, 0, LIMIT_M_AUTOMATIC };
            Expr* v = compute_limit(f_path, &sub);
            expr_free(point);
            expr_free(f_path);
            if (!v) ok = false;
            else {
                if (!common) common = v;
                else {
                    if (!expr_eq(v, common)) {
                        expr_free(common);
                        free(dirs);
                        for (size_t j = 0; j < n; j++) if (values[j]) expr_free(values[j]);
                        expr_free(v); expr_free(t_sym);
                        return mk_sym("Indeterminate");
                    }
                    expr_free(v);
                }
            }
        }
        for (size_t j = 0; j < n; j++) {
            if (values[j]) expr_free(values[j]);
            values[j] = NULL;
        }
    }
    free(dirs);
    expr_free(t_sym);
    if (!ok) { if (common) expr_free(common); return NULL; }
    return common;
}

/* One sample value for a polar angle: rational multiples of Pi and a
 * plain rational, all strictly inside (0, Pi/2) so they name valid
 * directions for both the origin (kind 0) and positive-orthant infinity
 * (kind 1) cases, while avoiding the axis degeneracies at 0 / Pi/2 where
 * Tan blows up. */
static Expr* mk_angle_value(int idx) {
    static const int tbl[][3] = {   /* {numerator, denominator, timesPi} */
        {1, 6, 1}, {1, 4, 1}, {1, 3, 1}, {1, 5, 1}, {2, 7, 1}, {2, 5, 0}
    };
    int n = (int)(sizeof(tbl) / sizeof(tbl[0]));
    const int* t = tbl[((idx % n) + n) % n];
    Expr* frac = mk_fn2("Times", mk_int(t[0]),
                        mk_fn2("Power", mk_int(t[1]), mk_int(-1)));
    return t[2] ? mk_times(frac, mk_sym("Pi")) : frac;
}
#define N_ANGLE_SAMPLES 6

/* Decide the joint limit from the polar radial limit `rlim` (a function of
 * the angle symbols, and possibly still of r_sym if the r-limit was only
 * partial):
 *   - reduces to an angle-free, r-free constant  -> that constant is the value
 *   - two sampled directions give different finite values -> Indeterminate
 *   - residual r, or fewer than two clean samples -> NULL (inconclusive)
 * This samples the clean 1-/2-parameter angular form directly, which is far
 * more reliable than re-deriving straight-line path limits from the original
 * (0/0-prone) integrand. Does not take ownership of `rlim`. */
static Expr* resolve_angular_limit(Expr* rlim, Expr** angle_syms,
                                   size_t n_angles, Expr* r_sym) {
    Expr* s = simp(expr_copy(rlim));
    if (is_divergent(s)) { expr_free(s); return NULL; }
    /* The r-limit must have eliminated r; if not, we cannot decide here. */
    if (expr_contains(s, r_sym)) { expr_free(s); return NULL; }
    bool angular = false;
    for (size_t j = 0; j < n_angles; j++)
        if (expr_contains(s, angle_syms[j])) { angular = true; break; }
    if (!angular) return s;              /* angle-free constant -> the limit */

    /* Genuinely angle-dependent: numerically confirm a disagreement. */
    Expr* common = NULL;
    int n_clean = 0;
    bool disagree = false;
    for (int k = 0; k < N_ANGLE_SAMPLES && !disagree; k++) {
        Expr* v = expr_copy(s);
        for (size_t j = 0; j < n_angles; j++) {
            Expr* a = mk_angle_value(k + (int)j);
            Expr* nv = subst_eval(v, angle_syms[j], a);
            expr_free(a); expr_free(v); v = nv;
        }
        bool clean = !is_divergent(v);
        for (size_t j = 0; j < n_angles && clean; j++)
            if (expr_contains(v, angle_syms[j])) clean = false;
        if (!clean) { expr_free(v); continue; }
        n_clean++;
        if (!common) { common = v; continue; }
        if (!expr_eq(v, common)) disagree = true;
        expr_free(v);
    }
    expr_free(s);
    if (disagree) {
        if (common) expr_free(common);
        return mk_sym("Indeterminate");
    }
    if (n_clean >= 2 && common) return common;   /* samples agree -> constant */
    if (common) expr_free(common);
    return NULL;                                 /* inconclusive */
}

/* Handle Limit[f, {x1,...,xn} -> {a1,...,an}] as a joint limit. */
static Expr* run_multivariate(Expr* f_in, Expr* vars, Expr* points) {
    if (!head_is(vars, SYM_List) || !head_is(points, SYM_List)) return NULL;
    size_t n = vars->data.function.arg_count;
    if (n != points->data.function.arg_count) return NULL;
    if (n < 2) return NULL;

    /* Limit is HoldAll, so any user-defined `f[x, y]` reaches us
     * unexpanded. Run one evaluation pass so subsequent structural
     * checks see the real expression tree (e.g. so the 0/0 scan can
     * find Power[x^2 + x^3, -1] inside an ArcTan). */
    Expr* f = simp(expr_copy(f_in));

    /* Simple substitution fast path: only trust it when no sub-expression
     * risks a 0/0 fold. Mathilda's arithmetic eagerly folds 0/0 to 0 at
     * Sin, ArcTan, and other non-rational heads; for those cases we MUST
     * go through the path-dependence analysis even if the top-level
     * Together denominator looks safe. Walking the tree for any
     * `Power[base, negative]` whose base vanishes at the joint point
     * catches the ArcTan[y^2/(x^2 + x^3)]-style 0/0. */
    Expr* tog = simp(mk_fn1("Together", expr_copy(f)));
    Expr* den = simp(mk_fn1("Denominator", expr_copy(tog)));
    Expr* den_at = expr_copy(den);
    for (size_t i = 0; i < n; i++) {
        Expr* next = subst_eval(den_at, vars->data.function.args[i],
                                        points->data.function.args[i]);
        expr_free(den_at);
        den_at = next;
    }
    bool den_bad = is_lit_zero(den_at) || is_divergent(den_at);
    expr_free(den); expr_free(tog); expr_free(den_at);
    /* Also refuse when any reciprocal sub-term's base vanishes. */
    bool inner_div_by_zero = false;
    if (!den_bad) {
        /* Recursive walk: find Power[b, k] with k negative OR 0 after
         * substitution, and check if b at the point is 0. */
        Expr* probe = expr_copy(f);
        for (size_t i = 0; i < n; i++) {
            Expr* next = subst_eval(probe, vars->data.function.args[i],
                                           points->data.function.args[i]);
            expr_free(probe);
            probe = next;
        }
        /* is_divergent on the *substituted* original catches any
         * ComplexInfinity / Indeterminate the evaluator produced; that's
         * a reliable signal of an internal 0/0 even when the outer form
         * folded to a plausible finite value. */
        if (is_divergent(probe)) inner_div_by_zero = true;
        else {
            /* Evaluate every reciprocal base at the point and see if any
             * goes to 0. */
            Expr* probe2 = expr_copy(f);
            /* lazy stack-based scan */
            Expr* stack[64]; int top = 0; stack[top++] = probe2;
            while (top > 0 && !inner_div_by_zero) {
                Expr* e = stack[--top];
                if (e->type == EXPR_FUNCTION) {
                    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
                        Expr* exp = e->data.function.args[1];
                        bool is_neg = false;
                        if (exp->type == EXPR_INTEGER && exp->data.integer < 0) is_neg = true;
                        else if (head_is(exp, SYM_Times) && exp->data.function.arg_count > 0 &&
                                 exp->data.function.args[0]->type == EXPR_INTEGER &&
                                 exp->data.function.args[0]->data.integer < 0) is_neg = true;
                        if (is_neg) {
                            Expr* b = expr_copy(e->data.function.args[0]);
                            for (size_t i = 0; i < n; i++) {
                                Expr* nb = subst_eval(b, vars->data.function.args[i],
                                                         points->data.function.args[i]);
                                expr_free(b); b = nb;
                            }
                            if (is_lit_zero(b)) inner_div_by_zero = true;
                            expr_free(b);
                        }
                    }
                    for (size_t i = 0; i < e->data.function.arg_count && top < 60; i++) {
                        stack[top++] = e->data.function.args[i];
                    }
                    if (top < 60) stack[top++] = e->data.function.head;
                }
            }
            expr_free(probe2);
        }
        expr_free(probe);
    }
    /* The simple-substitution shortcut is only safe at finite points. At
     * Infinity, Mathilda's arithmetic folds Infinity/Infinity into path-
     * dependent shortcut values (e.g. ArcTan[y/x] /. x,y -> Infinity
     * yields Pi/4 via Infinity/Infinity -> 1 -> ArcTan[1]) that hide the
     * genuine path-dependence. Gate the fast path to the origin case. */
    bool all_finite = true;
    for (size_t i = 0; i < n; i++) {
        Expr* p = points->data.function.args[i];
        if (is_infinity_sym(p) || is_neg_infinity(p) || is_complex_infinity(p)) {
            all_finite = false; break;
        }
    }
    if (all_finite && !den_bad && !inner_div_by_zero) {
        Expr* cur = expr_copy(f);
        for (size_t i = 0; i < n; i++) {
            Expr* s = subst_eval(cur, vars->data.function.args[i],
                                      points->data.function.args[i]);
            expr_free(cur);
            cur = s;
        }
        bool bad_residual = false;
        for (size_t i = 0; i < n; i++) {
            if (expr_contains(cur, vars->data.function.args[i])) {
                bad_residual = true; break;
            }
        }
        if (!is_divergent(cur) && !bad_residual) { expr_free(f); return cur; }
        expr_free(cur);
    }

    /* Path-dependence analysis via polar / spherical substitution.
     * We only handle the two canonical cases:
     *   - all points are 0 (joint limit at origin)
     *   - all points are +Infinity (joint limit at infinity along the
     *     positive orthant)
     * Mixed cases (x -> 0, y -> Infinity) stay unevaluated for now;
     * they're rare in practice and usually need a change of variables
     * from the user side anyway. */
    int kind;
    if (all_points_are(points, 0))      kind = 0;
    else if (all_points_are(points, 1)) kind = 1;
    else { expr_free(f); return NULL; }

    Expr** vars_arr = calloc(n, sizeof(Expr*));
    for (size_t i = 0; i < n; i++) vars_arr[i] = vars->data.function.args[i];

    Expr* r_sym = mk_sym("$LimitPolarR$");
    Expr* t_sym = mk_sym("$LimitPolarTheta$");
    Expr* result = NULL;

    if (n == 2) {
        Expr* rlim = polar_2d_limit(f, vars_arr, kind, r_sym, t_sym);
        if (rlim) {
            Expr* angs[1] = { t_sym };
            result = resolve_angular_limit(rlim, angs, 1, r_sym);
            expr_free(rlim);
        }
        /* Inconclusive polar analysis -> fall back to direction sampling. */
        if (!result) result = sample_joint_limit(f, vars_arr, n, kind);
    } else if (n == 3) {
        Expr* p_sym = mk_sym("$LimitPolarPhi$");
        Expr* rlim = polar_3d_limit(f, vars_arr, kind, r_sym, t_sym, p_sym);
        if (rlim) {
            Expr* angs[2] = { t_sym, p_sym };
            result = resolve_angular_limit(rlim, angs, 2, r_sym);
            expr_free(rlim);
        }
        if (!result) result = sample_joint_limit(f, vars_arr, n, kind);
        expr_free(p_sym);
    } else {
        /* n >= 4: skip polar, just sample. */
        result = sample_joint_limit(f, vars_arr, n, kind);
    }

    free(vars_arr);
    expr_free(r_sym); expr_free(t_sym);
    expr_free(f);
    return result;
}

/* ---------------------------------------------------------------------- */
/* Public entry point                                                      */
/* ---------------------------------------------------------------------- */
/* True iff `e` structurally contains the imaginary unit I -- as a bare
 * Symbol "I", as Complex[a, b] with b != 0, or nested inside any
 * function-call / Plus / Times children. Used by the branch-cut
 * post-pass to tell whether the principal-branch result actually
 * picked up an imaginary part. */
static bool contains_imaginary_unit(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == SYM_I;
    if (e->type == EXPR_FUNCTION) {
        if (head_is(e, SYM_Complex) && e->data.function.arg_count == 2) {
            Expr* im = e->data.function.args[1];
            if (im->type == EXPR_INTEGER) return im->data.integer != 0;
            if (im->type == EXPR_REAL)    return im->data.real != 0.0;
            return true;
        }
        if (contains_imaginary_unit(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (contains_imaginary_unit(e->data.function.args[i])) return true;
        }
    }
    return false;
}

/* Conjugate the imaginary part of a closed-form result. Mathilda's
 * `Conjugate` evaluator only folds for some forms (e.g. Complex[a, b]
 * literals) but leaves generic expressions like `I * Pi` unfolded, so
 * we implement this via ReplaceAll[I -> -I]. That substitution gives
 * exact complex-conjugation for any expression whose only non-real
 * content is the imaginary unit symbol, which covers Sqrt/Log branch
 * values at a negative real point. */
static Expr* conjugate_imaginary(Expr* e) {
    Expr* neg_i = mk_neg(mk_sym("I"));
    Expr* rule = mk_fn2("Rule", mk_sym("I"), neg_i);
    Expr* ra   = mk_fn2("ReplaceAll", expr_copy(e), rule);
    return simp(ra);
}

static Expr* builtin_limit_impl(Expr* res) {
    /* Contract: the evaluator frees `res` on a non-NULL return (see
     * src/eval.c); we must not free it ourselves. Return NULL to leave
     * the expression unevaluated. */
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    Expr* f   = res->data.function.args[0];
    Expr* spec= res->data.function.args[1];

    /* Thread over a top-level List in the first argument. Mathematica's
     * Limit is effectively listable on its first argument (with
     * HoldAll), so Limit[{a, b}, x -> c] maps to {Limit[a, x -> c],
     * Limit[b, x -> c]} with the options forwarded unchanged. */
    if (head_is(f, SYM_List)) {
        size_t k = f->data.function.arg_count;
        Expr** results = calloc(k, sizeof(Expr*));
        for (size_t i = 0; i < k; i++) {
            size_t nargs = res->data.function.arg_count;
            Expr** new_args = calloc(nargs, sizeof(Expr*));
            new_args[0] = expr_copy(f->data.function.args[i]);
            for (size_t j = 1; j < nargs; j++) {
                new_args[j] = expr_copy(res->data.function.args[j]);
            }
            Expr* call = expr_new_function(mk_sym("Limit"), new_args, nargs);
            free(new_args);
            results[i] = eval_and_free(call);
        }
        Expr* out = expr_new_function(mk_sym("List"), results, k);
        free(results);
        return out;
    }

    /* Collect option args (positions 2..argc-1). */
    Expr** opts = (argc > 2) ? &res->data.function.args[2] : NULL;
    size_t nopts = (argc > 2) ? argc - 2 : 0;

    Expr* dir_opt = find_direction_opt(opts, nopts);
    int dir = LIMIT_DIR_TWOSIDED;
    if (dir_opt && !parse_direction(dir_opt, &dir)) {
        /* Unknown direction -- keep unevaluated. */
        return NULL;
    }

    /* Method -> m selects the internal strategy. Automatic (the default)
     * runs the full cascade; a named method restricts the top-level call
     * to one strategy group and leaves Limit unevaluated if it fails. An
     * unrecognised value is reported and left unevaluated. */
    Expr* method_opt = find_method_opt(opts, nopts);
    int method = LIMIT_M_AUTOMATIC;
    if (method_opt && !parse_method(method_opt, &method)) {
        char* s = expr_to_string(method_opt);
        fprintf(stderr,
                "Limit::method: %s is not a recognised setting for Method.\n",
                s ? s : "?");
        free(s);
        return NULL;
    }

    /* Targeted Assumptions support: Limit[Power[B, var], var -> ±Infinity,
     * Assumptions -> Abs[B] R c]. The standard Limit machinery doesn't
     * carry assumption context, so we dispatch this specific (but useful)
     * shape here before the general pipeline. Falls through to the
     * standard machinery on no-match. */
    Expr* assumptions_opt = find_assumptions_opt(opts, nopts);
    if (assumptions_opt && head_is(spec, SYM_Rule) &&
        spec->data.function.arg_count == 2) {
        Expr* lvar  = spec->data.function.args[0];
        Expr* lpoint = spec->data.function.args[1];
        Expr* dispatched = limit_power_under_abs_assumption(f, lvar, lpoint,
                                                            assumptions_opt);
        if (dispatched) return dispatched;
    }

    /* --- Form A: Limit[f, x -> a]
     * --- Form C: Limit[f, {x1,...,xn} -> {a1,...,an}]   (multivariate) */
    if (head_is(spec, SYM_Rule) && spec->data.function.arg_count == 2) {
        if (head_is(spec->data.function.args[0], SYM_List) &&
            head_is(spec->data.function.args[1], SYM_List)) {
            return run_multivariate(f, spec->data.function.args[0],
                                       spec->data.function.args[1]);
        }
        Expr *var = NULL, *point = NULL;
        if (!split_rule(spec, &var, &point)) return NULL;
        /* Compute the base (principal-branch) limit first. The complex
         * directions (LIMIT_DIR_IMAGINARY, LIMIT_DIR_COMPLEX) are
         * routed through LIMIT_DIR_TWOSIDED for the analytic layers --
         * they mostly only affect the branch-cut post-pass below, not
         * the series / L'Hospital computation itself. */
        int inner_dir = dir;
        if (inner_dir == LIMIT_DIR_IMAGINARY || inner_dir == LIMIT_DIR_COMPLEX) {
            inner_dir = LIMIT_DIR_TWOSIDED;
        }
        LimitCtx ctx = { var, point, inner_dir, 0, method };
        Expr* base = compute_limit(f, &ctx);
        if (!base) return NULL;

        /* Branch-cut post-pass:
         *   Direction -> I        flips the imaginary part (the "other"
         *                         branch of Sqrt/Log at z = negative real)
         *   Direction -> Complexes returns Indeterminate when the base
         *                         result has a non-zero imaginary part
         *                         (the radial limit disagrees between
         *                         approach directions on either side of
         *                         the branch cut). Poles continue to
         *                         return ComplexInfinity via Layer 2. */
        if (dir == LIMIT_DIR_IMAGINARY && contains_imaginary_unit(base)) {
            Expr* flipped = conjugate_imaginary(base);
            expr_free(base);
            return flipped;
        }
        if (dir == LIMIT_DIR_COMPLEX) {
            /* Radial approach interpretation: any pole is ComplexInfinity
             * regardless of parity, and a branch-point value that picked
             * up an imaginary part is Indeterminate. Finite real results
             * pass through unchanged (they are the same from every
             * complex approach direction). */
            if (is_infinity_sym(base) || is_neg_infinity(base)) {
                expr_free(base);
                return mk_sym("ComplexInfinity");
            }
            if (contains_imaginary_unit(base) && !is_complex_infinity(base)) {
                expr_free(base);
                return mk_sym("Indeterminate");
            }
        }
        return base;
    }

    /* --- Form B: Limit[f, {x1 -> a1, ..., xn -> an}] iterated --- */
    if (head_is(spec, SYM_List)) {
        size_t n = spec->data.function.arg_count;
        if (n == 0) return NULL;
        for (size_t i = 0; i < n; i++) {
            if (!head_is(spec->data.function.args[i], SYM_Rule)) return NULL;
        }
        return run_iterated(f, spec, dir, method, 0);
    }

    return NULL;
}

/* Public entry point. Wraps the implementation with the arithmetic-warning
 * mute: Limit's internal probes (Together, Series, L'Hospital, polar
 * substitutions, direct sub-at-point attempts) routinely generate transient
 * Power::infy / Infinity::indet messages. Those are noise to the user --
 * the divergent sub-expression is expected and handled. Muting applies only
 * while Limit is running; nested Limit calls nest cleanly via the counter. */
Expr* builtin_limit(Expr* res) {
    /* Inexact coefficients break the symbolic limit machinery (which
     * leans on Together / Cancel / Series, all rational-coefficient
     * algorithms). Rationalise inputs, run, and numericalise the limit
     * value — done outside the mute_push so warnings still nest cleanly. */
    if (internal_args_contain_inexact(res)) {
        return internal_rationalize_then_numericalize(res, builtin_limit);
    }
    arith_warnings_mute_push();
    Expr* out = builtin_limit_impl(res);
    arith_warnings_mute_pop();
    return out;
}

/* ---------------------------------------------------------------------- */
/* Registration                                                            */
/* ---------------------------------------------------------------------- */
void limit_init(void) {
    symtab_add_builtin("Limit", builtin_limit);

    /* Limit does not hold its arguments in Mathematica -- Attributes[Limit]
     * is {Protected, ReadProtected}. The first argument f must be evaluated
     * so forms like Limit[%, x -> Infinity] (where % is Out[-1]) see the
     * actual expression, and the spec rule x -> a evaluates to Rule[x, a]
     * for a free symbol x. The internal layers evaluate/substitute as
     * needed, so they remain correct with pre-evaluated arguments. */
    symtab_get_def("Limit")->attributes |=
        ATTR_PROTECTED | ATTR_READPROTECTED;

    symtab_set_docstring("Limit",
        "Limit[f, x -> a]\n"
        "\tfinds the limit of f as x approaches a.\n"
        "Limit[f, {x1 -> a1, ..., xn -> an}]\n"
        "\titerated limit, applied rightmost-first.\n"
        "Limit[f, {x1, ..., xn} -> {a1, ..., an}]\n"
        "\tmultivariate (joint) limit.\n"
        "Limit[f, x -> a, Direction -> d]\n"
        "\tspecifies the direction of approach:\n"
        "\t  Reals or \"TwoSided\" -- default two-sided limit\n"
        "\t  \"FromAbove\" or -1   -- approach from above (x -> a^+)\n"
        "\t  \"FromBelow\" or +1   -- approach from below (x -> a^-)\n"
        "\t  Complexes           -- limit over all complex directions\n"
        "Limit[f, x -> a, Method -> m]\n"
        "\tselects the internal strategy:\n"
        "\t  Automatic          -- (default) try all strategies in order\n"
        "\t  \"Substitution\"     -- continuity, Abs kink, atom/one-sided probes\n"
        "\t  \"RationalFunction\" -- degree comparison for P(x)/Q(x)\n"
        "\t  \"Series\"           -- Taylor/Laurent/Puiseux leading term\n"
        "\t  \"LHospital\"        -- L'Hospital's rule for 0/0 and Inf/Inf\n"
        "\t  \"Asymptotic\"       -- dominant-term / log / exp reductions\n"
        "\t  \"Bounded\"          -- squeeze and bounded-oscillation Interval\n"
        "\tA named method leaves Limit unevaluated when it does not apply.\n"
        "\n"
        "May return a finite value, Infinity, -Infinity, ComplexInfinity,\n"
        "Indeterminate, Interval[{lo, hi}], or the original unevaluated\n"
        "expression when the limit cannot be determined.");
}

#ifndef ITER_H
#define ITER_H

#include <stdbool.h>
#include <math.h>     /* fabs: iter_real_slack */
#include "expr.h"
#include "symtab.h"   /* Rule, SymbolDef for the localization helpers */

Expr* builtin_do(Expr* res);
Expr* builtin_for(Expr* res);
Expr* builtin_while(Expr* res);
Expr* builtin_break(Expr* res);
Expr* builtin_continue(Expr* res);

void iter_init(void);

/*
 * ----------------------------------------------------------------------------
 *  Shared Wolfram iterator-spec helpers
 * ----------------------------------------------------------------------------
 *
 * A single structural parser + Block-style localization shared by every
 * iterator-driven builtin (Table, Do, Sum, ...).  This replaces the
 * near-identical spec-parsing blocks that previously lived inside each
 * builtin.  Each consumer decides for itself whether a parsed bound is
 * usable (numeric range vs. symbolic vs. Infinity), so the parser stays
 * policy-free: it only classifies and evaluates.
 */

typedef enum {
    ITER_KIND_COUNT,   /* bare n or {n}: repeat-count, var == NULL          */
    ITER_KIND_RANGE,   /* {i,imax} / {i,imin,imax} / {i,imin,imax,di}       */
    ITER_KIND_LIST     /* {i, {a,b,...}}: iterate over explicit list values */
} IterKind;

typedef struct {
    Expr*    var;    /* iterator symbol (owned copy); NULL for COUNT        */
    IterKind kind;
    Expr*    imin;   /* owned; RANGE only                                   */
    Expr*    imax;   /* owned; RANGE upper bound, or COUNT count            */
    Expr*    di;     /* owned; RANGE step (defaults to Integer 1)           */
    Expr*    list;   /* owned List[...] of values; LIST only                */
    struct LinkedList* ll;  /* LIST only, instead of `list`: an owned snapshot
                              * of a LinkedList bound ({i, ll}), iterated from
                              * its slots -- never turned into a List          */
} IterSpec;

/* The values of a LIST spec, whichever holds them: how many, and value i
 * (borrowed from the List, or a fresh box of a LinkedList slot when *owned;
 * free it after binding). */
int64_t iter_list_length(const IterSpec* s);
Expr*   iter_list_element(const IterSpec* s, int64_t i, bool* owned);

/*
 * Parse a held iterator spec into `out`.  Bound expressions are evaluate()d
 * so that symbolic bounds survive as symbols and constant arithmetic folds.
 * Returns true with `out` fully populated (caller must iter_spec_free), or
 * false for a malformed spec (variable slot not a symbol, or arity > 4), in
 * which case `out` is left freed/zeroed.
 */
bool iter_spec_parse(Expr* spec, IterSpec* out);

/* A RANGE spec whose bounds are symbolic but whose length is not --
 * {i, a, a + 2}, where (imax - imin)/di is the explicit number 2 -- is turned
 * into the LIST spec of its values imin + k di (a, 1 + a, 2 + a), so Do and
 * Table can run it.  Returns false (spec untouched) when the length is not an
 * explicit number or exceeds 10^7. */
bool iter_spec_symbolic_range_to_list(IterSpec* s);

/*
 * iter_spec_parse for a DISCRETE iterator (Table, Do, Sum, Product): as above,
 * then a RANGE whose bounds are real NumericQ but not explicit numbers is
 * normalised by iter_normalize_bounds, so {x, 0, 2 Pi, Pi/2} iterates the
 * lattice 0, Pi/2, ..., 2 Pi. The plotters keep plain iter_spec_parse: their
 * {x, 0, 2 Pi} is a continuous interval, and normalising it would clip the
 * plot at 6.
 */
bool iter_spec_parse_lattice(Expr* spec, IterSpec* out);

/* A Throw/Abort/Return the spec's evaluated bounds raised (Do[1, {i, Throw[3]}]),
 * detached from the spec -- the caller frees the spec and returns it, so it
 * reaches the enclosing Catch -- or NULL when the bounds are ordinary values. */
Expr* iter_spec_take_unwind(IterSpec* s);

/* Is the iterator variable an indexed one, v[k] (only the discrete iterators
 * accept it)? Such a variable is bound as the DownValue v[k] = value. */
bool iter_var_is_indexed(const Expr* var);

/* Bind the iterator variable to `value` for one step: an OwnValue of a symbol,
 * or the DownValue v[k] = value of an indexed variable.  Pair the loop with
 * iter_spec_shadow / iter_spec_restore, which save and restore either kind. */
void iter_bind(Expr* var, Expr* value);

/* Free every Expr* owned by `s` and zero it.  NULL-safe. */
void iter_spec_free(IterSpec* s);

/*
 * Resolve a COUNT/RANGE spec to machine doubles for in-loop iteration.
 * `allow_inf` permits an Infinity upper bound (Do/Sum: true; Table: false).
 * Returns false when a bound is non-numeric (and not an allowed Infinity),
 * or the step is zero.  LIST specs are not numeric and return false.
 * Expects a spec from iter_spec_parse_lattice: a non-explicit upper bound is
 * taken to be a lattice point and widened by half a step (see iter.c).
 */
bool iter_spec_resolve_numeric(const IterSpec* s, bool allow_inf,
                               double* min_val, double* max_val,
                               double* di_val, bool* is_real, bool* is_inf);

/*
 * Machine value of a real numeric bound: an explicit Integer / Real /
 * Rational, or any NumericQ expression whose machine value is a finite real
 * (Pi, 2 Pi, Sqrt[2], a BigInt). False for a free symbol, a complex value or
 * a non-finite one. Borrows `e`.
 */
bool iter_real_value(const Expr* e, double* out);

/*
 * Sum / Product with a single-element iterator {c} (a COUNT spec from
 * iter_spec_parse, whose imax holds the evaluated c), read as the reference
 * does:
 *   - c a real number (3, 5/2, Pi, 10^9): c copies of the summand, each
 *     evaluated separately -- *out is the range spec {dummy, 1, c} over a
 *     fresh private symbol `dummy_name`, so the caller's range engine
 *     enumerates (side effects run once per term) or takes a closed form;
 *   - c an unprotected symbol with no value (n): the indefinite form, so
 *     Sum[n^2, {n}] is Sum[n^2, n] -- *out is a copy of the symbol;
 *   - anything else ("s", -n, 3 + I, Infinity): `head`::vloc (or ::iterb for
 *     a non-real number), the call stays unevaluated.
 * Bare `Sum[f, 3]` (no braces) is not handled here and stays unevaluated.
 * `spec` is the held spec, `s` its parse; neither is consumed.
 */
typedef enum { ITER_COUNT_DECLINE, ITER_COUNT_RANGE, ITER_COUNT_INDEFINITE } IterCountKind;
IterCountKind iter_count_spec_for_sum(const Expr* spec, const IterSpec* s,
                                      const char* head, const char* dummy_name,
                                      Expr** out);

/*
 * NSum / NProduct with single-element iterators, by the same convention as
 * Sum / Product: each spec {c} among call's arguments [1, pos_end) whose c
 * evaluates to a real number becomes {dummy_i, 1, Floor[c]} over a fresh private
 * symbol `<dummy_prefix><i>`.  Returns the rewritten copy of `call` (caller
 * owns it), or NULL when no argument is a single-element spec.  A {c} whose c
 * is not real (a symbol: a numeric sum has no indefinite form) emits
 * `head`::itraw and sets *declined (the call then stays unevaluated).
 */
Expr* iter_numeric_count_specs(const Expr* call, size_t pos_end, const char* head,
                               const char* dummy_prefix, bool* declined);

/* True when the symbol named `name` (interned) occurs anywhere in `e`. */
bool iter_expr_mentions_symbol(const Expr* e, const char* name);

/*
 * Normalise a range {imin, imax, di} whose bounds are not all explicit numbers
 * but are all real numeric (Table[x, {x, 0, 2 Pi, 1.}], {x, Pi, 5},
 * {x, 0, 2 Pi, Pi/2}): replace *imax with the exact last point reached,
 * imin + n di with n = Floor[(imax - imin)/di], as Mathematica iterates.
 * Afterwards the bound lies on the lattice, so exact loops and closed forms
 * (Sum's F(imax + 1) - F(imin)) see the range they actually iterate -- a
 * closed form over the raw Pi bound would be a wrong answer. An empty range
 * gets imin - di. Returns true when *imax was replaced (the old one freed);
 * false, leaving everything untouched, for explicit or non-numeric bounds.
 */
bool iter_normalize_bounds(const Expr* imin, Expr** imax, const Expr* di);

/*
 * The lattice {imin, imin + di, ..., imin + n di} of a range whose bounds are
 * not numeric but whose length is: (imax - imin)/di evaluates to an explicit
 * real number (Range[a, a + 3], Table[i, {i, a, a + 3}], Range[0, 3 b, b]),
 * with n = Floor of it.  Returns a new List of the evaluated points, an empty
 * List for a negative length, or NULL when the helper does not apply (all
 * three bounds numeric, a zero step, a non-numeric length, or more than `cap`
 * points).  Iterating over this List is iterating the range exactly.
 */
Expr* iter_symbolic_lattice(const Expr* imin, const Expr* imax, const Expr* di, int64_t cap);

/*
 * Block-style localization of the iterator variable.  shadow() sets the
 * variable's current value aside (an OwnValue, or the v[k] DownValue of an
 * indexed variable) and returns a token; bind each iteration's value with
 * iter_bind(); restore() frees the iteration-created value and reinstates the
 * saved one.  The saved value lives on a binding frame (modular.h,
 * mth_bind_push), which a TimeConstrained timeout drains, so an interrupted
 * loop cannot leave its variable bound.  shadow() allocates nothing once the
 * frame stack has grown, and restores may come out of order.
 * `var` may be NULL (COUNT specs): both are no-ops then.
 */
typedef int IterShadow;
IterShadow iter_spec_shadow(Expr* var);
void       iter_spec_restore(Expr* var, IterShadow token);

/* curr + step for a machine-numeric iterator, without going through the
 * evaluator. Advancing by building Plus[curr, step] and calling evaluate() costs
 * three Expr allocations and a full evaluator dispatch per step, to compute i+1;
 * on `Do[Null, {i, 200000}]` that was 0.0251 s against 0.0034 s for the add.
 *
 * Returns NULL -- meaning "use the evaluator" -- for BigInt, Rational, MPFR or a
 * symbolic step, so those keep bit-identical values, and for an integer sum that
 * would overflow, leaving BigInt promotion to that path. Borrows both operands;
 * the result is owned by the caller. */
Expr* iter_step_add(const Expr* curr, const Expr* step);

/* Out-of-line slow path of iter_range_continue: an exact compare when a BigInt
 * is involved (curr_e overflowed int64 on the last advance, or a bound is big).
 * Only reached at the 2^63 boundary, so it stays off the hot loop. */
bool iter_range_continue_bigint(const Expr* curr_e, const Expr* imax_e, double di_val);

/* The rounding slack a real iterator's bound test allows, RELATIVE to the
 * iteration's scale: 1e-14 (|max| + |di|). It absorbs the drift of the
 * accumulated `val += di`. The historical absolute 1e-14 swamped a tiny range:
 * Range[1.*^-200, 2.*^-200, 1.*^-200] counted ~10^14 elements and Table over
 * it never finished (numeric-tower review E10). For a range of order 1 it is
 * the same slack as before. */
static inline double iter_real_slack(double max_val, double di_val) {
    return 1e-14 * (fabs(max_val) + fabs(di_val));
}

/* Loop-continuation test for an arithmetic-progression iterator, shared by
 * Table / Do / Sum / Product. Returns true while the loop should keep running.
 *
 * Real iteration compares the double running value `val` against `max_val` (with
 * a relative slack, iter_real_slack). EXACT integer iteration compares the exact
 * running value `curr_e` against the exact bound `imax_e` — an int64 compare, or
 * GMP when a BigInt is involved — because int64 values near 2^63 collapse to the
 * same double and the double test can no longer terminate the loop (issue #52).
 * `di_val`'s sign gives the direction; `is_inf` short-circuits to true.
 *
 * Kept inline (with only the BigInt arm out-of-line) so the common int64 loop
 * pays a single register compare per step, not a call — the microbenchmark PR
 * #50 tuned (Do[Null, {i, n}]) must not regress. */
static inline bool iter_range_continue(bool is_real, bool is_inf,
                                       const Expr* curr_e, const Expr* imax_e,
                                       double val, double max_val, double di_val) {
    if (is_inf) return true;
    if (!is_real && curr_e && imax_e) {
        if (curr_e->type == EXPR_INTEGER && imax_e->type == EXPR_INTEGER) {
            int64_t a = curr_e->data.integer, b = imax_e->data.integer;
            return (di_val > 0) ? (a <= b) : (a >= b);
        }
        if (expr_is_integer_like(curr_e) && expr_is_integer_like(imax_e))
            return iter_range_continue_bigint(curr_e, imax_e, di_val);
    }
    double slack = iter_real_slack(max_val, di_val);
    return (di_val > 0 && val <= max_val + slack)
        || (di_val < 0 && val >= max_val - slack);
}

/* The number of elements of an inexact range: imin + k di for
 * k = 0 .. Floor[(imax - imin)/di], as in the reference. The quotient is
 * computed once, with a relative slack for the rounding of imax - imin, of the
 * division and of di's own decimal reading (each ~2^-53 of
 * (|imin| + |imax|)/|di|; the slack is ~45 times that), so a bound the
 * progression reaches exactly in the reals is counted. Accumulating
 * `val += di` instead drifted by one ulp per step and outgrew a fixed slack:
 * Range[0., 100., 0.01] lost its last element (audit Range#1), and Table's
 * accumulated value drifted the same way (audit Table#5). The caller has
 * already rejected a zero or wrong-signed step and an over-large estimate. */
static inline size_t iter_real_count(double min_val, double max_val, double di_val) {
    double q = (max_val - min_val) / di_val;
    double tol = 1e-14 * (fabs(min_val) + fabs(max_val) + fabs(di_val)) / fabs(di_val);
    double k = floor(q + tol);
    return k < 0 ? 0 : (size_t)k + 1;
}

/* Element k of an inexact range: imin + k di, never past imax. Only the last
 * element can overshoot -- by the rounding the count's slack admitted -- and
 * it is then imax itself, as in the reference (Range[0., 0.3, 0.1] ends at
 * 0.3, not 0.30000000000000004). */
static inline double iter_real_at(double min_val, double max_val, double di_val, size_t k) {
    double v = min_val + (double)k * di_val;
    return (di_val > 0) ? (v > max_val ? max_val : v) : (v < max_val ? max_val : v);
}

#endif

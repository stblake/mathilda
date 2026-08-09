/* Nearest[list, x] -- the element(s) of `list` closest to the target `x`.
 *
 *   Nearest[{1, 5, 10}, 3]    -> {1, 5}
 *   Nearest[{10, 20, 30}, 100] -> {30}
 *   Nearest[{}, 3]             -> {}
 *
 * Distance is Abs[element - x], composed from the existing internal_subtract
 * and internal_abs rather than a bespoke helper, so a complex element uses its
 * modulus for free: Nearest[{3 + 4 I, 1}, 0] is {1} because the distances are
 * 5 and 1.
 *
 * ALL elements tied at the minimum distance are returned, in their original
 * order. That requirement is what fixes the algorithm. The obvious neighbour,
 * RankedMin (sort.c), selects with a quickselect whose comparator carries an
 * original-index tiebreak (sort.c:932) existing precisely to make ties
 * IMPOSSIBLE, so that exactly one element can win -- the opposite of what is
 * needed here. The shape used instead is MinimalBy's (sort.c:663-716): one pass
 * to find the minimum, a second to collect every distance equal to it. Input
 * order among ties then falls out of the ascending collect pass, so there is no
 * tie logic in this file at all.
 *
 * Nearest diverges from MinimalBy in exactly one respect, and deliberately.
 * Every distance must be a real number, or the whole call stays unevaluated.
 * MinimalBy[{1, a, 3}, Abs[# - 2] &] answers {1, 3}: expr_compare orders
 * symbols after all numbers (sort.c:379-380), so the symbolic element is never
 * minimal and silently vanishes from a result that still looks plausible.
 * Gating on the DISTANCE rather than on the element covers a symbolic element,
 * a symbolic target, and a non-real complex in a single check. A symbolic real
 * such as Pi is rejected too -- Abs[Pi - 3] stays as Abs[-3 + Pi] -- rather
 * than being numericalized the way RankedMin's ranked_numeric_key would.
 *
 * Cost: O(n) distance evaluations, O(n) comparisons, O(n) peak extra memory.
 * The two evaluate passes per element dominate, so this is an interpreter-speed
 * path and not a buffer path. There is no exact-hit short circuit: a later
 * element can tie at distance 0, and dropping it would break the tie contract.
 *
 * Only the two-argument form lives here. The n-nearest, radius, rule,
 * all-pairs, and NearestTo operator forms, and the DistanceFunction option, are
 * separate follow-ups. */

#include "list_common.h"
#include "internal.h"
#include "nearest.h"

#include <math.h>       /* isnan -- C99, needs no feature-test macro */

/* Abs[e - x], fully evaluated. Caller owns the result.
 *
 * internal_call_impl (internal.c:189-211) consumes the argument array, and when
 * the builtin declines it returns the UNEVALUATED node rather than NULL -- so a
 * symbolic operand comes back as Abs[...], which the caller's numeric gate then
 * rejects. Neither call can yield NULL, so neither needs a null check. The
 * composition and the eval_and_free wrapping follow comparisons.c:313-314,
 * which already builds a difference this way for the Equal zero test. */
static Expr* nearest_distance(Expr* e, Expr* x) {
    Expr* sub_args[2] = { expr_copy(e), expr_copy(x) };
    Expr* diff = eval_and_free(internal_subtract(sub_args, 2));
    Expr* abs_args[1] = { diff };          /* internal_abs takes ownership */
    return eval_and_free(internal_abs(abs_args, 1));
}

/* True for a real number this file can order exactly.
 *
 * Deliberately NOT list_common.h's is_real_numeric: that routes rationals
 * through is_rational (arithmetic.c:105-117), which requires int64 numerator
 * AND denominator, so an ordinary exact input like 1/10^25 -- a
 * Rational[1, BigInt] -- is rejected and the whole call declines. is_rational_like
 * (arithmetic.c:125) is the bigint-aware predicate. */
static bool nearest_is_real_number(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    if (e->type == EXPR_REAL) return !isnan(e->data.real);
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return !mpfr_nan_p(e->data.mpfr);
#endif
    return is_rational_like(e);
}

/* Sign of a real number: -1, 0, +1. Sets *ok false for anything this file
 * cannot read a sign from, which the caller turns into an unevaluated result.
 *
 * The *ok flag is the whole point. expr_numeric_sign (arithmetic.c) returns a
 * bare 0 both for "genuinely zero" and for "I do not recognise this", and it
 * recognises neither MPFR nor a bigint-component Rational -- so using it here
 * would silently report a tie for two distances it simply could not read. */
static int nearest_sign(Expr* e, bool* ok) {
    *ok = true;
    if (!e) { *ok = false; return 0; }
    if (e->type == EXPR_INTEGER)
        return (e->data.integer > 0) - (e->data.integer < 0);
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint);
    if (e->type == EXPR_REAL) {
        if (isnan(e->data.real)) { *ok = false; return 0; }
        return (e->data.real > 0.0) - (e->data.real < 0.0);
    }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        if (mpfr_nan_p(e->data.mpfr)) { *ok = false; return 0; }
        return mpfr_sgn(e->data.mpfr);
    }
#endif
    /* Rational[num, den], bigint components included. The denominator is
     * conventionally positive, but both signs are read so a hand-built
     * Rational with a negative denominator still orders correctly. */
    if (is_rational_like(e) && e->type == EXPR_FUNCTION) {
        bool ok_n = true, ok_d = true;
        int sn = nearest_sign(e->data.function.args[0], &ok_n);
        int sd = nearest_sign(e->data.function.args[1], &ok_d);
        if (!ok_n || !ok_d) { *ok = false; return 0; }
        return sn * sd;
    }
    *ok = false;
    return 0;
}

/* Order two distances by NUMERIC VALUE: -1, 0, +1. Sets *ok false when the
 * comparison cannot be decided.
 *
 * NOT expr_compare, which is a canonical total order and is wrong here in both
 * directions. It breaks a value tie between different ExprTypes on the type
 * enum (sort.c:376), so Nearest[{0, 2.0}, 1] would drop the 2.0 -- distances 1
 * and 1.0 are equal in value but Integer sorts before Real -- defeating the
 * all-ties contract this file exists to honour. And for atoms that are not both
 * integer-like it falls back to comparing get_numeric_value() doubles
 * (sort.c:372-377), so two exact rationals differing below double resolution
 * compare equal and a strictly-farther element is reported as tied.
 *
 * Subtracting instead is exact where it must be and inexact only where the
 * input already was: 1 - 1.0 is 0.0 (a real tie, as Mathematica has it), while
 * 1/3 - (1/3 + 1/10^18) is the exact Rational[-1, 10^18]. */
static int nearest_cmp(Expr* a, Expr* b, bool* ok) {
    *ok = true;
    /* Same-type fast paths for the two common cases, avoiding an allocation and
     * two evaluator passes per comparison. Both are exact for their type. */
    if (a->type == EXPR_INTEGER && b->type == EXPR_INTEGER)
        return (a->data.integer > b->data.integer) - (a->data.integer < b->data.integer);
    if (a->type == EXPR_REAL && b->type == EXPR_REAL) {
        if (isnan(a->data.real) || isnan(b->data.real)) { *ok = false; return 0; }
        return (a->data.real > b->data.real) - (a->data.real < b->data.real);
    }
    Expr* sub_args[2] = { expr_copy(a), expr_copy(b) };
    Expr* diff = eval_and_free(internal_subtract(sub_args, 2));
    int s = nearest_sign(diff, ok);
    expr_free(diff);
    return s;
}

Expr* builtin_nearest(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;

    Expr* list = res->data.function.args[0];
    Expr* x    = res->data.function.args[1];

    /* A PACKED list has already been materialised on the way in, because
     * Nearest is not on pack.c's AWARE list. A VISIBLE NDArray is not, and is
     * not a List either, so it lands here and stays unevaluated rather than
     * being silently truncated. Following RankedMin (sort.c:1014) rather than
     * MinimalBy, which accepts and preserves any head. */
    if (!is_listq(list)) return NULL;

    size_t n = list->data.function.arg_count;
    Expr** elem = list->data.function.args;

    /* Empty in, empty out -- checked before the gate, so Nearest[{}, a] is {}
     * and not unevaluated. Matches MaximalBy (sort.c:684). */
    if (n == 0) return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);

    Expr** dist = malloc(sizeof(Expr*) * n);
    if (!dist) return NULL;

    /* One pass to build the distances. A distance that is not a real number
     * means there is no definite answer: free what we built and decline. */
    for (size_t i = 0; i < n; i++) {
        dist[i] = nearest_distance(elem[i], x);
        if (!nearest_is_real_number(dist[i])) {
            for (size_t j = 0; j <= i; j++) expr_free(dist[j]);
            free(dist);
            return NULL;
        }
    }

    /* Two passes over nearest_cmp: find the minimum, then collect every
     * distance equal to it. An undecidable comparison declines the whole call
     * rather than guessing a tie. */
    bool ok = true;
    size_t best = 0;
    for (size_t i = 1; i < n && ok; i++)
        if (nearest_cmp(dist[i], dist[best], &ok) < 0 && ok) best = i;

    Expr** out = ok ? malloc(sizeof(Expr*) * n) : NULL;
    if (!out) {
        for (size_t i = 0; i < n; i++) expr_free(dist[i]);
        free(dist);
        return NULL;
    }

    /* Ascending index order is what preserves input order among ties. */
    size_t nout = 0;
    for (size_t i = 0; i < n && ok; i++)
        if (nearest_cmp(dist[i], dist[best], &ok) == 0 && ok)
            out[nout++] = expr_copy(elem[i]);

    for (size_t i = 0; i < n; i++) expr_free(dist[i]);
    free(dist);

    if (!ok) {
        for (size_t i = 0; i < nout; i++) expr_free(out[i]);
        free(out);
        return NULL;
    }

    /* The wrapper is always List. Stated explicitly because the input head is
     * necessarily List here, but the n-nearest follow-up keeps the same rule. */
    Expr* head = expr_new_symbol(SYM_List);
    Expr* result = expr_new_function(head, out, nout);
    if (!result) {                       /* OOM: expr_new_function took nothing */
        for (size_t i = 0; i < nout; i++) expr_free(out[i]);
        expr_free(head);
    }
    free(out);
    return result;
}

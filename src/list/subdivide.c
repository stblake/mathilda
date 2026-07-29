/* Subdivide — equally spaced points spanning an interval, endpoints included.
 *
 * Mathematica semantics:
 *
 *   Subdivide[n]             n + 1 points spanning 0 to 1:
 *                            {0, 1/n, 2/n, ..., 1}
 *   Subdivide[max, n]        n + 1 points spanning 0 to max
 *   Subdivide[min, max, n]   n + 1 points spanning min to max
 *
 * Every form yields n + 1 points, because n counts the *parts* the interval is
 * cut into, not the points produced. Point i (0-based) is
 *
 *     min + i (max - min) / n
 *
 * DESCENDING INTERVALS
 *
 * The formula above is applied verbatim for every point, with no special case.
 * When max < min the quantity (max - min) is negative, so the step is negative
 * and the points descend: Subdivide[3, 1, 4] -> {3, 5/2, 2, 3/2, 1}. Nothing
 * here branches on the sign of the step, compares min against max, sorts, or
 * takes an absolute value, so a descending call returns n + 1 points exactly
 * as an ascending one does. A degenerate interval (min == max) has step 0 and
 * simply repeats the endpoint: Subdivide[5, 5, 2] -> {5, 5, 5}.
 *
 * EXACTNESS
 *
 * Two independent guarantees:
 *
 *   1. The endpoints are *copied, never computed*. Element 0 is a copy of min
 *      and element n is a copy of max, so no arithmetic — and therefore no
 *      representation change — can touch them.
 *   2. Interior points are each computed directly from their index i, never as
 *      `previous + step`. With no accumulator there is nothing for error to
 *      accumulate in.
 *
 * Interior point i is formed as the single fraction
 *
 *     (min n + i (max - min)) / n
 *
 * and reduced once by make_rational(), which divides through by the gcd,
 * normalizes the sign onto the numerator, and returns a plain Integer when the
 * reduced denominator is 1. That is what lets whole-number points print as
 * integers alongside rationals in one list:
 * Subdivide[10, 4] -> {0, 5/2, 5, 15/2, 10}.
 *
 * VALIDITY
 *
 * n must be a positive machine integer; anything else (0, negative, Rational,
 * Real, symbolic, or a bigint demanding more elements than can be built)
 * leaves the expression unevaluated, per the builtin NULL convention.
 *
 * PERFORMANCE
 *
 * The int64 fast path below is chosen only when overflow is impossible by
 * construction (see SUBDIVIDE_MAX_ENDPOINT), so it needs no per-operation
 * overflow check. Everything it declines — bigint, rational, real, or symbolic
 * endpoints — falls back to building a Plus/Times/Power expression and letting
 * the core evaluator do the arithmetic, the same strategy src/list/rescale.c
 * uses. That keeps exactness and bigint correctness without duplicating the
 * evaluator's numeric tower here.
 */

#include "list_common.h"
#include "subdivide.h"

/* Upper bound on n. Matches the generation cap already used by Range (see
 * src/list/range.c): past this the result list is not a useful object, and
 * refusing beats attempting the allocation. */
#define SUBDIVIDE_MAX_N 1000000

/* Endpoint magnitude admitted by the int64 fast path. With |min|, |max| <= 2^31
 * and n <= 10^6:
 *
 *     |min * n|         <= 2^31 * 10^6  ~= 2.1e15
 *     |i * (max - min)| <= 10^6 * 2^32  ~= 4.3e15
 *
 * so the numerator stays under ~6.4e15, comfortably inside int64's ~9.2e18.
 * Overflow is thus excluded by precondition rather than by runtime checks. */
#define SUBDIVIDE_MAX_ENDPOINT ((int64_t)2147483648LL)

/* True when the int64 fast path can compute every interior point of this
 * interval without overflow. */
static bool subdivide_fits_int64(const Expr* min_e, const Expr* max_e) {
    if (min_e->type != EXPR_INTEGER || max_e->type != EXPR_INTEGER) return false;
    int64_t lo = min_e->data.integer;
    int64_t hi = max_e->data.integer;
    if (lo < -SUBDIVIDE_MAX_ENDPOINT || lo > SUBDIVIDE_MAX_ENDPOINT) return false;
    if (hi < -SUBDIVIDE_MAX_ENDPOINT || hi > SUBDIVIDE_MAX_ENDPOINT) return false;
    return true;
}

/* General interior point: min + i (max - min) / n, evaluated by the core
 * evaluator. Handles bigint, rational, real, and symbolic endpoints uniformly
 * and keeps exact input exact. `min_e` and `max_e` are borrowed. */
static Expr* subdivide_point_general(Expr* min_e, Expr* max_e,
                                     int64_t i, int64_t n) {
    /* max - min, as Plus[max, Times[-1, min]] */
    Expr* neg_min = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(-1), expr_copy(min_e) }, 2);
    Expr* span = expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy(max_e), neg_min }, 2);

    /* 1/n, as Power[n, -1] */
    Expr* inv_n = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_new_integer(n), expr_new_integer(-1) }, 2);

    /* i (max - min) / n */
    Expr* offset = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(i), span, inv_n }, 3);

    /* min + offset */
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy(min_e), offset }, 2);

    Expr* out = evaluate(sum);
    expr_free(sum);
    return out;
}

Expr* builtin_subdivide(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return NULL;

    /* Normalize all three surface forms to one (min, max, n) triple before a
     * single point is computed. The 1- and 2-argument forms differ from the
     * 3-argument form only in which endpoints default to 0 and 1. */
    Expr* min_e;
    Expr* max_e;
    Expr* n_e;
    Expr* implied_min = NULL;   /* owned here when the form implies min = 0 */
    Expr* implied_max = NULL;   /* owned here when the form implies max = 1 */

    if (argc == 1) {
        implied_min = expr_new_integer(0);
        implied_max = expr_new_integer(1);
        min_e = implied_min;
        max_e = implied_max;
        n_e   = res->data.function.args[0];
    } else if (argc == 2) {
        /* The implied lower endpoint takes on the exactness of the supplied
         * one, so an inexact interval does not carry a lone exact 0:
         * Subdivide[1.0, 4] -> {0., 0.25, 0.5, 0.75, 1.} rather than
         * {0, 0.25, ...}. */
        max_e = res->data.function.args[0];
        implied_min = (max_e->type == EXPR_REAL) ? expr_new_real(0.0)
                                                 : expr_new_integer(0);
        min_e = implied_min;
        n_e   = res->data.function.args[1];
    } else {
        min_e = res->data.function.args[0];
        max_e = res->data.function.args[1];
        n_e   = res->data.function.args[2];
    }

    /* n must be a positive machine integer. Zero, negatives, rationals, reals,
     * symbols, and bigints all leave the call unevaluated. */
    if (n_e->type != EXPR_INTEGER ||
        n_e->data.integer <= 0 ||
        n_e->data.integer > SUBDIVIDE_MAX_N) {
        expr_free(implied_min);
        expr_free(implied_max);
        return NULL;
    }
    int64_t n = n_e->data.integer;

    size_t count = (size_t)n + 1;
    Expr** points = malloc(sizeof(Expr*) * count);
    if (!points) {
        expr_free(implied_min);
        expr_free(implied_max);
        return NULL;
    }

    /* Endpoints verbatim — copied, never computed. */
    points[0]     = expr_copy(min_e);
    points[count - 1] = expr_copy(max_e);

    /* Interior points, each derived directly from its own index. */
    if (subdivide_fits_int64(min_e, max_e)) {
        int64_t lo   = min_e->data.integer;
        int64_t span = max_e->data.integer - lo;
        int64_t base = lo * n;
        for (int64_t i = 1; i < n; i++) {
            points[i] = make_rational(base + i * span, n);
        }
    } else {
        for (int64_t i = 1; i < n; i++) {
            points[i] = subdivide_point_general(min_e, max_e, i, n);
        }
    }

    Expr* out = expr_new_function(expr_new_symbol(SYM_List), points, count);
    free(points);
    expr_free(implied_min);
    expr_free(implied_max);
    return out;
}

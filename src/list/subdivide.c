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
#include "pack.h"

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

/* EITHER endpoint being a machine Real makes EVERY point one, the other
 * endpoint included -- see common.h's note on machine-real contagion. This is
 * Mathematica's rule and it used to be got wrong here:
 *
 *     Subdivide[0, 1., 4]     {0., 0.25, 0.5, 0.75, 1.}   (was {0, 0.25, ...})
 *     Subdivide[3., 1, 4]     {3., 2.5, 2., 1.5, 1.}      (was {..., 1})
 *     Subdivide[1/2, 1., 4]   {0.5, 0.625, 0.75, 0.875, 1.}
 *
 * The old answer kept the exact endpoint because endpoints are copied rather
 * than computed -- true, and beside the point: WHICH exact value gets copied is
 * decided by the interval's exactness, not by the copying.
 *
 * Only machine Real is contagious. An MPFR endpoint keeps the exact one
 * (`Subdivide[0, 1.\`30, 4]` starts at `0`), and a symbolic endpoint keeps the
 * general path, so both decline here. Writes the two endpoints as doubles on
 * success. */
static bool subdivide_fits_real(const Expr* min_e, const Expr* max_e,
                                double* lo, double* hi) {
    if (min_e->type != EXPR_REAL && max_e->type != EXPR_REAL) return false;
    return common_machine_real_value(min_e, lo) &&
           common_machine_real_value(max_e, hi);
}

/* Interior point i for two Real endpoints, in doubles: min + i * step, with
 * step = (max - min) / n.
 *
 * THIS IS A DELIBERATE CHANGE OF THE LAST BIT, and the reason to write it down
 * is that the previous value was not a choice. subdivide_point_general builds
 * Times[i, span, Power[n, -1]] and hands it to the evaluator; Times is
 * Orderless, so its factors are sorted CANONICALLY -- by value -- and then
 * folded left to right, turning exact into Real at whichever factor the Real
 * span happens to sort after. Subdivide[0., 1., 10] came out as 0.3 for i = 3
 * because 1/10 sorted first and (1/10 * 1.) * 3 rounds that way; a different
 * interval folds in a different order. The old result was an artifact of
 * argument ordering inside a Times node, not a documented rounding rule, and
 * it is not reproducible without re-running the evaluator -- which is the cost
 * this path exists to remove.
 *
 * min + i * step is the rule BOTH reference systems use: Mathematica gives
 * 0.30000000000000004 for that element, and NumPy's linspace computes
 * `start + arange(num) * ((stop - start) / div)` and gives the same. So the
 * change moves toward both, and the packed and unpacked paths here compute it
 * identically, which is what N2 actually requires.
 *
 * The ENDPOINTS are not computed by either path: element 0 is a copy of min and
 * element n a copy of max. That is the exactness guarantee at the top of this
 * file, not an optimisation -- min + n * step is not max in floating point. */
static double subdivide_point_real(double lo, double step, int64_t i) {
    return lo + (double)i * step;
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

    double real_lo = 0.0, real_hi = 0.0;
    bool is_real = subdivide_fits_real(min_e, max_e, &real_lo, &real_hi);
    double real_step = is_real ? (real_hi - real_lo) / (double)n : 0.0;

    /* THE BUFFER PATH. A machine-real interval makes every element a machine
     * real, so the whole result is one float64 array and no Expr node need be
     * built at all. Without it, Subdivide[0., 1., 999999] built 10^6
     * Plus/Times/Power trees and called the evaluator on each: 1.87 s against
     * NumPy's 1.11 ms for np.linspace (HPC plan item C.5).
     *
     * ndbuild_open_f64 declining -- packing off, or a count under the
     * threshold -- falls through to the List path below, which computes the
     * same doubles directly rather than through the evaluator, so a short
     * Subdivide gets most of the win too. */
    if (is_real) {
        double* buf = NULL;
        Expr* packed = ndbuild_open_f64((int64_t)count, &buf);
        if (packed) {
            buf[0]         = real_lo;             /* endpoints not computed */
            buf[count - 1] = real_hi;
            for (int64_t i = 1; i < n; i++)
                buf[i] = subdivide_point_real(real_lo, real_step, i);
            expr_free(implied_min);
            expr_free(implied_max);
            return packed;
        }
    }

    Expr** points = malloc(sizeof(Expr*) * count);
    if (!points) {
        expr_free(implied_min);
        expr_free(implied_max);
        return NULL;
    }

    /* Endpoints not computed -- taken from the input, at the interval's own
     * exactness. A machine-real interval carries an exact endpoint across as a
     * Real (Subdivide[0, 1., 4] starts at 0., not 0); everything else copies. */
    if (is_real) {
        points[0]         = expr_new_real(real_lo);
        points[count - 1] = expr_new_real(real_hi);
    } else {
        points[0]         = expr_copy(min_e);
        points[count - 1] = expr_copy(max_e);
    }

    /* Interior points, each derived directly from its own index. */
    if (subdivide_fits_int64(min_e, max_e)) {   /* disjoint from is_real */
        int64_t lo   = min_e->data.integer;
        int64_t span = max_e->data.integer - lo;
        int64_t base = lo * n;
        for (int64_t i = 1; i < n; i++) {
            points[i] = make_rational(base + i * span, n);
        }
    } else if (is_real) {
        /* Same doubles as the buffer above -- reached when the count is under
         * the packing threshold, or packing is off. */
        for (int64_t i = 1; i < n; i++) {
            points[i] = expr_new_real(subdivide_point_real(real_lo, real_step, i));
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

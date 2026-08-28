#include "list_common.h"
#include "pad_schemes.h"

/* =====================================================================
 * Shared 1-D padding-scheme engine (see pad_schemes.h).
 *
 * The value-dependent schemes build symbolic arithmetic trees
 * (Plus/Times/Subtract) from copies of the fiber elements and reduce them
 * with the evaluator, so integers, reals, bignums, symbols, and whole
 * sub-arrays (via Listable Plus/Times) all combine correctly.
 * ===================================================================== */

/* Floored modulo: result has the sign of the divisor m (> 0 here). */
static int64_t floormod(int64_t a, int64_t m) {
    int64_t r = a % m;
    if (r < 0) r += m;
    return r;
}

/* ---- small symbolic-arithmetic builders (each consumes its owned args) ---- */

/* -x  (x consumed) -> Times[-1, x], evaluated. */
static Expr* ps_neg(Expr* x) {
    Expr** a = malloc(sizeof(Expr*) * 2);
    a[0] = expr_new_integer(-1);
    a[1] = x;
    Expr* t = expr_new_function(expr_new_symbol(SYM_Times), a, 2);
    free(a);
    return eval_and_free(t);
}

/* Expand[x]  (x consumed) -> distribute products over sums, evaluated.  The
 * extrapolation combos build Times[k, Plus[...]] terms, which the evaluator
 * leaves undistributed (2 (2a-b), not 4a-2b); Wolfram's ArrayPad shows the
 * expanded polynomial, so match it. */
static Expr* ps_expand(Expr* x) {
    Expr** a = malloc(sizeof(Expr*) * 1);
    a[0] = x;
    Expr* e = expr_new_function(expr_new_symbol(intern_symbol("Expand")), a, 1);
    free(a);
    return eval_and_free(e);
}

/* 2*edge - inner  (edge borrowed, inner consumed), evaluated. */
static Expr* ps_two_edge_minus(const Expr* edge, Expr* inner) {
    Expr** ma = malloc(sizeof(Expr*) * 2);
    ma[0] = expr_new_integer(2);
    ma[1] = expr_copy((Expr*)edge);
    Expr* two_edge = expr_new_function(expr_new_symbol(SYM_Times), ma, 2);
    free(ma);
    Expr** sa = malloc(sizeof(Expr*) * 2);
    sa[0] = two_edge;
    sa[1] = inner;
    Expr* sub = expr_new_function(expr_new_symbol(SYM_Subtract), sa, 2);
    free(sa);
    return eval_and_free(sub);
}

/* ---- scheme classification --------------------------------------------- */

static bool name_is(const Expr* s, const char* name) {
    if (s->type == EXPR_STRING) return strcmp(s->data.string, name) == 0;
    if (s->type == EXPR_SYMBOL) return strcmp(s->data.symbol.name, name) == 0;
    return false;
}

PadScheme pad_scheme_classify(const Expr* s) {
    if (!s) return PS_CONSTANT;                 /* default fill 0 */
    if (s->type == EXPR_STRING || s->type == EXPR_SYMBOL) {
        /* A bare symbol that is NOT a known scheme name is still a valid
         * constant fill (e.g. ArrayPad[{1,2,3},2,x]).  Only STRINGS (and the
         * scheme-named symbols) select a named scheme. */
        static const struct { const char* n; PadScheme s; } names[] = {
            {"Fixed", PS_FIXED}, {"Periodic", PS_PERIODIC},
            {"Reflected", PS_REFLECTED}, {"Reversed", PS_REVERSED},
            {"ReversedNegation", PS_REVERSEDNEG},
            {"ReflectedDifferences", PS_REFLECTEDDIFF},
            {"ReversedDifferences", PS_REVERSEDDIFF},
            {"Extrapolated", PS_EXTRAPOLATED},
        };
        for (size_t i = 0; i < sizeof(names) / sizeof(names[0]); i++)
            if (name_is(s, names[i].n)) return names[i].s;
        /* Unknown symbol -> treat as a constant fill; unknown string -> invalid. */
        return (s->type == EXPR_SYMBOL) ? PS_CONSTANT : PS_INVALID;
    }
    if (s->type == EXPR_FUNCTION && s->data.function.head->type == EXPR_SYMBOL
        && s->data.function.head->data.symbol.name == SYM_List) {
        return PS_CYCLIC;
    }
    return PS_CONSTANT;                          /* any other scalar */
}

bool pad_scheme_value_dependent(PadScheme sc) {
    return sc != PS_CONSTANT && sc != PS_CYCLIC && sc != PS_INVALID;
}

bool pad_scheme_needs_two(PadScheme sc) {
    return sc == PS_REFLECTED || sc == PS_REFLECTEDDIFF || sc == PS_REVERSEDDIFF;
}

/* ---- direct-index schemes: value at signed offset t --------------------- *
 * For the non-arithmetic value schemes, a pad element at offset t (t < 0 or
 * t >= L; the caller never asks for 0 <= t < L here) is a copy (possibly
 * negated) of some original element.  Returns an owned Expr*. */
static Expr* ps_direct_at(const Expr** v, size_t L, PadScheme sc, int64_t t) {
    int64_t Li = (int64_t)L;
    switch (sc) {
    case PS_FIXED:
        return expr_copy((Expr*)v[t < 0 ? 0 : Li - 1]);
    case PS_PERIODIC:
        return expr_copy((Expr*)v[floormod(t, Li)]);
    case PS_REVERSED: {
        int64_t p = floormod(t, 2 * Li);
        return expr_copy((Expr*)v[p < Li ? p : 2 * Li - 1 - p]);
    }
    case PS_REFLECTED: {
        int64_t p = floormod(t, 2 * Li - 2);
        return expr_copy((Expr*)v[p < Li ? p : 2 * Li - 2 - p]);
    }
    case PS_REVERSEDNEG: {
        int64_t p = floormod(t, 2 * Li);
        if (p < Li) return expr_copy((Expr*)v[p]);
        return ps_neg(expr_copy((Expr*)v[2 * Li - 1 - p]));
    }
    default:
        return NULL;
    }
}

/* ---- antisymmetric fold schemes (ReflectedDifferences/ReversedDifferences)
 * w(t) for any t, folded into [0, L-1] applying `2*edge - x` at each boundary
 * reflection.  `repeat` selects the Reversed (boundary duplicated) variant.
 * Requires L >= 2 (guaranteed by the caller's mindimsize guard). */
static Expr* ps_fold_diff(const Expr** v, size_t L, int64_t t, bool repeat) {
    int64_t Li = (int64_t)L;
    if (t >= 0 && t < Li) return expr_copy((Expr*)v[t]);
    if (t >= Li) {
        int64_t tp = repeat ? (2 * Li - 1 - t) : (2 * (Li - 1) - t);
        return ps_two_edge_minus(v[Li - 1], ps_fold_diff(v, L, tp, repeat));
    }
    /* t < 0 */
    int64_t tp = repeat ? (-1 - t) : (-t);
    return ps_two_edge_minus(v[0], ps_fold_diff(v, L, tp, repeat));
}

/* ---- Extrapolated: one extrapolated value from a nearest-first window ---- *
 * window[0..d] are the d+1 neighbours ordered nearest-first.  Result is
 * Sum_{i=1..d+1} (-1)^(i+1) C(d+1,i) window[i-1].  Returns NULL if a binomial
 * coefficient would overflow int64 (very high order). */
static Expr* ps_extrap_combo(const Expr** window, int d) {
    int m = d + 1;                              /* number of terms */
    int64_t c = 1;                              /* C(m, 0) = 1 */
    Expr** terms = malloc(sizeof(Expr*) * (size_t)m);
    for (int i = 1; i <= m; i++) {
        /* C(m,i) = C(m,i-1) * (m-i+1) / i, exact */
        if (c > INT64_MAX / (m - i + 1)) { /* overflow guard */
            for (int j = 0; j < i - 1; j++) expr_free(terms[j]);
            free(terms);
            return NULL;
        }
        c = c * (m - i + 1) / i;
        int64_t coeff = (i % 2 == 1) ? c : -c;  /* (-1)^(i+1) */
        Expr** ta = malloc(sizeof(Expr*) * 2);
        ta[0] = expr_new_integer(coeff);
        ta[1] = expr_copy((Expr*)window[i - 1]);
        terms[i - 1] = expr_new_function(expr_new_symbol(SYM_Times), ta, 2);
        free(ta);
    }
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)m);
    free(terms);
    return ps_expand(eval_and_free(sum));
}

/* ---- the public engine -------------------------------------------------- */

bool pad_scheme_extend(const Expr** v, size_t L, int64_t front, int64_t back,
                       const Expr* scheme, PadScheme sc, int64_t order,
                       Expr*** out, size_t* out_len) {
    if (sc == PS_INVALID) return false;
    if (front < 0) front = 0;
    if (back < 0) back = 0;
    size_t N = (size_t)front + L + (size_t)back;
    if (pad_scheme_value_dependent(sc) && (front > 0 || back > 0) && L == 0)
        return false;                           /* nothing to derive from */
    if (pad_scheme_needs_two(sc) && (front > 0 || back > 0) && L < 2)
        return false;                           /* mindimsize: caller reports */

    Expr** res = malloc(sizeof(Expr*) * (N == 0 ? 1 : N));

    /* original region */
    for (size_t j = 0; j < L; j++) res[front + j] = expr_copy((Expr*)v[j]);

    if (sc == PS_EXTRAPOLATED) {
        int d = 0;
        if (order < 0) d = (int)(L - 1);        /* Infinity */
        else d = (order < (int64_t)L - 1) ? (int)order : (int)(L - 1);
        if (d < 0) d = 0;
        if (d >= 63) {                          /* window/binomial would overflow */
            for (size_t q = front; q < front + L; q++) expr_free(res[q]);
            free(res);
            return false;
        }
        /* back: pos increasing, window = out[pos-1 .. pos-1-d] nearest-first */
        for (size_t pos = front + L; pos < N; pos++) {
            const Expr* window[64];
            for (int j = 0; j <= d; j++) window[j] = res[pos - 1 - (size_t)j];
            Expr* nv = ps_extrap_combo(window, d);
            if (!nv) { for (size_t q = front; q < pos; q++) expr_free(res[q]); free(res); return false; }
            res[pos] = nv;
        }
        /* front: pos decreasing, window = out[pos+1 .. pos+1+d] nearest-first */
        for (int64_t pos = front - 1; pos >= 0; pos--) {
            const Expr* window[64];
            for (int j = 0; j <= d; j++) window[j] = res[(size_t)pos + 1 + (size_t)j];
            Expr* nv = ps_extrap_combo(window, d);
            if (!nv) { for (size_t q = (size_t)pos + 1; q < N; q++) expr_free(res[q]); free(res); return false; }
            res[(size_t)pos] = nv;
        }
    } else if (sc == PS_REFLECTEDDIFF || sc == PS_REVERSEDDIFF) {
        bool repeat = (sc == PS_REVERSEDDIFF);
        for (int64_t i = 0; i < front; i++)
            res[i] = ps_fold_diff(v, L, i - front, repeat);
        for (size_t i = front + L; i < N; i++)
            res[i] = ps_fold_diff(v, L, (int64_t)i - front, repeat);
    } else if (pad_scheme_value_dependent(sc)) {
        /* direct-index value schemes */
        for (int64_t i = 0; i < front; i++)
            res[i] = ps_direct_at(v, L, sc, i - front);
        for (size_t i = front + L; i < N; i++)
            res[i] = ps_direct_at(v, L, sc, (int64_t)i - front);
    } else {
        /* PS_CONSTANT / PS_CYCLIC: fill from the block, not from v. */
        for (size_t i = 0; i < N; i++) {
            if (i >= (size_t)front && i < (size_t)front + L) continue; /* original */
            int64_t t = (int64_t)i - front;
            int64_t coord = (t >= (int64_t)L) ? (t - (int64_t)L) : t;   /* back starts at 0 */
            if (sc == PS_CYCLIC) {
                int64_t len = (int64_t)scheme->data.function.arg_count;
                if (len == 0) { res[i] = expr_new_integer(0); }
                else res[i] = expr_copy(scheme->data.function.args[floormod(coord, len)]);
            } else {
                res[i] = scheme ? expr_copy((Expr*)scheme) : expr_new_integer(0);
            }
        }
    }

    *out = res;
    *out_len = N;
    return true;
}

/* Distance functions on numeric vectors.
 *
 *   EuclideanDistance[u, v]         Sqrt[Sum Abs[u_i - v_i]^2]
 *   SquaredEuclideanDistance[u, v]  Sum Abs[u_i - v_i]^2
 *   ManhattanDistance[u, v]         Sum Abs[u_i - v_i]
 *   CosineDistance[u, v]            1 - (u . Conjugate[v]) / (Norm[u] Norm[v])
 *
 * These existed only as *names* before: FindClusters' DistanceFunction option
 * validated the symbol and then ignored it, which is why its documentation could
 * truthfully say all four choices give the same partition -- in 1D they do, and
 * nothing was ever computed. `Names["*Distance*"]` returned only GraphDistance.
 * n-dimensional clustering needs the real thing, and the functions are worth
 * having on their own.
 *
 * Every one is composed from internal_subtract / internal_abs / internal_power /
 * internal_plus / internal_sqrt rather than reimplementing arithmetic, following
 * nearest.c. Three properties fall out of that and are not accidents:
 *
 *   1. EXACT INPUT STAYS EXACT where the result is rational.
 *      SquaredEuclideanDistance[{1, 2}, {4, 6}] is 25, not 25., and
 *      SquaredEuclideanDistance[{1/3, 0}, {0, 1/7}] is 58/441. This is what lets
 *      the clustering methods that merely RANK distances stay exact in n
 *      dimensions: squared Euclidean is monotone in Euclidean, so ranking on the
 *      square gives the same order without ever taking an irrational root.
 *   2. COMPLEX VECTORS WORK, because Abs of a complex number is its modulus.
 *      Following Mathematica, the definition is Abs-then-square rather than
 *      square-then-sum, which differ for complex components.
 *   3. SYMBOLIC INPUT SURVIVES rather than being rejected. ManhattanDistance[
 *      {a}, {b}] is Abs[a - b]; Mathematica answers the same way. Callers that
 *      need a number (FindClusters) gate on the result themselves.
 *
 * Shape rules: both arguments scalars, or both Lists of the same length. A
 * length mismatch, or a component that is itself a List, leaves the call
 * unevaluated -- these are vector functions, and a silently threaded answer for
 * matrix input would be worse than no answer. */

#include "list_common.h"
#include "internal.h"
#include "distance.h"

/* Shape check, in two strictnesses.
 *
 * `dist_shape_general` accepts both scalars, or two compound expressions sharing
 * ONE head with equal arity and no list-valued component. Sharing the head is
 * what lets RGBColor[r, g, b] and friends be treated as points without any
 * colour-specific code: a colour is an EXPR_FUNCTION carrying numeric arguments,
 * which is the same shape as a vector.
 *
 * `dist_shape` is the stricter rule the public builtins use -- scalars or Lists
 * only. Their documented domain is numeric vectors, and answering
 * EuclideanDistance[RGBColor[...], RGBColor[...]] with a number would be
 * asserting a colour-space convention that has not been checked against
 * Mathematica. The clustering path uses the general form deliberately.
 *
 * `*n` is the component count (1 for the scalar case). */
static bool dist_shape_general(Expr* u, Expr* v, size_t* n) {
    bool cu = u->type == EXPR_FUNCTION, cv = v->type == EXPR_FUNCTION;
    if (cu != cv) return false;
    if (!cu) {
        *n = 1;
        return true;
    }
    if (!expr_eq(u->data.function.head, v->data.function.head)) return false;
    size_t nu = u->data.function.arg_count;
    if (nu != v->data.function.arg_count) return false;
    if (nu == 0) return false;                 /* no distance between points of no dimension */
    for (size_t i = 0; i < nu; i++)
        if (is_listq(u->data.function.args[i]) || is_listq(v->data.function.args[i]))
            return false;                      /* matrix input: decline, don't thread */
    *n = nu;
    return true;
}

static bool dist_shape(Expr* u, Expr* v) {
    bool lu = is_listq(u), lv = is_listq(v);
    if (lu != lv) return false;
    if (!lu) return u->type != EXPR_FUNCTION;  /* a scalar, not some other head */
    return true;
}

/* Component i, borrowed. Treats a scalar as a 1-vector.
 *
 * Keyed on EXPR_FUNCTION rather than on List, so RGBColor[r, g, b] yields its
 * three channels. Testing is_listq here would hand back the whole colour as if
 * it were a scalar and silently compute nonsense. */
static Expr* dist_part(Expr* e, size_t i) {
    return (e->type == EXPR_FUNCTION) ? e->data.function.args[i] : e;
}

/* Absolute value. Takes ownership of `e` and returns a new expression.
 *
 * NOT simply internal_abs, because Abs declines on a Rational whose numerator or
 * denominator is a bigint: `Abs[-1/10^20]` stays unevaluated (as does `Sign` of
 * it), while `Abs[-1/7]` is fine. Composing internal_abs blindly therefore left
 * every distance between exact vectors at that scale unevaluated, and
 * FindClusters declined data it should have clustered -- exactly the
 * high-precision case the exact path exists for.
 *
 * For a real argument the sign is decided by list_numeric_sign, which reads
 * bigint rational components correctly, and the value is negated if needed.
 * Complex and symbolic arguments still go through internal_abs, which is what
 * supplies the modulus. This routes around the Abs gap rather than fixing it;
 * the underlying bug is in Abs itself and is worth fixing separately. */
static Expr* dist_abs(Expr* e) {
    if (list_real_number_q(e)) {
        bool ok = true;
        int s = list_numeric_sign(e, &ok);
        if (ok) {
            if (s >= 0) return e;
            Expr* neg[2] = { expr_new_integer(0), e };
            return eval_and_free(internal_subtract(neg, 2));
        }
    }
    Expr* ab[1] = { e };
    return eval_and_free(internal_abs(ab, 1));
}

/* Sum_i Abs[u_i - v_i]^p, for p == 1 or p == 2. Caller owns the result.
 *
 * `v == NULL` means "no subtraction": the sum becomes Sum_i Abs[u_i]^p, which is
 * what the cosine norms need. Sharing the loop keeps one definition of the
 * Abs-then-power convention that makes complex components use their modulus.
 *
 * Neither internal_subtract nor internal_abs can return NULL: internal_call_impl
 * hands back the unevaluated node when a builtin declines, so a symbolic operand
 * comes back as Abs[...] and flows into the sum. */
static Expr* dist_sum(Expr* u, Expr* v, int p, size_t n) {
    Expr** terms = malloc(sizeof(Expr*) * n);
    if (!terms) return NULL;

    for (size_t i = 0; i < n; i++) {
        Expr* d;
        if (v) {
            Expr* sub[2] = { expr_copy(dist_part(u, i)), expr_copy(dist_part(v, i)) };
            d = eval_and_free(internal_subtract(sub, 2));
        } else {
            d = expr_copy(dist_part(u, i));
        }
        /* Squaring a real number makes the absolute value redundant, so skip it:
         * one less node, and it cannot meet the Abs gap described above. */
        Expr* t = (p == 2 && list_real_number_q(d)) ? d : dist_abs(d);
        if (p == 2) {
            Expr* pw[2] = { t, expr_new_integer(2) };
            t = eval_and_free(internal_power(pw, 2));
        }
        if (!t) {                              /* OOM mid-way: unwind what we hold */
            for (size_t j = 0; j < i; j++) expr_free(terms[j]);
            free(terms);
            return NULL;
        }
        terms[i] = t;
    }

    /* internal_plus consumes the terms; with one term it is the identity, which
     * is what makes the scalar case fall out with no special path. */
    Expr* sum = eval_and_free(internal_plus(terms, n));
    free(terms);
    return sum;
}

/* Sum_i u_i Conjugate[v_i] -- the cosine numerator. Conjugate is what makes the
 * definition correct for complex vectors and is a no-op on reals; Mathematica
 * writes it the same way. Not internal_dot, which is undefined on two scalars. */
static Expr* dist_dot_conj(Expr* u, Expr* v, size_t n) {
    Expr** terms = malloc(sizeof(Expr*) * n);
    if (!terms) return NULL;

    for (size_t i = 0; i < n; i++) {
        Expr* cj[1] = { expr_copy(dist_part(v, i)) };
        Expr* c = eval_and_free(internal_conjugate(cj, 1));
        Expr* mul[2] = { expr_copy(dist_part(u, i)), c };
        Expr* t = eval_and_free(internal_times(mul, 2));
        if (!t) {
            for (size_t j = 0; j < i; j++) expr_free(terms[j]);
            free(terms);
            return NULL;
        }
        terms[i] = t;
    }
    Expr* sum = eval_and_free(internal_plus(terms, n));
    free(terms);
    return sum;
}

/* Sqrt[Sum Abs[u_i]^2] -- the Euclidean norm of one point. */
static Expr* dist_norm(Expr* u, size_t n) {
    Expr* s = dist_sum(u, NULL, 2, n);
    if (!s) return NULL;
    Expr* sq[1] = { s };
    return eval_and_free(internal_sqrt(sq, 1));
}

/* True when `e` is a real number we can see is exactly zero. Deliberately
 * conservative: a symbolic norm answers false, so the symbolic form is built
 * instead of a shortcut being taken on a value we cannot evaluate. */
static bool dist_is_zero(Expr* e) {
    if (!list_real_number_q(e)) return false;
    bool ok = true;
    int s = list_numeric_sign(e, &ok);
    return ok && s == 0;
}

/* Exact pairwise distances for other modules -- FindClusters ranks n-dimensional
 * points with these. Shared rather than reimplemented so there is one definition
 * of each distance in the system. Both return NULL when the shapes disagree, so
 * a caller can treat that as "not a pair of comparable points".
 *
 * Squared Euclidean is the one FindClusters ranks on: it is rational for
 * rational input and monotone in the true distance, so it orders points
 * identically without introducing a root. */
Expr* distance_squared_euclidean(Expr* u, Expr* v) {
    size_t n = 0;
    if (!dist_shape_general(u, v, &n)) return NULL;
    return dist_sum(u, v, 2, n);
}

Expr* distance_manhattan(Expr* u, Expr* v) {
    size_t n = 0;
    if (!dist_shape_general(u, v, &n)) return NULL;
    return dist_sum(u, v, 1, n);
}

/* 1 - (u . Conjugate[v]) / (Norm[u] Norm[v]).
 *
 * Ranges over [0, 2]: 0 for parallel, 1 for orthogonal, 2 for antiparallel.
 * Unlike the Euclidean family this is NOT a metric (it violates the triangle
 * inequality) and there is no squared form that ranks identically, so callers
 * use it directly.
 *
 * A zero vector on either side gives 0, following Mathematica --
 * CosineDistance[{0, 0}, {1, 2}] is 0 there. That is a convention, not a
 * derivation: the quotient is 0/0, so without the special case the result would
 * be Indeterminate. Checked against wolframscript for {0,0} against {0,0},
 * {1,2} and the reverse. */
Expr* distance_cosine(Expr* u, Expr* v) {
    size_t n = 0;
    if (!dist_shape_general(u, v, &n)) return NULL;

    Expr* nu = dist_norm(u, n);
    Expr* nv = dist_norm(v, n);
    if (!nu || !nv) { expr_free(nu); expr_free(nv); return NULL; }

    if (dist_is_zero(nu) || dist_is_zero(nv)) {
        expr_free(nu); expr_free(nv);
        return expr_new_integer(0);
    }

    Expr* num = dist_dot_conj(u, v, n);
    if (!num) { expr_free(nu); expr_free(nv); return NULL; }

    Expr* dm[2] = { nu, nv };
    Expr* den = eval_and_free(internal_times(dm, 2));
    Expr* dv[2] = { num, den };
    Expr* q = eval_and_free(internal_divide(dv, 2));
    if (!q) return NULL;

    Expr* sub[2] = { expr_new_integer(1), q };
    return eval_and_free(internal_subtract(sub, 2));
}

/* ------------------------------------------------------------------------- */
/* Sequence distances                                                         */
/* ------------------------------------------------------------------------- */

/* Levenshtein distance: the fewest single-element insertions, deletions and
 * substitutions turning one sequence into the other.
 *
 * Two rows of the DP table rather than the full matrix, so O(min(m, n)) memory
 * and O(m*n) time. Returns -1 on allocation failure.
 *
 * Comparison is by expr_eq on elements, which makes the same routine serve
 * strings (compared character by character) and lists of arbitrary expressions,
 * matching Mathematica: EditDistance[{1, 2, 3}, {1, 3}] is 1 there. */
static long dist_levenshtein(Expr** a, size_t m, Expr** b, size_t n) {
    if (m == 0) return (long)n;
    if (n == 0) return (long)m;

    size_t* prev = malloc(sizeof(size_t) * (n + 1));
    size_t* cur  = malloc(sizeof(size_t) * (n + 1));
    if (!prev || !cur) { free(prev); free(cur); return -1; }

    for (size_t j = 0; j <= n; j++) prev[j] = j;
    for (size_t i = 1; i <= m; i++) {
        cur[0] = i;
        for (size_t j = 1; j <= n; j++) {
            size_t sub = prev[j - 1] + (expr_eq(a[i - 1], b[j - 1]) ? 0 : 1);
            size_t del = prev[j] + 1;
            size_t ins = cur[j - 1] + 1;
            size_t best = sub < del ? sub : del;
            cur[j] = best < ins ? best : ins;
        }
        size_t* t = prev; prev = cur; cur = t;
    }
    long r = (long)prev[n];
    free(prev);
    free(cur);
    return r;
}

/* Explode a String into one Expr per character, or borrow a List's elements.
 * `*owned` says whether the caller must free the array's elements.
 *
 * Bytes, not code points: a multi-byte UTF-8 character counts as several
 * elements. That matches the ASCII and DNA cases these functions are used for
 * and is noted in the docstring rather than silently assumed. */
static Expr** dist_seq(Expr* e, size_t* n, bool* owned) {
    if (is_listq(e)) {
        *n = e->data.function.arg_count;
        *owned = false;
        return e->data.function.args;
    }
    if (e->type != EXPR_STRING) return NULL;

    const char* s = e->data.string;
    size_t len = strlen(s);
    Expr** out = malloc(sizeof(Expr*) * (len ? len : 1));
    if (!out) return NULL;
    for (size_t i = 0; i < len; i++) {
        char buf[2] = { s[i], '\0' };
        out[i] = expr_new_string(buf);
        if (!out[i]) {
            for (size_t j = 0; j < i; j++) expr_free(out[j]);
            free(out);
            return NULL;
        }
    }
    *n = len;
    *owned = true;
    return out;
}

static void dist_seq_free(Expr** a, size_t n, bool owned) {
    if (!owned) return;
    for (size_t i = 0; i < n; i++) expr_free(a[i]);
    free(a);
}

/* Both arguments Strings, or both Lists. Mixed shapes decline. */
static bool dist_seq_pair(Expr* u, Expr* v) {
    if (u->type == EXPR_STRING && v->type == EXPR_STRING) return true;
    return is_listq(u) && is_listq(v);
}

Expr* distance_edit(Expr* u, Expr* v) {
    if (!dist_seq_pair(u, v)) return NULL;

    size_t m = 0, n = 0;
    bool om = false, on = false;
    Expr** a = dist_seq(u, &m, &om);
    Expr** b = dist_seq(v, &n, &on);
    long r = (a && b) ? dist_levenshtein(a, m, b, n) : -1;
    dist_seq_free(a, m, om);
    dist_seq_free(b, n, on);
    return (r < 0) ? NULL : expr_new_integer(r);
}

/* Hamming distance: how many positions differ. Requires equal lengths --
 * Mathematica raises ::idim and leaves the call unevaluated otherwise, so
 * declining is the faithful behaviour. */
Expr* distance_hamming(Expr* u, Expr* v) {
    if (!dist_seq_pair(u, v)) return NULL;

    size_t m = 0, n = 0;
    bool om = false, on = false;
    Expr** a = dist_seq(u, &m, &om);
    Expr** b = dist_seq(v, &n, &on);

    Expr* out = NULL;
    if (a && b && m == n) {
        long c = 0;
        for (size_t i = 0; i < n; i++) if (!expr_eq(a[i], b[i])) c++;
        out = expr_new_integer(c);
    }
    dist_seq_free(a, m, om);
    dist_seq_free(b, n, on);
    return out;
}

Expr* builtin_edit_distance(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    return distance_edit(res->data.function.args[0], res->data.function.args[1]);
}

Expr* builtin_hamming_distance(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    return distance_hamming(res->data.function.args[0], res->data.function.args[1]);
}

/* The three builtins differ only in p and whether a root is taken. */
static Expr* dist_builtin(Expr* res, int p, bool root) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;

    Expr* u = res->data.function.args[0];
    Expr* v = res->data.function.args[1];

    size_t n = 0;
    if (!dist_shape(u, v) || !dist_shape_general(u, v, &n)) return NULL;

    Expr* sum = dist_sum(u, v, p, n);
    if (!sum) return NULL;
    if (!root) return sum;

    Expr* sq[1] = { sum };                     /* internal_sqrt takes ownership */
    return eval_and_free(internal_sqrt(sq, 1));
}

Expr* builtin_euclidean_distance(Expr* res) {
    return dist_builtin(res, 2, true);
}

Expr* builtin_squared_euclidean_distance(Expr* res) {
    return dist_builtin(res, 2, false);
}

Expr* builtin_manhattan_distance(Expr* res) {
    return dist_builtin(res, 1, false);
}

Expr* builtin_cosine_distance(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* u = res->data.function.args[0];
    Expr* v = res->data.function.args[1];
    if (!dist_shape(u, v)) return NULL;
    return distance_cosine(u, v);
}

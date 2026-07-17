#include "expand.h"
#include "eval.h"
#include "symtab.h"
#include "arithmetic.h"
#include "match.h"
#include "sym_names.h"
#include "flint_bridge.h"
#include <string.h>
#include <stdlib.h>

static bool expr_contains_patt(Expr* e, Expr* patt) {
    if (!patt) return true; // NULL pattern matches everything
    
    MatchEnv* env = env_new();
    bool is_match = match(e, patt, env);
    env_free(env);
    if (is_match) return true;
    
    if (e->type == EXPR_FUNCTION) {
        if (expr_contains_patt(e->data.function.head, patt)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (expr_contains_patt(e->data.function.args[i], patt)) return true;
        }
    }
    return false;
}

static Expr* multiply_two(Expr* a, Expr* b) {
    bool a_is_plus = (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL && a->data.function.head->data.symbol.name == SYM_Plus);
    bool b_is_plus = (b->type == EXPR_FUNCTION && b->data.function.head->type == EXPR_SYMBOL && b->data.function.head->data.symbol.name == SYM_Plus);

    if (a_is_plus && b_is_plus) {
        size_t count = a->data.function.arg_count * b->data.function.arg_count;
        Expr** args = malloc(sizeof(Expr*) * count);
        size_t k = 0;
        for (size_t i = 0; i < a->data.function.arg_count; i++) {
            for (size_t j = 0; j < b->data.function.arg_count; j++) {
                Expr* t_args[2] = { expr_copy(a->data.function.args[i]), expr_copy(b->data.function.args[j]) };
                args[k++] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), t_args, 2));
            }
        }
        Expr* res = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), args, count));
        free(args);
        return res;
    } else if (a_is_plus) {
        size_t count = a->data.function.arg_count;
        Expr** args = malloc(sizeof(Expr*) * count);
        for (size_t i = 0; i < count; i++) {
            Expr* t_args[2] = { expr_copy(a->data.function.args[i]), expr_copy(b) };
            args[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), t_args, 2));
        }
        Expr* res = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), args, count));
        free(args);
        return res;
    } else if (b_is_plus) {
        size_t count = b->data.function.arg_count;
        Expr** args = malloc(sizeof(Expr*) * count);
        for (size_t i = 0; i < count; i++) {
            Expr* t_args[2] = { expr_copy(a), expr_copy(b->data.function.args[i]) };
            args[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), t_args, 2));
        }
        Expr* res = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), args, count));
        free(args);
        return res;
    } else {
        Expr* t_args[2] = { expr_copy(a), expr_copy(b) };
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), t_args, 2));
    }
}

static Expr* multiply_all(Expr** args, size_t start, size_t end) {
    if (start == end) return expr_copy(args[start]);
    if (start + 1 == end) return multiply_two(args[start], args[end]);
    size_t mid = start + (end - start) / 2;
    Expr* left = multiply_all(args, start, mid);
    Expr* right = multiply_all(args, mid + 1, end);
    Expr* res = multiply_two(left, right);
    expr_free(left);
    expr_free(right);
    return res;
}

/* Multiply the already-expanded active factors of a Times.  When `can_flint`
 * and at least two factors are polynomials over Q, their sub-product is handed
 * to FLINT (packed fmpq_mpoly multiply) in one shot and only the remaining
 * non-polynomial factors are multiplied with the generic tree multiplier —
 * so a product like Log[x] (1+y)^45 (1+z)^45 (1+w)^45 no longer distributes the
 * three dense polynomials term-by-term.  The factors are borrowed (freed by the
 * caller); the returned Expr is fresh. */
static Expr* multiply_active(Expr** active, size_t na, bool can_flint) {
    if (na == 0) return expr_new_integer(1);
    if (na == 1) return expr_copy(active[0]);

    if (can_flint) {
        Expr** poly = malloc(sizeof(Expr*) * na);
        Expr** rest = malloc(sizeof(Expr*) * na);
        size_t np = 0, nr = 0;
        for (size_t i = 0; i < na; i++) {
            if (flint_is_polynomial_over_q(active[i])) poly[np++] = active[i];
            else                                       rest[nr++] = active[i];
        }
        if (np >= 2) {
            /* Raw (unevaluated) Times of the polynomial factors: to_mpoly reads
             * the product structure directly, so evaluating it here — which
             * would invoke the generic multiply we are trying to avoid — is
             * neither needed nor wanted. */
            Expr** cp = malloc(sizeof(Expr*) * np);
            for (size_t i = 0; i < np; i++) cp[i] = expr_copy(poly[i]);
            Expr* tmp = expr_new_function(expr_new_symbol(SYM_Times), cp, np);
            free(cp);
            Expr* P = flint_expand_polynomial(tmp);
            expr_free(tmp);
            if (P) {
                Expr* res;
                if (nr == 0) {
                    res = P;
                } else {
                    Expr** all = malloc(sizeof(Expr*) * (nr + 1));
                    all[0] = P;
                    for (size_t i = 0; i < nr; i++) all[i + 1] = rest[i];
                    res = multiply_all(all, 0, nr);
                    free(all);
                    expr_free(P);
                }
                free(poly); free(rest);
                return res;
            }
        }
        free(poly); free(rest);
    }
    return multiply_all(active, 0, na - 1);
}

static Expr* power_expand(Expr* base, int64_t exp) {
    if (exp == 0) return expr_new_integer(1);
    if (exp == 1) return expr_copy(base);
    
    Expr* res = expr_new_integer(1);
    Expr* b = expr_copy(base);
    int64_t e = exp;
    
    while (e > 0) {
        if (e % 2 == 1) {
            Expr* next_res = multiply_two(res, b);
            expr_free(res);
            res = next_res;
        }
        e /= 2;
        if (e > 0) {
            Expr* next_b = multiply_two(b, b);
            expr_free(b);
            b = next_b;
        }
    }
    expr_free(b);
    return res;
}

/* Two ceilings on the estimated monomial count of an expanded integer power.
 * Neither is a mathematical cap on degree or exponent: univariate and
 * low-dimensional expansions of ANY degree pass both, because the estimate
 * below is exact for them.  They differ only in what happens to a genuinely
 * high-dimensional multinomial that would exhaust memory:
 *
 *   EXPAND_FACTOR_CAP   — internal callers (Together, Apart, Cancel, the
 *     integrators, linear algebra, ...) reach `expr_expand` as a helper and
 *     need a *usable* expression back.  Above this modest ceiling the power is
 *     simply left factored, which is a valid, correct input for those
 *     algorithms and avoids needless combinatorial work.
 *
 *   EXPAND_OVERFLOW_CAP — the user-facing `Expand` never declines: it expands
 *     everything that plausibly fits in memory (each result term costs a few
 *     hundred bytes, so tens of millions of terms approach a multi-GB tree).
 *     Only beyond this genuine memory ceiling does it yield `Overflow[]`
 *     instead of an answer.  Heuristic and machine-independent; tune if needed. */
#define EXPAND_FACTOR_CAP   2.0e5
#define EXPAND_OVERFLOW_CAP 5.0e7

/* An upper bound on the size of an expanded (sub)expression, as a Newton-box
 * estimate: `terms` bounds the monomial count, and `deg` records, per distinct
 * kernel (indeterminate), the maximum total degree that kernel reaches.  Any
 * atom that is not a rational coefficient — a symbol, or an opaque head such as
 * Sin[x] or x^(1/2) — is treated as a kernel of degree 1, so the estimate is a
 * true upper bound for ANY expression, not just polynomials over Q. */
typedef struct {
    double     terms;
    struct KDeg { Expr* kernel; long deg; } *deg;
    size_t     n, cap;
} SizeEst;

static void se_free(SizeEst* s) { free(s->deg); s->deg = NULL; s->n = s->cap = 0; }

/* Record degree `d` for `kernel`, combining with any existing entry via `f`
 * (fmax for Plus alternatives, addition for Times/Power products). */
static void se_put(SizeEst* s, Expr* kernel, long d, long (*combine)(long, long)) {
    for (size_t i = 0; i < s->n; i++) {
        if (expr_eq(s->deg[i].kernel, kernel)) { s->deg[i].deg = combine(s->deg[i].deg, d); return; }
    }
    if (s->n == s->cap) { s->cap = s->cap ? s->cap * 2 : 8; s->deg = realloc(s->deg, s->cap * sizeof(*s->deg)); }
    s->deg[s->n].kernel = kernel;
    s->deg[s->n].deg = d;
    s->n++;
}
static long comb_max(long a, long b) { return a > b ? a : b; }
static long comb_add(long a, long b) { return a + b; }

/* Bounding-box monomial count prod_kernel (deg + 1), clamped at `cap`. */
static double se_box(const SizeEst* s, double cap) {
    double b = 1.0;
    for (size_t i = 0; i < s->n && b <= cap; i++) b *= (double)s->deg[i].deg + 1.0;
    return b;
}

static SizeEst estimate_terms(const Expr* e, double cap);

/* Term bound for base^n given the base's estimate `c`: the tighter of the
 * multinomial count C(n + m - 1, m - 1) (m = base term count, exact when the
 * base terms are independent) and the bounding box (exact for a true
 * polynomial).  Fills `r`'s degree map with the base degrees scaled by n. */
static double se_power(const SizeEst* c, int64_t n, SizeEst* r, double cap) {
    double multinom = 1.0;                        /* C(n + m - 1, m - 1) */
    for (double i = 1.0; i <= c->terms - 1.0 && multinom <= cap; i += 1.0)
        multinom *= (double)((double)n + i) / i;
    for (size_t j = 0; j < c->n; j++)
        se_put(r, c->deg[j].kernel, (long)n * c->deg[j].deg, comb_max);
    double box = se_box(r, cap);
    return multinom < box ? multinom : box;
}

/* Recursively bound the number of monomials in the fully expanded `e`.  Always
 * returns a valid upper bound: collection during real expansion only reduces
 * the count.  The `cap` clamps intermediate products so the estimate is
 * overflow-safe even for astronomically large results. */
static SizeEst estimate_terms(const Expr* e, double cap) {
    SizeEst r = { 1.0, NULL, 0, 0 };
    if (!e) return r;

    switch (e->type) {
        case EXPR_INTEGER: case EXPR_REAL: case EXPR_BIGINT:
            return r;                                  /* numeric coefficient */
        case EXPR_SYMBOL:
            se_put(&r, (Expr*)e, 1, comb_max);
            return r;
        case EXPR_FUNCTION:
            break;
        default:
            se_put(&r, (Expr*)e, 1, comb_max);         /* opaque atom */
            return r;
    }

    const char* h = (e->data.function.head->type == EXPR_SYMBOL)
        ? e->data.function.head->data.symbol.name : "";
    size_t n = e->data.function.arg_count;

    if (strcmp(h, "Rational") == 0 || strcmp(h, "Complex") == 0)
        return r;                                      /* numeric coefficient */

    if (strcmp(h, "Plus") == 0) {
        double sum = 0.0;
        for (size_t i = 0; i < n; i++) {
            SizeEst c = estimate_terms(e->data.function.args[i], cap);
            sum += c.terms; if (sum > cap) sum = cap + 1.0;
            for (size_t j = 0; j < c.n; j++) se_put(&r, c.deg[j].kernel, c.deg[j].deg, comb_max);
            se_free(&c);
        }
        double box = se_box(&r, cap);
        r.terms = sum < box ? sum : box;
        return r;
    }

    if (strcmp(h, "Times") == 0) {
        double prod = 1.0;
        for (size_t i = 0; i < n; i++) {
            SizeEst c = estimate_terms(e->data.function.args[i], cap);
            prod *= c.terms; if (prod > cap) prod = cap + 1.0;
            for (size_t j = 0; j < c.n; j++) se_put(&r, c.deg[j].kernel, c.deg[j].deg, comb_add);
            se_free(&c);
        }
        double box = se_box(&r, cap);
        r.terms = prod < box ? prod : box;
        return r;
    }

    if (strcmp(h, "Power") == 0 && n == 2) {
        const Expr* ex = e->data.function.args[1];
        if (ex->type == EXPR_INTEGER && ex->data.integer > 0) {
            SizeEst c = estimate_terms(e->data.function.args[0], cap);
            r.terms = se_power(&c, ex->data.integer, &r, cap);
            se_free(&c);
            return r;
        }
        if (ex->type == EXPR_INTEGER && ex->data.integer == 0)
            return r;                                  /* x^0 -> 1 */
        /* negative / fractional / symbolic exponent: opaque kernel */
        se_put(&r, (Expr*)e, 1, comb_max);
        return r;
    }

    /* any other head (Sin, Log, ...): opaque kernel */
    se_put(&r, (Expr*)e, 1, comb_max);
    return r;
}

/* Upper bound on the expanded monomial count of `e`, clamped at `cap`. */
static double estimate_expand_terms(const Expr* e, double cap) {
    SizeEst s = estimate_terms(e, cap);
    double t = s.terms;
    se_free(&s);
    return t;
}

/* Core expander.  `overflow_mode` selects the behavior for a power whose
 * expansion is estimated too large: user-facing Expand (true) uses the memory
 * ceiling and yields Overflow[]; internal helpers (false) use the modest factor
 * ceiling and leave the power factored.  The flag is threaded through every
 * recursive call so an oversized subexpression is handled consistently; because
 * Plus/Times absorb Overflow[], an oversized part collapses the whole tree to
 * Overflow[] in user-facing mode. */
static Expr* expr_expand_impl(Expr* e, Expr* patt, bool overflow_mode) {
    if (!e) return NULL;
    if (patt && !expr_contains_patt(e, patt)) return expr_copy(e);

    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    const char* head = e->data.function.head->type == EXPR_SYMBOL ? e->data.function.head->data.symbol.name : "";

    /* Fast path: for a plain (single-argument) Expand of an arithmetic head,
     * size-check the result and hand the whole expansion to FLINT, which
     * distributes and collects in packed fmpq_mpoly arithmetic — vastly faster
     * than the generic Expr multiply on dense, high-degree inputs.  The estimate
     * is an upper bound on the monomial count:
     *   - above the ceiling, user-facing Expand yields Overflow[] (it never
     *     silently declines); an internal helper instead falls through to the
     *     classical path, which leaves an oversized power factored — a valid
     *     input for its caller.
     *   - within the ceiling, try FLINT; if it declines (not a polynomial over
     *     Q — a transcendental head, symbolic exponent, inexact coefficient),
     *     fall through to the classical distributor, which still expands the
     *     polynomial parts of a mixed expression via this same fast path on
     *     recursion. */
    double cap = overflow_mode ? EXPAND_OVERFLOW_CAP : EXPAND_FACTOR_CAP;
    if (patt == NULL &&
        (strcmp(head, "Plus") == 0 || strcmp(head, "Times") == 0 || strcmp(head, "Power") == 0)) {
        double est = estimate_expand_terms(e, cap);
        if (est > cap) {
            if (overflow_mode)
                return expr_new_function(expr_new_symbol(SYM_Overflow), NULL, 0);
            /* internal mode: fall through to leave oversized powers factored */
        } else if (flint_bridge_available()) {
            Expr* fast = flint_expand_polynomial(e);
            if (fast) return fast;
        }
    }

    // Thread over lists, equations, inequalities, logic. Inequality has
    // operator-symbol slots at odd indices that must be passed through.
    if (strcmp(head, "List") == 0 || strcmp(head, "Equal") == 0 || strcmp(head, "Less") == 0 ||
        strcmp(head, "LessEqual") == 0 || strcmp(head, "Greater") == 0 || strcmp(head, "GreaterEqual") == 0 ||
        strcmp(head, "Inequality") == 0 ||
        strcmp(head, "And") == 0 || strcmp(head, "Or") == 0 || strcmp(head, "Not") == 0) {
        bool is_ineq = (strcmp(head, "Inequality") == 0);
        Expr** args = malloc(sizeof(Expr*) * e->data.function.arg_count);
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (is_ineq && (i & 1u) == 1) {
                args[i] = expr_copy(e->data.function.args[i]);
            } else {
                args[i] = expr_expand_impl(e->data.function.args[i], patt, overflow_mode);
            }
        }
        Expr* res = eval_and_free(expr_new_function(expr_copy(e->data.function.head), args, e->data.function.arg_count));
        free(args);
        return res;
    }

    if (strcmp(head, "Plus") == 0) {
        Expr** args = malloc(sizeof(Expr*) * e->data.function.arg_count);
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            args[i] = expr_expand_impl(e->data.function.args[i], patt, overflow_mode);
        }
        Expr* res = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), args, e->data.function.arg_count));
        free(args);
        return res;
    }

    if (strcmp(head, "Times") == 0) {
        size_t count = e->data.function.arg_count;
        if (count == 0) return expr_new_integer(1);

        /* Two-arg Expand[expr, patt] leaves parts free of patt unexpanded.
         * Partition the factors: pattern-free factors are carried as an
         * atomic coefficient (never distributed into, never expanded), and
         * only pattern-containing factors are expanded and multiplied out.
         * With patt == NULL every factor is active (expr_contains_patt returns
         * true), recovering plain single-arg Expand. */
        Expr** active = malloc(sizeof(Expr*) * count);
        Expr** freef  = malloc(sizeof(Expr*) * count);
        size_t na = 0, nf = 0;
        for (size_t i = 0; i < count; i++) {
            Expr* arg = e->data.function.args[i];
            if (patt && !expr_contains_patt(arg, patt)) {
                freef[nf++] = expr_copy(arg);
            } else {
                active[na++] = expr_expand_impl(arg, patt, overflow_mode);
            }
        }

        /* FLINT the polynomial sub-product only in the plain single-argument
         * case (patt == NULL), where the whole-expression size was already
         * bounded by the top-level fast path; the two-argument selective Expand
         * keeps the generic multiplier. */
        Expr* active_result = multiply_active(active, na,
                                              patt == NULL && flint_bridge_available());
        for (size_t i = 0; i < na; i++) expr_free(active[i]);
        free(active);

        if (nf == 0) {
            free(freef);
            return active_result;
        }

        /* Build the unexpanded free coefficient (a single Times of the
         * pattern-free factors). */
        Expr* free_coeff;
        if (nf == 1) {
            free_coeff = freef[0]; /* take ownership */
        } else {
            free_coeff = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), freef, nf));
        }
        free(freef);

        /* Distribute the free coefficient over each summand of the expanded
         * active part, keeping the coefficient itself unexpanded. */
        Expr* ret;
        if (active_result->type == EXPR_FUNCTION
            && active_result->data.function.head->type == EXPR_SYMBOL
            && active_result->data.function.head->data.symbol.name == SYM_Plus) {
            size_t k = active_result->data.function.arg_count;
            Expr** terms = malloc(sizeof(Expr*) * k);
            for (size_t i = 0; i < k; i++) {
                Expr* t_args[2] = { expr_copy(free_coeff),
                                    expr_copy(active_result->data.function.args[i]) };
                terms[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), t_args, 2));
            }
            ret = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms, k));
            free(terms);
        } else {
            Expr* t_args[2] = { expr_copy(free_coeff), expr_copy(active_result) };
            ret = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), t_args, 2));
        }
        expr_free(free_coeff);
        expr_free(active_result);
        return ret;
    }

    if (strcmp(head, "Power") == 0 && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp = e->data.function.args[1];
        if (exp->type == EXPR_INTEGER && exp->data.integer > 0) {
            int64_t n = exp->data.integer;
            Expr* exp_base = expr_expand_impl(base, patt, overflow_mode);
            /* Reached only when the top-level FLINT fast path did not handle
             * this power (FLINT absent, a non-polynomial base, or a
             * pattern-restricted Expand).  Gate on the estimated result size and
             * fall back to the generic distributor. */
            double pcap = overflow_mode ? EXPAND_OVERFLOW_CAP : EXPAND_FACTOR_CAP;
            SizeEst be = estimate_terms(exp_base, pcap);
            SizeEst pe = { 0.0, NULL, 0, 0 };
            double est = se_power(&be, n, &pe, pcap);
            se_free(&be); se_free(&pe);
            if (est <= pcap) {
                Expr* res = power_expand(exp_base, n);
                expr_free(exp_base);
                return res;
            }
            expr_free(exp_base);
            if (overflow_mode)
                return expr_new_function(expr_new_symbol(SYM_Overflow), NULL, 0);
            return expr_copy(e); /* internal: leave factored */
        }
        // For negative integer or non-integer power, we still don't go into subexpressions
        return expr_copy(e);
    }

    // Default: do not go into subexpressions
    return expr_copy(e);
}

/* Internal expander (public API).  Leaves an oversized power factored rather
 * than yielding Overflow[], so callers such as Together/Cancel/Apart and the
 * integrators always receive a usable expression. */
Expr* expr_expand_patt(Expr* e, Expr* patt) {
    return expr_expand_impl(e, patt, /*overflow_mode=*/false);
}

Expr* expr_expand(Expr* e) {
    return expr_expand_impl(e, NULL, /*overflow_mode=*/false);
}

/* User-facing Expand[]: never silently declines — an expansion too large to fit
 * in memory yields Overflow[]. */
Expr* builtin_expand(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 2) return NULL;
    Expr* patt = NULL;
    if (res->data.function.arg_count == 2) patt = res->data.function.args[1];
    Expr* ret = expr_expand_impl(res->data.function.args[0], patt, /*overflow_mode=*/true);
    return ret;
}

/* Returns true when e has the form Power[base, k] with k a negative integer. */
static bool is_negative_int_power(Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    if (e->data.function.head->data.symbol.name != SYM_Power) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* exp = e->data.function.args[1];
    return exp->type == EXPR_INTEGER && exp->data.integer < 0;
}

/* Threading head test: ExpandNumerator/ExpandDenominator descend into List,
 * Equal, Unequal, Less, LessEqual, Greater, GreaterEqual, Inequality, And,
 * Or, Not, and Plus. (Plus is handled because the operations apply
 * per-summand.) */
static bool is_thread_head(const char* head) {
    return strcmp(head, "List") == 0 || strcmp(head, "Equal") == 0 ||
           strcmp(head, "Unequal") == 0 || strcmp(head, "Less") == 0 ||
           strcmp(head, "LessEqual") == 0 || strcmp(head, "Greater") == 0 ||
           strcmp(head, "GreaterEqual") == 0 || strcmp(head, "Inequality") == 0 ||
           strcmp(head, "And") == 0 ||
           strcmp(head, "Or") == 0 || strcmp(head, "Not") == 0 ||
           strcmp(head, "Plus") == 0;
}

Expr* expr_expand_numerator(Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    const char* head = (e->data.function.head->type == EXPR_SYMBOL)
        ? e->data.function.head->data.symbol.name : "";

    if (is_thread_head(head)) {
        bool is_ineq = (strcmp(head, "Inequality") == 0);
        size_t n = e->data.function.arg_count;
        Expr** args = malloc(sizeof(Expr*) * (n > 0 ? n : 1));
        for (size_t i = 0; i < n; i++) {
            if (is_ineq && (i & 1u) == 1)
                args[i] = expr_copy(e->data.function.args[i]);
            else
                args[i] = expr_expand_numerator(e->data.function.args[i]);
        }
        Expr* ret = eval_and_free(expr_new_function(expr_copy(e->data.function.head), args, n));
        free(args);
        return ret;
    }

    if (strcmp(head, "Times") == 0) {
        size_t n = e->data.function.arg_count;
        Expr** num_args = malloc(sizeof(Expr*) * (n > 0 ? n : 1));
        Expr** den_args = malloc(sizeof(Expr*) * (n > 0 ? n : 1));
        size_t nc = 0, dc = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* arg = e->data.function.args[i];
            if (is_negative_int_power(arg)) {
                den_args[dc++] = expr_copy(arg);
            } else {
                num_args[nc++] = expr_copy(arg);
            }
        }

        Expr* num;
        if (nc == 0) {
            num = expr_new_integer(1);
        } else if (nc == 1) {
            num = num_args[0]; /* takes ownership */
        } else {
            num = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), num_args, nc));
        }
        free(num_args);

        Expr* expanded_num = expr_expand(num);
        expr_free(num);

        if (dc == 0) {
            free(den_args);
            return expanded_num;
        }

        Expr** result_args = malloc(sizeof(Expr*) * (dc + 1));
        result_args[0] = expanded_num;
        for (size_t i = 0; i < dc; i++) result_args[i + 1] = den_args[i];
        free(den_args);
        Expr* ret = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), result_args, dc + 1));
        free(result_args);
        return ret;
    }

    if (strcmp(head, "Power") == 0 && e->data.function.arg_count == 2) {
        Expr* exp = e->data.function.args[1];
        if (exp->type == EXPR_INTEGER && exp->data.integer < 0) {
            /* Pure denominator: ExpandNumerator leaves it unchanged. */
            return expr_copy(e);
        }
        /* Positive integer power (or symbolic): try to expand at the top level. */
        return expr_expand(e);
    }

    return expr_copy(e);
}

Expr* expr_expand_denominator(Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    const char* head = (e->data.function.head->type == EXPR_SYMBOL)
        ? e->data.function.head->data.symbol.name : "";

    if (is_thread_head(head)) {
        bool is_ineq = (strcmp(head, "Inequality") == 0);
        size_t n = e->data.function.arg_count;
        Expr** args = malloc(sizeof(Expr*) * (n > 0 ? n : 1));
        for (size_t i = 0; i < n; i++) {
            if (is_ineq && (i & 1u) == 1)
                args[i] = expr_copy(e->data.function.args[i]);
            else
                args[i] = expr_expand_denominator(e->data.function.args[i]);
        }
        Expr* ret = eval_and_free(expr_new_function(expr_copy(e->data.function.head), args, n));
        free(args);
        return ret;
    }

    if (strcmp(head, "Times") == 0) {
        size_t n = e->data.function.arg_count;
        Expr** num_args = malloc(sizeof(Expr*) * (n > 0 ? n : 1));
        Expr** den_pos = malloc(sizeof(Expr*) * (n > 0 ? n : 1));
        size_t nc = 0, dc = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* arg = e->data.function.args[i];
            if (is_negative_int_power(arg)) {
                Expr* base = arg->data.function.args[0];
                int64_t k = -arg->data.function.args[1]->data.integer; /* k > 0 */
                den_pos[dc++] = eval_and_free(expr_new_function(
                    expr_new_symbol(SYM_Power),
                    (Expr*[]){expr_copy(base), expr_new_integer(k)}, 2));
            } else {
                num_args[nc++] = expr_copy(arg);
            }
        }

        if (dc == 0) {
            for (size_t i = 0; i < nc; i++) expr_free(num_args[i]);
            free(num_args);
            free(den_pos);
            return expr_copy(e);
        }

        Expr* den_product;
        if (dc == 1) {
            den_product = den_pos[0];
        } else {
            den_product = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), den_pos, dc));
        }
        free(den_pos);

        Expr* expanded_den = expr_expand(den_product);
        expr_free(den_product);

        Expr* den_inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){expanded_den, expr_new_integer(-1)}, 2));

        Expr** result_args = malloc(sizeof(Expr*) * (nc + 1));
        for (size_t i = 0; i < nc; i++) result_args[i] = num_args[i];
        result_args[nc] = den_inv;
        free(num_args);
        Expr* ret = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), result_args, nc + 1));
        free(result_args);
        return ret;
    }

    if (strcmp(head, "Power") == 0 && e->data.function.arg_count == 2) {
        Expr* exp = e->data.function.args[1];
        if (exp->type == EXPR_INTEGER && exp->data.integer < 0) {
            int64_t k = -exp->data.integer;
            Expr* base = e->data.function.args[0];
            Expr* pos = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                (Expr*[]){expr_copy(base), expr_new_integer(k)}, 2));
            Expr* expanded = expr_expand(pos);
            expr_free(pos);
            Expr* ret = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                (Expr*[]){expanded, expr_new_integer(-1)}, 2));
            return ret;
        }
        return expr_copy(e);
    }

    return expr_copy(e);
}

Expr* builtin_expand_numerator(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    return expr_expand_numerator(res->data.function.args[0]);
}

Expr* builtin_expand_denominator(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    return expr_expand_denominator(res->data.function.args[0]);
}

void expand_init(void) {
    symtab_add_builtin("Expand", builtin_expand);
    symtab_get_def("Expand")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("ExpandNumerator", builtin_expand_numerator);
    symtab_get_def("ExpandNumerator")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("ExpandDenominator", builtin_expand_denominator);
    symtab_get_def("ExpandDenominator")->attributes |= ATTR_PROTECTED;
}

/* qafactor.c — Trager algebraic-factoring helpers (Phase G).
 *
 * Phase G3 lives here: norm of f ∈ Q(α)[x] via Resultant_y(P_α, g(x, y)).
 *
 * Strategy: lift QAExt and QAUPoly to picocas Expr*'s in two free
 * symbols (x_name, y_name), then call `internal_resultant` to drop into
 * the Sylvester-matrix routine in poly.c.  This keeps qafactor.c a
 * thin glue layer — all multivariate arithmetic is delegated. */

#include "qafactor.h"
#include "expr.h"
#include "expand.h"
#include "internal.h"
#include "eval.h"
#include "deriv.h"
#include "poly.h"
#include "arithmetic.h"
#include "sym_names.h"
#include <gmp.h>
#include <stdlib.h>
#include <string.h>

/* Forward decls for builtins we call via internal_call_impl. */
extern Expr* internal_call_impl(const char* name,
                                Expr* (*builtin_func)(Expr*),
                                Expr** args, size_t count);

/* ============================== mpq → Expr ============================== */

/* Convert an mpq_t to a picocas Expr.  Returns Integer when the
 * denominator is 1, otherwise Rational[num, den] (num, den as
 * normalised Integer / BigInt nodes).  Handles arbitrary precision. */
static Expr* mpq_to_expr(const mpq_t q) {
    mpz_t num, den;
    mpz_init(num);
    mpz_init(den);
    mpq_get_num(num, q);
    mpq_get_den(den, q);

    Expr* result;
    if (mpz_cmp_ui(den, 1) == 0) {
        result = expr_bigint_normalize(expr_new_bigint_from_mpz(num));
    } else {
        Expr* n = expr_bigint_normalize(expr_new_bigint_from_mpz(num));
        Expr* d = expr_bigint_normalize(expr_new_bigint_from_mpz(den));
        Expr* args[2] = { n, d };
        result = expr_new_function(expr_new_symbol("Rational"), args, 2);
    }

    mpz_clear(num);
    mpz_clear(den);
    return result;
}

/* ====================== Expr-tree construction helpers ================== */

/* Build var^k (k >= 0).  k == 0 → Integer 1, k == 1 → Symbol[var]. */
static Expr* var_power(const char* var, long k) {
    if (k == 0) return expr_new_integer(1);
    if (k == 1) return expr_new_symbol(var);
    Expr* args[2] = { expr_new_symbol(var), expr_new_integer((int64_t)k) };
    return expr_new_function(expr_new_symbol("Power"), args, 2);
}

/* Multiply `coef` (an Expr*, ownership taken) by var^k, with sensible
 * elision when coef is ±1 / k is 0. */
static Expr* term_coef_times_var_power(Expr* coef, const char* var, long k) {
    if (k == 0) return coef;

    Expr* vp = var_power(var, k);

    if (coef->type == EXPR_INTEGER && coef->data.integer == 1) {
        expr_free(coef);
        return vp;
    }
    if (coef->type == EXPR_INTEGER && coef->data.integer == -1) {
        expr_free(coef);
        Expr* args[2] = { expr_new_integer(-1), vp };
        return expr_new_function(expr_new_symbol("Times"), args, 2);
    }

    Expr* args[2] = { coef, vp };
    return expr_new_function(expr_new_symbol("Times"), args, 2);
}

/* Build a polynomial expression from an mpq_t-array of coefficients
 * coefs[0..n], where coefs[i] is the coefficient of var^i. */
static Expr* mpq_array_to_poly_expr(const mpq_t* coefs,
                                    size_t n_plus_1,
                                    const char* var) {
    /* First pass: count non-zero terms so we can size the Plus arg list. */
    size_t n_terms = 0;
    for (size_t i = 0; i < n_plus_1; i++) {
        if (mpq_sgn(coefs[i]) != 0) n_terms++;
    }
    if (n_terms == 0) return expr_new_integer(0);

    Expr** terms = (Expr**)malloc(sizeof(Expr*) * n_terms);
    size_t out = 0;
    for (size_t i = 0; i < n_plus_1; i++) {
        if (mpq_sgn(coefs[i]) == 0) continue;
        Expr* coef = mpq_to_expr(coefs[i]);
        terms[out++] = term_coef_times_var_power(coef, var, (long)i);
    }

    if (n_terms == 1) {
        Expr* r = terms[0];
        free(terms);
        return r;
    }

    Expr* result = expr_new_function(expr_new_symbol("Plus"), terms, n_terms);
    free(terms);
    return result;
}

/* ============================== Public API ============================== */

Expr* qaext_to_expr(const QAExt* ext, const char* y_name) {
    /* P_α(y) = sum_{i=0..deg} ext->coef[i] * y^i */
    return mpq_array_to_poly_expr(ext->coef, ext->deg + 1, y_name);
}

Expr* qaupoly_to_expr(const QAUPoly* f,
                      const char* x_name,
                      const char* y_name) {
    if (f->deg < 0) return expr_new_integer(0);

    /* Build f(x, α) as a sum of x-monomials whose coefficients are
     * polynomials-in-y derived from each QANum's α-coefficients. */
    size_t n_x_terms = 0;
    for (int i = 0; i <= f->deg; i++) {
        if (!qa_is_zero(f->c[i])) n_x_terms++;
    }
    if (n_x_terms == 0) return expr_new_integer(0);

    Expr** x_terms = (Expr**)malloc(sizeof(Expr*) * n_x_terms);
    size_t out = 0;
    for (int i = 0; i <= f->deg; i++) {
        if (qa_is_zero(f->c[i])) continue;
        /* Coefficient is a polynomial-in-y of α-degree < ext->deg. */
        Expr* y_poly = mpq_array_to_poly_expr(f->c[i]->coef,
                                              f->ext->deg,
                                              y_name);
        x_terms[out++] = term_coef_times_var_power(y_poly, x_name, (long)i);
    }

    Expr* g;
    if (n_x_terms == 1) {
        g = x_terms[0];
    } else {
        g = expr_new_function(expr_new_symbol("Plus"), x_terms, n_x_terms);
    }
    free(x_terms);
    return g;
}

Expr* qaupoly_norm(const QAUPoly* f,
                   const char* x_name,
                   const char* y_name) {
    if (f->deg < 0) return NULL;

    /* Build P_α(y) and g(x, y). */
    Expr* p_alpha = qaext_to_expr(f->ext, y_name);
    Expr* g       = qaupoly_to_expr(f, x_name, y_name);
    Expr* y_var   = expr_new_symbol(y_name);

    /* Resultant_y(P_α, g) — internal_resultant takes ownership of args
     * (they get attached to a temporary Resultant[...] node which is
     * freed inside internal_call_impl regardless of outcome). */
    Expr* args[3] = { p_alpha, g, y_var };
    Expr* r = internal_resultant(args, 3);
    if (!r) return NULL;

    /* The resultant routine already calls Expand internally on the
     * Sylvester determinant, but the fast paths (Times / Power
     * factorisations of P or Q) skip that — re-expand to guarantee a
     * canonical sum-of-monomials form. */
    Expr* expanded = expr_expand(r);
    expr_free(r);
    return expanded;
}

/* ============================== Phase G4 ============================== */

/* Lift a rational/integer Expr coefficient to mpq_t.  Returns false
 * for non-numeric or non-Q-rational inputs (e.g. symbolic).  Handles
 * Integer, BigInt, and Rational[n, d] with small-int n/d. */
static bool expr_coef_to_mpq(const Expr* e, mpq_t out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) {
        mpq_set_si(out, (long)e->data.integer, 1);
        return true;
    }
    if (e->type == EXPR_BIGINT) {
        mpq_set_z(out, e->data.bigint);
        return true;
    }
    int64_t n, d;
    if (is_rational((Expr*)e, &n, &d)) {
        mpq_set_si(out, (long)n, (unsigned long)d);
        mpq_canonicalize(out);
        return true;
    }
    return false;
}

/* Build a QAUPoly over Q(α) from a Q[x] Expr polynomial (no α
 * dependency).  The α-component of every coefficient is zero, so each
 * QANum sits in the Q-rational subspace of Q(α).  Returns NULL if any
 * coefficient is not Q-rational. */
static QAUPoly* expr_q_poly_to_qaupoly(Expr* poly,
                                       const char* x_name,
                                       const QAExt* ext) {
    Expr* x_var = expr_new_symbol(x_name);
    int deg = get_degree_poly(poly, x_var);
    if (deg < 0) {
        expr_free(x_var);
        return qaupoly_zero(ext);
    }
    QAUPoly* p = qaupoly_new(ext, deg + 1);
    bool ok = true;
    for (int i = 0; i <= deg; i++) {
        Expr* c = get_coeff(poly, x_var, i);
        mpq_t q;
        mpq_init(q);
        if (!expr_coef_to_mpq(c, q)) {
            ok = false;
            mpq_clear(q);
            expr_free(c);
            break;
        }
        QANum* qn = qa_from_mpq(ext, q);
        qaupoly_setcoef(p, i, qn);
        qa_free(qn);
        mpq_clear(q);
        expr_free(c);
    }
    expr_free(x_var);
    if (!ok) {
        qaupoly_free(p);
        return NULL;
    }
    qaupoly_normalize(p);
    return p;
}

/* True if R(x) ∈ Q[x] is squarefree (gcd(R, R') has degree 0). */
static bool expr_qx_is_squarefree(Expr* R, const char* x_name) {
    Expr* d_args[2] = { expr_copy(R), expr_new_symbol(x_name) };
    Expr* dR = internal_call_impl("D", builtin_d, d_args, 2);
    /* internal_call_impl returns the Expr unevaluated if builtin_d
     * yielded NULL; force a full evaluate so the derivative is in
     * canonical Plus-of-monomials form. */
    Expr* dR_eval = evaluate(dR);
    expr_free(dR);

    /* PolynomialGCD is variadic over polynomials (no variable arg);
     * passing a third symbol "x" silently coerces it into a
     * three-way gcd(R, R', x) which is always 1. */
    Expr* g_args[2] = { expr_copy(R), dR_eval };
    Expr* g = internal_polynomialgcd(g_args, 2);
    Expr* g_eval = evaluate(g);
    expr_free(g);

    Expr* x_check = expr_new_symbol(x_name);
    int deg = get_degree_poly(g_eval, x_check);
    expr_free(x_check);
    expr_free(g_eval);
    return (deg <= 0);
}

/* Build the QANum c = -s · α  in Q(α). */
static QANum* qa_neg_s_alpha(const QAExt* ext, int s) {
    QANum* alpha = qa_alpha(ext);
    QANum* c = qa_scale_si(alpha, -(long)s, 1);
    qa_free(alpha);
    return c;
}

void qa_sqfr_norm_result_free(QASqfrNormResult* r) {
    if (!r) return;
    if (r->g) qaupoly_free(r->g);
    if (r->R) expr_free(r->R);
    r->g = NULL;
    r->R = NULL;
    r->s = -1;
}

QASqfrNormResult qa_sqfr_norm(const QAUPoly* f,
                              const char* x_name,
                              const char* y_name) {
    QASqfrNormResult result = { -1, NULL, NULL };
    if (!f || f->deg < 0) return result;

    /* s = 0 fast path: try the un-shifted norm first. */
    QAUPoly* g = qaupoly_copy(f);
    Expr* R = qaupoly_norm(g, x_name, y_name);
    if (!R) {
        qaupoly_free(g);
        return result;
    }
    if (expr_qx_is_squarefree(R, x_name)) {
        result.s = 0;
        result.g = g;
        result.R = R;
        return result;
    }
    expr_free(R);
    qaupoly_free(g);

    /* Walk s = 1, 2, ...  At each step build g = f(x − sα) via
     * qaupoly_shift, recompute R, retest. */
    for (int s = 1; s <= QA_SQFR_NORM_MAX_TRIES; s++) {
        QANum* c = qa_neg_s_alpha(f->ext, s);
        QAUPoly* g_s = qaupoly_shift(f, c);
        qa_free(c);
        if (!g_s) continue;

        Expr* R_s = qaupoly_norm(g_s, x_name, y_name);
        if (!R_s) {
            qaupoly_free(g_s);
            continue;
        }
        if (expr_qx_is_squarefree(R_s, x_name)) {
            result.s = s;
            result.g = g_s;
            result.R = R_s;
            return result;
        }
        expr_free(R_s);
        qaupoly_free(g_s);
    }
    /* Never converged: result.s stays −1, callers treat as failure. */
    return result;
}

/* ====================== alg_factor: result walker ===================== */

/* Append to a growing (capacity-doubling) array of Q[x] factor Exprs.
 * Each entry is a borrowed reference — we expr_copy when extracting
 * from the Factor result so the caller (qa_alg_factor) owns them. */
typedef struct {
    Expr** v;
    size_t n, cap;
} ExprList;

static void exprlist_init(ExprList* L) { L->v = NULL; L->n = 0; L->cap = 0; }
static void exprlist_push(ExprList* L, Expr* e) {
    if (L->n + 1 > L->cap) {
        L->cap = L->cap ? L->cap * 2 : 8;
        L->v = (Expr**)realloc(L->v, sizeof(Expr*) * L->cap);
    }
    L->v[L->n++] = e;
}
static void exprlist_free_contents(ExprList* L) {
    for (size_t i = 0; i < L->n; i++) expr_free(L->v[i]);
    free(L->v);
    L->v = NULL; L->n = 0; L->cap = 0;
}

/* Walk a Factor[] result and collect its non-trivial polynomial
 * factors as freshly-copied Exprs.  Numeric scalars are skipped (they
 * vanish into the leading-coefficient unit).  Power[h, k] is treated
 * as h appearing k times — but for squarefree-norm input, k == 1
 * in practice. */
static void collect_factors(Expr* factored, const char* x_name, ExprList* out) {
    if (!factored) return;

    /* Times[a, b, c, ...]: recurse on each arg. */
    if (factored->type == EXPR_FUNCTION
        && factored->data.function.head
        && factored->data.function.head->type == EXPR_SYMBOL
        && factored->data.function.head->data.symbol == SYM_Times) {
        for (size_t i = 0; i < factored->data.function.arg_count; i++) {
            collect_factors(factored->data.function.args[i], x_name, out);
        }
        return;
    }

    /* Power[h, k]: emit h k times. */
    if (factored->type == EXPR_FUNCTION
        && factored->data.function.head
        && factored->data.function.head->type == EXPR_SYMBOL
        && factored->data.function.head->data.symbol == SYM_Power
        && factored->data.function.arg_count == 2) {
        Expr* base = factored->data.function.args[0];
        Expr* exp  = factored->data.function.args[1];
        if (exp->type == EXPR_INTEGER && exp->data.integer >= 1) {
            for (int64_t k = 0; k < exp->data.integer; k++) {
                collect_factors(base, x_name, out);
            }
            return;
        }
    }

    /* Otherwise: a bare polynomial or scalar.  Emit only if it has
     * positive degree in x (i.e. is a non-trivial polynomial). */
    Expr* x_var = expr_new_symbol(x_name);
    int deg = get_degree_poly(factored, x_var);
    expr_free(x_var);
    if (deg > 0) {
        exprlist_push(out, expr_copy(factored));
    }
}

QAUPoly** qa_alg_factor(const QAUPoly* f,
                        const char* x_name,
                        const char* y_name,
                        int* n_out) {
    if (n_out) *n_out = 0;
    if (!f || f->deg < 0) return NULL;

    /* Squarefree-norm: locate the smallest shift s for which the
     * norm of g(x) = f(x − sα) is squarefree over Q. */
    QASqfrNormResult sqn = qa_sqfr_norm(f, x_name, y_name);
    if (sqn.s < 0) return NULL;

    /* Factor R(x) over Q via picocas's existing Factor builtin
     * (single-arg variadic — Factor[poly], not Factor[poly, var]). */
    Expr* fac_args[1] = { expr_copy(sqn.R) };
    Expr* factored = internal_factor(fac_args, 1);
    /* Re-evaluate to catch any stray un-canonical Times shapes. */
    Expr* factored_eval = evaluate(factored);
    expr_free(factored);

    ExprList qfacs;
    exprlist_init(&qfacs);
    collect_factors(factored_eval, x_name, &qfacs);
    expr_free(factored_eval);

    /* Trager: if R has a single Q-irreducible factor, f is irreducible
     * over Q(α) — return [monic copy of f]. */
    if (qfacs.n <= 1) {
        exprlist_free_contents(&qfacs);
        QAUPoly** out = (QAUPoly**)malloc(sizeof(QAUPoly*));
        QAUPoly* monic = qaupoly_make_monic(f);
        if (!monic) {
            free(out);
            qa_sqfr_norm_result_free(&sqn);
            return NULL;
        }
        out[0] = monic;
        if (n_out) *n_out = 1;
        qa_sqfr_norm_result_free(&sqn);
        return out;
    }

    /* Trager loop: lift each Q-rational factor of R(x) to a Q(α)-factor
     * of g via gcd, divide it out of g, and undo the shift. */
    QAUPoly** out = (QAUPoly**)malloc(sizeof(QAUPoly*) * qfacs.n);
    int n_factors = 0;
    QAUPoly* g = qaupoly_copy(sqn.g);
    bool failed = false;

    for (size_t i = 0; i < qfacs.n; i++) {
        QAUPoly* h_q = expr_q_poly_to_qaupoly(qfacs.v[i], x_name, f->ext);
        if (!h_q) { failed = true; break; }

        /* gcd_{Q(α)[x]}( h(x), g(x, α) ) — the Trager lift. */
        QAUPoly* lift = qaupoly_gcd(h_q, g);
        qaupoly_free(h_q);
        if (!lift || lift->deg <= 0) {
            /* Defensive: in theory each Q-factor must lift to a non-
             * trivial Q(α)-factor of g.  If not, abort and report
             * failure — the caller can fall back to leaving f
             * unfactored. */
            if (lift) qaupoly_free(lift);
            failed = true;
            break;
        }

        /* Divide the lift out of g (exact division in Q(α)[x]). */
        QAUPoly *q_div, *r_div;
        if (!qaupoly_divrem(g, lift, &q_div, &r_div) || !qaupoly_is_zero(r_div)) {
            qaupoly_free(lift);
            if (q_div) qaupoly_free(q_div);
            if (r_div) qaupoly_free(r_div);
            failed = true;
            break;
        }
        qaupoly_free(r_div);
        qaupoly_free(g);
        g = q_div;

        /* Undo the linear shift: replace x → x + sα in lift. */
        if (sqn.s != 0) {
            QANum* alpha = qa_alpha(f->ext);
            QANum* pos_s_alpha = qa_scale_si(alpha, (long)sqn.s, 1);
            qa_free(alpha);
            QAUPoly* unshifted = qaupoly_shift(lift, pos_s_alpha);
            qa_free(pos_s_alpha);
            qaupoly_free(lift);
            lift = unshifted;
        }

        out[n_factors++] = lift;
    }

    qaupoly_free(g);
    exprlist_free_contents(&qfacs);
    qa_sqfr_norm_result_free(&sqn);

    if (failed) {
        for (int i = 0; i < n_factors; i++) qaupoly_free(out[i]);
        free(out);
        return NULL;
    }
    if (n_out) *n_out = n_factors;
    return out;
}

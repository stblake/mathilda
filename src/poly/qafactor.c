/* qafactor.c — Trager algebraic-factoring helpers (Phase G).
 *
 * Phase G3 lives here: norm of f ∈ Q(α)[x] via Resultant_y(P_α, g(x, y)).
 *
 * Strategy: lift QAExt and QAUPoly to Mathilda Expr*'s in two free
 * symbols (x_name, y_name), then call `internal_resultant` to drop into
 * the Sylvester-matrix routine in poly.c.  This keeps qafactor.c a
 * thin glue layer — all multivariate arithmetic is delegated. */

#include "qafactor.h"
#include "rootofunity.h"
#include "expr.h"
#include "expand.h"
#include "internal.h"
#include "eval.h"
#include "deriv.h"
#include "poly.h"
#include "flint_bridge.h"
#include "arithmetic.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "core.h"
#include <gmp.h>
#include <stdlib.h>
#include <string.h>

/* Forward decls for builtins we call via internal_call_impl. */
extern Expr* internal_call_impl(const char* name,
                                Expr* (*builtin_func)(Expr*),
                                Expr** args, size_t count);

/* Forward decl: nested-radical extension recogniser (Phase G8).
 * Defined at the bottom of this file; called from qa_resolve_extension
 * when the surface form is `Sqrt[base]` / `base^(1/n)` with non-integer
 * `base`.  Returns a fresh QAExt (caller owns) and writes
 * `*render_out` on success; NULL on unrecognised shapes. */
static QAExt* qa_resolve_nested_radical(const Expr* alpha_expr,
                                        Expr** render_out);

/* ============================== mpq → Expr ============================== */

/* Convert an mpq_t to a Mathilda Expr.  Returns Integer when the
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
        result = expr_new_function(expr_new_symbol(SYM_Rational), args, 2);
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
    return expr_new_function(expr_new_symbol(SYM_Power), args, 2);
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
        return expr_new_function(expr_new_symbol(SYM_Times), args, 2);
    }

    Expr* args[2] = { coef, vp };
    return expr_new_function(expr_new_symbol(SYM_Times), args, 2);
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

    Expr* result = expr_new_function(expr_new_symbol(SYM_Plus), terms, n_terms);
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
        g = expr_new_function(expr_new_symbol(SYM_Plus), x_terms, n_x_terms);
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

    /* Fast path: gcd(R, R') over Q via FLINT (fmpq_mpoly_gcd).  R is the
     * univariate rational norm polynomial, so this is squarely in FLINT's
     * wheelhouse and collapses the squarefree-shift search that otherwise
     * dominates Trager's runtime (the classical content/pseudo-remainder
     * PolynomialGCD below is ~88% of Factor[…, Extension -> …] on higher
     * norm degrees).  flint_multivariate_gcd returns NULL for anything it
     * cannot represent (numeric, non-polynomial) or when USE_FLINT is off,
     * in which case we fall through to the classical path unchanged.  Only
     * the degree of the gcd matters here, so the FLINT-monic normalisation
     * of its output is irrelevant. */
    Expr* fg = flint_multivariate_gcd(R, dR_eval);
    if (fg) {
        expr_free(dR_eval);
        Expr* x_fcheck = expr_new_symbol(x_name);
        int fdeg = get_degree_poly(fg, x_fcheck);
        expr_free(x_fcheck);
        expr_free(fg);
        return (fdeg <= 0);
    }

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
        && factored->data.function.head->data.symbol.name == SYM_Times) {
        for (size_t i = 0; i < factored->data.function.arg_count; i++) {
            collect_factors(factored->data.function.args[i], x_name, out);
        }
        return;
    }

    /* Power[h, k]: emit h k times. */
    if (factored->type == EXPR_FUNCTION
        && factored->data.function.head
        && factored->data.function.head->type == EXPR_SYMBOL
        && factored->data.function.head->data.symbol.name == SYM_Power
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

    /* Factor R(x) over Q via Mathilda's existing Factor builtin
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

/* ============================== Phase G5 ============================== */

/* Internal symbol used as the α-placeholder when round-tripping through
 * Mathilda's polynomial machinery (Coefficient, Expand).  Picked to be
 * unlikely to collide with user input; in the rare case the user
 * actually defines this symbol the worst outcome is a NULL return from
 * qa_factor_with_extension, which the caller treats as "leave
 * unfactored". */
#define QA_ALPHA_INTERNAL "$qa$alpha$"

/* Recursive structural substitution: replace every sub-expression
 * structurally equal to `target` with a fresh copy of `replacement`.
 * Returns a freshly-allocated tree; caller owns. */
static Expr* expr_subst(const Expr* e,
                        const Expr* target,
                        const Expr* replacement) {
    if (!e) return NULL;
    if (expr_eq((Expr*)e, (Expr*)target)) {
        return expr_copy((Expr*)replacement);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    Expr* new_head = expr_subst(e->data.function.head, target, replacement);
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = expr_subst(e->data.function.args[i], target, replacement);
    }
    Expr* result = expr_new_function(new_head, new_args, n);
    free(new_args);
    return result;
}

/* Recognise the imaginary unit, in either of its Mathilda forms:
 * the bare symbol I or the canonical Complex[0, 1]. */
static bool expr_is_imaginary_unit(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_I) return true;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Complex
        && e->data.function.arg_count == 2) {
        Expr* re = e->data.function.args[0];
        Expr* im = e->data.function.args[1];
        if (re->type == EXPR_INTEGER && re->data.integer == 0
            && im->type == EXPR_INTEGER && im->data.integer == 1) return true;
    }
    return false;
}

/* Recognise Sqrt[c] with integer c, returning c via *out_c.  Mathilda
 * canonicalises Sqrt[c] as Power[c, Rational[1, 2]] so we accept that
 * form too. */
static bool expr_is_sqrt_int(const Expr* e, long* out_c) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (!e->data.function.head
        || e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.arg_count < 1) return false;
    const char* head = e->data.function.head->data.symbol.name;
    if (head == SYM_Sqrt && e->data.function.arg_count == 1) {
        Expr* c = e->data.function.args[0];
        if (c->type != EXPR_INTEGER) return false;
        *out_c = (long)c->data.integer;
        return true;
    }
    if (head == SYM_Power && e->data.function.arg_count == 2) {
        Expr* c   = e->data.function.args[0];
        Expr* exp = e->data.function.args[1];
        int64_t p, q;
        if (c->type != EXPR_INTEGER) return false;
        if (!is_rational(exp, &p, &q) || p != 1 || q != 2) return false;
        *out_c = (long)c->data.integer;
        return true;
    }
    return false;
}

QAExt* qa_resolve_extension(const Expr* alpha_expr, Expr** render_out) {
    if (!alpha_expr || !render_out) return NULL;

    /* I  →  P_α(y) = y² + 1, render = I.  Mathilda auto-evaluates the
     * literal `I` into `Complex[0, 1]` and `Sqrt[-1]` likewise, so by
     * the time we see the option both forms have collapsed to the
     * canonical Complex[0, 1].  We accept the bare symbol I as a
     * defensive convenience. */
    if (expr_is_imaginary_unit(alpha_expr)) {
        QAExt* ext = qaext_new(2);
        qaext_set_coef_si(ext, 0,  1, 1);   /* +1 */
        qaext_set_coef_si(ext, 1,  0, 1);
        qaext_set_coef_si(ext, 2,  1, 1);   /* y² */
        *render_out = expr_new_symbol(SYM_I);
        return ext;
    }

    /* Root of unity (-1)^(p/q)  →  cyclotomic field Q(ζ_{2q}).
     *
     * Mathilda canonicalises every root of unity as Power[-1, Rational[p,q]]
     * (or Complex[0, ±1]).  The natural generator is ζ_{2q} = (-1)^(1/q),
     * whose minimal polynomial is the cyclotomic polynomial Φ_{2q}; the
     * user's α = (-1)^(p/q) = ζ_{2q}^p lives in Q(ζ_{2q}).  Representing it
     * in the fixed Q-basis modulo Φ_{2q} makes the identities ζ^{2q}=1 and
     * Σ ζ^k = 0 structural, so cyclotomic Together/Cancel/Simplify no longer
     * blow up the rational-over-Q path.  render = the natural generator
     * Power[-1, Rational[1, q]].  q == 1 (±1, rational) is not an extension. */
    {
        int64_t ru_p, ru_q;
        if (expr_is_root_of_unity_pow(alpha_expr, &ru_p, &ru_q) && ru_q >= 2) {
            QAExt* ext = qaext_cyclotomic(2ul * (unsigned long)ru_q);
            if (ext) {
                Expr* rat = expr_new_function(expr_new_symbol(SYM_Rational),
                    (Expr*[]){ expr_new_integer(1), expr_new_integer(ru_q) }, 2);
                *render_out = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_new_integer(-1), rat }, 2);
                return ext;
            }
        }
    }

    /* Sqrt[c]  →  P_α(y) = y² − c, render = Sqrt[c]. */
    {
        long c;
        if (expr_is_sqrt_int(alpha_expr, &c)) {
            QAExt* ext = qaext_sqrt_si(c);
            *render_out = expr_copy((Expr*)alpha_expr);
            return ext;
        }
    }

    /* Times[I, Sqrt[c]]  =  Sqrt[-c]  →  P_α(y) = y² + c.  Mathilda
     * auto-evaluates Sqrt[-c] for negative argument into I·Sqrt[c],
     * so `Extension -> Sqrt[-3]` arrives here as Times[Complex[0,1],
     * Sqrt[3]].  We accept either ordering of the two factors. */
    if (alpha_expr->type == EXPR_FUNCTION
        && alpha_expr->data.function.head
        && alpha_expr->data.function.head->type == EXPR_SYMBOL
        && alpha_expr->data.function.head->data.symbol.name == SYM_Times
        && alpha_expr->data.function.arg_count == 2) {
        Expr* a = alpha_expr->data.function.args[0];
        Expr* b = alpha_expr->data.function.args[1];
        long c = 0;
        Expr* rad = NULL;
        if (expr_is_imaginary_unit(a) && expr_is_sqrt_int(b, &c)) rad = b;
        else if (expr_is_imaginary_unit(b) && expr_is_sqrt_int(a, &c)) rad = a;
        if (rad && c > 0) {
            QAExt* ext = qaext_sqrt_si(-c);   /* P_α(y) = y² + c */
            *render_out = expr_copy((Expr*)alpha_expr);
            return ext;
        }
    }

    /* Power[c, p/q]  →  natural generator β = c^(1/q_red) with minimal
     * polynomial P_β(y) = y^q_red - c, where q_red = q / gcd(|p|, q).
     *
     * For integer c, |c| >= 2: we accept ANY integer p (not just p == 1)
     * by reducing p/q to lowest terms and using the natural generator
     * β = c^(1/q_red).  The user's α = c^(p/q) = β^p_red lives in Q(β).
     *
     * `*render_out` is set to the natural form Power[c, 1/q_red]
     * (preserving the user's `Sqrt[c]` shape when q_red == 2 and the
     * user typed Sqrt[]).  Substitution in qa_expr_to_qaupoly_with_alpha
     * matches by structural equality on this natural form, with a
     * preprocessing step (`expand_radicals_to_atomic_poly`) rewriting
     * any non-natural `Power[c, p_user/q_user]` in the input into a
     * polynomial in the alpha-symbol.
     *
     * Non-integer base falls through to G8 only when p == 1 — the G8
     * `z^n - base` minpoly construction is correct only for
     * α = base^(1/n).  General p with non-integer base is deferred. */
    if (alpha_expr->type == EXPR_FUNCTION
        && alpha_expr->data.function.head
        && alpha_expr->data.function.head->type == EXPR_SYMBOL
        && alpha_expr->data.function.head->data.symbol.name == SYM_Power
        && alpha_expr->data.function.arg_count == 2) {
        Expr* c   = alpha_expr->data.function.args[0];
        Expr* exp = alpha_expr->data.function.args[1];
        int64_t p, q;
        if (c->type == EXPR_INTEGER
            && is_rational(exp, &p, &q) && q != 1) {
            /* Reduce p/q to lowest terms and normalise sign. */
            if (q < 0) { p = -p; q = -q; }
            int64_t ap = p < 0 ? -p : p;
            int64_t g  = 1;
            {
                int64_t a = ap, b = q;
                while (b) { int64_t r = a % b; a = b; b = r; }
                g = a;
            }
            int64_t p_red = (g > 1) ? (p / g) : p;
            int64_t q_red = (g > 1) ? (q / g) : q;
            (void)p_red;  /* used by lift preprocessing, not here */

            int64_t cv = c->data.integer;
            if (q_red >= 2 && cv != 0 && cv != 1 && cv != -1) {
                QAExt* ext = qaext_root_si((long)cv, (unsigned)q_red);
                /* Byte-for-byte compatible with the pre-change branch
                 * when p == 1 && q_red == q: just copy alpha_expr. */
                if (p == 1 && q_red == q) {
                    *render_out = expr_copy((Expr*)alpha_expr);
                } else {
                    /* Synthesise the natural form Power[c, Rational[1, q_red]]. */
                    Expr* rat_args[2] = { expr_new_integer(1),
                                          expr_new_integer(q_red) };
                    Expr* exp_natural = expr_new_function(
                        expr_new_symbol(SYM_Rational), rat_args, 2);
                    Expr* pow_args[2] = { expr_new_integer(cv), exp_natural };
                    *render_out = expr_new_function(
                        expr_new_symbol(SYM_Power), pow_args, 2);
                }
                return ext;
            }
        }
        /* Phase G8: `base^(1/n)` with non-integer base — try to resolve
         * `base` as an element of a sub-tower built from the atomic
         * radicals it contains.  Restricted to p == 1: the `z^n - base`
         * minpoly construction in G8 is only correct for α = base^(1/n). */
        if (is_rational(exp, &p, &q) && p == 1 && q >= 2) {
            QAExt* ext = qa_resolve_nested_radical(alpha_expr, render_out);
            if (ext) return ext;
        }
        return NULL;
    }

    return NULL;
}

/* Build α^j ∈ Q(α) (j ≥ 0).  α^0 = 1, α^1 = α, then iterated multiply. */
static QANum* qa_alpha_power(const QAExt* ext, int j) {
    if (j == 0) return qa_one(ext);
    QANum* alpha = qa_alpha(ext);
    if (j == 1) return alpha;
    QANum* acc = qa_copy(alpha);
    for (int k = 1; k < j; k++) {
        QANum* tmp = qa_mul(acc, alpha);
        qa_free(acc);
        acc = tmp;
    }
    qa_free(alpha);
    return acc;
}

/* gcd64: positive int64 gcd via Euclid. */
static int64_t gcd64(int64_t a, int64_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int64_t r = a % b; a = b; b = r; }
    return a;
}

/* expand_radicals_to_atomic_poly: walk `poly` and rewrite every
 * `Power[c_base, p_user/q_user]` whose reduced denominator divides
 * `q_natural` into the polynomial-in-`alpha_sym_name` form
 * representing α^k, where α = c_base^(1/q_natural) and
 * k = p_red * (q_natural / q_red).  The polynomial form is reduced
 * modulo `α^q_natural - c_base` via `qa_alpha_power_signed`, producing
 * an Expr that contains `alpha_sym_name` only with integer exponents
 * (e.g. `Power[α, 2]`) — safe for the downstream coefficient-peeling
 * step in `qa_expr_to_qaupoly_with_alpha` which uses `Coefficient`.
 *
 * Raw structural rewrite: does NOT call `evaluate()` to avoid having
 * Mathilda's Times canonicaliser route the result back through the
 * radical-absorption logic that produces `Power[c, p/q]` forms in the
 * first place.
 *
 * `ext` must be a `qaext_root_si(c_base, q_natural)` extension.
 *
 * Caller owns the returned Expr. */
static Expr* expand_radicals_to_atomic_poly(const Expr* e,
                                            int64_t c_base,
                                            int64_t q_natural,
                                            const char* alpha_sym_name,
                                            const QAExt* ext) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    /* Match Power[c_base, p/q] with q != 1.  When the reduced
     * denominator divides q_natural, rewrite. */
    if (e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Power
        && e->data.function.arg_count == 2) {
        const Expr* base  = e->data.function.args[0];
        const Expr* exp_e = e->data.function.args[1];
        int64_t p, q;
        if (base->type == EXPR_INTEGER
            && base->data.integer == c_base
            && is_rational((Expr*)exp_e, &p, &q) && q != 1) {
            if (q < 0) { p = -p; q = -q; }
            int64_t ap = p < 0 ? -p : p;
            int64_t g  = gcd64(ap, q);
            int64_t p_red = g > 1 ? p / g : p;
            int64_t q_red = g > 1 ? q / g : q;

            if (q_red >= 2 && (q_natural % q_red) == 0) {
                int64_t scale = q_natural / q_red;
                /* Overflow guard: p_red * scale must fit in long. */
                if (scale != 0
                    && (p_red > INT64_MAX / scale || p_red < INT64_MIN / scale)) {
                    /* Skip rewrite; fall through to default recursion. */
                } else {
                    long k = (long)(p_red * scale);
                    QANum* qn = qa_alpha_power_signed(ext, k);
                    if (qn) {
                        Expr* out = mpq_array_to_poly_expr(qn->coef,
                                                           qn->ext->deg,
                                                           alpha_sym_name);
                        qa_free(qn);
                        return out;
                    }
                    /* qa_alpha_power_signed failed (inversion on
                     * reducible extension); fall through to default. */
                }
            }
        }
    }

    /* Also handle `Sqrt[c_base]` shape when q_natural is even.
     * Mathilda keeps Sqrt[c] as a special head (not Power[c, 1/2]) in
     * some surface forms. */
    if (e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Sqrt
        && e->data.function.arg_count == 1) {
        const Expr* base = e->data.function.args[0];
        if (base->type == EXPR_INTEGER && base->data.integer == c_base
            && (q_natural % 2) == 0) {
            long k = (long)(q_natural / 2);
            QANum* qn = qa_alpha_power_signed(ext, k);
            if (qn) {
                Expr* out = mpq_array_to_poly_expr(qn->coef, qn->ext->deg,
                                                   alpha_sym_name);
                qa_free(qn);
                return out;
            }
        }
    }

    /* Default: recurse into head + args, preserving structure. */
    size_t count = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (count ? count : 1));
    for (size_t i = 0; i < count; i++) {
        new_args[i] = expand_radicals_to_atomic_poly(
            e->data.function.args[i], c_base, q_natural, alpha_sym_name, ext);
    }
    Expr* new_head = expand_radicals_to_atomic_poly(
        e->data.function.head, c_base, q_natural, alpha_sym_name, ext);
    Expr* out = expr_new_function(new_head, new_args, count);
    free(new_args);
    return out;
}

/* expand_roots_of_unity_to_atomic_poly: the cyclotomic analogue of
 * expand_radicals_to_atomic_poly.  Walk `e` and rewrite every
 * root-of-unity atom (-1)^(p/q) (or Complex[0, ±1]) whose denominator q
 * divides `q_natural` into the polynomial-in-`alpha_sym_name` form
 * representing α^k, where α = ζ_{2·q_natural} = (-1)^(1/q_natural) and
 * k = p · (q_natural / q).  The power is reduced modulo the cyclotomic
 * minimal polynomial Φ_{2·q_natural} via `qa_alpha_power_signed`, so the
 * result contains `alpha_sym_name` only with integer exponents < φ(2q_nat)
 * — ready for the coefficient-peeling step in qa_expr_to_qaupoly_with_alpha.
 *
 * `ext` must be the matching `qaext_cyclotomic(2·q_natural)` extension.
 * Caller owns the returned Expr. */
static Expr* expand_roots_of_unity_to_atomic_poly(const Expr* e,
                                                  int64_t q_natural,
                                                  const char* alpha_sym_name,
                                                  const QAExt* ext) {
    if (!e) return NULL;

    int64_t p, q;
    if (expr_is_root_of_unity_pow(e, &p, &q) && q >= 1
        && q_natural > 0 && (q_natural % q) == 0) {
        int64_t scale = q_natural / q;
        if (scale != 0
            && !(p > INT64_MAX / scale || p < INT64_MIN / scale)) {
            long k = (long)(p * scale);
            QANum* qn = qa_alpha_power_signed(ext, k);
            if (qn) {
                Expr* out = mpq_array_to_poly_expr(qn->coef, qn->ext->deg,
                                                   alpha_sym_name);
                qa_free(qn);
                return out;
            }
        }
    }

    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    size_t count = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (count ? count : 1));
    for (size_t i = 0; i < count; i++) {
        new_args[i] = expand_roots_of_unity_to_atomic_poly(
            e->data.function.args[i], q_natural, alpha_sym_name, ext);
    }
    Expr* new_head = expand_roots_of_unity_to_atomic_poly(
        e->data.function.head, q_natural, alpha_sym_name, ext);
    Expr* out = expr_new_function(new_head, new_args, count);
    free(new_args);
    return out;
}

/* Lift an Expr polynomial in `var` (whose coefficients may be
 * polynomials in `alpha_render`) to a QAUPoly over `ext`.  Returns
 * NULL if any coefficient is not in Q(α) (e.g. contains a free symbol
 * other than `var` and `alpha_render`).
 *
 * Strategy: substitute alpha_render → an internal placeholder symbol,
 * then use Coefficient[..] to peel off (var^i)(α^j) terms.
 *
 * When `c_base > 0`, runs `expand_radicals_to_atomic_poly` first to
 * rewrite any `Power[c_base, p/q]` in the input (including non-natural
 * forms like `Power[2, -2/3]`) into a polynomial in `QA_ALPHA_INTERNAL`.
 * Pass `c_base = 0` to skip this preprocessing (used when the
 * extension isn't a `Power[c, 1/q]` form — e.g. towers, nested
 * radicals). */
static QAUPoly* qa_expr_to_qaupoly_with_alpha(const Expr* poly,
                                              const Expr* var,
                                              const Expr* alpha_render,
                                              const QAExt* ext) {
    Expr* alpha_sym = expr_new_symbol(QA_ALPHA_INTERNAL);

    /* Stage 0 — when the input contains `Power[c, p/q]` with `p != 1`
     * but `q_red` divides our `ext`'s natural denominator, rewrite to
     * an explicit polynomial in alpha_sym so the structural
     * `expr_subst` and `get_degree_poly` steps below see the radical.
     *
     * Detect (c_base, q_natural) from `ext` when `ext` has the
     * `y^q - c` shape: coef[0] = -c, coef[i] = 0 for 0 < i < deg,
     * coef[deg] = 1, with c an integer (mpq denominator == 1). */
    int64_t c_base_inferred = 0;
    int64_t q_natural_inferred = 0;
    if (ext && ext->deg >= 2) {
        bool is_yn_minus_c = true;
        /* coef[deg] must be 1. */
        if (mpz_cmp_ui(mpq_numref(ext->coef[ext->deg]), 1) != 0
            || mpz_cmp_ui(mpq_denref(ext->coef[ext->deg]), 1) != 0) {
            is_yn_minus_c = false;
        }
        for (size_t i = 1; is_yn_minus_c && i < ext->deg; i++) {
            if (mpq_sgn(ext->coef[i]) != 0) is_yn_minus_c = false;
        }
        if (is_yn_minus_c) {
            /* coef[0] should be -c with c a non-zero integer. */
            if (mpz_cmp_ui(mpq_denref(ext->coef[0]), 1) == 0
                && mpq_sgn(ext->coef[0]) != 0
                && mpz_fits_slong_p(mpq_numref(ext->coef[0]))) {
                long neg_c = mpz_get_si(mpq_numref(ext->coef[0]));
                c_base_inferred = -neg_c;
                q_natural_inferred = (int64_t)ext->deg;
            }
        }
    }

    /* Cyclotomic generator: `alpha_render` is the natural root of unity
     * (-1)^(1/q_nat).  Rewrite every root-of-unity atom in the input into a
     * polynomial in QA_ALPHA_INTERNAL modulo Φ_{2 q_nat} (the y^q−c shape
     * detection above does not match a cyclotomic minimal polynomial, so
     * c_base_inferred is 0 here). */
    int64_t ru_p, ru_q;
    bool is_cyclotomic = alpha_render
        && expr_is_root_of_unity_pow(alpha_render, &ru_p, &ru_q) && ru_q >= 2;

    Expr* poly_expanded;
    if (is_cyclotomic) {
        poly_expanded = expand_roots_of_unity_to_atomic_poly(
            poly, ru_q, QA_ALPHA_INTERNAL, ext);
    } else if (c_base_inferred != 0
               && c_base_inferred != 1 && c_base_inferred != -1) {
        poly_expanded = expand_radicals_to_atomic_poly(
            poly, c_base_inferred, q_natural_inferred, QA_ALPHA_INTERNAL, ext);
    } else {
        poly_expanded = expr_copy((Expr*)poly);
    }

    /* Stage 1 — substitute α's surface form with an opaque symbol.
     * After Stage 0, any natural-form `Power[c, 1/q_natural]` has
     * already been rewritten to alpha_sym, so this is a no-op for
     * radical generators.  It remains needed for non-radical
     * extensions (Sqrt[I] / nested radicals) where `c_base_inferred`
     * was 0. */
    Expr* poly_sub = alpha_render
        ? expr_subst(poly_expanded, alpha_render, alpha_sym)
        : expr_copy(poly_expanded);
    expr_free(poly_expanded);

    /* Stage 2 — Expand so that nested Power / Times terms collapse to
     * a canonical Plus-of-monomials form.  Without this,
     * get_degree_poly on (a*x + b)^2 returns 0 (it doesn't see x). */
    Expr* poly_x = expr_expand(poly_sub);
    expr_free(poly_sub);

    int xdeg = get_degree_poly(poly_x, (Expr*)var);
    if (xdeg < 0) {
        expr_free(poly_x);
        expr_free(alpha_sym);
        return qaupoly_zero(ext);
    }

    QAUPoly* p = qaupoly_new(ext, xdeg + 1);
    bool ok = true;

    for (int i = 0; i <= xdeg && ok; i++) {
        Expr* c_in_x = get_coeff(poly_x, (Expr*)var, i);
        Expr* c_expanded = expr_expand(c_in_x);
        expr_free(c_in_x);

        int adeg = get_degree_poly(c_expanded, alpha_sym);
        if (adeg < 0) {
            /* Zero coefficient — leave c[i] = 0 in QAUPoly. */
            expr_free(c_expanded);
            continue;
        }
        /* NOTE: `adeg` may be >= ext->deg.  The internal α-symbol always
         * denotes the extension generator α (Stages 0/1 only ever introduce
         * it as α), so any α-power is legitimately reducible modulo the
         * minimal polynomial.  For radical generators the evaluator usually
         * pre-reduces (Sqrt[2]^3 → 2 Sqrt[2]), but the opaque internal symbol
         * carries no such auto-rule — and, critically, cyclotomic generators
         * surface high powers after Stage-2 Expand of products/powers of the
         * rewritten root-of-unity atoms (e.g. ((-1)^(2/3)+...)^3 expands to
         * α^2, α^3, ...).  qa_alpha_power(ext, j) below reduces α^j mod P_α
         * for ANY j via field multiplication, so we fold every power down
         * rather than refusing.  A genuinely out-of-field atom is NOT an
         * α-power: it stays a Power[base, p/q] node, contributes nothing to
         * the α-degree, and is rejected by expr_coef_to_mpq() instead. */

        QANum* qn = qa_zero(ext);
        for (int j = 0; j <= adeg; j++) {
            Expr* aj = get_coeff(c_expanded, alpha_sym, j);
            mpq_t q;
            mpq_init(q);
            if (!expr_coef_to_mpq(aj, q)) {
                mpq_clear(q);
                expr_free(aj);
                qa_free(qn);
                qn = NULL;
                ok = false;
                break;
            }
            QANum* aj_pow  = qa_alpha_power(ext, j);
            QANum* term    = qa_scale_mpq(aj_pow, q);
            qa_free(aj_pow);
            QANum* sum     = qa_add(qn, term);
            qa_free(qn);
            qa_free(term);
            qn = sum;
            mpq_clear(q);
            expr_free(aj);
        }
        expr_free(c_expanded);

        if (ok) {
            qaupoly_setcoef(p, i, qn);
            qa_free(qn);
        }
    }
    expr_free(poly_x);
    expr_free(alpha_sym);

    if (!ok) {
        qaupoly_free(p);
        return NULL;
    }
    qaupoly_normalize(p);
    return p;
}

/* Public wrapper for the internal `_with_alpha` variant, exposed via
 * qafactor.h so that other modules (poly.c, rat.c, parfrac.c) can lift
 * polynomials into Q(α)[x] without duplicating the alpha-substitution +
 * coefficient-peeling logic. */
QAUPoly* qa_expr_to_qaupoly(const Expr* poly,
                            const Expr* var,
                            const Expr* alpha_render,
                            const QAExt* ext) {
    return qa_expr_to_qaupoly_with_alpha(poly, var, alpha_render, ext);
}

Expr* qaupoly_to_expr_alpha(const QAUPoly* f,
                            const char* x_name,
                            const Expr* alpha_render) {
    /* Round-trip through qaupoly_to_expr with an internal y-symbol,
     * then ReplaceAll y → alpha_render and evaluate.  Evaluation
     * collapses Sqrt[c]^k / c^(k/n) etc. into canonical form.
     *
     * Compound `alpha_render` (Phase G6 — γ = α_1 + s_2 α_2 + ... is a
     * Plus expression) needs an extra Expand pass: Mathilda's evaluator
     * does not auto-distribute (a + b)^k, so γ^k would otherwise leak
     * through unsimplified.  Expand turns it into a sum of products
     * of α_i powers, after which auto-evaluation collapses Sqrt[c]^2
     * → c and friends. */
    Expr* in_y = qaupoly_to_expr(f, x_name, QA_ALPHA_INTERNAL);
    Expr* y_sym = expr_new_symbol(QA_ALPHA_INTERNAL);
    Expr* in_alpha = expr_subst(in_y, y_sym, alpha_render);
    expr_free(in_y);
    expr_free(y_sym);
    Expr* expanded = expr_expand(in_alpha);
    expr_free(in_alpha);
    Expr* canon = evaluate(expanded);
    expr_free(expanded);
    return canon;
}

/* Trager core, factored out so both G5 (single α) and G6 (tower) can
 * share it.  Lifts `poly` → QAUPoly[x] over `ext`, runs the squarefree
 * pre-pass + qa_alg_factor + multiplicity trial-division, and renders
 * each factor back as a Mathilda Expr using `alpha_render_output` as the
 * surface form for the algebraic generator.
 *
 * `poly` may already contain `alpha_render_input` literally; that
 * occurrence is structurally substituted with the internal placeholder
 * symbol (QA_ALPHA_INTERNAL) before lifting.  Pass NULL to skip the
 * substitution (when `poly` already uses the internal symbol).
 *
 * `ext` and the two render Exprs are borrowed; the caller frees them. */
static Expr* qa_factor_inner(const Expr* poly,
                             const QAExt* ext,
                             const Expr* alpha_render_input,
                             const Expr* alpha_render_output,
                             const Expr* var) {
    if (!poly || !ext || !alpha_render_output || !var) return NULL;
    if (var->type != EXPR_SYMBOL) return NULL;

    const char* x_name = var->data.symbol.name;
    const char* y_name = QA_ALPHA_INTERNAL;

    /* 1) Lift input poly → QAUPoly. */
    QAUPoly* f = qa_expr_to_qaupoly_with_alpha(poly, var, alpha_render_input, ext);
    if (!f) return NULL;

    /* Degenerate cases: zero or constant poly. */
    if (f->deg <= 0) {
        Expr* result = qaupoly_to_expr_alpha(f, x_name, alpha_render_output);
        qaupoly_free(f);
        return result;
    }

    /* Trager's algorithm requires squarefree input.  Reduce by
     * dividing out gcd(f, f'). */
    QAUPoly* f_for_alg;
    {
        QAUPoly* fp = qaupoly_new(f->ext, f->deg > 0 ? f->deg : 1);
        for (int i = 1; i <= f->deg; i++) {
            QANum* scaled = qa_scale_si(f->c[i], (long)i, 1);
            qaupoly_setcoef(fp, i - 1, scaled);
            qa_free(scaled);
        }
        qaupoly_normalize(fp);
        QAUPoly* g = qaupoly_gcd(f, fp);
        qaupoly_free(fp);
        if (!g || g->deg <= 0) {
            if (g) qaupoly_free(g);
            f_for_alg = qaupoly_copy(f);
        } else {
            QAUPoly *q_div, *r_div;
            if (!qaupoly_divrem(f, g, &q_div, &r_div) || !qaupoly_is_zero(r_div)) {
                if (q_div) qaupoly_free(q_div);
                if (r_div) qaupoly_free(r_div);
                qaupoly_free(g);
                qaupoly_free(f);
                return NULL;
            }
            qaupoly_free(r_div);
            qaupoly_free(g);
            f_for_alg = q_div;
        }
    }

    /* 2) Factor the squarefree part. */
    int n_factors = 0;
    QAUPoly** factors = qa_alg_factor(f_for_alg, x_name, y_name, &n_factors);
    qaupoly_free(f_for_alg);

    if (!factors || n_factors <= 0) {
        if (factors) free(factors);
        qaupoly_free(f);
        return NULL;
    }

    /* 3) Multiplicity recovery: trial-divide `f` by each squarefree
     * factor as many times as it goes. */
    int* mult = (int*)calloc((size_t)n_factors, sizeof(int));
    QAUPoly* remaining = qaupoly_copy(f);
    for (int i = 0; i < n_factors; i++) {
        while (remaining->deg >= factors[i]->deg) {
            QAUPoly *q_div, *r_div;
            if (!qaupoly_divrem(remaining, factors[i], &q_div, &r_div)) {
                if (q_div) qaupoly_free(q_div);
                if (r_div) qaupoly_free(r_div);
                break;
            }
            if (!qaupoly_is_zero(r_div)) {
                qaupoly_free(q_div);
                qaupoly_free(r_div);
                break;
            }
            qaupoly_free(r_div);
            qaupoly_free(remaining);
            remaining = q_div;
            mult[i]++;
        }
    }
    QANum* leading = remaining->deg < 0 ? qa_one(ext) : qa_copy(remaining->c[0]);
    qaupoly_free(remaining);

    /* 4) Render. */
    size_t n_terms = 0;
    for (int i = 0; i < n_factors; i++) if (mult[i] > 0) n_terms++;
    bool unit_is_one = qa_is_one(leading);
    if (!unit_is_one) n_terms++;

    Expr** terms = (Expr**)malloc(sizeof(Expr*) * (n_terms > 0 ? n_terms : 1));
    size_t out = 0;

    if (!unit_is_one) {
        QAUPoly* unit_poly = qaupoly_from_qa(leading);
        terms[out++] = qaupoly_to_expr_alpha(unit_poly, x_name, alpha_render_output);
        qaupoly_free(unit_poly);
    }
    qa_free(leading);

    for (int i = 0; i < n_factors; i++) {
        if (mult[i] <= 0) continue;
        Expr* h = qaupoly_to_expr_alpha(factors[i], x_name, alpha_render_output);
        if (mult[i] == 1) {
            terms[out++] = h;
        } else {
            Expr* p_args[2] = { h, expr_new_integer((int64_t)mult[i]) };
            terms[out++] = expr_new_function(expr_new_symbol(SYM_Power),
                                             p_args, 2);
        }
    }

    Expr* result;
    if (n_terms == 0)      result = expr_new_integer(1);
    else if (n_terms == 1) result = terms[0];
    else                   result = expr_new_function(expr_new_symbol(SYM_Times),
                                                      terms, n_terms);
    free(terms);
    free(mult);

    Expr* canon = evaluate(result);
    expr_free(result);

    for (int i = 0; i < n_factors; i++) qaupoly_free(factors[i]);
    free(factors);
    qaupoly_free(f);
    return canon;
}

Expr* qa_factor_with_extension(const Expr* poly,
                               const Expr* alpha_expr,
                               const Expr* var) {
    if (!poly || !alpha_expr || !var) return NULL;

    Expr* alpha_render = NULL;
    QAExt* ext = qa_resolve_extension(alpha_expr, &alpha_render);
    if (!ext) return NULL;

    /* G5: input is in user surface form; let qa_factor_inner handle the
     * α_user → QA_ALPHA_INTERNAL substitution by passing the surface
     * render as `alpha_render_input`. */
    Expr* result = qa_factor_inner(poly, ext, alpha_render, alpha_render, var);

    qaext_free(ext);
    expr_free(alpha_render);
    return result;
}

/* ============================== Phase G6 ============================== */
/* Tower of extensions Q(α_1, ..., α_n) → primitive element γ.
 *
 * Iterative construction (Trager §3):
 *   - Resolve α_1 individually: ext_1 has min poly P_1(y), γ_1 = α_1.
 *   - For i = 2..n:
 *       (a) Resolve α_i individually: it has min poly Q_i(z) ∈ Q[z].
 *       (b) Find smallest s ∈ [0, QA_SQFR_NORM_MAX_TRIES] such that
 *
 *               R_i(w) = Res_z(Q_i(z), P_{i-1}(w − s·z))
 *
 *           is squarefree over Q.  R_i is the min poly of the new
 *           primitive element γ_i = γ_{i-1} + s·α_i.
 *       (c) Recover α_i ∈ Q(γ_i) via gcd_z(Q_i(z), P_{i-1}(γ_i − s·z))
 *           in Q(γ_i)[z].  By construction the gcd has degree 1
 *           (z − α_i), so α_i = −(constant term of monic gcd).
 *       (d) Re-embed each previous α_j ∈ Q(γ_{i-1}) into Q(γ_i): we
 *           have γ_{i-1} = γ_i − s·α_i ∈ Q(γ_i), so each α_j (a
 *           polynomial-in-γ_{i-1} with rational coefficients) is
 *           Horner-evaluated at the new base.
 *
 * Symbol layout: γ uses QA_ALPHA_INTERNAL (also reused as the algebraic
 * placeholder during input-poly lifting; safe because both are
 * construction-time only). The resultant elimination variable z uses
 * a separate internal symbol QA_Z_INTERNAL. */

#define QA_Z_INTERNAL "$qa$z$"

/* ===== Tower-construction helpers ===== */

/* Reify a QAExt's minimal polynomial as a Mathilda Expr in `y_name`. */
static Expr* qaext_min_poly_expr(const QAExt* ext, const char* y_name) {
    return mpq_array_to_poly_expr(ext->coef, ext->deg + 1, y_name);
}

/* Build (y_name − s·z_name) as an Expr. */
static Expr* expr_y_minus_s_z(const char* y_name, const char* z_name, int s) {
    if (s == 0) return expr_new_symbol(y_name);
    Expr* sz_args[2] = { expr_new_integer((int64_t)s),
                         expr_new_symbol(z_name) };
    Expr* sz = expr_new_function(expr_new_symbol(SYM_Times), sz_args, 2);
    Expr* neg_args[2] = { expr_new_integer(-1), sz };
    Expr* neg_sz = expr_new_function(expr_new_symbol(SYM_Times), neg_args, 2);
    Expr* sum_args[2] = { expr_new_symbol(y_name), neg_sz };
    return expr_new_function(expr_new_symbol(SYM_Plus), sum_args, 2);
}

/* Substitute y_name → (y_name − s·z_name) in P, then expand. */
static Expr* expr_p_at_y_minus_sz(const Expr* P, const char* y_name,
                                  const char* z_name, int s) {
    Expr* y_sym = expr_new_symbol(y_name);
    Expr* repl  = expr_y_minus_s_z(y_name, z_name, s);
    Expr* sub   = expr_subst(P, y_sym, repl);
    expr_free(y_sym);
    expr_free(repl);
    Expr* expanded = expr_expand(sub);
    expr_free(sub);
    return expanded;
}

/* Compute Res_z(Q, P) ∈ Q[w] via Mathilda's internal_resultant. */
static Expr* expr_resultant_z(const Expr* P, const Expr* Q,
                              const char* z_name) {
    Expr* args[3] = { expr_copy((Expr*)P),
                      expr_copy((Expr*)Q),
                      expr_new_symbol(z_name) };
    Expr* r = internal_resultant(args, 3);
    if (!r) return NULL;
    Expr* expanded = expr_expand(r);
    expr_free(r);
    return expanded;
}

/* Walk a Q[w] Expr and convert to a fresh monic-mpq array of length
 * deg+1 (caller mpq_clears each entry and frees the array).  Returns
 * the original (pre-monic) leading coefficient via *lc_out (caller
 * mpq_inits + mpq_clears) when non-NULL.  Returns deg, or -1 on
 * non-rational input or deg < 1. */
static int q_poly_to_mpq_array_monic(const Expr* P, const char* w_name,
                                     mpq_t** out) {
    Expr* w_var = expr_new_symbol(w_name);
    int deg = get_degree_poly((Expr*)P, w_var);
    if (deg < 1) {
        expr_free(w_var);
        return -1;
    }
    mpq_t* coef = (mpq_t*)malloc(sizeof(mpq_t) * (deg + 1));
    for (int i = 0; i <= deg; i++) mpq_init(coef[i]);

    bool ok = true;
    for (int i = 0; i <= deg; i++) {
        Expr* c = get_coeff((Expr*)P, w_var, i);
        if (!expr_coef_to_mpq(c, coef[i])) ok = false;
        expr_free(c);
        if (!ok) break;
    }
    expr_free(w_var);
    if (!ok) {
        for (int i = 0; i <= deg; i++) mpq_clear(coef[i]);
        free(coef);
        return -1;
    }
    if (mpq_sgn(coef[deg]) == 0) {
        for (int i = 0; i <= deg; i++) mpq_clear(coef[i]);
        free(coef);
        return -1;
    }
    if (mpq_cmp_ui(coef[deg], 1, 1) != 0) {
        mpq_t lc; mpq_init(lc); mpq_set(lc, coef[deg]);
        for (int i = 0; i <= deg; i++) mpq_div(coef[i], coef[i], lc);
        mpq_clear(lc);
    }
    *out = coef;
    return deg;
}

/* Convert a monic Q[w] Expr to a QAExt. */
static QAExt* qaext_from_q_expr(const Expr* P, const char* w_name) {
    mpq_t* coef = NULL;
    int deg = q_poly_to_mpq_array_monic(P, w_name, &coef);
    if (deg < 1) return NULL;
    QAExt* ext = qaext_new((size_t)deg);
    for (int i = 0; i <= deg; i++) {
        qaext_set_coef_mpq(ext, (size_t)i, coef[i]);
        mpq_clear(coef[i]);
    }
    free(coef);
    return ext;
}

/* Find smallest s such that R(w) = Res_z(Q(z), P(w − s·z)) is squarefree.
 * Returns the resulting Q[w] Expr (caller owns) and writes `*s_out`;
 * NULL on non-convergence within QA_SQFR_NORM_MAX_TRIES. */
static Expr* find_primitive_shift(const Expr* P_old_in_w,
                                  const Expr* Q_new_in_z,
                                  const char* w_name, const char* z_name,
                                  int* s_out) {
    for (int s = 0; s <= QA_SQFR_NORM_MAX_TRIES; s++) {
        Expr* P_shifted = expr_p_at_y_minus_sz(P_old_in_w, w_name, z_name, s);
        Expr* R = expr_resultant_z(P_shifted, Q_new_in_z, z_name);
        expr_free(P_shifted);
        if (!R) continue;
        if (expr_qx_is_squarefree(R, w_name)) {
            *s_out = s;
            return R;
        }
        expr_free(R);
    }
    return NULL;
}

/* Recover α_new ∈ Q(γ_new) via gcd_z(Q_new(z), P_old(γ_new − s·z))
 * in Q(γ_new)[z].  By Trager that gcd is z − α_new (linear, monic). */
static QANum* qa_recover_alpha_in_new(const QAExt* new_ext,
                                      const Expr* P_old_in_w,
                                      const Expr* Q_new_in_z,
                                      const char* w_name, const char* z_name,
                                      int s) {
    Expr* P_at_g_minus_sz = expr_p_at_y_minus_sz(P_old_in_w, w_name, z_name, s);

    Expr* z_var     = expr_new_symbol(z_name);
    Expr* gamma_sym = expr_new_symbol(w_name);

    QAUPoly* P_lifted = qa_expr_to_qaupoly_with_alpha(
        P_at_g_minus_sz, z_var, gamma_sym, new_ext);
    QAUPoly* Q_lifted = qa_expr_to_qaupoly_with_alpha(
        Q_new_in_z, z_var, gamma_sym, new_ext);

    expr_free(P_at_g_minus_sz);
    expr_free(z_var);
    expr_free(gamma_sym);

    if (!P_lifted || !Q_lifted) {
        if (P_lifted) qaupoly_free(P_lifted);
        if (Q_lifted) qaupoly_free(Q_lifted);
        return NULL;
    }
    QAUPoly* g = qaupoly_gcd(P_lifted, Q_lifted);
    qaupoly_free(P_lifted);
    qaupoly_free(Q_lifted);
    if (!g || g->deg != 1) {
        if (g) qaupoly_free(g);
        return NULL;
    }
    QAUPoly* monic = qaupoly_make_monic(g);
    qaupoly_free(g);
    if (!monic) return NULL;
    /* α_new = -(constant term) of (z − α_new). */
    QANum* alpha = qa_neg(monic->c[0]);
    qaupoly_free(monic);
    return alpha;
}

/* Re-embed an old QANum (polynomial-in-γ_old with rational coeffs)
 * into Q(γ_new) by Horner-evaluating at `gamma_old_in_new` (= γ_new
 * − s·α_new ∈ Q(γ_new)). */
static QANum* qa_reembed_old(const QANum* old, const QAExt* new_ext,
                             const QANum* gamma_old_in_new) {
    int top = -1;
    for (int i = (int)old->ext->deg - 1; i >= 0; i--) {
        if (mpq_sgn(old->coef[i]) != 0) { top = i; break; }
    }
    if (top < 0) return qa_zero(new_ext);

    QANum* acc = qa_from_mpq(new_ext, old->coef[top]);
    for (int i = top - 1; i >= 0; i--) {
        QANum* m    = qa_mul(acc, gamma_old_in_new);
        QANum* term = qa_from_mpq(new_ext, old->coef[i]);
        QANum* sum  = qa_add(m, term);
        qa_free(m);
        qa_free(term);
        qa_free(acc);
        acc = sum;
    }
    return acc;
}

/* Render a QANum as Σ_k coef[k] · γ^k where γ is `gamma_name`. */
static Expr* qanum_to_expr_in_gamma_sym(const QANum* qn, const char* gamma_name) {
    return mpq_array_to_poly_expr(qn->coef, qn->ext->deg, gamma_name);
}

/* Extend a tower in place by absorbing one more α.  Returns false on
 * any failure (caller is responsible for freeing the tower). */
static bool qa_tower_extend(QATower* t, const Expr* alpha_new) {
    Expr* alpha_render_new = NULL;
    QAExt* alpha_ext = qa_resolve_extension(alpha_new, &alpha_render_new);
    if (!alpha_ext) return false;

    const char* w_name = QA_ALPHA_INTERNAL;
    const char* z_name = QA_Z_INTERNAL;

    Expr* P_old_in_w = qaext_min_poly_expr(t->ext, w_name);
    Expr* Q_new_in_z = qaext_min_poly_expr(alpha_ext, z_name);

    int s = -1;
    Expr* R = find_primitive_shift(P_old_in_w, Q_new_in_z, w_name, z_name, &s);
    if (!R) {
        expr_free(P_old_in_w);
        expr_free(Q_new_in_z);
        qaext_free(alpha_ext);
        expr_free(alpha_render_new);
        return false;
    }

    QAExt* new_ext = qaext_from_q_expr(R, w_name);
    expr_free(R);
    if (!new_ext) {
        expr_free(P_old_in_w);
        expr_free(Q_new_in_z);
        qaext_free(alpha_ext);
        expr_free(alpha_render_new);
        return false;
    }

    QANum* alpha_new_in_new = qa_recover_alpha_in_new(
        new_ext, P_old_in_w, Q_new_in_z, w_name, z_name, s);
    expr_free(P_old_in_w);
    expr_free(Q_new_in_z);
    qaext_free(alpha_ext);
    if (!alpha_new_in_new) {
        qaext_free(new_ext);
        expr_free(alpha_render_new);
        return false;
    }

    /* γ_old ∈ Q(γ_new) = γ_new − s·α_new. */
    QANum* gamma_new_qa     = qa_alpha(new_ext);
    QANum* s_alpha_new      = qa_scale_si(alpha_new_in_new, (long)s, 1);
    QANum* gamma_old_in_new = qa_sub(gamma_new_qa, s_alpha_new);
    qa_free(gamma_new_qa);
    qa_free(s_alpha_new);

    QANum** new_alphas = (QANum**)malloc(sizeof(QANum*) * (t->n + 1));
    bool ok = true;
    for (int i = 0; i < t->n; i++) {
        QANum* re = qa_reembed_old(t->alphas[i], new_ext, gamma_old_in_new);
        if (!re) { ok = false; break; }
        new_alphas[i] = re;
    }
    qa_free(gamma_old_in_new);
    if (!ok) {
        for (int i = 0; i < t->n; i++) if (new_alphas[i]) qa_free(new_alphas[i]);
        free(new_alphas);
        qa_free(alpha_new_in_new);
        qaext_free(new_ext);
        expr_free(alpha_render_new);
        return false;
    }
    new_alphas[t->n] = alpha_new_in_new;

    Expr** new_renders = (Expr**)malloc(sizeof(Expr*) * (t->n + 1));
    for (int i = 0; i < t->n; i++) new_renders[i] = t->alpha_renders[i];
    new_renders[t->n] = alpha_render_new;

    Expr* new_gamma_render;
    if (s == 0) {
        /* γ_new = γ_old + 0 — render is unchanged. */
        new_gamma_render = t->gamma_render;
    } else {
        Expr* term;
        if (s == 1) {
            term = expr_copy(alpha_render_new);
        } else {
            Expr* args[2] = { expr_new_integer((int64_t)s),
                              expr_copy(alpha_render_new) };
            term = expr_new_function(expr_new_symbol(SYM_Times), args, 2);
        }
        Expr* args[2] = { t->gamma_render, term };
        new_gamma_render = expr_new_function(expr_new_symbol(SYM_Plus), args, 2);
        /* t->gamma_render is consumed into the new tree — null out to
         * avoid double-free below. */
    }

    /* Replace tower fields. Old t->alpha_renders[i] elements are now
     * owned by new_renders[i]; we just free the outer array.  Old
     * t->alphas[i] are replaced (we owned them), so qa_free each. */
    for (int i = 0; i < t->n; i++) qa_free(t->alphas[i]);
    free(t->alphas);
    free(t->alpha_renders);
    qaext_free(t->ext);

    t->ext = new_ext;
    t->n  += 1;
    t->alphas = new_alphas;
    t->alpha_renders = new_renders;
    t->gamma_render = new_gamma_render;
    return true;
}

void qa_tower_free(QATower* t) {
    if (!t) return;
    if (t->alphas) {
        for (int i = 0; i < t->n; i++) if (t->alphas[i]) qa_free(t->alphas[i]);
        free(t->alphas);
    }
    if (t->alpha_renders) {
        for (int i = 0; i < t->n; i++)
            if (t->alpha_renders[i]) expr_free(t->alpha_renders[i]);
        free(t->alpha_renders);
    }
    if (t->gamma_render) expr_free(t->gamma_render);
    if (t->ext) qaext_free(t->ext);
    free(t);
}

QATower* qa_resolve_extension_tower(Expr* const* alpha_exprs, int n) {
    if (!alpha_exprs || n < 1) return NULL;

    Expr* render_1 = NULL;
    QAExt* ext_1 = qa_resolve_extension(alpha_exprs[0], &render_1);
    if (!ext_1) return NULL;

    QATower* t = (QATower*)calloc(1, sizeof(QATower));
    if (!t) {
        qaext_free(ext_1);
        expr_free(render_1);
        return NULL;
    }
    t->ext = ext_1;
    t->n   = 1;
    t->alphas = (QANum**)malloc(sizeof(QANum*));
    t->alphas[0] = qa_alpha(ext_1);
    t->alpha_renders = (Expr**)malloc(sizeof(Expr*));
    t->alpha_renders[0] = render_1;
    t->gamma_render = expr_copy(render_1);

    for (int i = 1; i < n; i++) {
        if (!qa_tower_extend(t, alpha_exprs[i])) {
            qa_tower_free(t);
            return NULL;
        }
    }
    return t;
}

/* Forward decls: defined below (shared with the Cancel/Together tower path). */
static Expr* decompose_redundant_sqrts(const Expr* e, const QATower* t);
static Expr* tower_substitute_alphas(const Expr* e, const QATower* t);

Expr* qa_factor_with_extension_tower(const Expr* poly,
                                     Expr* const* alpha_exprs,
                                     int n_alphas,
                                     const Expr* var) {
    if (!poly || !alpha_exprs || n_alphas < 1 || !var) return NULL;
    if (var->type != EXPR_SYMBOL) return NULL;

    /* Single-element list reduces to G5. */
    if (n_alphas == 1) {
        return qa_factor_with_extension(poly, alpha_exprs[0], var);
    }

    QATower* t = qa_resolve_extension_tower(alpha_exprs, n_alphas);
    if (!t) return NULL;

    /* Decompose composite radicals (e.g. Sqrt[6] -> Times[Sqrt[2], Sqrt[3]])
     * into the tower's prime-Sqrt generators BEFORE substitution, exactly as
     * qa_cancel_with_tower does. Without this, a cross term like Sqrt[6] that
     * arises from expanding (x - Sqrt[2])(x - Sqrt[3]) is neither Sqrt[2] nor
     * Sqrt[3], so the per-generator substitution below misses it, leaving a
     * stray radical that qa_factor_inner cannot lift -> the factor silently
     * fails. */
    Expr* poly_decomposed = decompose_redundant_sqrts(poly, t);

    /* Substitute each tower generator α_i with its polynomial-in-γ form
     * (a polynomial in QA_ALPHA_INTERNAL). Reuse tower_substitute_alphas —
     * the same routine the Cancel/Together tower path uses — which matches
     * radicals by base value via expand_radicals_to_atomic_poly, so it
     * catches the Sqrt-headed factors produced by decompose_redundant_sqrts
     * (a bare structural expr_subst against the Power[c,1/2]-form renders
     * would miss them). After this pass `poly_internal`'s only extension
     * symbol is QA_ALPHA_INTERNAL, so qa_factor_inner lifts it directly. */
    Expr* poly_internal = tower_substitute_alphas(poly_decomposed, t);
    expr_free(poly_decomposed);

    Expr* result = qa_factor_inner(
        poly_internal, t->ext,
        /* alpha_render_input = */ NULL,    /* poly already uses internal sym */
        /* alpha_render_output = */ t->gamma_render,
        var);

    expr_free(poly_internal);
    qa_tower_free(t);
    return result;
}

/* Walk `e` and rewrite every `Sqrt[u]` / `Power[u, p/q]` whose base `u`
 * is a non-integer expression mentioning the integer base `c_base` by
 * canonicalising `u` via `Together[u, Extension -> Power[c_base, 1/q_natural]]`.
 * The result is then used as the radicand.
 *
 * This is the *input-side* counterpart to the Phase F canonicalise-post
 * step (which canonicalised the generator-set surface forms): it makes
 * mathematically-equal but structurally-distinct nested radicals in the
 * INPUT collapse to a single canonical form, so the downstream tower
 * substitution (`expr_subst(input, t->alpha_renders[i], ...)`) finds
 * all instances.
 *
 * Without this pass, `Together[Sqrt[A] - Sqrt[B], Extension -> Automatic]`
 * with A = B mathematically (but A != B structurally) misses one of the
 * Sqrts at substitution time. */
static Expr* canonicalise_nested_radicands(const Expr* e,
                                           int64_t c_base, int64_t q_natural) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    /* Detect Sqrt[u] / Power[u, p/q] with non-integer base `u`. */
    bool is_radical = false;
    const Expr* radicand = NULL;
    const Expr* exp_e = NULL;
    bool is_sqrt_head = false;
    if (e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Sqrt
        && e->data.function.arg_count == 1) {
        radicand = e->data.function.args[0];
        is_radical = true;
        is_sqrt_head = true;
    } else if (e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Power
        && e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp  = e->data.function.args[1];
        int64_t p, q;
        if (is_rational((Expr*)exp, &p, &q) && q != 1
            && base->type != EXPR_INTEGER) {
            radicand = base;
            exp_e = exp;
            is_radical = true;
        }
    }

    if (is_radical && radicand) {
        /* Canonicalise the radicand via Together-with-Extension. */
        Expr* base_e = expr_new_integer(c_base);
        Expr* rat_args[2] = { expr_new_integer(1),
                              expr_new_integer(q_natural) };
        Expr* exp_natural = expr_new_function(expr_new_symbol(SYM_Rational),
                                              rat_args, 2);
        Expr* pow_args[2] = { base_e, exp_natural };
        Expr* alpha = expr_new_function(expr_new_symbol(SYM_Power),
                                        pow_args, 2);
        Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
            (Expr*[]){expr_new_symbol(SYM_Extension), alpha}, 2);
        Expr* tg_call = expr_new_function(expr_new_symbol(SYM_Together),
            (Expr*[]){expr_copy((Expr*)radicand), rule}, 2);
        Expr* canon = evaluate(tg_call);
        expr_free(tg_call);

        if (canon && !expr_eq(canon, (Expr*)radicand)) {
            /* Rebuild the radical with the canonical radicand. */
            Expr* new_radical;
            if (is_sqrt_head) {
                new_radical = expr_new_function(expr_new_symbol(SYM_Sqrt),
                                                (Expr*[]){canon}, 1);
            } else {
                new_radical = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){canon, expr_copy((Expr*)exp_e)}, 2);
            }
            /* Evaluate so the new radical is canonical. */
            return eval_and_free(new_radical);
        }
        if (canon) expr_free(canon);
        /* Fall through to recursive walk if no change. */
    }

    /* Default: recurse into head + args, preserving structure. */
    size_t count = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (count ? count : 1));
    for (size_t i = 0; i < count; i++) {
        new_args[i] = canonicalise_nested_radicands(
            e->data.function.args[i], c_base, q_natural);
    }
    Expr* new_head = canonicalise_nested_radicands(
        e->data.function.head, c_base, q_natural);
    Expr* out = eval_and_free(expr_new_function(new_head, new_args, count));
    free(new_args);
    return out;
}

/* Layer-2 input-side decomposition helper. Walks `e` and rewrites every
 * Sqrt[c] / Power[c, 1/2] (c composite squarefree, all prime factors in
 * `prime_set`) into a raw `Times[Sqrt[p_1], ..., Sqrt[p_k]]` constructed
 * via expr_new_function — NOT through evaluate(), so Mathilda's Times
 * canonicaliser does not immediately re-combine the product back to
 * Sqrt[c]. After this pass, qa_cancel_with_tower's main substitution
 * loop replaces each Sqrt[p_i] with its γ-poly form, after which the
 * Times factors are γ-polynomials with no Sqrt nodes left to recombine
 * — and the downstream Together / lift sees a clean polynomial in γ.
 *
 * Conservative: returns expr_copy(e) when no rewrite fires, so callers
 * always get a freshly-allocated owned tree. */
static Expr* decompose_redundant_sqrts_walk(const Expr* e,
                                            const int64_t* prime_set,
                                            size_t n_primes) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    /* Detect Sqrt[c] or Power[c, 1/2] with c a composite squarefree
     * integer whose prime factors are all in prime_set. */
    if (e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        const Expr* base_e = NULL;
        bool is_sqrt = false;
        if (h == SYM_Sqrt && e->data.function.arg_count == 1) {
            base_e = e->data.function.args[0];
            is_sqrt = true;
        } else if (h == SYM_Power && e->data.function.arg_count == 2) {
            int64_t p, q;
            if (is_rational(e->data.function.args[1], &p, &q)
                && p == 1 && q == 2) {
                base_e = e->data.function.args[0];
                is_sqrt = true;
            }
        }
        if (is_sqrt && base_e && base_e->type == EXPR_INTEGER) {
            int64_t c = base_e->data.integer;
            if (c >= 2) {
                /* Trial-divide c, checking squarefree AND all primes
                 * are in prime_set. Accumulate factors. */
                int64_t factors[64];
                int n_factors = 0;
                int64_t r = c;
                bool ok = true;
                for (int64_t p = 2; p * p <= r && ok; p++) {
                    if (r % p == 0) {
                        int cnt = 0;
                        while (r % p == 0) { r /= p; cnt++; }
                        if (cnt != 1) { ok = false; break; }
                        bool found = false;
                        for (size_t j = 0; j < n_primes; j++) {
                            if (prime_set[j] == p) { found = true; break; }
                        }
                        if (!found) { ok = false; break; }
                        if (n_factors < 64) factors[n_factors++] = p;
                        else { ok = false; break; }
                    }
                }
                if (ok && r > 1) {
                    bool found = false;
                    for (size_t j = 0; j < n_primes; j++) {
                        if (prime_set[j] == r) { found = true; break; }
                    }
                    if (!found) ok = false;
                    else if (n_factors < 64) factors[n_factors++] = r;
                    else ok = false;
                }
                /* Only rewrite genuine composites (n_factors >= 2). */
                if (ok && n_factors >= 2) {
                    Expr** factor_exprs =
                        (Expr**)malloc(sizeof(Expr*) * n_factors);
                    for (int i = 0; i < n_factors; i++) {
                        Expr* b = expr_new_integer(factors[i]);
                        Expr* args_one[1] = {b};
                        factor_exprs[i] = expr_new_function(
                            expr_new_symbol(SYM_Sqrt), args_one, 1);
                    }
                    Expr* out = expr_new_function(
                        expr_new_symbol(SYM_Times), factor_exprs, n_factors);
                    free(factor_exprs);
                    return out;
                }
            }
        }
    }

    /* Generic recursion: rebuild with rewritten children. */
    Expr* new_head = decompose_redundant_sqrts_walk(
        e->data.function.head, prime_set, n_primes);
    size_t nc = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (nc > 0 ? nc : 1));
    for (size_t i = 0; i < nc; i++) {
        new_args[i] = decompose_redundant_sqrts_walk(
            e->data.function.args[i], prime_set, n_primes);
    }
    Expr* out = expr_new_function(new_head, new_args, nc);
    free(new_args);
    return out;
}

/* Public-to-module wrapper. Extracts the prime-Sqrt set from the tower
 * and runs the walker. Returns expr_copy(e) when the tower has < 2
 * prime-Sqrt generators (the precondition for any composite Sqrt to be
 * coalescable). */
static Expr* decompose_redundant_sqrts(const Expr* e, const QATower* t) {
    if (!t || t->n < 2) return expr_copy((Expr*)e);
    int64_t primes[QA_AUTODETECT_MAX_GENS];
    size_t n_primes = 0;
    for (int i = 0; i < t->n; i++) {
        const Expr* a = t->alpha_renders[i];
        if (!a || a->type != EXPR_FUNCTION) continue;
        if (!a->data.function.head
            || a->data.function.head->type != EXPR_SYMBOL) continue;
        const char* h = a->data.function.head->data.symbol.name;
        const Expr* base = NULL;
        bool is_sqrt = false;
        if (h == SYM_Sqrt && a->data.function.arg_count == 1) {
            base = a->data.function.args[0];
            is_sqrt = true;
        } else if (h == SYM_Power && a->data.function.arg_count == 2) {
            int64_t p, q;
            if (is_rational(a->data.function.args[1], &p, &q)
                && p == 1 && q == 2) {
                base = a->data.function.args[0];
                is_sqrt = true;
            }
        }
        if (!is_sqrt || !base || base->type != EXPR_INTEGER) continue;
        int64_t b = base->data.integer;
        if (b < 2) continue;
        /* Trial-divide b for primality. */
        bool is_prime = true;
        for (int64_t p = 2; p * p <= b; p++) {
            if (b % p == 0) { is_prime = false; break; }
        }
        if (is_prime && n_primes < QA_AUTODETECT_MAX_GENS) {
            primes[n_primes++] = b;
        }
    }
    if (n_primes < 2) return expr_copy((Expr*)e);
    return decompose_redundant_sqrts_walk(e, primes, n_primes);
}

/* Count Power[c, p/q] (q ≥ 2) occurrences in `e` — both atomic Sqrt[c]
 * and Power[c, p/q] forms. Used by the Layer-1 pre-substitution gate
 * in qa_cancel_with_tower to predict the post-sub leaf-count blow-up. */
static void count_radical_power_occurrences(const Expr* e, int* count) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_Sqrt && e->data.function.arg_count == 1) {
            (*count)++;
        } else if (h == SYM_Power && e->data.function.arg_count == 2) {
            int64_t p, q;
            if (is_rational(e->data.function.args[1], &p, &q) && q >= 2) {
                (*count)++;
            }
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        count_radical_power_occurrences(e->data.function.args[i], count);
    }
}

Expr* qa_cancel_with_tower(const Expr* arg, const QATower* t) {
    if (!arg || !t || t->n < 1) return NULL;

    /* 0a. Expand Power[Plus[..], integer] forms in the input so the
     * algebraic atoms (Sqrt[c], c^(1/q)) become free at the surface
     * level for the Step 1 substitution loop.  Without this an input
     * like `(Sqrt[2]+Sqrt[3])^2 - 5 - 2 Sqrt[6]` keeps the Power[Plus,2]
     * opaque all the way through Step 5 (PolynomialRemainder /
     * qa_expr_to_qaupoly cannot lift Power[Plus[γ_poly], 2] when the
     * Plus is not a polynomial in QA_ALPHA_INTERNAL after surface
     * substitution).  Expand fully distributes the binomial, after
     * which Step 0b (decompose_redundant_sqrts) sees the cross term
     * `2 Sqrt[6]` and rewrites it to `2 Sqrt[2] Sqrt[3]` for the
     * smaller tower.  Expand on rational expressions is monotone or a
     * no-op in leaf count for the inputs we hit here; the downstream
     * leaf-count gates still apply if a pathological Expand inflates
     * the surface form. */
    Expr* expanded_call = expr_new_function(
        expr_new_symbol(SYM_Expand),
        (Expr*[]){expr_copy((Expr*)arg)}, 1);
    Expr* arg_expanded = evaluate(expanded_call);
    expr_free(expanded_call);
    if (!arg_expanded) arg_expanded = expr_copy((Expr*)arg);

    /* 0b. Canonicalise nested radicands in the input against each
     * integer-base generator.  This makes mathematically-equal but
     * structurally-distinct `Sqrt[u]` / `Power[u, p/q]` sub-expressions
     * collapse to a single canonical form, so the substitution loop in
     * step 1 finds every instance.
     *
     * Layer-2 step: rewrite composite Sqrt[c] occurrences (e.g. Sqrt[15]
     * when the tower holds Sqrt[3] and Sqrt[5] as primes) into raw
     * Times[Sqrt[p_1], ..., Sqrt[p_k]] before the substitution loop.
     * This lets a smaller tower (degree 4 over Q(Sqrt[3], Sqrt[5]) vs
     * degree 8 over Q(Sqrt[3], Sqrt[5], Sqrt[15])) handle the input
     * correctly: the per-prime substitution replaces each Sqrt[p_i]
     * with its γ-poly, so the Times factors become γ-polys (no Sqrt
     * left for Mathilda's Times canonicaliser to re-combine into Sqrt[c]
     * when Together is evaluated in Step 4). */
    Expr* arg_internal = decompose_redundant_sqrts(arg_expanded, t);
    expr_free(arg_expanded);
    for (int i = 0; i < t->n; i++) {
        const Expr* a_r = t->alpha_renders[i];
        int64_t c_i = 0, q_i = 0;
        if (a_r && a_r->type == EXPR_FUNCTION
            && a_r->data.function.head
            && a_r->data.function.head->type == EXPR_SYMBOL
            && a_r->data.function.head->data.symbol.name == SYM_Power
            && a_r->data.function.arg_count == 2
            && a_r->data.function.args[0]->type == EXPR_INTEGER) {
            int64_t pe, qe;
            if (is_rational(a_r->data.function.args[1], &pe, &qe)
                && pe == 1 && qe >= 2) {
                c_i = a_r->data.function.args[0]->data.integer;
                q_i = qe;
            }
        }
        if (c_i == 0 && a_r && a_r->type == EXPR_FUNCTION
            && a_r->data.function.head
            && a_r->data.function.head->type == EXPR_SYMBOL
            && a_r->data.function.head->data.symbol.name == SYM_Sqrt
            && a_r->data.function.arg_count == 1
            && a_r->data.function.args[0]->type == EXPR_INTEGER) {
            c_i = a_r->data.function.args[0]->data.integer;
            q_i = 2;
        }
        if (c_i != 0 && c_i != 1 && c_i != -1 && q_i >= 2) {
            Expr* canonicalised = canonicalise_nested_radicands(
                arg_internal, c_i, q_i);
            if (canonicalised) {
                expr_free(arg_internal);
                arg_internal = canonicalised;
            }
        }
    }

    /* Layer-1 pre-substitution gate: the substitution loop below replaces
     * each `Power[c, p/q]` occurrence with a γ-polynomial of degree
     * `deg(γ) - 1`. With k radical occurrences and γ of degree d, the
     * blow-up is at minimum k·(d-1) extra leaves (each polynomial term
     * is `coef · γ^j`, ≥ 2 leaves).  A strict lower bound on lc_out is
     * `lc_in + k·(d-1)·2`.  When even that bound trips the post-sub
     * gate, the substitution + Together + qaupoly-lift + rejection
     * cycle is guaranteed wasted work (~30-100 ms for high-degree γ).
     * Bail before paying for it.
     *
     * Conservative direction: predictor never returns false-positive
     * (never bails when post-sub gate would have accepted).  Cases
     * where the post-sub form happens to be smaller than the lower
     * bound (rare cancellations during substitution) still get a
     * chance through the existing post-sub gate. */
    {
        size_t deg_gamma = (t->ext && t->ext->deg > 0) ? t->ext->deg : 1;
        if (deg_gamma >= 2) {
            int radical_occ = 0;
            count_radical_power_occurrences(arg_internal, &radical_occ);
            int64_t lc_in_pred = leaf_count_internal((Expr*)arg, true);
            int64_t lc_pred = lc_in_pred
                            + (int64_t)radical_occ * (int64_t)(deg_gamma - 1) * 2;
            if (lc_in_pred > 0
                && ((lc_pred > 100 && lc_pred > 2 * lc_in_pred)
                    || lc_pred > 200)) {
                expr_free(arg_internal);
                return NULL;
            }
        }
    }

    /* 1. Substitute each user-side α_i with its polynomial-in-γ form
     * (a polynomial in QA_ALPHA_INTERNAL).
     *
     * For each `Power[c, 1/q]`-shape α_i: first rewrite every
     * `Power[c, p/q]` form in the input to a polynomial in a per-α_i
     * temp symbol via `expand_radicals_to_atomic_poly`, then substitute
     * the temp symbol with the γ-poly form.  This catches the non-
     * natural `Power[c, p/q]` shapes (e.g. `Power[2, -2/3]`) that
     * Mathilda's Times canonicaliser produces, which a bare
     * `expr_subst(input, t->alpha_renders[i], ...)` would miss. */
    for (int i = 0; i < t->n; i++) {
        /* Parse alpha_renders[i] as `Power[c_i, 1/q_i]` or `Sqrt[c_i]`. */
        int64_t c_i = 0, q_i = 0;
        const Expr* a_r = t->alpha_renders[i];
        if (a_r && a_r->type == EXPR_FUNCTION
            && a_r->data.function.head
            && a_r->data.function.head->type == EXPR_SYMBOL
            && a_r->data.function.head->data.symbol.name == SYM_Power
            && a_r->data.function.arg_count == 2
            && a_r->data.function.args[0]->type == EXPR_INTEGER) {
            int64_t pe, qe;
            if (is_rational(a_r->data.function.args[1], &pe, &qe)
                && pe == 1 && qe >= 2) {
                c_i = a_r->data.function.args[0]->data.integer;
                q_i = qe;
            }
        }
        if (c_i == 0 && a_r && a_r->type == EXPR_FUNCTION
            && a_r->data.function.head
            && a_r->data.function.head->type == EXPR_SYMBOL
            && a_r->data.function.head->data.symbol.name == SYM_Sqrt
            && a_r->data.function.arg_count == 1
            && a_r->data.function.args[0]->type == EXPR_INTEGER) {
            c_i = a_r->data.function.args[0]->data.integer;
            q_i = 2;
        }

        Expr* alpha_in_gamma = qanum_to_expr_in_gamma_sym(
            t->alphas[i], QA_ALPHA_INTERNAL);

        if (c_i != 0 && c_i != 1 && c_i != -1 && q_i >= 2) {
            char tmp_name[64];
            snprintf(tmp_name, sizeof(tmp_name), "$qa$twrtmp$%d$", i);
            const char* tmp_iname = intern_symbol(tmp_name);
            QAExt* atom_ext = qaext_root_si((long)c_i, (unsigned)q_i);
            Expr* expanded = expand_radicals_to_atomic_poly(
                arg_internal, c_i, q_i, tmp_iname, atom_ext);
            qaext_free(atom_ext);

            Expr* tmp_sym_expr = expr_new_symbol(tmp_iname);
            Expr* substituted = expr_subst(expanded, tmp_sym_expr,
                                           alpha_in_gamma);
            expr_free(tmp_sym_expr);
            expr_free(expanded);
            expr_free(arg_internal);
            arg_internal = substituted;
        } else {
            /* Non-Power-shape α_i (e.g. Sqrt[non-integer], I): fall
             * back to the original structural substitution. */
            Expr* old = arg_internal;
            arg_internal = expr_subst(old, t->alpha_renders[i],
                                      alpha_in_gamma);
            expr_free(old);
        }
        expr_free(alpha_in_gamma);
    }

    /* Step 1 substitution can blow up the leaf count by O(deg(γ)) on
     * every Power[c, p/q] occurrence in the input.  The downstream
     * Together (no-extension) call drops into multivariate polynomial
     * GCD over Q[γ_internal, x, ...] which is intractable on the
     * resulting expressions.  Gate: if the post-substitution leaf
     * count exceeds an empirically-tuned threshold relative to the
     * input, bail out so the caller falls back to the non-tower path.
     * Closing this gap requires a multi-generator analogue of
     * together_recursive_ext that threads the tower's polynomial
     * relations through PolynomialLCM/PolynomialQuotient without
     * pre-expanding to γ-polynomial form; deferred. */
    {
        int64_t lc_in  = leaf_count_internal((Expr*)arg, true);
        int64_t lc_out = leaf_count_internal(arg_internal, true);
        if (lc_in > 0
            && ((lc_out > 100 && lc_out > 2 * lc_in)
                || lc_out > 200)) {
            expr_free(arg_internal);
            return NULL;
        }
    }

    /* 2. Combine into a single fraction via the standard (no-extension)
     * Together.  After step 1 the only "algebraic symbol" in arg_internal
     * is QA_ALPHA_INTERNAL — treated as a free polynomial symbol by the
     * standard path, which is exactly what we want before lifting. */
    Expr* together_call = expr_new_function(
        expr_new_symbol(SYM_Together),
        (Expr*[]){arg_internal}, 1);
    Expr* arg_combined = evaluate(together_call);
    expr_free(together_call);

    /* 3. Extract numerator and denominator. */
    Expr* num_call = expr_new_function(expr_new_symbol(SYM_Numerator),
        (Expr*[]){expr_copy(arg_combined)}, 1);
    Expr* num_expr = evaluate(num_call);
    expr_free(num_call);
    Expr* den_call = expr_new_function(expr_new_symbol(SYM_Denominator),
        (Expr*[]){expr_copy(arg_combined)}, 1);
    Expr* den_expr = evaluate(den_call);
    expr_free(den_call);
    expr_free(arg_combined);

    /* If denominator is 1, no cancellation needed; just rebuild the
     * numerator with γ → user-surface substitution.  qaupoly_to_expr_alpha
     * needs a poly path so we go through the standard lift just to
     * substitute QA_ALPHA_INTERNAL back to gamma_render. */

    /* 4. Pick the polynomial variable: the first free symbol that is
     * NOT QA_ALPHA_INTERNAL.  If none, this is a constant in Q(γ); we
     * still need to render back, but qaupoly arithmetic needs a
     * polynomial variable so bail in that case. */
    Expr* var = NULL;
    {
        size_t vc = 0, vcap = 8;
        Expr** vars = (Expr**)malloc(sizeof(Expr*) * vcap);
        collect_variables(num_expr, &vars, &vc, &vcap);
        collect_variables(den_expr, &vars, &vc, &vcap);
        for (size_t i = 0; i < vc; i++) {
            if (vars[i]->type == EXPR_SYMBOL && vars[i]->data.symbol.name
                && strcmp(vars[i]->data.symbol.name, QA_ALPHA_INTERNAL) == 0) continue;
            if (!var) var = expr_copy(vars[i]);
        }
        for (size_t i = 0; i < vc; i++) expr_free(vars[i]);
        free(vars);
    }

    if (!var) {
        /* Constant in Q(γ): reduce num and den modulo γ's minimal
         * polynomial (so both become polynomials in γ_internal of
         * degree < ext->deg), lift to QANums, divide via qa_div, then
         * render the resulting QANum's coefs in t->gamma_render.
         *
         * This is the proper Q(γ)-arithmetic path: it shrinks
         * `(Sqrt[2]+Sqrt[3])(Sqrt[2]-Sqrt[3])` to its canonical QANum
         * value `-1` instead of dropping the polynomial-in-γ form back
         * verbatim (which after substitution renders as the unwieldy
         * γ·(γ^3 - 10γ) surface form). */
        Expr* gamma_minpoly = qaext_to_expr(t->ext, QA_ALPHA_INTERNAL);
        Expr* gamma_sym_e   = expr_new_symbol(QA_ALPHA_INTERNAL);

        Expr* num_reduced_call = expr_new_function(
            expr_new_symbol(SYM_PolynomialRemainder),
            (Expr*[]){expr_copy(num_expr), expr_copy(gamma_minpoly),
                      expr_copy(gamma_sym_e)}, 3);
        Expr* num_reduced = evaluate(num_reduced_call);
        expr_free(num_reduced_call);

        Expr* den_reduced_call = expr_new_function(
            expr_new_symbol(SYM_PolynomialRemainder),
            (Expr*[]){expr_copy(den_expr), expr_copy(gamma_minpoly),
                      expr_copy(gamma_sym_e)}, 3);
        Expr* den_reduced = evaluate(den_reduced_call);
        expr_free(den_reduced_call);

        expr_free(gamma_minpoly);
        expr_free(num_expr); expr_free(den_expr);

        /* Lift to QAUPoly with a dummy polynomial variable.  Both
         * inputs are now polynomials in γ_internal of degree < ext->deg
         * with no other free variables, so the lift produces a degree-0
         * QAUPoly whose c[0] is the canonical QANum. */
        const char* dummy_name = "$qa$twrdummy$";
        Expr* dummy_var = expr_new_symbol(dummy_name);
        QAUPoly* num_qp = qa_expr_to_qaupoly(num_reduced, dummy_var,
                                             gamma_sym_e, t->ext);
        QAUPoly* den_qp = qa_expr_to_qaupoly(den_reduced, dummy_var,
                                             gamma_sym_e, t->ext);
        expr_free(dummy_var);
        expr_free(gamma_sym_e);
        expr_free(num_reduced); expr_free(den_reduced);

        Expr* candidate = NULL;
        /* Zero numerator: result is 0 directly. */
        if (num_qp && qaupoly_is_zero(num_qp) && den_qp
            && !qaupoly_is_zero(den_qp)) {
            candidate = expr_new_integer(0);
        }
        if (!candidate && num_qp && den_qp && !qaupoly_is_zero(den_qp)
            && num_qp->deg == 0 && den_qp->deg == 0) {
            QANum* result_qn = qa_div(num_qp->c[0], den_qp->c[0]);
            if (result_qn) {
                /* Render: polynomial in γ_internal via mpq_array_to_poly_expr,
                 * then substitute γ_internal → t->gamma_render. */
                Expr* gamma_poly = mpq_array_to_poly_expr(
                    result_qn->coef, t->ext->deg, QA_ALPHA_INTERNAL);
                qa_free(result_qn);

                Expr* gs = expr_new_symbol(QA_ALPHA_INTERNAL);
                candidate = expr_subst(gamma_poly, gs, t->gamma_render);
                expr_free(gs);
                expr_free(gamma_poly);
                /* Final canonicalisation: Expand then Together (no
                 * extension) to reduce γ-polynomial form to canonical
                 * linear-basis when γ is a sum of atomic radicals.
                 *
                 * Without Expand, `Together[1/(Sqrt[2]+Sqrt[3]) +
                 * 1/(Sqrt[2]-Sqrt[3]), Extension -> Automatic]` returns
                 * the unsimplified γ-polynomial form
                 * `-(Sqrt[2]+Sqrt[3])^3 + 9(Sqrt[2]+Sqrt[3])` instead
                 * of the canonical `-2 Sqrt[2]`.
                 *
                 * Without Together-no-extension, Expand of an
                 * arithmetically-zero γ-polynomial like
                 * `(Sqrt[2]+Sqrt[3])^2 - 5 - 2 Sqrt[6]` can produce a
                 * spurious nonzero linear combination of `Sqrt[6]` and
                 * `1/Sqrt[6]` terms (the qa_div lift doesn't fully kill
                 * the algebraic relation between independent autodetect
                 * generators).  Running `Together` over `Q` after
                 * Expand collapses these to 0 via Sqrt-base polynomial
                 * GCD.  The downstream leaf-count gate keeps the
                 * candidate only if it beats the input. */
                Expr* expand_call = expr_new_function(
                    expr_new_symbol(SYM_Expand),
                    (Expr*[]){candidate}, 1);
                Expr* expanded = evaluate(expand_call);
                expr_free(expand_call);

                Expr* tog_call = expr_new_function(
                    expr_new_symbol(SYM_Together),
                    (Expr*[]){expanded}, 1);
                candidate = evaluate(tog_call);
                expr_free(tog_call);
            }
        }
        if (num_qp) qaupoly_free(num_qp);
        if (den_qp) qaupoly_free(den_qp);

        if (!candidate) return NULL;

        int64_t in_size  = leaf_count_internal((Expr*)arg, true);
        int64_t out_size = leaf_count_internal(candidate, true);
        if (out_size > in_size) {
            expr_free(candidate);
            return NULL;
        }
        return candidate;
    }

    /* 5. Lift to QAUPoly over Q(γ).  alpha_render = NULL signals that
     * the polynomial already uses QA_ALPHA_INTERNAL as the α symbol. */
    QAUPoly* num = qa_expr_to_qaupoly(num_expr, var, NULL, t->ext);
    QAUPoly* den = qa_expr_to_qaupoly(den_expr, var, NULL, t->ext);
    expr_free(num_expr); expr_free(den_expr);

    if (!num || !den || qaupoly_is_zero(den)) {
        if (num) qaupoly_free(num);
        if (den) qaupoly_free(den);
        expr_free(var);
        return NULL;
    }

    /* 6. g = gcd(num, den); num /= g, den /= g. */
    QAUPoly* g = qaupoly_gcd(num, den);
    if (!g || qaupoly_is_zero(g)) {
        qaupoly_free(num); qaupoly_free(den);
        if (g) qaupoly_free(g);
        expr_free(var);
        return NULL;
    }
    QAUPoly *q_num = NULL, *r_num = NULL;
    QAUPoly *q_den = NULL, *r_den = NULL;
    bool ok = qaupoly_divrem(num, g, &q_num, &r_num)
           && qaupoly_divrem(den, g, &q_den, &r_den);
    qaupoly_free(num); qaupoly_free(den); qaupoly_free(g);

    Expr* result = NULL;
    if (ok && q_num && q_den) {
        Expr* num_out = qaupoly_to_expr_alpha(q_num, var->data.symbol.name,
                                              t->gamma_render);
        Expr* den_out = qaupoly_to_expr_alpha(q_den, var->data.symbol.name,
                                              t->gamma_render);
        Expr* den_inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){den_out, expr_new_integer(-1)}, 2));
        result = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){num_out, den_inv}, 2));
    }

    if (q_num) qaupoly_free(q_num);
    if (r_num) qaupoly_free(r_num);
    if (q_den) qaupoly_free(q_den);
    if (r_den) qaupoly_free(r_den);
    expr_free(var);

    /* Complexity gate.  When the compositum's primitive element γ is a
     * sum of generators (e.g. γ = 2^(1/3) + Sqrt[α/6]), the Q(γ)
     * representation of the cancelled fraction expresses every term as
     * a polynomial in γ of high degree — typically much larger than
     * the user's input.  Return NULL to let the caller fall back to the
     * no-tower path when the tower-cancelled form is bigger than the
     * input by LeafCount.  Uses heads=true so Power/Times nodes count
     * (the canonical Mathilda complexity proxy used by simp_search). */
    if (result) {
        int64_t in_size  = leaf_count_internal((Expr*)arg,    true);
        int64_t out_size = leaf_count_internal((Expr*)result, true);
        if (out_size > in_size) {
            expr_free(result);
            return NULL;
        }
    }
    return result;
}

/* Phase C helper: substitute every tower generator α_i in `e` with its
 * polynomial-in-γ form (a polynomial in QA_ALPHA_INTERNAL).  Mirrors
 * the substitution loop in `qa_cancel_with_tower` Step 1, including the
 * `expand_radicals_to_atomic_poly` preprocessing for integer-base
 * generators so non-natural `Power[c, -2/3]`-shape radicals are caught
 * alongside the canonical `Power[c, 1/q]` form.
 *
 * Returns a fresh Expr; caller owns. */
static Expr* tower_substitute_alphas(const Expr* e, const QATower* t) {
    if (!e || !t) return NULL;
    Expr* arg_internal = expr_copy((Expr*)e);
    for (int i = 0; i < t->n; i++) {
        int64_t c_i = 0, q_i = 0;
        const Expr* a_r = t->alpha_renders[i];
        if (a_r && a_r->type == EXPR_FUNCTION
            && a_r->data.function.head
            && a_r->data.function.head->type == EXPR_SYMBOL
            && a_r->data.function.head->data.symbol.name == SYM_Power
            && a_r->data.function.arg_count == 2
            && a_r->data.function.args[0]->type == EXPR_INTEGER) {
            int64_t pe, qe;
            if (is_rational(a_r->data.function.args[1], &pe, &qe)
                && pe == 1 && qe >= 2) {
                c_i = a_r->data.function.args[0]->data.integer;
                q_i = qe;
            }
        }
        if (c_i == 0 && a_r && a_r->type == EXPR_FUNCTION
            && a_r->data.function.head
            && a_r->data.function.head->type == EXPR_SYMBOL
            && a_r->data.function.head->data.symbol.name == SYM_Sqrt
            && a_r->data.function.arg_count == 1
            && a_r->data.function.args[0]->type == EXPR_INTEGER) {
            c_i = a_r->data.function.args[0]->data.integer;
            q_i = 2;
        }

        Expr* alpha_in_gamma = qanum_to_expr_in_gamma_sym(
            t->alphas[i], QA_ALPHA_INTERNAL);

        if (c_i != 0 && c_i != 1 && c_i != -1 && q_i >= 2) {
            char tmp_name[64];
            snprintf(tmp_name, sizeof(tmp_name), "$qa$twrtmp$%d$", i);
            const char* tmp_iname = intern_symbol(tmp_name);
            QAExt* atom_ext = qaext_root_si((long)c_i, (unsigned)q_i);
            Expr* expanded = expand_radicals_to_atomic_poly(
                arg_internal, c_i, q_i, tmp_iname, atom_ext);
            qaext_free(atom_ext);

            Expr* tmp_sym_expr = expr_new_symbol(tmp_iname);
            Expr* substituted = expr_subst(expanded, tmp_sym_expr,
                                           alpha_in_gamma);
            expr_free(tmp_sym_expr);
            expr_free(expanded);
            expr_free(arg_internal);
            arg_internal = substituted;
        } else {
            /* Non-Power-shape α_i (Sqrt[non-integer], I, nested radical):
             * fall back to the bare structural substitution. */
            Expr* old = arg_internal;
            arg_internal = expr_subst(old, t->alpha_renders[i],
                                      alpha_in_gamma);
            expr_free(old);
        }
        expr_free(alpha_in_gamma);
    }
    return arg_internal;
}

/* Common preamble for the tower-based GCD/LCM helpers below: substitute
 * each α_i in every input, collect the free polynomial variable, and
 * lift the substituted inputs to QAUPoly[x] over the compositum.
 *
 * On success populates *poly_var_out (borrowed pointer into vars_out[]),
 * vars_out[] / *vc_out (owned by caller), and ps_out[] (owned by caller).
 * On failure frees everything it allocated and returns false.  The
 * subst_inputs array itself is always returned to the caller — caller
 * must free its contents.  Designed for the common scaffolding both
 * GCD and LCM need before they diverge on the per-poly arithmetic. */
static bool tower_lift_inputs(Expr* const* argv, size_t argc,
                              const QATower* t,
                              Expr*** subst_out, size_t* subst_n_out,
                              Expr*** vars_out, size_t* vc_out,
                              const Expr** poly_var_out,
                              QAUPoly*** ps_out) {
    *subst_out = NULL; *subst_n_out = 0;
    *vars_out = NULL;  *vc_out = 0;
    *poly_var_out = NULL;
    *ps_out = NULL;

    /* 1. Substitute each α_i in each input. */
    Expr** subst = (Expr**)malloc(sizeof(Expr*) * argc);
    if (!subst) return false;
    for (size_t i = 0; i < argc; i++) subst[i] = NULL;
    for (size_t i = 0; i < argc; i++) {
        subst[i] = tower_substitute_alphas(argv[i], t);
        if (!subst[i]) {
            for (size_t j = 0; j < argc; j++) if (subst[j]) expr_free(subst[j]);
            free(subst);
            return false;
        }
    }

    /* 2. Collect free vars; require exactly one besides QA_ALPHA_INTERNAL. */
    Expr* alpha_internal = expr_new_symbol(QA_ALPHA_INTERNAL);
    size_t vc = 0, vcap = 8;
    Expr** vars = (Expr**)malloc(sizeof(Expr*) * vcap);
    for (size_t i = 0; i < argc; i++) {
        collect_variables(subst[i], &vars, &vc, &vcap);
    }
    Expr* poly_var = NULL;
    size_t live = 0;
    for (size_t i = 0; i < vc; i++) {
        if (expr_eq(vars[i], alpha_internal)) continue;
        poly_var = vars[i];
        live++;
    }
    expr_free(alpha_internal);
    if (live != 1 || !poly_var || poly_var->type != EXPR_SYMBOL) {
        for (size_t i = 0; i < vc; i++) expr_free(vars[i]);
        free(vars);
        for (size_t i = 0; i < argc; i++) expr_free(subst[i]);
        free(subst);
        return false;
    }

    /* 3. Lift each substituted input to QAUPoly[x] over the compositum. */
    QAUPoly** ps = (QAUPoly**)malloc(sizeof(QAUPoly*) * argc);
    for (size_t i = 0; i < argc; i++) ps[i] = NULL;
    for (size_t i = 0; i < argc; i++) {
        ps[i] = qa_expr_to_qaupoly(subst[i], poly_var, NULL, t->ext);
        if (!ps[i]) {
            for (size_t j = 0; j < argc; j++) if (ps[j]) qaupoly_free(ps[j]);
            free(ps);
            for (size_t j = 0; j < vc; j++) expr_free(vars[j]);
            free(vars);
            for (size_t j = 0; j < argc; j++) expr_free(subst[j]);
            free(subst);
            return false;
        }
    }

    *subst_out = subst;     *subst_n_out = argc;
    *vars_out = vars;       *vc_out = vc;
    *poly_var_out = poly_var;
    *ps_out = ps;
    return true;
}

static void tower_lift_cleanup(Expr** subst, size_t subst_n,
                               Expr** vars, size_t vc,
                               QAUPoly** ps, size_t ps_n) {
    if (ps) {
        for (size_t i = 0; i < ps_n; i++) if (ps[i]) qaupoly_free(ps[i]);
        free(ps);
    }
    if (vars) {
        for (size_t i = 0; i < vc; i++) expr_free(vars[i]);
        free(vars);
    }
    if (subst) {
        for (size_t i = 0; i < subst_n; i++) if (subst[i]) expr_free(subst[i]);
        free(subst);
    }
}

Expr* qa_polynomialgcd_with_tower(Expr* const* argv, size_t argc,
                                  const QATower* t) {
    if (argc < 1 || !t || t->n < 1) return NULL;
    /* Single-generator case: callers should route through
     * polynomialgcd_with_extension; signal NULL so the wrapper falls
     * through to the single-α path. */
    if (t->n == 1) return NULL;

    Expr** subst = NULL; size_t subst_n = 0;
    Expr** vars  = NULL; size_t vc = 0;
    const Expr* poly_var = NULL;
    QAUPoly** ps = NULL;
    if (!tower_lift_inputs(argv, argc, t,
                           &subst, &subst_n, &vars, &vc, &poly_var, &ps)) {
        return NULL;
    }

    Expr* result = NULL;
    QAUPoly* g = qaupoly_copy(ps[0]);
    for (size_t i = 1; i < argc && g; i++) {
        QAUPoly* next = qaupoly_gcd(g, ps[i]);
        qaupoly_free(g);
        g = next;
    }
    if (g && !qaupoly_is_zero(g)) {
        result = qaupoly_to_expr_alpha(g, poly_var->data.symbol.name,
                                       t->gamma_render);
    }
    if (g) qaupoly_free(g);

    tower_lift_cleanup(subst, subst_n, vars, vc, ps, argc);
    /* Strip any uniform rational content factor left by the single-variable
     * Q(γ) lift so the GCD is monic in all polynomial variables. */
    if (result) result = qa_make_poly_numerically_monic(result, false);
    return result;
}

Expr* qa_polynomiallcm_with_tower(Expr* const* argv, size_t argc,
                                  const QATower* t) {
    if (argc < 1 || !t || t->n < 1) return NULL;
    if (t->n == 1) return NULL;

    Expr** subst = NULL; size_t subst_n = 0;
    Expr** vars  = NULL; size_t vc = 0;
    const Expr* poly_var = NULL;
    QAUPoly** ps = NULL;
    if (!tower_lift_inputs(argv, argc, t,
                           &subst, &subst_n, &vars, &vc, &poly_var, &ps)) {
        return NULL;
    }

    Expr* result = NULL;
    /* lcm(a, b) = a * b / gcd(a, b); fold left-to-right. */
    QAUPoly* L = qaupoly_copy(ps[0]);
    for (size_t i = 1; i < argc && L; i++) {
        QAUPoly* gp = qaupoly_gcd(L, ps[i]);
        QAUPoly* prod = qaupoly_mul(L, ps[i]);
        QAUPoly *quot = NULL, *rem = NULL;
        if (gp && qaupoly_divrem(prod, gp, &quot, &rem)) {
            qaupoly_free(L);
            L = quot;
            qaupoly_free(rem);
        } else {
            if (quot) qaupoly_free(quot);
            if (rem)  qaupoly_free(rem);
            qaupoly_free(L);
            L = NULL;
        }
        if (gp) qaupoly_free(gp);
        qaupoly_free(prod);
    }
    if (L && !qaupoly_is_zero(L)) {
        QAUPoly* monic = qaupoly_make_monic(L);
        if (monic) {
            result = qaupoly_to_expr_alpha(monic, poly_var->data.symbol.name,
                                           t->gamma_render);
            qaupoly_free(monic);
        }
    }
    if (L) qaupoly_free(L);

    tower_lift_cleanup(subst, subst_n, vars, vc, ps, argc);
    /* qaupoly_make_monic normalises the leading coefficient over Q(γ) only in
     * the single lifted variable; a second polynomial variable folded into the
     * coefficient field (or a γ-denominator) can still leave a uniform rational
     * content factor.  Strip it so the LCM is monic in all polynomial
     * variables. */
    if (result) result = qa_make_poly_numerically_monic(result, false);
    return result;
}

/* Numeric prefactor of a polynomial term: an OWNED copy of its leading
 * Integer/Rational/Real/Bigint coefficient, or Integer 1 when the term carries
 * no numeric factor (a bare monomial, or a monomial with an algebraic
 * coefficient such as Sqrt[2] x). */
static Expr* qa_term_numeric_coeff(const Expr* t) {
    if (!t) return expr_new_integer(1);
    if (t->type == EXPR_INTEGER || t->type == EXPR_BIGINT
        || t->type == EXPR_REAL || is_rational(t, NULL, NULL))
        return expr_copy((Expr*)t);
    if (t->type == EXPR_FUNCTION && t->data.function.head
        && t->data.function.head->type == EXPR_SYMBOL
        && t->data.function.head->data.symbol.name == SYM_Times
        && t->data.function.arg_count >= 1) {
        const Expr* first = t->data.function.args[0];
        if (first->type == EXPR_INTEGER || first->type == EXPR_BIGINT
            || first->type == EXPR_REAL || is_rational(first, NULL, NULL))
            return expr_copy((Expr*)first);
    }
    return expr_new_integer(1);
}

/* Render a multivariate polynomial `p` monic in its polynomial variables by
 * dividing out the (rational/real) leading coefficient of its lex-leading
 * monomial.  `p` is consumed; returns a new Expr (which may be structurally
 * `p` unchanged when no numeric normalisation applies).
 *
 * Why this exists: an extension GCD/LCM assembled over a primitive-element
 * tower can carry a spurious rational content factor.  The primitive element
 * γ has γ-polynomial images of the generators with rational denominators
 * (e.g. Sqrt[2] = (γ^3 − 9γ)/2 for γ = Sqrt[2] + Sqrt[3]), and the intermediate
 * Q[γ, x, …] GCD/LCM content normalisation leaks that denominator into the
 * result as a uniform scalar.  Wolfram normalises PolynomialGCD / PolynomialLCM
 * so the leading polynomial term is monic; this restores that convention and
 * matches the sibling passing tests (`x + y`, `2^(1/3) + x`, …).
 *
 * Only a *numeric* (Integer/Rational/Real) leading coefficient is divided out.
 * When the leading monomial's coefficient is algebraic (e.g. Sqrt[2] + Sqrt[3])
 * the result is a legitimate associate over Q(γ) that must be preserved — see
 * the documented Phase-D unit-factor boundary tests.  When `only_noninteger`
 * is true, a pure *integer* leading coefficient is also left intact (the
 * integer-domain primitive-part convention, e.g. PolynomialGCD[2x−3, 4x²−9] =
 * 2x−3); only genuine fractions / inexact leads are cleared.  When false, any
 * numeric lead ≠ ±1 is divided out (the tower path, whose inputs are monic in
 * the polynomial variables so the result should be monic too). */
Expr* qa_make_poly_numerically_monic(Expr* p, bool only_noninteger) {
    if (!p) return NULL;
    /* Fully expand + evaluate so the polynomial is a flat, canonically ordered
     * sum of monomials (tower renderers hand back un-evaluated Power/Times
     * trees; Coefficient/Exponent are unreliable across algebraic-number
     * coefficients, so we read the leading monomial off the canonical order
     * directly instead). */
    p = eval_and_free(expr_new_function(expr_new_symbol(SYM_Expand),
        (Expr*[]){ p }, 1));
    if (!p) return NULL;

    /* Mathilda's canonical Plus ordering places the leading (highest-degree)
     * monomial last; a non-Plus is its own leading term. */
    const Expr* lead = p;
    if (p->type == EXPR_FUNCTION && p->data.function.head
        && p->data.function.head->type == EXPR_SYMBOL
        && p->data.function.head->data.symbol.name == SYM_Plus
        && p->data.function.arg_count >= 1) {
        lead = p->data.function.args[p->data.function.arg_count - 1];
    }
    Expr* lc = qa_term_numeric_coeff(lead);

    /* Decide whether `lc` is a numeric factor worth dividing out. */
    int divide = 0;
    int64_t rp, rq;
    if (lc->type == EXPR_INTEGER) {
        divide = (!only_noninteger) && lc->data.integer != 0
                 && lc->data.integer != 1 && lc->data.integer != -1;
    } else if (lc->type == EXPR_BIGINT) {
        divide = !only_noninteger;
    } else if (lc->type == EXPR_REAL) {
        divide = (lc->data.real != 0.0 && lc->data.real != 1.0
                  && lc->data.real != -1.0);
    } else if (is_rational(lc, &rp, &rq)) {
        divide = 1; /* genuine fraction (denominator != 1) */
    }
    if (!divide) { expr_free(lc); return p; }

    Expr* inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ lc, expr_new_integer(-1) }, 2));
    Expr* scaled = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ p, inv }, 2));
    Expr* expanded = eval_and_free(expr_new_function(expr_new_symbol(SYM_Expand),
        (Expr*[]){ scaled }, 1));
    return expanded;
}

/* Phase D core: multivariate tower-based GCD/LCM via substitute-back.
 *
 * Substitutes each α_i in every input with its polynomial-in-γ form
 * (using QA_ALPHA_INTERNAL as the γ symbol), then calls the
 * no-extension multivariate `internal_polynomialgcd` or
 * `internal_polynomiallcm` on the substituted forms, treating
 * QA_ALPHA_INTERNAL as just another polynomial variable.  Finally
 * substitutes QA_ALPHA_INTERNAL → t->gamma_render and runs
 * Expand + evaluate to canonicalise.
 *
 * Correctness caveat: this computes the GCD over `Q[γ, x, y, ...]`
 * (γ treated as polynomial variable) rather than the canonical GCD
 * over `Q(γ)[x, y, ...]` (γ algebraic, modulo its minimal polynomial).
 * The two coincide exactly when the Q[γ,...]-GCD has γ-degree
 * strictly less than `deg(γ_min)` AND the inputs are γ-primitive.
 * When they differ:
 *   - GCD returned is a Q(γ)-divisor of both inputs, but possibly not
 *     the maximal one (no "extra" common factor unlocked by γ_min).
 *   - LCM is correspondingly a Q(γ)-multiple of both inputs, but
 *     possibly larger than the canonical Q(γ)-LCM.
 * Both results remain mathematically equivalent; downstream Cancel /
 * Together passes can reduce further when needed.  For the practical
 * inputs Phase D targets (inputs linear or low-degree in γ — the
 * Cardano sub-problem shape `PolynomialLCM[18, R]` with R linear in
 * the radicals), the two GCD definitions agree. */
static Expr* tower_multivar_combine(Expr* const* argv, size_t argc,
                                    const QATower* t,
                                    const char* op_head) {
    if (!argv || argc < 1 || !t || t->n < 1 || !op_head) return NULL;

    /* 1. Substitute each α_i → γ-poly in each input. */
    Expr** subst = (Expr**)malloc(sizeof(Expr*) * argc);
    if (!subst) return NULL;
    for (size_t i = 0; i < argc; i++) subst[i] = NULL;
    for (size_t i = 0; i < argc; i++) {
        subst[i] = tower_substitute_alphas(argv[i], t);
        if (!subst[i]) {
            for (size_t j = 0; j < argc; j++) if (subst[j]) expr_free(subst[j]);
            free(subst);
            return NULL;
        }
    }

    /* 2. Build the call: head[subst[0], subst[1], ...] and evaluate.
     * The substituted forms own their Exprs; the new function takes
     * ownership of those copies. */
    Expr** call_args = (Expr**)malloc(sizeof(Expr*) * argc);
    for (size_t i = 0; i < argc; i++) call_args[i] = expr_copy(subst[i]);
    Expr* call = expr_new_function(expr_new_symbol(op_head), call_args, argc);
    free(call_args);

    Expr* in_gamma = evaluate(call);
    expr_free(call);

    for (size_t i = 0; i < argc; i++) expr_free(subst[i]);
    free(subst);

    if (!in_gamma) return NULL;

    /* 3. Guard: if evaluate returned the original head unchanged (i.e.
     * no-ext path also failed), don't render — fall back to caller. */
    if (in_gamma->type == EXPR_FUNCTION
        && in_gamma->data.function.head
        && in_gamma->data.function.head->type == EXPR_SYMBOL
        && in_gamma->data.function.head->data.symbol.name
        && strcmp(in_gamma->data.function.head->data.symbol.name, op_head) == 0) {
        expr_free(in_gamma);
        return NULL;
    }

    /* 4. Substitute γ_internal → γ_render, Expand to distribute
     * Power[Plus[...], n] terms (γ_render is typically a Plus for
     * multi-α towers), and evaluate to canonicalise (Sqrt[c]^2 → c,
     * etc.). */
    Expr* gamma_sym = expr_new_symbol(QA_ALPHA_INTERNAL);
    Expr* in_alpha = expr_subst(in_gamma, gamma_sym, t->gamma_render);
    expr_free(gamma_sym);
    expr_free(in_gamma);

    Expr* expanded = expr_expand(in_alpha);
    expr_free(in_alpha);

    Expr* canon = evaluate(expanded);
    expr_free(expanded);

    /* Normalise away the spurious rational content factor the primitive-element
     * tower arithmetic can leave behind, so the result is monic in the
     * polynomial variables (matching Wolfram and the sibling passing tests).
     * only_noninteger = false: over a field extension the result should be
     * monic even when the stray factor is a plain integer. */
    canon = qa_make_poly_numerically_monic(canon, false);
    return canon;
}

Expr* qa_polynomialgcd_with_tower_multivar(Expr* const* argv, size_t argc,
                                           const QATower* t) {
    return tower_multivar_combine(argv, argc, t, "PolynomialGCD");
}

Expr* qa_polynomiallcm_with_tower_multivar(Expr* const* argv, size_t argc,
                                           const QATower* t) {
    return tower_multivar_combine(argv, argc, t, "PolynomialLCM");
}

bool qa_tower_has_nested_radical(const QATower* t) {
    if (!t) return false;
    for (int i = 0; i < t->n; i++) {
        const Expr* a = t->alpha_renders[i];
        if (!a || a->type != EXPR_FUNCTION) continue;
        if (!a->data.function.head ||
            a->data.function.head->type != EXPR_SYMBOL) continue;
        const char* h = a->data.function.head->data.symbol.name;
        const Expr* base = NULL;
        if (h == SYM_Sqrt && a->data.function.arg_count == 1) {
            base = a->data.function.args[0];
        } else if (h == SYM_Power && a->data.function.arg_count == 2) {
            base = a->data.function.args[0];
        }
        if (base && base->type != EXPR_INTEGER
                 && base->type != EXPR_BIGINT) {
            return true;
        }
    }
    return false;
}

/* Lightweight free-polynomial-variable detector: returns true on first
 * non-numeric, non-radical leaf. Mirrors collect_variables's traversal
 * rules (Plus/Times/List/integer-Power recurse into args; everything
 * else is a "variable" if it isn't a number or radical). Skips
 * allocation by short-circuiting on first find. */
static bool expr_has_free_var(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT ||
        e->type == EXPR_REAL || e->type == EXPR_STRING) return false;
    if (e->type == EXPR_SYMBOL) {
        const char* s = e->data.symbol.name;
        if (!s) return false;
        /* Mathematical constants are not free variables. */
        if (s == SYM_Pi || s == SYM_E || s == SYM_I ||
            s == SYM_EulerGamma || s == SYM_Catalan ||
            s == SYM_Degree || s == SYM_GoldenRatio ||
            s == SYM_Infinity || s == SYM_True || s == SYM_False) return false;
        return true;
    }
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        /* Sqrt[c] / Power[c, p/q] with integer-only base is a "radical
         * constant", not a variable. Walk into the base to confirm. */
        if (h == SYM_Sqrt && e->data.function.arg_count == 1) {
            return expr_has_free_var(e->data.function.args[0]);
        }
        if (h == SYM_Power && e->data.function.arg_count == 2) {
            int64_t p, q;
            if (is_rational(e->data.function.args[1], &p, &q) && q >= 2) {
                return expr_has_free_var(e->data.function.args[0]);
            }
            /* Integer-exponent Power: recurse into base only (matches
             * collect_variables). */
            if (e->data.function.args[1]->type == EXPR_INTEGER) {
                return expr_has_free_var(e->data.function.args[0]);
            }
        }
        if (h == SYM_Plus || h == SYM_Times || h == SYM_List ||
            h == SYM_Rational || h == SYM_Complex) {
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (expr_has_free_var(e->data.function.args[i])) return true;
            }
            return false;
        }
    }
    /* Other function heads (Sin, Log, user-defined, ...): treat as
     * opaque variable atoms — they participate in polynomial structure
     * as opaque generators. */
    return true;
}

/* Lightweight radical detector: any Sqrt[*] or Power[*, p/q] (q ≥ 2). */
static bool expr_contains_any_radical(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_Sqrt && e->data.function.arg_count == 1) return true;
        if (h == SYM_Power && e->data.function.arg_count == 2) {
            int64_t p, q;
            if (is_rational(e->data.function.args[1], &p, &q) && q >= 2)
                return true;
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (expr_contains_any_radical(e->data.function.args[i])) return true;
    }
    return false;
}

static bool expr_has_nested_radical_radicand_walk(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        const Expr* base = NULL;
        bool is_radical_node = false;
        if (h == SYM_Sqrt && e->data.function.arg_count == 1) {
            base = e->data.function.args[0];
            is_radical_node = true;
        } else if (h == SYM_Power && e->data.function.arg_count == 2) {
            int64_t p, q;
            if (is_rational(e->data.function.args[1], &p, &q) && q >= 2) {
                base = e->data.function.args[0];
                is_radical_node = true;
            }
        }
        if (is_radical_node && base
            && base->type != EXPR_INTEGER
            && base->type != EXPR_BIGINT
            && expr_contains_any_radical(base)) {
            return true;
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (expr_has_nested_radical_radicand_walk(e->data.function.args[i]))
            return true;
    }
    return false;
}

bool expr_has_nested_radical_radicand(const Expr* e) {
    /* Only nested-radical inputs with NO free polynomial variable can be
     * safely bailed.  With free vars present (the canonical case is
     * `D[Integrate[a x/(x^3+2), x], x]` with `a`, `x`), the tower path
     * runs qaupoly_gcd over Q(γ)[x] and is the only thing that recovers
     * a closed form — even when the substituted form transiently
     * inflates, the GCD-based cancellation collapses it back.
     *
     * For variable-free nested-radical inputs (the user's
     * `1/(5+2 Sqrt[6])^(3/2) + (Sqrt[2]-Sqrt[3])^3` shape) qaupoly_gcd
     * has nothing to bite on — both num and den are degree-0 in the
     * dummy variable — so the tower path's substitution + Together
     * + lift + GCD-on-constants pipeline does no useful work, just
     * costs ~150 ms before the leaf-count gate rejects. */
    if (!expr_has_nested_radical_radicand_walk(e)) return false;
    return !expr_has_free_var(e);
}

/* ============================== Phase G8 ============================== */
/* Nested radical generators: `Sqrt[base]` / `base^(1/n)` where `base`
 * is a polynomial expression in atomic radicals (Sqrt[c], c^(1/m), I).
 * See plans/FACTOR_PLAN.md §14 / Phase G8 for the design. */

/* Atomic-algebraic recogniser: returns true iff `e` is one of the four
 * surface forms the original `qa_resolve_extension` accepts:
 *   - the imaginary unit (bare I or Complex[0, 1])
 *   - Sqrt[c] / Power[c, 1/2] with integer c
 *   - I·Sqrt[c]   (auto-evaluated form of Sqrt[-c])
 *   - c^(1/n) with integer c, integer n ≥ 2
 *
 * This is a strict-shape predicate: it does NOT recurse into nested
 * radicals.  By construction, anything passing this test resolves
 * through the original (non-G8) branches of qa_resolve_extension. */
static bool expr_is_atomic_algebraic(const Expr* e) {
    if (!e) return false;

    if (expr_is_imaginary_unit(e)) return true;

    long c;
    if (expr_is_sqrt_int(e, &c)) return true;

    /* I·Sqrt[c]  =  Sqrt[-c]  (Mathilda auto-evaluation). */
    if (e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Times
        && e->data.function.arg_count == 2) {
        Expr* a = e->data.function.args[0];
        Expr* b = e->data.function.args[1];
        long cc = 0;
        if ((expr_is_imaginary_unit(a) && expr_is_sqrt_int(b, &cc)) ||
            (expr_is_imaginary_unit(b) && expr_is_sqrt_int(a, &cc))) {
            return cc > 0;
        }
    }

    /* Power[c, p/q] with integer c, |c| >= 2, after gcd-reducing p/q to
     * lowest terms p_red/q_red, q_red >= 2.  This includes the original
     * c^(1/n) case (p == 1) and now also c^(p/q) for any integer p
     * (e.g. Power[2, -2/3], Power[2, 2/3], Power[6, 5/3]).  Excludes
     * |c| < 2 (trivial roots of unity / 0) and the Sqrt[c] form
     * already matched above. */
    if (e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Power
        && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp  = e->data.function.args[1];
        int64_t p, q;
        if (base->type == EXPR_INTEGER
            && is_rational(exp, &p, &q) && q != 1) {
            if (q < 0) { p = -p; q = -q; }
            int64_t ap = p < 0 ? -p : p;
            int64_t a = ap, b = q;
            while (b) { int64_t r = a % b; a = b; b = r; }
            int64_t q_red = a > 1 ? q / a : q;
            int64_t cv = base->data.integer;
            if (q_red >= 2 && cv != 0 && cv != 1 && cv != -1) return true;
        }
    }

    return false;
}

/* Push `e` onto a growing dedup'd Expr-list (deduplicated by expr_eq).
 * The list owns no references — entries point into the original tree
 * the caller is walking; callers must not free entries individually. */
static void exprlist_push_unique_borrowed(const Expr* e,
                                          const Expr*** list,
                                          int* n, int* cap) {
    for (int i = 0; i < *n; i++) {
        if (expr_eq((Expr*)(*list)[i], (Expr*)e)) return;
    }
    if (*n + 1 > *cap) {
        *cap = *cap ? (*cap) * 2 : 4;
        *list = (const Expr**)realloc(*list, sizeof(Expr*) * (*cap));
    }
    (*list)[(*n)++] = e;
}

/* Walk `e` and collect every subtree recognised by
 * `expr_is_atomic_algebraic` into `*atoms`.  Recursion stops at atomic
 * subtrees (so nested-inside-atomic content is not enumerated).
 *
 * Returns false if `e` contains any non-atomic, non-arithmetic
 * subexpression — i.e. an unrecognised radical-like Power node, an
 * unbound symbol, etc.  In that case the caller should reject the
 * input as out-of-MVP-scope.
 *
 * Accepted polynomial structure (the "skeleton" between atomics):
 * Plus, Times, Power[base, integer-exp], integer/rational/bigint
 * constants. */
static bool expr_collect_atomic_algebraics(const Expr* e,
                                           const Expr*** atoms,
                                           int* n_atoms, int* cap) {
    if (!e) return false;

    /* Constants are fine. */
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    int64_t p, q;
    if (is_rational((Expr*)e, &p, &q)) return true;

    /* An atomic algebraic generator: record and stop. */
    if (expr_is_atomic_algebraic(e)) {
        exprlist_push_unique_borrowed(e, atoms, n_atoms, cap);
        return true;
    }

    if (e->type != EXPR_FUNCTION
        || !e->data.function.head
        || e->data.function.head->type != EXPR_SYMBOL) return false;

    const char* head = e->data.function.head->data.symbol.name;

    /* Plus / Times: recurse on every argument. */
    if (head == SYM_Plus || head == SYM_Times) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (!expr_collect_atomic_algebraics(
                    e->data.function.args[i], atoms, n_atoms, cap)) {
                return false;
            }
        }
        return true;
    }

    /* Power[base, integer]: recurse into base. */
    if (head == SYM_Power && e->data.function.arg_count == 2) {
        Expr* exp = e->data.function.args[1];
        if (exp->type == EXPR_INTEGER) {
            return expr_collect_atomic_algebraics(
                e->data.function.args[0], atoms, n_atoms, cap);
        }
        /* Power[base, non-integer]: only allowed if it's atomic itself
         * (handled above).  Anything else (nested radical inside the
         * skeleton, fractional Power of a non-rational base) is out
         * of MVP scope. */
        return false;
    }

    /* Anything else (Cos, Log, free symbol, ...) is unrecognised. */
    return false;
}

/* Build z^n  −  base_internal  as an Expr in (alpha_internal, z_internal).
 * Used as the second resultant argument when computing the absolute
 * minimal polynomial of α = base^(1/n) over Q. */
static Expr* expr_zn_minus_base(const Expr* base_internal,
                                int n, const char* z_name) {
    /* z^n */
    Expr* zpow_args[2] = { expr_new_symbol(z_name),
                           expr_new_integer((int64_t)n) };
    Expr* zpow = expr_new_function(expr_new_symbol(SYM_Power),
                                   zpow_args, 2);
    /* −base_internal */
    Expr* neg_args[2] = { expr_new_integer(-1),
                          expr_copy((Expr*)base_internal) };
    Expr* neg_base = expr_new_function(expr_new_symbol(SYM_Times),
                                       neg_args, 2);
    /* z^n + (−base_internal) */
    Expr* sum_args[2] = { zpow, neg_base };
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus),
                                  sum_args, 2);
    Expr* expanded = expr_expand(sum);
    expr_free(sum);
    return expanded;
}

/* Reduce R(z) ∈ Q[z] to its squarefree part: R / gcd(R, R').  Returns
 * a new Expr (caller owns).  Returns a copy of R unchanged when R is
 * already squarefree. */
static Expr* expr_qx_squarefree_part(const Expr* R, const char* z_name) {
    /* dR/dz */
    Expr* d_args[2] = { expr_copy((Expr*)R), expr_new_symbol(z_name) };
    Expr* dR = internal_call_impl("D", builtin_d, d_args, 2);
    Expr* dR_eval = evaluate(dR);
    expr_free(dR);

    /* gcd(R, R') */
    Expr* g_args[2] = { expr_copy((Expr*)R), dR_eval };
    Expr* g = internal_polynomialgcd(g_args, 2);
    Expr* g_eval = evaluate(g);
    expr_free(g);

    /* If gcd is degree 0 (a nonzero constant), R is already squarefree. */
    Expr* z_sym = expr_new_symbol(z_name);
    int gdeg = get_degree_poly(g_eval, z_sym);
    if (gdeg <= 0) {
        expr_free(z_sym);
        expr_free(g_eval);
        return expr_copy((Expr*)R);
    }

    /* R / gcd(R, R') via PolynomialQuotient. */
    Expr* q_args[3] = { expr_copy((Expr*)R), g_eval, z_sym };
    Expr* qexpr = internal_polynomialquotient(q_args, 3);
    Expr* qeval = evaluate(qexpr);
    expr_free(qexpr);
    return qeval;
}

/* True if R ∈ Q[z] is irreducible over Q.  Conservative: factors via
 * Mathilda's `Factor`; returns true iff the Q-factorisation contains a
 * single non-constant factor (with multiplicity 1).  Returns false
 * (out parameter unchanged) on any analysis failure — caller treats
 * that as "not provably irreducible". */
static bool expr_qx_is_irreducible(const Expr* R, const char* z_name) {
    Expr* fac_args[1] = { expr_copy((Expr*)R) };
    Expr* factored = internal_factor(fac_args, 1);
    if (!factored) return false;
    Expr* fe = evaluate(factored);
    expr_free(factored);

    int n_nontrivial = 0;
    Expr* z_sym = expr_new_symbol(z_name);

    /* Walk the Factor result.  Times[poly1, poly2, ...] / Power[poly, k]
     * / a single poly all need handling.  We count distinct non-constant
     * factors, treating any Power exponent > 1 as multiplicity > 1
     * (which means the Q-factorisation is NOT irreducible). */
    if (fe->type == EXPR_FUNCTION
        && fe->data.function.head
        && fe->data.function.head->type == EXPR_SYMBOL
        && fe->data.function.head->data.symbol.name == SYM_Times) {
        for (size_t i = 0; i < fe->data.function.arg_count; i++) {
            Expr* a = fe->data.function.args[i];
            if (a->type == EXPR_FUNCTION
                && a->data.function.head
                && a->data.function.head->type == EXPR_SYMBOL
                && a->data.function.head->data.symbol.name == SYM_Power
                && a->data.function.arg_count == 2) {
                Expr* k = a->data.function.args[1];
                if (k->type == EXPR_INTEGER && k->data.integer > 1) {
                    expr_free(z_sym);
                    expr_free(fe);
                    return false;  /* multiplicity > 1 ⇒ not irreducible */
                }
                int d = get_degree_poly(a->data.function.args[0], z_sym);
                if (d > 0) n_nontrivial++;
            } else {
                int d = get_degree_poly(a, z_sym);
                if (d > 0) n_nontrivial++;
            }
        }
    } else if (fe->type == EXPR_FUNCTION
               && fe->data.function.head
               && fe->data.function.head->type == EXPR_SYMBOL
               && fe->data.function.head->data.symbol.name == SYM_Power
               && fe->data.function.arg_count == 2) {
        Expr* k = fe->data.function.args[1];
        if (k->type == EXPR_INTEGER && k->data.integer > 1) {
            expr_free(z_sym);
            expr_free(fe);
            return false;
        }
        int d = get_degree_poly(fe->data.function.args[0], z_sym);
        if (d > 0) n_nontrivial++;
    } else {
        int d = get_degree_poly(fe, z_sym);
        if (d > 0) n_nontrivial++;
    }

    expr_free(z_sym);
    expr_free(fe);
    return n_nontrivial == 1;
}

/* Phase G8 entry: resolve `Sqrt[base]` / `base^(1/n)` with non-rational
 * `base` to a (QAExt, render) pair via the algorithm in FACTOR_PLAN
 * §14 Phase G8. */
static QAExt* qa_resolve_nested_radical(const Expr* alpha_expr,
                                        Expr** render_out) {
    if (!alpha_expr || alpha_expr->type != EXPR_FUNCTION) return NULL;
    if (!alpha_expr->data.function.head
        || alpha_expr->data.function.head->type != EXPR_SYMBOL
        || alpha_expr->data.function.head->data.symbol.name != SYM_Power
        || alpha_expr->data.function.arg_count != 2) return NULL;

    Expr* base = alpha_expr->data.function.args[0];
    Expr* exp  = alpha_expr->data.function.args[1];

    int64_t p, q;
    if (!is_rational(exp, &p, &q) || p != 1 || q < 2) return NULL;
    int n = (int)q;

    /* 1) Collect atomic radicals from `base`.  Reject if `base`
     * contains any unrecognised structure. */
    const Expr** atoms = NULL;
    int n_atoms = 0;
    int cap = 0;
    if (!expr_collect_atomic_algebraics(base, &atoms, &n_atoms, &cap)) {
        free(atoms);
        return NULL;
    }
    if (n_atoms == 0) {
        /* `base` is purely Q-rational — caller's c^(1/n) branch should
         * have caught this when c was integer.  Reject non-integer
         * rational `base` for now (not in MVP). */
        free(atoms);
        return NULL;
    }

    /* 2) Build sub-tower Q(γ_sub) = Q(β_1, ..., β_k). */
    Expr** atom_args = (Expr**)malloc(sizeof(Expr*) * n_atoms);
    for (int i = 0; i < n_atoms; i++) {
        atom_args[i] = (Expr*)atoms[i];   /* qa_resolve_extension_tower
                                           * borrows; we own the tree
                                           * the atoms point into. */
    }
    QATower* sub = qa_resolve_extension_tower(atom_args, n_atoms);
    free(atom_args);
    free(atoms);
    if (!sub) return NULL;

    /* 3) Substitute each atom in `base` with its Q(γ_sub) representation
     * (a polynomial in QA_ALPHA_INTERNAL).
     *
     * For each natural-form atom `Power[c_i, 1/q_i]`, we first rewrite
     * any `Power[c_i, p/q]` in `base` into a polynomial in a per-atom
     * temp symbol via `expand_radicals_to_atomic_poly`, then substitute
     * that temp symbol with the gamma-representation of the atom.  The
     * two-step approach is needed because Mathilda's Times canonicaliser
     * may have rewritten the user's `Power[c_i, 1/q_i]` to non-natural
     * forms like `Power[c_i, -2/3]`, and structural substitution by
     * `sub->alpha_renders[i]` would miss those. */
    Expr* base_internal = expr_copy((Expr*)base);
    for (int i = 0; i < sub->n; i++) {
        /* Try to parse `sub->alpha_renders[i]` as `Power[c_i, 1/q_i]`.
         * When it has this shape, run the expand-radicals-to-atomic
         * preprocessing.  Otherwise (Sqrt, I, nested-of-nested), fall
         * back to the original structural substitution. */
        int64_t c_i = 0, q_i = 0;
        const Expr* a_r = sub->alpha_renders[i];
        if (a_r && a_r->type == EXPR_FUNCTION
            && a_r->data.function.head
            && a_r->data.function.head->type == EXPR_SYMBOL
            && a_r->data.function.head->data.symbol.name == SYM_Power
            && a_r->data.function.arg_count == 2
            && a_r->data.function.args[0]->type == EXPR_INTEGER) {
            int64_t pe, qe;
            if (is_rational(a_r->data.function.args[1], &pe, &qe)
                && pe == 1 && qe >= 2) {
                c_i = a_r->data.function.args[0]->data.integer;
                q_i = qe;
            }
        }
        /* Also accept Sqrt[c] as Power[c, 1/2]. */
        if (c_i == 0 && a_r && a_r->type == EXPR_FUNCTION
            && a_r->data.function.head
            && a_r->data.function.head->type == EXPR_SYMBOL
            && a_r->data.function.head->data.symbol.name == SYM_Sqrt
            && a_r->data.function.arg_count == 1
            && a_r->data.function.args[0]->type == EXPR_INTEGER) {
            c_i = a_r->data.function.args[0]->data.integer;
            q_i = 2;
        }

        Expr* atom_in_gamma_int = qanum_to_expr_in_gamma_sym(
            sub->alphas[i], QA_ALPHA_INTERNAL);

        if (c_i != 0 && c_i != 1 && c_i != -1 && q_i >= 2) {
            /* Per-atom temp symbol for the expansion result.  Then
             * substitute temp_sym → atom_in_gamma_int.  We use a fixed
             * sentinel per loop iteration; intern_symbol ensures
             * pointer-equality matches. */
            char tmp_name[64];
            snprintf(tmp_name, sizeof(tmp_name), "$qa$g8tmp$%d$", i);
            const char* tmp_iname = intern_symbol(tmp_name);

            QAExt* atom_ext = qaext_root_si((long)c_i, (unsigned)q_i);
            Expr* expanded = expand_radicals_to_atomic_poly(
                base_internal, c_i, q_i, tmp_iname, atom_ext);
            qaext_free(atom_ext);

            Expr* tmp_sym_expr = expr_new_symbol(tmp_iname);
            Expr* substituted = expr_subst(expanded, tmp_sym_expr,
                                           atom_in_gamma_int);
            expr_free(tmp_sym_expr);
            expr_free(expanded);
            expr_free(base_internal);
            base_internal = substituted;
        } else {
            /* Non-Power-shape atom (e.g. Sqrt[non-integer], I): fall
             * back to the original structural substitution. */
            Expr* old = base_internal;
            base_internal = expr_subst(old, sub->alpha_renders[i],
                                       atom_in_gamma_int);
            expr_free(old);
        }
        expr_free(atom_in_gamma_int);
    }
    Expr* base_canon = evaluate(expr_expand(base_internal));
    expr_free(base_internal);
    base_internal = base_canon;

    /* 4) Compute Res_w(P_sub(w), z^n − base(w))  ∈  Q[z]. */
    const char* w_name = QA_ALPHA_INTERNAL;
    const char* z_name = QA_Z_INTERNAL;

    Expr* P_sub_in_w = qaext_min_poly_expr(sub->ext, w_name);
    Expr* zn_minus_base = expr_zn_minus_base(base_internal, n, z_name);

    Expr* res_args[3] = { P_sub_in_w,
                          zn_minus_base,
                          expr_new_symbol(w_name) };
    Expr* min_poly_raw = internal_resultant(res_args, 3);
    expr_free(base_internal);

    if (!min_poly_raw) {
        qa_tower_free(sub);
        return NULL;
    }
    Expr* min_poly = expr_expand(min_poly_raw);
    expr_free(min_poly_raw);

    /* 5) Take the squarefree part — drops any spurious multiplicity
     * coming from algebraic dependencies among γ_sub conjugates. */
    Expr* min_poly_sf = expr_qx_squarefree_part(min_poly, z_name);
    expr_free(min_poly);

    /* 6) MVP guard: require the squarefree part to be Q-irreducible.
     * If not, α actually lies in a smaller subextension or has
     * conjugates that split — out of MVP scope for now. */
    if (!expr_qx_is_irreducible(min_poly_sf, z_name)) {
        expr_free(min_poly_sf);
        qa_tower_free(sub);
        return NULL;
    }

    /* 7) Build the QAExt from the irreducible min poly. */
    QAExt* ext = qaext_from_q_expr(min_poly_sf, z_name);
    expr_free(min_poly_sf);
    qa_tower_free(sub);

    if (!ext) return NULL;
    *render_out = expr_copy((Expr*)alpha_expr);
    return ext;
}

/* ============================== Phase G9 ============================== */
/* Automatic algebraic-extension detection — extension_autodetect.        */

typedef enum {
    GEN_INT_BASE,    /* Power[integer, 1/q_lcm], represented compactly */
    GEN_NESTED,     /* Power[non-integer-base, 1/q]: stored by borrowed Expr */
    GEN_ABSORBED    /* Nested generator whose canonical surface simplified
                     * to a Times of existing-generator powers via the Power
                     * canonicaliser (e.g. Sqrt[1/3/2^(2/3)] -> 1/(2^(1/3)
                     * Sqrt[3])). Marker only; dropped by post-canon dedup. */
} AutodetectGenKind;

typedef struct {
    AutodetectGenKind kind;
    /* GEN_INT_BASE: */
    int64_t base;        /* integer base of the radical, |base| >= 2 */
    int64_t q_lcm;       /* LCM of exponent denominators seen for `base` */
    /* GEN_NESTED: */
    const Expr* render;  /* surface form Power[base, 1/q] — borrowed pointer
                          * into the input expression; owner of input owns. */
} AutodetectGen;

static int64_t autodetect_gcd_i64(int64_t a, int64_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int64_t r = a % b; a = b; b = r; }
    return a;
}

static int64_t autodetect_lcm_i64(int64_t a, int64_t b) {
    if (a == 0 || b == 0) return 0;
    int64_t g = autodetect_gcd_i64(a, b);
    return (a / g) * b;
}

/* Merge an integer-base (base, q) into the generator set.  Returns false
 * (with *bail = true) if the set is full and `base` isn't already
 * present.  q must satisfy q >= 2. */
static bool autodetect_add_int(AutodetectGen* gens, size_t* n, size_t max,
                               int64_t base, int64_t q, bool* bail) {
    for (size_t i = 0; i < *n; i++) {
        if (gens[i].kind == GEN_INT_BASE && gens[i].base == base) {
            gens[i].q_lcm = autodetect_lcm_i64(gens[i].q_lcm, q);
            return true;
        }
    }
    if (*n >= max) { *bail = true; return false; }
    gens[*n].kind = GEN_INT_BASE;
    gens[*n].base  = base;
    gens[*n].q_lcm = q;
    gens[*n].render = NULL;
    (*n)++;
    return true;
}

/* Merge a nested-radical surface form into the generator set.  Dedup is
 * by structural Expr equality on the full Power[base, 1/q] form. */
static bool autodetect_add_nested(AutodetectGen* gens, size_t* n, size_t max,
                                  const Expr* surface, bool* bail) {
    for (size_t i = 0; i < *n; i++) {
        if (gens[i].kind == GEN_NESTED && gens[i].render
            && expr_eq((Expr*)gens[i].render, (Expr*)surface)) {
            return true;
        }
    }
    if (*n >= max) { *bail = true; return false; }
    gens[*n].kind = GEN_NESTED;
    gens[*n].base = 0;
    gens[*n].q_lcm = 0;
    gens[*n].render = surface;
    (*n)++;
    return true;
}

/* INT-only sub-walker: recurse into `e` collecting integer-base
 * generators ONLY for `Power[c, p/q]` forms that G8's
 * `expr_is_atomic_algebraic` would NOT accept — i.e., p != 1 (or
 * negative p).  The atomic `c^(1/n)` and `Sqrt[c]` forms are handled
 * by G8's own sub-tower construction on the parent GEN_NESTED, so
 * surfacing them again at the top level creates a redundant compositum
 * (e.g. Sqrt[1+Sqrt[2]] would become a degree-8 tower with both
 * Sqrt[2] and Sqrt[1+Sqrt[2]] when the correct degree is 4).
 *
 * The point of this walker is to surface integer bases needed by
 * `autodetect_canonicalise_radicand`, which would otherwise have
 * nothing to match against. */
static void autodetect_walk_intbase_only(const Expr* e, AutodetectGen* gens,
                                         size_t* n, size_t max, bool* bail) {
    if (!e || *bail) return;
    if (e->type != EXPR_FUNCTION) return;
    if (!e->data.function.head
        || e->data.function.head->type != EXPR_SYMBOL) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            autodetect_walk_intbase_only(e->data.function.args[i],
                                         gens, n, max, bail);
            if (*bail) return;
        }
        return;
    }
    const char* head = e->data.function.head->data.symbol.name;
    if (head == SYM_Sqrt && e->data.function.arg_count == 1) {
        /* Sqrt[c] with integer c is atomic for G8 — skip. */
        return;
    }
    if (head == SYM_Power && e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp_e = e->data.function.args[1];
        int64_t p, q;
        if (is_rational((Expr*)exp_e, &p, &q) && q != 1) {
            if (q < 0) { p = -p; q = -q; }
            int64_t ap = p < 0 ? -p : p;
            int64_t g  = autodetect_gcd_i64(ap, q);
            if (g > 1) { p /= g; q /= g; }
            /* Only surface for non-atomic exponent form (p != 1).
             * Atomic c^(1/q) is handled by G8 — surfacing would
             * double-count. */
            if (q >= 2 && p != 1 && base->type == EXPR_INTEGER) {
                int64_t c = base->data.integer;
                if (c != 0 && c != 1 && c != -1) {
                    if (!autodetect_add_int(gens, n, max, c, q, bail)) return;
                }
            }
            return;
        }
    }
    /* Recurse into args (head walk is unusual). */
    autodetect_walk_intbase_only(e->data.function.head, gens, n, max, bail);
    for (size_t i = 0; i < e->data.function.arg_count && !*bail; i++) {
        autodetect_walk_intbase_only(e->data.function.args[i],
                                     gens, n, max, bail);
    }
}

/* Walker: scan `e` for `Power[base, p/q]` / `Sqrt[base]` sub-expressions
 * and accumulate generators into `gens[]`.
 *
 * Integer bases are merged by base under the LCM of exponent denominators
 * (e.g. 2^(1/3) and 2^(1/2) -> one generator with q_lcm = 6).
 *
 * Non-integer bases (rational, polynomial, nested-radical) are surfaced
 * as GEN_NESTED generators, dispatched at tower-build time to G8
 * (`qa_resolve_nested_radical`).  The walker does NOT recurse into a
 * nested-radical's base — the G8 machinery builds its own sub-tower from
 * the base's atomic radicals.  However, integer-base generators inside
 * the base are surfaced via `autodetect_walk_intbase_only` so the
 * canonicalise_post pass can use them.
 *
 * Sets `*bail` on generator-array overflow. */
static void autodetect_walk(const Expr* e, AutodetectGen* gens,
                            size_t* n, size_t max, bool* bail) {
    if (!e || *bail) return;
    if (e->type != EXPR_FUNCTION) return;
    if (!e->data.function.head
        || e->data.function.head->type != EXPR_SYMBOL) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            autodetect_walk(e->data.function.args[i], gens, n, max, bail);
            if (*bail) return;
        }
        return;
    }

    const char* head = e->data.function.head->data.symbol.name;
    const Expr* base = NULL;
    int64_t p = 0, q = 1;
    bool is_radical = false;

    if (head == SYM_Sqrt && e->data.function.arg_count == 1) {
        base = e->data.function.args[0];
        p = 1; q = 2;
        is_radical = true;
    } else if (head == SYM_Power && e->data.function.arg_count == 2) {
        const Expr* exp_e = e->data.function.args[1];
        if (is_rational((Expr*)exp_e, &p, &q) && q != 1) {
            base = e->data.function.args[0];
            is_radical = true;
        }
    }

    if (is_radical) {
        /* Normalise sign of denominator and reduce p/q to lowest terms. */
        if (q < 0) { p = -p; q = -q; }
        int64_t ap = p < 0 ? -p : p;
        int64_t g  = autodetect_gcd_i64(ap, q);
        if (g > 1) { p /= g; q /= g; }

        if (q >= 2) {
            if (base && base->type == EXPR_INTEGER) {
                /* Integer base: merge into the (base, q_lcm) generator set.
                 * Base -1 is a root of unity (-1)^(1/q) = ζ_{2q}, resolved
                 * by qa_resolve_extension's cyclotomic branch; merging by
                 * the LCM of the q's gives a single ζ_{2·lcm} generator. */
                int64_t c = base->data.integer;
                if (c != 0 && c != 1) {
                    if (!autodetect_add_int(gens, n, max, c, q, bail)) return;
                }
                /* Integer base has no sub-radicals to chase. */
                return;
            } else {
                /* Non-integer base: surface the whole `e` as a nested-
                 * radical generator.  G8 handles polynomial-in-atomic-
                 * radicals bases at tower-build time.
                 *
                 * Phase F canonicalise-dedup: ALSO scan the base for
                 * integer-base generators (e.g. 2^(1/3) inside
                 * Sqrt[2^(1/3)/6]) so the canonicalise_post pass can use
                 * them to normalise radicands.  Inner nested radicals
                 * (Sqrt[Sqrt[2]]-style) are NOT re-surfaced here — G8
                 * builds its own sub-tower from them.  We use a
                 * dedicated INT-only sub-walker to avoid the double-
                 * count regression. */
                if (!autodetect_add_nested(gens, n, max, e, bail)) return;
                autodetect_walk_intbase_only(base, gens, n, max, bail);
                return;
            }
        }
        return;
    }

    /* Generic recursion: walk head (unusual) and every argument. */
    autodetect_walk(e->data.function.head, gens, n, max, bail);
    for (size_t i = 0; i < e->data.function.arg_count && !*bail; i++) {
        autodetect_walk(e->data.function.args[i], gens, n, max, bail);
    }
}

/* Floor-divide `p` by `q` (q > 0); result is the floor of p/q. */
static int64_t autodetect_floor_div(int64_t p, int64_t q) {
    int64_t a = p / q;
    int64_t b = p - a * q;
    if (b < 0) { a -= 1; b += q; }
    return a;
}

/* Recursively rewrite every `Power[c, p/q]` (c integer with |c| >= 2,
 * q >= 2) into `c^a · Power[c, b/q]` where p = a·q + b, 0 ≤ b < q.
 * The transformation is needed before passing canonicalised radicands
 * to G8 (`qa_resolve_nested_radical`), whose
 * `expr_is_atomic_algebraic` recogniser only accepts the pure
 * `Power[c, 1/n]` form, not `Power[c, p/q]` with p > 1 or p < 0. */
static Expr* autodetect_rewrite_to_atomic(const Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    if (e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Power
        && e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp_e = e->data.function.args[1];
        int64_t p, q;
        if (base->type == EXPR_INTEGER
            && is_rational((Expr*)exp_e, &p, &q) && q >= 2) {
            int64_t c = base->data.integer;
            if (c == 0 || c == 1 || c == -1) {
                return expr_copy((Expr*)e);
            }
            int64_t a = autodetect_floor_div(p, q);
            int64_t b = p - a * q;
            if (b == 0) {
                /* Power[c, a]: pure integer exponent. */
                return eval_and_free(expr_new_function(
                    expr_new_symbol(SYM_Power),
                    (Expr*[]){expr_new_integer(c), expr_new_integer(a)}, 2));
            }
            /* Power[c, a] * Power[c, b/q]  with 1 <= b < q. */
            Expr* int_part = (a == 0)
                ? expr_new_integer(1)
                : eval_and_free(expr_new_function(
                    expr_new_symbol(SYM_Power),
                    (Expr*[]){expr_new_integer(c), expr_new_integer(a)}, 2));
            Expr* rat_args[2] = { expr_new_integer(b), expr_new_integer(q) };
            Expr* exp_atomic = expr_new_function(expr_new_symbol(SYM_Rational),
                                                 rat_args, 2);
            Expr* atomic_part = expr_new_function(expr_new_symbol(SYM_Power),
                (Expr*[]){expr_new_integer(c), exp_atomic}, 2);
            return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){int_part, atomic_part}, 2));
        }
    }

    size_t count = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (count ? count : 1));
    for (size_t i = 0; i < count; i++) {
        new_args[i] = autodetect_rewrite_to_atomic(e->data.function.args[i]);
    }
    Expr* new_head = autodetect_rewrite_to_atomic(e->data.function.head);
    Expr* out = eval_and_free(expr_new_function(new_head, new_args, count));
    free(new_args);
    return out;
}

/* canonicalise_radicand: run `Together[u, Extension -> α]` against the
 * given integer-base generator α (`Power[base, 1/q_lcm]`) so that u is
 * rewritten to its canonical c0 + c1·α + c2·α² + ... form modulo the
 * algebraic relation α^q = base.  Then apply `autodetect_rewrite_to_atomic`
 * so any `Power[c, p/q]` with |p| > 1 is split into a pure `c^a · Power[c,
 * b/q]` form acceptable to G8's atomic-radical recogniser.
 *
 * Used by the post-walker canonical-dedup pass on nested radicals:
 * `Sqrt[-1/(9·2^(2/3)) + (2/9)·2^(1/3)]` and `Sqrt[1/(3·2^(2/3))]` both
 * reduce to `Sqrt[2^(1/3)/6]` after canonicalisation against
 * α = 2^(1/3).
 *
 * Returns a freshly-allocated Expr (caller frees via expr_free), or
 * NULL on canonicalisation failure. */
static Expr* autodetect_canonicalise_radicand(const Expr* u,
                                              int64_t base, int64_t q_lcm) {
    Expr* base_e = expr_new_integer(base);
    Expr* rat_args[2] = { expr_new_integer(1), expr_new_integer(q_lcm) };
    Expr* exp_e = expr_new_function(expr_new_symbol(SYM_Rational), rat_args, 2);
    Expr* pow_args[2] = { base_e, exp_e };
    Expr* alpha = expr_new_function(expr_new_symbol(SYM_Power), pow_args, 2);

    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){expr_new_symbol(SYM_Extension), alpha}, 2);
    Expr* tg_call = expr_new_function(expr_new_symbol(SYM_Together),
        (Expr*[]){expr_copy((Expr*)u), rule}, 2);
    Expr* canon = evaluate(tg_call);
    expr_free(tg_call);
    if (!canon) return NULL;

    Expr* atomic = autodetect_rewrite_to_atomic(canon);
    expr_free(canon);
    return atomic;
}

/* True when `e` is a single atomic algebraic surface, i.e. `Sqrt[X]` or
 * `Power[X, p/q]` with reduced denominator q >= 2.  Used to detect
 * generators that the Power canonicaliser absorbed into existing
 * integer-base generators (e.g. `Sqrt[1/3/2^(2/3)]` -> `1/(2^(1/3)
 * Sqrt[3])`, a Times of inverse powers — no longer a single
 * algebraic generator). */
static bool autodetect_is_atomic_radical_surface(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (!e->data.function.head
        || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (h == SYM_Sqrt && e->data.function.arg_count == 1) return true;
    if (h == SYM_Power && e->data.function.arg_count == 2) {
        int64_t p, q;
        if (is_rational((Expr*)e->data.function.args[1], &p, &q)) {
            if (q < 0) { p = -p; q = -q; }
            int64_t ap = p < 0 ? -p : p;
            int64_t a = ap, b = q;
            while (b) { int64_t r = a % b; a = b; b = r; }
            int64_t q_red = a > 1 ? q / a : q;
            if (q_red >= 2) return true;
        }
    }
    return false;
}

/* For each GEN_NESTED generator: try to canonicalise its radicand
 * against each GEN_INT_BASE generator in the set, replacing the borrowed
 * render with a freshly-allocated canonical surface form when the
 * canonicalisation produces something structurally different.  Owned
 * canonical renders are written to `owned[]` so callers can free them.
 *
 * If the canonical surface is no longer a single atomic radical (e.g.
 * the Power canonicaliser simplified it to a Times of existing-generator
 * powers), the generator has been ABSORBED into the integer-base
 * generators we already have.  In that case: harvest any newly-exposed
 * components by re-walking the canonical surface (this can merge with or
 * widen existing INT_BASE entries via autodetect_add_int), and mark the
 * absorbed entry for removal in the dedup step.
 *
 * After canonicalisation, re-deduplicate the GEN_NESTED entries:
 * structurally-equal renders (now reflecting canonical forms) collapse
 * to a single generator.  Absorbed entries are also dropped.  Re-dedup
 * is in-place by sweeping the array and compacting. */
static void autodetect_canonicalise_post(AutodetectGen* gens, size_t* n,
                                         Expr** owned) {
    /* Step 1: canonicalise each nested radical against every integer-base
     * generator we have, AND normalise the surface form to use a positive
     * 1/q exponent.  Mathilda's Power canonicaliser rewrites `1/Sqrt[X]`
     * to `Power[X, -1/2]`, so the same algebraic generator can appear
     * in the input with negative or positive exponent — both generate
     * the same field Q(Sqrt[X]) and should dedup. */
    for (size_t i = 0; i < *n; i++) {
        if (gens[i].kind != GEN_NESTED) continue;
        const Expr* surface = gens[i].render;
        if (!surface || surface->type != EXPR_FUNCTION
            || surface->data.function.arg_count == 0) continue;

        /* Parse the surface form: either `Sqrt[radicand]` (1 arg) or
         * `Power[radicand, p/q]` (2 args, rational exp). */
        const Expr* radicand = NULL;
        int64_t q_natural = 0;
        if (surface->data.function.head
            && surface->data.function.head->type == EXPR_SYMBOL
            && surface->data.function.head->data.symbol.name == SYM_Sqrt
            && surface->data.function.arg_count == 1) {
            radicand = surface->data.function.args[0];
            q_natural = 2;
        } else if (surface->data.function.head
            && surface->data.function.head->type == EXPR_SYMBOL
            && surface->data.function.head->data.symbol.name == SYM_Power
            && surface->data.function.arg_count == 2) {
            const Expr* e_x = surface->data.function.args[1];
            int64_t pe, qe;
            if (is_rational((Expr*)e_x, &pe, &qe) && qe != 1) {
                if (qe < 0) { pe = -pe; qe = -qe; }
                int64_t ape = pe < 0 ? -pe : pe;
                int64_t a = ape, b = qe;
                while (b) { int64_t r = a % b; a = b; b = r; }
                int64_t q_red = a > 1 ? qe / a : qe;
                if (q_red >= 2) {
                    radicand = surface->data.function.args[0];
                    q_natural = q_red;
                }
            }
        }
        if (!radicand) continue;

        /* Try canonicalising the radicand against integer-base generators. */
        Expr* canon_radicand = NULL;
        for (size_t j = 0; j < *n; j++) {
            if (gens[j].kind != GEN_INT_BASE) continue;
            Expr* c = autodetect_canonicalise_radicand(
                radicand, gens[j].base, gens[j].q_lcm);
            if (c) { canon_radicand = c; break; }
        }
        /* If no canonicalisation happened, use a copy of the original
         * radicand so the rebuilt surface form has consistent positive
         * exponent. */
        if (!canon_radicand) canon_radicand = expr_copy((Expr*)radicand);

        /* Always rebuild surface form with positive 1/q_natural exponent.
         * This collapses Power[X, -1/q] and Sqrt[X] (= Power[X, 1/q] for
         * q=2) to the same canonical form when X is the same. */
        Expr* new_surface;
        if (q_natural == 2) {
            new_surface = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Sqrt),
                (Expr*[]){canon_radicand}, 1));
        } else {
            Expr* rat_args[2] = { expr_new_integer(1),
                                  expr_new_integer(q_natural) };
            Expr* exp_pos = expr_new_function(
                expr_new_symbol(SYM_Rational), rat_args, 2);
            new_surface = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Power),
                (Expr*[]){canon_radicand, exp_pos}, 2));
        }
        /* Only adopt the new render if it actually differs from the
         * old one — preserves the existing-test invariant when the
         * walker already produced a positive-exponent natural form. */
        if (!expr_eq(new_surface, (Expr*)surface)) {
            owned[i] = new_surface;
            gens[i].render = new_surface;

            /* If the canonical surface is no longer a single atomic
             * radical, the generator was absorbed into the existing
             * integer-base set.  Re-walk the simplified surface to
             * harvest any newly-exposed integer-base components (the
             * walk dedups against existing entries via
             * autodetect_add_int), then mark this entry as absorbed
             * so step 2's compaction drops it. */
            if (!autodetect_is_atomic_radical_surface(new_surface)) {
                bool harvest_bail = false;
                autodetect_walk(new_surface, gens, n,
                                QA_AUTODETECT_MAX_GENS, &harvest_bail);
                gens[i].kind = GEN_ABSORBED;
                /* harvest_bail (cap full) is non-fatal here: the
                 * generator is still dropped, we just won't surface
                 * any additional generators beyond the cap.  Future
                 * iterations of the outer loop may pick up any newly-
                 * appended GEN_NESTED entries. */
            }
        } else {
            expr_free(new_surface);
        }
    }

    /* Step 2: re-dedup + drop absorbed entries.  Sweep i from 0 upward;
     * for each i, decide whether to drop or keep.  Drop conditions:
     *   - kind == GEN_ABSORBED (canonical surface no longer atomic).
     *   - GEN_NESTED with a duplicate render among earlier entries. */
    size_t out = 0;
    for (size_t i = 0; i < *n; i++) {
        bool drop = false;
        if (gens[i].kind == GEN_ABSORBED) {
            drop = true;
        } else if (gens[i].kind == GEN_NESTED) {
            for (size_t j = 0; j < out; j++) {
                if (gens[j].kind == GEN_NESTED && gens[j].render && gens[i].render
                    && expr_eq((Expr*)gens[j].render, (Expr*)gens[i].render)) {
                    drop = true; break;
                }
            }
        }
        if (drop) {
            /* Free the owned canonical surface, if any, since we're
             * dropping this entry. */
            if (owned[i]) { expr_free(owned[i]); owned[i] = NULL; }
        } else {
            if (out != i) {
                gens[out] = gens[i];
                owned[out] = owned[i];
                owned[i] = NULL;
            }
            out++;
        }
    }
    *n = out;
}

/* Layer 2: redundant-Sqrt coalesce.
 *
 * extension_autodetect collects every integer-base Sqrt[c] it sees as a
 * separate GEN_INT_BASE generator. For c = p·q (squarefree, p,q prime),
 * `Sqrt[c]` is algebraically dependent on Sqrt[p] and Sqrt[q]:
 * `Sqrt[p·q] = Sqrt[p]·Sqrt[q]`. If both Sqrt[p] and Sqrt[q] are already
 * generators, including Sqrt[c] inflates the tower's primitive-element
 * degree (it would be deg-4 over Q(Sqrt[p],Sqrt[q]) but the tower-builder
 * doesn't know about the dependency and produces deg-8). The qa_cancel_
 * with_tower work cost is super-linear in deg(γ), so the inflation hurts.
 *
 * This pass marks every GEN_INT_BASE Sqrt[c] (q_lcm == 2, c > 1) as
 * GEN_ABSORBED when c factors into squarefree primes p_i, all already
 * present as GEN_INT_BASE Sqrt[p_i]. The downstream dedup compaction
 * drops the absorbed entries; qa_cancel_with_tower's caller-side
 * decomposition (decompose_redundant_sqrts) rewrites input occurrences
 * of Sqrt[c] as `Times[Sqrt[p_1], ..., Sqrt[p_k]]` so the substitution
 * loop still handles them correctly.
 *
 * Conservative scope: only squarefree positive integer bases with q=2.
 * Higher-q (cube roots etc.) and negative bases are left alone — their
 * dependence patterns are more involved and rarer in practice. */
static void autodetect_coalesce_int_sqrts(AutodetectGen* gens, size_t* n,
                                          Expr** owned) {
    /* Collect prime bases p where Sqrt[p] (q_lcm == 2) is a generator. */
    int64_t primes_in_set[QA_AUTODETECT_MAX_GENS];
    size_t  n_primes = 0;
    for (size_t i = 0; i < *n; i++) {
        if (gens[i].kind != GEN_INT_BASE) continue;
        if (gens[i].q_lcm != 2) continue;
        int64_t b = gens[i].base;
        if (b < 2) continue;
        /* Primality test by trial division. */
        bool is_prime = true;
        for (int64_t p = 2; p * p <= b; p++) {
            if (b % p == 0) { is_prime = false; break; }
        }
        if (is_prime) primes_in_set[n_primes++] = b;
    }
    if (n_primes < 2) return;  /* Need at least two primes to coalesce. */

    /* Mark non-prime GEN_INT_BASE Sqrt[c] absorbed when c factors into
     * squarefree primes all in primes_in_set. */
    for (size_t i = 0; i < *n; i++) {
        if (gens[i].kind != GEN_INT_BASE) continue;
        if (gens[i].q_lcm != 2) continue;
        int64_t c = gens[i].base;
        if (c < 2) continue;
        /* Skip primes (they're the irreducible generators themselves). */
        bool is_prime = true;
        for (int64_t p = 2; p * p <= c; p++) {
            if (c % p == 0) { is_prime = false; break; }
        }
        if (is_prime) continue;
        /* Trial-divide c, requiring each prime factor to appear exactly
         * once (squarefree) AND be in primes_in_set. */
        int64_t r = c;
        bool ok = true;
        for (int64_t p = 2; p * p <= r && ok; p++) {
            if (r % p == 0) {
                int cnt = 0;
                while (r % p == 0) { r /= p; cnt++; }
                if (cnt != 1) { ok = false; break; }
                bool found = false;
                for (size_t j = 0; j < n_primes; j++) {
                    if (primes_in_set[j] == p) { found = true; break; }
                }
                if (!found) { ok = false; break; }
            }
        }
        if (ok && r > 1) {
            bool found = false;
            for (size_t j = 0; j < n_primes; j++) {
                if (primes_in_set[j] == r) { found = true; break; }
            }
            if (!found) ok = false;
        }
        if (ok) {
            gens[i].kind = GEN_ABSORBED;
            /* No owned[] entry to free (GEN_INT_BASE doesn't allocate). */
            (void)owned;
        }
    }

    /* Compact: sweep + drop GEN_ABSORBED. */
    size_t out = 0;
    for (size_t i = 0; i < *n; i++) {
        if (gens[i].kind == GEN_ABSORBED) {
            if (owned[i]) { expr_free(owned[i]); owned[i] = NULL; }
            continue;
        }
        if (out != i) {
            gens[out] = gens[i];
            owned[out] = owned[i];
            owned[i] = NULL;
        }
        out++;
    }
    *n = out;
}

/* Materialise an AutodetectGen[] into a tower by building a surface-form
 * generator Expr for each entry and dispatching to
 * qa_resolve_extension_tower.  Used by both the single-Expr and multi-Expr
 * public entry points.  Returns NULL when n == 0 or the tower build fails. */
static QATower* autodetect_build_tower(const AutodetectGen* gens, size_t n) {
    if (n == 0) return NULL;

    Expr** alpha_exprs = (Expr**)malloc(sizeof(Expr*) * n);
    if (!alpha_exprs) return NULL;
    for (size_t i = 0; i < n; i++) {
        if (gens[i].kind == GEN_INT_BASE) {
            /* `Power[c, Rational[1, q]]` — qa_resolve_extension's
             * `c^(1/n)` branch accepts this directly. */
            Expr* base_e = expr_new_integer(gens[i].base);
            Expr* rat_args[2] = { expr_new_integer(1),
                                  expr_new_integer(gens[i].q_lcm) };
            Expr* exp_e = expr_new_function(expr_new_symbol(SYM_Rational),
                                            rat_args, 2);
            Expr* pow_args[2] = { base_e, exp_e };
            alpha_exprs[i] = expr_new_function(expr_new_symbol(SYM_Power),
                                               pow_args, 2);
        } else {
            /* GEN_NESTED: surface form is borrowed from the input.  Deep-
             * copy so qa_resolve_extension_tower can claim it (it borrows
             * its arguments but the tower's stored render is a copy). */
            alpha_exprs[i] = expr_copy((Expr*)gens[i].render);
        }
    }

    QATower* t = qa_resolve_extension_tower(alpha_exprs, (int)n);

    for (size_t i = 0; i < n; i++) expr_free(alpha_exprs[i]);
    free(alpha_exprs);
    return t;
}

QATower* extension_autodetect(const Expr* e) {
    if (!e) return NULL;

    AutodetectGen gens[QA_AUTODETECT_MAX_GENS];
    Expr*         owned[QA_AUTODETECT_MAX_GENS] = {NULL};
    size_t n = 0;
    bool bail = false;
    autodetect_walk(e, gens, &n, QA_AUTODETECT_MAX_GENS, &bail);
    if (bail) {
        for (size_t i = 0; i < QA_AUTODETECT_MAX_GENS; i++)
            if (owned[i]) expr_free(owned[i]);
        return NULL;
    }
    autodetect_canonicalise_post(gens, &n, owned);
    autodetect_coalesce_int_sqrts(gens, &n, owned);
    QATower* t = autodetect_build_tower(gens, n);
    for (size_t i = 0; i < QA_AUTODETECT_MAX_GENS; i++)
        if (owned[i]) expr_free(owned[i]);
    return t;
}

QATower* extension_autodetect_args(struct Expr* const* args, size_t argc) {
    if (!args || argc == 0) return NULL;

    AutodetectGen gens[QA_AUTODETECT_MAX_GENS];
    Expr*         owned[QA_AUTODETECT_MAX_GENS] = {NULL};
    size_t n = 0;
    bool bail = false;
    for (size_t i = 0; i < argc && !bail; i++) {
        autodetect_walk(args[i], gens, &n, QA_AUTODETECT_MAX_GENS, &bail);
    }
    if (bail) {
        for (size_t i = 0; i < QA_AUTODETECT_MAX_GENS; i++)
            if (owned[i]) expr_free(owned[i]);
        return NULL;
    }
    autodetect_canonicalise_post(gens, &n, owned);
    autodetect_coalesce_int_sqrts(gens, &n, owned);
    QATower* t = autodetect_build_tower(gens, n);
    for (size_t i = 0; i < QA_AUTODETECT_MAX_GENS; i++)
        if (owned[i]) expr_free(owned[i]);
    return t;
}

/* ============================== Phase E ============================== */
/* Single-generator polynomial-radicand Cancel/Together path.
 *
 * Handles inputs containing exactly one distinct radical of the form
 * `Sqrt[poly]` / `Power[poly, 1/q]` (q >= 2), where `poly` is a polynomial
 * expression with free symbols (e.g. `Sqrt[p+q]`, `Power[1+x^2, 1/3]`).
 *
 * These inputs are rejected by extension_autodetect because
 * qa_resolve_nested_radical's expr_collect_atomic_algebraics rejects free
 * symbols in the radicand.  Rather than enlarge the QAExt machinery to
 * support function-field coefficients, this path implements the simpler
 * substitute-then-reduce identity for the one-radical case:
 *
 *   1. Substitute every occurrence of α with a fresh symbol S.
 *   2. Run Together (no extension) on the substituted input — combines
 *      sum-of-fractions over Q[params, vars, S].
 *   3. Extract Numerator and Denominator.
 *   4. For each, run PolynomialRemainder[..., S^q - radicand, S]:
 *      reduces every S^k with k >= q via S^q = radicand.
 *   5. Compute PolynomialGCD[num_reduced, den_reduced] (over Q[params,
 *      vars, S]) and divide both sides by it to cancel any common
 *      factor surfaced by the reduction.
 *   6. Substitute S -> α back and emit num / den.
 *
 * Limitations (documented in changelog):
 *   - Single distinct radical only.  Multi-radical inputs like the
 *     Cardano identity (which involves three radicals whose product
 *     simplification requires conjugate-pair reasoning) need a Groebner-
 *     basis treatment that this path does not attempt.
 *   - No denominator rationalisation: the output may still carry α in
 *     the denominator.  `simp_rationalize_denom` (in src/simp/) does that
 *     for the integer-base analogue via PolynomialExtendedGCD; this path
 *     stops at the substitute + reduce step.
 *   - `Power[poly, p/q]` with reduced p != 1 is rejected (only the
 *     natural-form generator base^(1/q) is recognised).
 *
 * Returns NULL when detection fails, no actual simplification occurs,
 * or any intermediate evaluation step fails — caller falls back to the
 * non-extension path. */

/* Test if `e` is purely an integer / bigint / rational literal (no
 * symbols, no functions).  Used to filter integer-base radicals from the
 * polynomial-radicand detection — integer bases are handled by the main
 * extension_autodetect path. */
static bool polyrad_radicand_is_const(const Expr* e) {
    if (!e) return true;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    int64_t p, q;
    if (is_rational((Expr*)e, &p, &q)) return true;
    return false;
}

/* Test if `e` contains a free polynomial variable.  Returns true iff a
 * symbol other than the mathematical constants (Pi, E, I, ...) appears
 * somewhere in the tree.  Sqrt[c] / Power[c, p/q] with integer-only `c`
 * is treated as a "radical constant" (not a variable) — matching the
 * `expr_has_free_var` predicate defined above.  Used by Phase E to
 * reject radicals whose radicand is purely algebraic-over-Q (which
 * should flow through the standard `extension_autodetect` / G8 nested-
 * radical path, not Phase E). */
static bool polyrad_radicand_has_free_var(const Expr* e) {
    return expr_has_free_var(e);
}

/* Extract the radicand and the (reduced) exponent denominator q from a
 * radical surface form `e`.  Accepted: `Sqrt[base]` -> (base, 2),
 * `Power[base, Rational[1, q]]` -> (base, q), more generally
 * `Power[base, p/q]` -> (base, q_red) iff p/q reduces to 1/q_red.
 *
 * Returns false (with *out_base / *out_q untouched) for other shapes. */
static bool polyrad_parse_radical(const Expr* e,
                                  const Expr** out_base,
                                  int64_t* out_q) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (!e->data.function.head
        || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (h == SYM_Sqrt && e->data.function.arg_count == 1) {
        *out_base = e->data.function.args[0];
        *out_q = 2;
        return true;
    }
    if (h == SYM_Power && e->data.function.arg_count == 2) {
        const Expr* exp = e->data.function.args[1];
        int64_t p, q;
        if (!is_rational((Expr*)exp, &p, &q) || q == 1) return false;
        if (q < 0) { p = -p; q = -q; }
        int64_t ap = p < 0 ? -p : p;
        int64_t a = ap, b = q;
        while (b) { int64_t r = a % b; a = b; b = r; }
        int64_t g = a;
        int64_t p_red = g > 1 ? p / g : p;
        int64_t q_red = g > 1 ? q / g : q;
        /* Only natural-form 1/q is recognised — Power[poly, p/q] with
         * reduced p != 1 would need additional expand-radicals-to-atomic
         * machinery on the substitution side. */
        if (p_red != 1 || q_red < 2) return false;
        *out_base = e->data.function.args[0];
        *out_q = q_red;
        return true;
    }
    return false;
}

/* Walk `e` collecting distinct radical surface forms (by structural
 * equality).  Each entry is a borrowed pointer into the input tree.
 * Sets *overflow = true when more than `cap` distinct entries appear
 * (caller treats as detection-failed). */
static void polyrad_collect_distinct(const Expr* e,
                                     const Expr*** list,
                                     int* n, int* cap,
                                     bool* overflow) {
    if (!e || *overflow) return;
    if (e->type != EXPR_FUNCTION) return;
    const Expr* radicand = NULL;
    int64_t q = 0;
    if (polyrad_parse_radical(e, &radicand, &q)
        && !polyrad_radicand_is_const(radicand)
        && polyrad_radicand_has_free_var(radicand)) {
        /* Dedup by structural equality on the whole radical surface. */
        for (int i = 0; i < *n; i++) {
            if (expr_eq((Expr*)(*list)[i], (Expr*)e)) goto recurse;
        }
        if (*n + 1 > *cap) {
            int new_cap = (*cap == 0) ? 4 : (*cap) * 2;
            const Expr** new_list = (const Expr**)realloc(
                *list, sizeof(Expr*) * (size_t)new_cap);
            if (!new_list) { *overflow = true; return; }
            *list = new_list;
            *cap = new_cap;
        }
        (*list)[(*n)++] = e;
        /* Threshold: more than 1 distinct radical -> bail. We could
         * extend this with a triangular-ideal reduction later. */
        if (*n > 1) { *overflow = true; return; }
    }
recurse:
    if (e->data.function.head)
        polyrad_collect_distinct(e->data.function.head, list, n, cap, overflow);
    for (size_t i = 0; i < e->data.function.arg_count && !*overflow; i++) {
        polyrad_collect_distinct(e->data.function.args[i], list, n, cap, overflow);
    }
}

/* Leaf-count helper for the post-substitution change detector.  Local
 * minimal version that matches `leaf_count_internal`'s definition
 * sufficiently for the "did this simplify?" check: every integer,
 * symbol, and float counts as 1; function nodes count their head + each
 * arg leaf. */
static int64_t polyrad_leaf_count(const Expr* e) {
    if (!e) return 0;
    if (e->type != EXPR_FUNCTION) return 1;
    int64_t c = polyrad_leaf_count(e->data.function.head);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        c += polyrad_leaf_count(e->data.function.args[i]);
    }
    return c;
}

/* Sign of a polynomial term's numeric coefficient: -1, 0, or +1. A bare
 * number reports its own sign; a Times reports the sign of its (canonically
 * leading) numeric factor; anything else has an implicit +1 coefficient. */
static int polyrad_term_sign(const Expr* t) {
    if (!t) return 1;
    if (t->type == EXPR_INTEGER)
        return t->data.integer < 0 ? -1 : (t->data.integer > 0 ? 1 : 0);
    if (t->type == EXPR_REAL)
        return t->data.real < 0.0 ? -1 : (t->data.real > 0.0 ? 1 : 0);
    if (t->type == EXPR_BIGINT)
        return mpz_sgn(t->data.bigint);
    { int64_t np, nq; if (is_rational(t, &np, &nq)) return np < 0 ? -1 : 1; }
    if (t->type == EXPR_FUNCTION && t->data.function.head
        && t->data.function.head->type == EXPR_SYMBOL
        && t->data.function.head->data.symbol.name == SYM_Times
        && t->data.function.arg_count >= 1) {
        return polyrad_term_sign(t->data.function.args[0]);
    }
    return 1;
}

/* True when the leading (highest-degree) term of `den` carries a negative
 * coefficient.  Mathilda's canonical Plus ordering places the leading
 * monomial last, so the last summand is inspected; a non-Plus is treated as
 * its own leading term. */
static bool polyrad_leading_coeff_is_negative(const Expr* den) {
    if (!den) return false;
    const Expr* lead = den;
    if (den->type == EXPR_FUNCTION && den->data.function.head
        && den->data.function.head->type == EXPR_SYMBOL
        && den->data.function.head->data.symbol.name == SYM_Plus
        && den->data.function.arg_count >= 1) {
        lead = den->data.function.args[den->data.function.arg_count - 1];
    }
    return polyrad_term_sign(lead) < 0;
}

/* Canonicalise the sign of a rational result so its denominator's leading term
 * is positive (Wolfram's Together convention): if the denominator leads with a
 * negative coefficient, negate numerator and denominator together (value-
 * preserving).  A non-fraction, or one whose denominator already leads
 * positive, is returned unchanged.  `e` is consumed.
 *
 * The algebraic Together engine (flint_algebraic_field_together) can emit the
 * sign-flipped associate -(2 x)/((p+q)^(2/3) - x^2) where the equivalent
 * Sqrt-radicand case yields the canonical (2 x)/(x^2 - (p+q)^(2/3)); this makes
 * the two consistent. */
Expr* qa_normalize_fraction_sign(Expr* e) {
    if (!e) return NULL;
    Expr* den = eval_and_free(expr_new_function(expr_new_symbol(SYM_Denominator),
        (Expr*[]){ expr_copy(e) }, 1));
    bool neg = polyrad_leading_coeff_is_negative(den);
    /* Only canonicalise the sign of a genuine *polynomial* denominator.  A
     * purely numeric denominator (e.g. 1 - 2^(2/3), a single negative constant)
     * has no leading monomial to make positive; flipping it would gratuitously
     * disagree with the explicit Extension -> alpha form (2/(1 - 2^(2/3))). */
    bool has_var = false;
    if (neg && den) {
        Expr* vars = eval_and_free(expr_new_function(expr_new_symbol(SYM_Variables),
            (Expr*[]){ expr_copy(den) }, 1));
        has_var = vars && vars->type == EXPR_FUNCTION && vars->data.function.head
            && vars->data.function.head->type == EXPR_SYMBOL
            && vars->data.function.head->data.symbol.name == SYM_List
            && vars->data.function.arg_count > 0;
        if (vars) expr_free(vars);
    }
    if (den) expr_free(den);
    if (!neg || !has_var) return e;

    Expr* num = eval_and_free(expr_new_function(expr_new_symbol(SYM_Numerator),
        (Expr*[]){ expr_copy(e) }, 1));
    Expr* den2 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Denominator),
        (Expr*[]){ expr_copy(e) }, 1));
    expr_free(e);
    /* Expand the negation so -1 distributes into the Plus (an un-distributed
     * Times[-1, Plus[...]] would re-surface the sign through Power[...,-1] and
     * defeat the normalisation). */
    Expr* nneg = eval_and_free(expr_new_function(expr_new_symbol(SYM_Expand),
        (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_new_integer(-1), num }, 2) }, 1));
    Expr* dneg = eval_and_free(expr_new_function(expr_new_symbol(SYM_Expand),
        (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_new_integer(-1), den2 }, 2) }, 1));
    Expr* inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ dneg, expr_new_integer(-1) }, 2));
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ nneg, inv }, 2));
}

/* Public entry point.  Returns the simplified Expr (caller owns) or
 * NULL when this path can't be applied / doesn't help.  See the comment
 * block above for algorithm details. */
static Expr* qa_cancel_with_poly_radical_impl(const Expr* arg) {
    if (!arg) return NULL;

    /* Step 1: detect. */
    const Expr** rads = NULL;
    int n = 0, cap = 0;
    bool overflow = false;
    polyrad_collect_distinct(arg, &rads, &n, &cap, &overflow);
    if (overflow || n != 1 || !rads) {
        free(rads);
        return NULL;
    }
    const Expr* alpha_render = rads[0];
    const Expr* radicand = NULL;
    int64_t q = 0;
    if (!polyrad_parse_radical(alpha_render, &radicand, &q) || q < 2) {
        free(rads);
        return NULL;
    }

    /* Step 2: substitute alpha -> fresh symbol S. */
    const char* gen_name = intern_symbol("$qa$polyrad$");
    Expr* gen_sym = expr_new_symbol(gen_name);
    Expr* substituted = expr_subst(arg, alpha_render, gen_sym);

    /* Step 3: Together (no extension). */
    Expr* together_call = expr_new_function(
        expr_new_symbol(SYM_Together),
        (Expr*[]){substituted}, 1);
    Expr* combined = evaluate(together_call);
    expr_free(together_call);
    if (!combined) {
        expr_free(gen_sym);
        free(rads);
        return NULL;
    }

    /* Step 4: extract num and den. */
    Expr* num_call = expr_new_function(expr_new_symbol(SYM_Numerator),
        (Expr*[]){expr_copy(combined)}, 1);
    Expr* num_expr = evaluate(num_call);
    expr_free(num_call);
    Expr* den_call = expr_new_function(expr_new_symbol(SYM_Denominator),
        (Expr*[]){expr_copy(combined)}, 1);
    Expr* den_expr = evaluate(den_call);
    expr_free(den_call);
    expr_free(combined);

    if (!num_expr || !den_expr) {
        if (num_expr) expr_free(num_expr);
        if (den_expr) expr_free(den_expr);
        expr_free(gen_sym);
        free(rads);
        return NULL;
    }

    /* Step 5: build the relation polynomial S^q - radicand and reduce
     * num and den modulo it w.r.t. S. */
    Expr* gen_pow_args[2] = { expr_new_symbol(gen_name),
                               expr_new_integer(q) };
    Expr* gen_pow = expr_new_function(expr_new_symbol(SYM_Power),
                                      gen_pow_args, 2);
    Expr* gen_pow_e = evaluate(gen_pow);
    expr_free(gen_pow);
    Expr* neg_radicand_args[2] = { expr_new_integer(-1),
                                    expr_copy((Expr*)radicand) };
    Expr* neg_radicand = expr_new_function(expr_new_symbol(SYM_Times),
                                            neg_radicand_args, 2);
    Expr* relation_args[2] = { gen_pow_e, neg_radicand };
    Expr* relation_call = expr_new_function(expr_new_symbol(SYM_Plus),
                                             relation_args, 2);
    Expr* relation = evaluate(relation_call);
    expr_free(relation_call);

    Expr* num_rem_args[3] = { num_expr,
                               expr_copy(relation),
                               expr_new_symbol(gen_name) };
    Expr* num_rem_call = expr_new_function(
        expr_new_symbol(SYM_PolynomialRemainder), num_rem_args, 3);
    Expr* num_reduced = evaluate(num_rem_call);
    expr_free(num_rem_call);

    Expr* den_rem_args[3] = { den_expr,
                               expr_copy(relation),
                               expr_new_symbol(gen_name) };
    Expr* den_rem_call = expr_new_function(
        expr_new_symbol(SYM_PolynomialRemainder), den_rem_args, 3);
    Expr* den_reduced = evaluate(den_rem_call);
    expr_free(den_rem_call);

    if (!num_reduced || !den_reduced) {
        if (num_reduced) expr_free(num_reduced);
        if (den_reduced) expr_free(den_reduced);
        expr_free(relation);
        expr_free(gen_sym);
        free(rads);
        return NULL;
    }

    /* Step 5.5 (rationalisation): if the post-reduction denominator
     * still contains S (i.e. its S-degree is > 0), invert it modulo
     * (S^q - radicand) via PolynomialExtendedGCD.  When the gcd comes
     * out S-free, the substitution N/D -> N*u (mod relation) clears S
     * from the denominator, producing the canonical
     * `(N reduced)/gcd_e` form.  When the gcd has S (meaning D shares
     * a factor with S^q - radicand — the expression is undefined at
     * the corresponding root), skip rationalisation and keep the un-
     * rationalised result. */
    {
        Expr* sym_e = expr_new_symbol(gen_name);
        int den_deg_in_s = get_degree_poly(den_reduced, sym_e);
        expr_free(sym_e);
        if (den_deg_in_s > 0) {
            Expr* xgcd_args[3] = { expr_copy(den_reduced),
                                    expr_copy(relation),
                                    expr_new_symbol(gen_name) };
            Expr* xgcd_call = expr_new_function(
                expr_new_symbol(SYM_PolynomialExtendedGCD), xgcd_args, 3);
            Expr* xgcd_result = evaluate(xgcd_call);
            expr_free(xgcd_call);

            if (xgcd_result
                && xgcd_result->type == EXPR_FUNCTION
                && xgcd_result->data.function.head
                && xgcd_result->data.function.head->type == EXPR_SYMBOL
                && xgcd_result->data.function.head->data.symbol.name == SYM_List
                && xgcd_result->data.function.arg_count == 2) {
                Expr* gcd_e = xgcd_result->data.function.args[0];
                Expr* coeffs = xgcd_result->data.function.args[1];
                if (coeffs && coeffs->type == EXPR_FUNCTION
                    && coeffs->data.function.head
                    && coeffs->data.function.head->type == EXPR_SYMBOL
                    && coeffs->data.function.head->data.symbol.name == SYM_List
                    && coeffs->data.function.arg_count >= 1) {
                    Expr* u_in_s = coeffs->data.function.args[0];

                    /* gcd must be S-free for the rationalisation to
                     * actually clear S from the denominator. */
                    Expr* sym_c = expr_new_symbol(gen_name);
                    int gcd_deg_in_s = get_degree_poly(gcd_e, sym_c);
                    expr_free(sym_c);
                    if (gcd_deg_in_s == 0) {
                        /* prod = num * u; Together; then PolynomialRemainder
                         * the numerator mod relation. */
                        Expr* prod_args[2] = { expr_copy(num_reduced),
                                                expr_copy(u_in_s) };
                        Expr* prod_call = expr_new_function(
                            expr_new_symbol(SYM_Times), prod_args, 2);
                        Expr* prod = evaluate(prod_call);
                        expr_free(prod_call);

                        Expr* prod_tg_call = expr_new_function(
                            expr_new_symbol(SYM_Together),
                            (Expr*[]){prod}, 1);
                        Expr* prod_combined = evaluate(prod_tg_call);
                        expr_free(prod_tg_call);

                        Expr* prod_num_call = expr_new_function(
                            expr_new_symbol(SYM_Numerator),
                            (Expr*[]){expr_copy(prod_combined)}, 1);
                        Expr* prod_num = evaluate(prod_num_call);
                        expr_free(prod_num_call);
                        Expr* prod_den_call = expr_new_function(
                            expr_new_symbol(SYM_Denominator),
                            (Expr*[]){expr_copy(prod_combined)}, 1);
                        Expr* prod_den = evaluate(prod_den_call);
                        expr_free(prod_den_call);
                        expr_free(prod_combined);

                        Expr* nrem_args[3] = { prod_num,
                                                expr_copy(relation),
                                                expr_new_symbol(gen_name) };
                        Expr* nrem_call = expr_new_function(
                            expr_new_symbol(SYM_PolynomialRemainder),
                            nrem_args, 3);
                        Expr* num_new = evaluate(nrem_call);
                        expr_free(nrem_call);

                        /* new_den = gcd_e * prod_den, simplified via
                         * Together so any rational parts canonicalise. */
                        Expr* dprod_args[2] = { expr_copy(gcd_e),
                                                 prod_den };
                        Expr* dprod_call = expr_new_function(
                            expr_new_symbol(SYM_Times), dprod_args, 2);
                        Expr* new_den_raw = evaluate(dprod_call);
                        expr_free(dprod_call);
                        Expr* den_tg_call = expr_new_function(
                            expr_new_symbol(SYM_Together),
                            (Expr*[]){new_den_raw}, 1);
                        Expr* den_new = evaluate(den_tg_call);
                        expr_free(den_tg_call);

                        if (num_new && den_new) {
                            /* Only accept the rationalised candidate if
                             * it didn't inflate the combined size.  For
                             * q=2 the conjugate trick reliably shrinks
                             * the structure; for q >= 3 the norm-style
                             * cofactor can produce a result that is
                             * larger in total leaf count than the
                             * un-rationalised form (q=2: a + b S -> a^2 -
                             * b^2 r drops S without bloating; q=3:
                             * a + b S + c S^2 -> norm with up to 4
                             * monomials in (a, b, c, r)).  Keep whichever
                             * is smaller. */
                            int64_t lc_unrat = polyrad_leaf_count(num_reduced)
                                             + polyrad_leaf_count(den_reduced);
                            int64_t lc_rat = polyrad_leaf_count(num_new)
                                           + polyrad_leaf_count(den_new);
                            if (lc_rat <= lc_unrat) {
                                expr_free(num_reduced);
                                expr_free(den_reduced);
                                num_reduced = num_new;
                                den_reduced = den_new;
                            } else {
                                expr_free(num_new);
                                expr_free(den_new);
                            }
                        } else {
                            if (num_new) expr_free(num_new);
                            if (den_new) expr_free(den_new);
                        }
                    }
                }
            }
            if (xgcd_result) expr_free(xgcd_result);
        }
    }

    expr_free(relation);

    /* Step 6: cancel common factor in Q[params, vars, S]. */
    {
        Expr* gcd_args[3] = { expr_copy(num_reduced),
                               expr_copy(den_reduced),
                               expr_new_symbol(gen_name) };
        Expr* gcd_call = expr_new_function(
            expr_new_symbol(SYM_PolynomialGCD), gcd_args, 3);
        Expr* gcd_e = evaluate(gcd_call);
        expr_free(gcd_call);
        /* If gcd is a non-trivial polynomial in S (or params/vars),
         * divide both sides by it.  When gcd is a unit (Integer or
         * Rational != 0), skip — division by a constant would just
         * scale num and den identically. */
        bool is_unit = false;
        if (gcd_e) {
            int64_t gp, gq;
            if (gcd_e->type == EXPR_INTEGER
                || gcd_e->type == EXPR_BIGINT) is_unit = true;
            else if (is_rational(gcd_e, &gp, &gq)) is_unit = true;
        }
        if (gcd_e && !is_unit) {
            Expr* nq_args[3] = { expr_copy(num_reduced),
                                  expr_copy(gcd_e),
                                  expr_new_symbol(gen_name) };
            Expr* nq_call = expr_new_function(
                expr_new_symbol(SYM_PolynomialQuotient), nq_args, 3);
            Expr* num_div = evaluate(nq_call);
            expr_free(nq_call);

            Expr* dq_args[3] = { expr_copy(den_reduced),
                                  expr_copy(gcd_e),
                                  expr_new_symbol(gen_name) };
            Expr* dq_call = expr_new_function(
                expr_new_symbol(SYM_PolynomialQuotient), dq_args, 3);
            Expr* den_div = evaluate(dq_call);
            expr_free(dq_call);

            if (num_div && den_div) {
                expr_free(num_reduced);
                expr_free(den_reduced);
                num_reduced = num_div;
                den_reduced = den_div;
            } else {
                if (num_div) expr_free(num_div);
                if (den_div) expr_free(den_div);
            }
        }
        if (gcd_e) expr_free(gcd_e);
    }

    /* Step 7: substitute S -> alpha_render back. */
    Expr* num_back = expr_subst(num_reduced, gen_sym, alpha_render);
    Expr* den_back = expr_subst(den_reduced, gen_sym, alpha_render);
    expr_free(num_reduced);
    expr_free(den_reduced);
    expr_free(gen_sym);

    /* Build num / den (with den==1 shortcut). */
    Expr* result;
    bool den_is_one = (den_back->type == EXPR_INTEGER
                      && den_back->data.integer == 1);
    if (den_is_one) {
        expr_free(den_back);
        result = evaluate(num_back);
        expr_free(num_back);
    } else {
        Expr* inv_args[2] = { den_back, expr_new_integer(-1) };
        Expr* inv = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Power), inv_args, 2));
        Expr* times_args[2] = { num_back, inv };
        result = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Times), times_args, 2));
    }

    free(rads);

    /* Step 8: change-detector.  If the result is structurally identical
     * to the input (no reduction fired), return NULL so the caller
     * falls back to the no-extension path — avoids wrapping the input
     * in a meaningless Times/Power layer. */
    if (result && expr_eq(result, (Expr*)arg)) {
        expr_free(result);
        return NULL;
    }
    /* Cheap "did not shrink" guard: bail if the result's leaf count is
     * larger than the input's by more than the typical Together
     * combine-then-reduce slack.  This keeps non-helpful expansions
     * from being returned. */
    if (result) {
        int64_t lc_in  = polyrad_leaf_count(arg);
        int64_t lc_out = polyrad_leaf_count(result);
        if (lc_in > 0 && lc_out > 2 * lc_in && lc_out > 40) {
            expr_free(result);
            return NULL;
        }
    }
    return result;
}

struct Expr* qa_cancel_with_poly_radical(const struct Expr* arg) {
    return qa_cancel_with_poly_radical_impl((const Expr*)arg);
}

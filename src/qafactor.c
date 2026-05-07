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

/* ============================== Phase G5 ============================== */

/* Internal symbol used as the α-placeholder when round-tripping through
 * picocas's polynomial machinery (Coefficient, Expand).  Picked to be
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

/* Recognise the imaginary unit, in either of its picocas forms:
 * the bare symbol I or the canonical Complex[0, 1]. */
static bool expr_is_imaginary_unit(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL && e->data.symbol == SYM_I) return true;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Complex
        && e->data.function.arg_count == 2) {
        Expr* re = e->data.function.args[0];
        Expr* im = e->data.function.args[1];
        if (re->type == EXPR_INTEGER && re->data.integer == 0
            && im->type == EXPR_INTEGER && im->data.integer == 1) return true;
    }
    return false;
}

/* Recognise Sqrt[c] with integer c, returning c via *out_c.  picocas
 * canonicalises Sqrt[c] as Power[c, Rational[1, 2]] so we accept that
 * form too. */
static bool expr_is_sqrt_int(const Expr* e, long* out_c) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (!e->data.function.head
        || e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.arg_count < 1) return false;
    const char* head = e->data.function.head->data.symbol;
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

    /* I  →  P_α(y) = y² + 1, render = I.  picocas auto-evaluates the
     * literal `I` into `Complex[0, 1]` and `Sqrt[-1]` likewise, so by
     * the time we see the option both forms have collapsed to the
     * canonical Complex[0, 1].  We accept the bare symbol I as a
     * defensive convenience. */
    if (expr_is_imaginary_unit(alpha_expr)) {
        QAExt* ext = qaext_new(2);
        qaext_set_coef_si(ext, 0,  1, 1);   /* +1 */
        qaext_set_coef_si(ext, 1,  0, 1);
        qaext_set_coef_si(ext, 2,  1, 1);   /* y² */
        *render_out = expr_new_symbol("I");
        return ext;
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

    /* Times[I, Sqrt[c]]  =  Sqrt[-c]  →  P_α(y) = y² + c.  picocas
     * auto-evaluates Sqrt[-c] for negative argument into I·Sqrt[c],
     * so `Extension -> Sqrt[-3]` arrives here as Times[Complex[0,1],
     * Sqrt[3]].  We accept either ordering of the two factors. */
    if (alpha_expr->type == EXPR_FUNCTION
        && alpha_expr->data.function.head
        && alpha_expr->data.function.head->type == EXPR_SYMBOL
        && alpha_expr->data.function.head->data.symbol == SYM_Times
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

    /* c^(1/n)  →  P_α(y) = yⁿ − c, render = c^(1/n). */
    if (alpha_expr->type == EXPR_FUNCTION
        && alpha_expr->data.function.head
        && alpha_expr->data.function.head->type == EXPR_SYMBOL
        && alpha_expr->data.function.head->data.symbol == SYM_Power
        && alpha_expr->data.function.arg_count == 2) {
        Expr* c   = alpha_expr->data.function.args[0];
        Expr* exp = alpha_expr->data.function.args[1];
        int64_t p, q;
        if (c->type == EXPR_INTEGER
            && is_rational(exp, &p, &q)
            && p == 1 && q >= 2) {
            QAExt* ext = qaext_root_si((long)c->data.integer, (unsigned)q);
            *render_out = expr_copy((Expr*)alpha_expr);
            return ext;
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

/* Lift an Expr polynomial in `var` (whose coefficients may be
 * polynomials in `alpha_render`) to a QAUPoly over `ext`.  Returns
 * NULL if any coefficient is not in Q(α) (e.g. contains a free symbol
 * other than `var` and `alpha_render`).
 *
 * Strategy: substitute alpha_render → an internal placeholder symbol,
 * then use Coefficient[..] to peel off (var^i)(α^j) terms. */
static QAUPoly* qa_expr_to_qaupoly_with_alpha(const Expr* poly,
                                              const Expr* var,
                                              const Expr* alpha_render,
                                              const QAExt* ext) {
    Expr* alpha_sym = expr_new_symbol(QA_ALPHA_INTERNAL);

    /* Stage 1 — substitute α's surface form with an opaque symbol. */
    Expr* poly_sub = alpha_render
        ? expr_subst(poly, alpha_render, alpha_sym)
        : expr_copy((Expr*)poly);

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
        if ((size_t)adeg >= ext->deg) {
            /* α-degree exceeds the extension's; refuse.  In practice
             * picocas auto-reduces e.g. Sqrt[2]^3 → 2 Sqrt[2] before
             * we even see the input, so this branch only fires when
             * the input is genuinely outside Q(α). */
            expr_free(c_expanded);
            ok = false;
            break;
        }

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

Expr* qaupoly_to_expr_alpha(const QAUPoly* f,
                            const char* x_name,
                            const Expr* alpha_render) {
    /* Round-trip through qaupoly_to_expr with an internal y-symbol,
     * then ReplaceAll y → alpha_render and evaluate.  Evaluation
     * collapses Sqrt[c]^k / c^(k/n) etc. into canonical form.
     *
     * Compound `alpha_render` (Phase G6 — γ = α_1 + s_2 α_2 + ... is a
     * Plus expression) needs an extra Expand pass: picocas's evaluator
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
 * each factor back as a picocas Expr using `alpha_render_output` as the
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

    const char* x_name = var->data.symbol;
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
            terms[out++] = expr_new_function(expr_new_symbol("Power"),
                                             p_args, 2);
        }
    }

    Expr* result;
    if (n_terms == 0)      result = expr_new_integer(1);
    else if (n_terms == 1) result = terms[0];
    else                   result = expr_new_function(expr_new_symbol("Times"),
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

/* Reify a QAExt's minimal polynomial as a picocas Expr in `y_name`. */
static Expr* qaext_min_poly_expr(const QAExt* ext, const char* y_name) {
    return mpq_array_to_poly_expr(ext->coef, ext->deg + 1, y_name);
}

/* Build (y_name − s·z_name) as an Expr. */
static Expr* expr_y_minus_s_z(const char* y_name, const char* z_name, int s) {
    if (s == 0) return expr_new_symbol(y_name);
    Expr* sz_args[2] = { expr_new_integer((int64_t)s),
                         expr_new_symbol(z_name) };
    Expr* sz = expr_new_function(expr_new_symbol("Times"), sz_args, 2);
    Expr* neg_args[2] = { expr_new_integer(-1), sz };
    Expr* neg_sz = expr_new_function(expr_new_symbol("Times"), neg_args, 2);
    Expr* sum_args[2] = { expr_new_symbol(y_name), neg_sz };
    return expr_new_function(expr_new_symbol("Plus"), sum_args, 2);
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

/* Compute Res_z(Q, P) ∈ Q[w] via picocas's internal_resultant. */
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
            term = expr_new_function(expr_new_symbol("Times"), args, 2);
        }
        Expr* args[2] = { t->gamma_render, term };
        new_gamma_render = expr_new_function(expr_new_symbol("Plus"), args, 2);
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

    /* Substitute each user-side α_i (its surface form) with its
     * polynomial-in-γ_internal Expr in the input.  After this pass,
     * `poly_internal`'s only "extension symbol" is QA_ALPHA_INTERNAL,
     * so qa_factor_inner can lift it directly without further surface-
     * level substitution. */
    Expr* poly_internal = expr_copy((Expr*)poly);
    for (int i = 0; i < t->n; i++) {
        Expr* alpha_in_gamma_int = qanum_to_expr_in_gamma_sym(
            t->alphas[i], QA_ALPHA_INTERNAL);
        Expr* old = poly_internal;
        poly_internal = expr_subst(old, t->alpha_renders[i], alpha_in_gamma_int);
        expr_free(old);
        expr_free(alpha_in_gamma_int);
    }

    Expr* result = qa_factor_inner(
        poly_internal, t->ext,
        /* alpha_render_input = */ NULL,    /* poly already uses internal sym */
        /* alpha_render_output = */ t->gamma_render,
        var);

    expr_free(poly_internal);
    qa_tower_free(t);
    return result;
}

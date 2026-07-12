/* risch_hypertangent.c — the hypertangent case (Bronstein §5.10).
 *
 * See risch_hypertangent.h.  IntegrateHypertangentPolynomial reduces p modulo
 * derivatives (PolynomialReduce, since a hypertangent is a nonlinear monomial)
 * to a remainder of degree <= 1, then reads off the log-term coefficient
 * c = coeff(r, t) / (2a), where a = Dt/(t^2+1).
 */

#include "risch_hypertangent.h"

#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_intern.h"
#include "risch_field.h"

#include <stdlib.h>

/* ------------------------------------------------------------------ */
/* Evaluation helpers.                                                 */
/* ------------------------------------------------------------------ */

static Expr* rh_cp(const Expr* e) { return expr_copy((Expr*)e); }
static Expr* rh_eval_adopt(Expr* call) { Expr* r = evaluate(call); expr_free(call); return r; }
static Expr* rh_fn(const char* head, Expr** args, size_t n) {
    return rh_eval_adopt(expr_new_function(expr_new_symbol(head), args, n));
}
static Expr* rh_call1(const char* head, Expr* a) { return rh_fn(head, (Expr*[]){ a }, 1); }
static Expr* rh_call3(const char* head, Expr* a, Expr* b, Expr* c) {
    return rh_fn(head, (Expr*[]){ a, b, c }, 3);
}
static Expr* rh_times(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol("Times"), (Expr*[]){ a, b }, 2);
}
static Expr* rh_pow(Expr* a, long n) {
    return expr_new_function(expr_new_symbol("Power"), (Expr*[]){ a, expr_new_integer(n) }, 2);
}

/* ------------------------------------------------------------------ */
/* IntegrateHypertangentPolynomial (§5.10, p.167).                     */
/* ------------------------------------------------------------------ */

bool risch_integrate_hypertangent_poly(const Expr* p, const Expr* t,
                                       const RischDeriv* d, Expr** qo, Expr** co) {
    const Expr* Dt = risch_deriv_lookup(d, t);
    if (!Dt) return false;
    /* a = Dt/(t^2+1) must lie in k (t is hypertangent). */
    Expr* tsq1 = rh_eval_adopt(expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ rh_pow(rh_cp(t), 2), expr_new_integer(1) }, 2));
    Expr* a = rh_call1("Cancel", rh_times(rh_cp(Dt), rh_pow(rh_cp(tsq1), -1)));
    expr_free(tsq1);
    /* a must be free of t (element of the base field). */
    Expr* fq = rh_fn("FreeQ", (Expr*[]){ rh_cp(a), rh_cp(t) }, 2);
    bool a_in_k = fq && fq->type == EXPR_SYMBOL && fq->data.symbol == intern_symbol("True");
    expr_free(fq);
    if (!a_in_k) { expr_free(a); return false; }

    Expr* q = NULL; Expr* r = NULL;
    if (!risch_field_polynomial_reduce(p, t, d, &q, &r)) { expr_free(a); return false; }

    /* c = coeff(r, t) / (2 a).  (deg_t(r) <= 1 after the reduction.) */
    Expr* r1 = rh_call3("Coefficient", r, rh_cp(t), expr_new_integer(1));  /* adopts r */
    Expr* twoa = rh_eval_adopt(rh_times(expr_new_integer(2), a));          /* adopts a */
    Expr* c = rh_call1("Cancel", rh_times(r1, rh_pow(twoa, -1)));

    *qo = q;
    *co = c;
    return true;
}

/* ------------------------------------------------------------------ */
/* Builtin.                                                            */
/* ------------------------------------------------------------------ */

/* Risch`IntegrateHypertangentPolynomial[p, t, deriv] -> {q, c}. */
static Expr* builtin_risch_int_hypertan_poly(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    const Expr* p = res->data.function.args[0];
    const Expr* t = res->data.function.args[1];
    if (t->type != EXPR_SYMBOL) return NULL;
    RischDeriv d;
    if (!risch_deriv_from_rules(res->data.function.args[2], &d)) return NULL;
    Expr* q = NULL; Expr* c = NULL;
    bool ok = risch_integrate_hypertangent_poly(p, t, &d, &q, &c);
    risch_deriv_free(&d);
    if (!ok) return NULL;
    return expr_new_function(expr_new_symbol("List"), (Expr*[]){ q, c }, 2);
}

/* ------------------------------------------------------------------ */
/* Registration.                                                       */
/* ------------------------------------------------------------------ */

static void rh_install(const char* name, Expr* (*fn)(Expr*), const char* doc) {
    symtab_add_builtin(name, fn);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    if (doc) symtab_set_docstring(name, doc);
}

void integrate_risch_hypertangent_init(void) {
    rh_install("Risch`IntegrateHypertangentPolynomial", builtin_risch_int_hypertan_poly,
        "Risch`IntegrateHypertangentPolynomial[p, t, deriv] integrates the\n"
        "polynomial p over a hypertangent monomial t (Dt = a (t^2+1)): returns\n"
        "{q, c} with p - D[q] - c D[t^2+1]/(t^2+1) in the base field, so that the\n"
        "integral is q + c Log[t^2+1] plus a base-field integral (Bronstein 5.10).\n"
        "If D[c] is nonzero the integral is not elementary.");
}

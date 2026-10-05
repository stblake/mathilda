/* integrate_intrep.c
 *
 * Definite integration by recognising a classical INTEGRAL REPRESENTATION of a
 * special function.  See integrate_intrep.h.  Each family is matched
 * structurally, its convergence gate is proved from the caller's Assumptions,
 * and the closed form is emitted directly -- correct-by-construction, no
 * NIntegrate crosscheck.  All families live on the half line [0, Infinity).
 */

#include "integrate_intrep.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"

#include <string.h>
#include <stdbool.h>

/* ---- construction / evaluation helpers ---------------------------------- */

static Expr* cp(const Expr* e) { return expr_copy((Expr*)e); }
static Expr* mk_int(long v) { return expr_new_integer((int64_t)v); }
static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }

static Expr* mk_fn(const char* head, Expr** args, size_t n) {
    return expr_new_function(expr_new_symbol(head), args, n);
}
static Expr* mk_fn1(const char* h, Expr* a) { Expr* v[1]={a}; return mk_fn(h,v,1); }
static Expr* mk_fn2(const char* h, Expr* a, Expr* b) { Expr* v[2]={a,b}; return mk_fn(h,v,2); }
static Expr* mk_fn3(const char* h, Expr* a, Expr* b, Expr* c) { Expr* v[3]={a,b,c}; return mk_fn(h,v,3); }

static Expr* Times_(Expr* a, Expr* b) { return mk_fn2("Times", a, b); }
static Expr* Plus_(Expr* a, Expr* b)  { return mk_fn2("Plus", a, b); }
static Expr* Pow_(Expr* a, Expr* b)   { return mk_fn2("Power", a, b); }
static Expr* mk_rat(long p, long q)   { return mk_fn2("Rational", mk_int(p), mk_int(q)); }

static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);   /* evaluate does not free its input */
    expr_free(call);
    return r;
}
static Expr* ev1(const char* name, Expr* a) { return eval_take(mk_fn1(name, a)); }
static Expr* simp(Expr* e) { return ev1("Simplify", e); }
static Expr* simp2(Expr* e, Expr* as) {
    if (!as) return simp(e);
    return eval_take(mk_fn2("Simplify", e, cp(as)));
}

static bool head_name_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}
static bool sym_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}
static bool is_symbol(const Expr* e, const Expr* x) {
    return e && e->type == EXPR_SYMBOL && x->type == EXPR_SYMBOL &&
           e->data.symbol.name == x->data.symbol.name;
}
static bool contains_symbol(const Expr* e, const Expr* x) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return is_symbol(e, x);
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_symbol(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_symbol(e->data.function.args[i], x)) return true;
    return false;
}
static bool free_of_x(const Expr* e, const Expr* x) { return !contains_symbol(e, x); }

static bool is_zero_expr(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
}
/* b == +Infinity (bare symbol Infinity, or DirectedInfinity[1]). */
static bool is_pos_inf(const Expr* b) {
    if (sym_is(b, "Infinity")) return true;
    if (head_name_is(b, "DirectedInfinity") && b->data.function.arg_count == 1) {
        Expr* d = b->data.function.args[0];
        return d->type == EXPR_INTEGER && d->data.integer == 1;
    }
    return false;
}

/* Prove a predicate under `as` via Simplify: +1 True, -1 False, 0 undecided. */
static int prove(Expr* pred, Expr* as) {
    Expr* s = simp2(pred, as);
    int r = 0;
    if (sym_is(s, "True")) r = 1;
    else if (sym_is(s, "False")) r = -1;
    expr_free(s);
    return r;
}
static bool prove_pos(Expr* e, Expr* as) {   /* e > 0 ? */
    int r = prove(mk_fn2("Greater", cp(e), mk_int(0)), as);
    return r == 1;
}
static bool prove_real(Expr* e, Expr* as) {  /* e in Reals ? (a literal real passes) */
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL ||
        head_name_is(e, "Rational")) return true;
    int r = prove(mk_fn2("Element", cp(e), mk_sym("Reals")), as);
    return r == 1;
}

/* Coefficient[e, x, k] (owned). */
static Expr* coeff(const Expr* e, const Expr* x, long k) {
    return eval_take(mk_fn3("Coefficient", cp(e), cp(x), mk_int(k)));
}
/* True iff Simplify[e] === 0. */
static bool is_zero_simp(Expr* e) {   /* consumes e */
    Expr* s = simp(e);
    bool z = (s->type == EXPR_INTEGER && s->data.integer == 0);
    expr_free(s);
    return z;
}

/* Collect the top-level product factors of f into fv[] (borrowed pointers),
 * returns the count (<= cap). A non-Times f is a single factor. */
static size_t collect_factors(Expr* f, Expr** fv, size_t cap) {
    size_t n = 0;
    if (head_name_is(f, "Times")) {
        for (size_t i = 0; i < f->data.function.arg_count && n < cap; i++)
            fv[n++] = f->data.function.args[i];
    } else if (n < cap) {
        fv[n++] = f;
    }
    return n;
}

/* The exponent g of a factor E^g (Power[E, g]); NULL if not such a factor. */
static Expr* exp_arg(Expr* fac) {
    if (head_name_is(fac, "Power") && fac->data.function.arg_count == 2 &&
        sym_is(fac->data.function.args[0], "E"))
        return fac->data.function.args[1];
    return NULL;
}

/* ---- (1) Laplace-Bessel: E^(-p x) BesselJ[nu, q x] ---------------------- */
static Expr* rec_laplace_bessel(Expr* f, Expr* x, Expr* as) {
    Expr* fv[16]; size_t n = collect_factors(f, fv, 16);
    Expr* C = mk_int(1);
    Expr* g = NULL;        /* exponent of E^g */
    Expr* nu = NULL, *barg = NULL;
    for (size_t i = 0; i < n; i++) {
        Expr* fac = fv[i];
        if (free_of_x(fac, x)) { C = Times_(C, cp(fac)); continue; }
        Expr* e = exp_arg(fac);
        if (e && !g) { g = e; continue; }                 /* E^g */
        if (head_name_is(fac, "BesselJ") && fac->data.function.arg_count == 2 &&
            !nu) { nu = fac->data.function.args[0]; barg = fac->data.function.args[1]; continue; }
        expr_free(C); return NULL;                          /* unrecognised x-factor */
    }
    if (!g || !nu) { expr_free(C); return NULL; }
    /* g must be -p x (linear, no constant); barg must be q x. */
    Expr* g1 = coeff(g, x, 1), *g0 = coeff(g, x, 0);
    Expr* q  = coeff(barg, x, 1), *b0 = coeff(barg, x, 0);
    bool shape = is_zero_simp(cp(g0)) && is_zero_simp(cp(b0)) &&
                 is_zero_simp(Plus_(cp(g), Times_(mk_int(-1), Times_(cp(g1), cp(x))))) &&
                 is_zero_simp(Plus_(cp(barg), Times_(mk_int(-1), Times_(cp(q), cp(x)))));
    expr_free(g0); expr_free(b0);
    Expr* p = simp(Times_(mk_int(-1), cp(g1)));    /* p = -Coefficient[g,x,1] */
    expr_free(g1);
    if (!shape || !prove_pos(p, as)) {
        expr_free(C); expr_free(p); expr_free(q); return NULL;
    }
    /* value = C (Sqrt[q^2+p^2]-p)^nu / (q^nu Sqrt[q^2+p^2]); nu=0 -> C/Sqrt[q^2+p^2]
     * (q sign irrelevant). For nu != 0 the q^nu factor needs q > 0. */
    Expr* root = mk_fn1("Sqrt", Plus_(Pow_(cp(q), mk_int(2)), Pow_(cp(p), mk_int(2))));
    Expr* val;
    bool nu_zero = (nu->type == EXPR_INTEGER && nu->data.integer == 0);
    if (nu_zero) {
        val = Times_(C, Pow_(root, mk_int(-1)));
    } else if (prove_pos(q, as)) {
        val = Times_(C, Times_(
                  Pow_(Plus_(cp(root), Times_(mk_int(-1), cp(p))), cp(nu)),
                  Times_(Pow_(cp(q), Times_(mk_int(-1), cp(nu))), Pow_(root, mk_int(-1)))));
    } else { expr_free(C); expr_free(root); expr_free(p); expr_free(q); return NULL; }
    expr_free(p); expr_free(q);
    return simp2(val, as);
}

/* ---- (2) Bessel-K (cosh): E^(-A Cosh[x]) Cosh[n x] ---------------------- */
static Expr* rec_besselk_cosh(Expr* f, Expr* x, Expr* as) {
    Expr* fv[16]; size_t n = collect_factors(f, fv, 16);
    Expr* C = mk_int(1);
    Expr* g = NULL, *narg = NULL;
    for (size_t i = 0; i < n; i++) {
        Expr* fac = fv[i];
        if (free_of_x(fac, x)) { C = Times_(C, cp(fac)); continue; }
        Expr* e = exp_arg(fac);
        if (e && !g) { g = e; continue; }
        if (head_name_is(fac, "Cosh") && fac->data.function.arg_count == 1 && !narg) {
            narg = fac->data.function.args[0]; continue;
        }
        expr_free(C); return NULL;
    }
    if (!g) { expr_free(C); return NULL; }
    /* g = -A Cosh[x]: A = -g/Cosh[x], must be x-free.  The outer Cosh[n x] factor
     * may be ABSENT (n = 0: Exp[-A Cosh[x]] alone -> BesselK[0, A]). */
    Expr* A = simp(Times_(mk_int(-1), Times_(cp(g), Pow_(mk_fn1("Cosh", cp(x)), mk_int(-1)))));
    Expr* nn;
    bool shape;
    if (narg) {
        Expr* nc = coeff(narg, x, 1), *n0 = coeff(narg, x, 0);  /* narg = n x */
        shape = A && free_of_x(A, x) && is_zero_simp(cp(n0)) &&
                is_zero_simp(Plus_(cp(narg), Times_(mk_int(-1), Times_(cp(nc), cp(x)))));
        expr_free(n0);
        nn = nc;
    } else {
        shape = A && free_of_x(A, x);
        nn = mk_int(0);
    }
    if (!shape || !prove_pos(A, as)) {
        expr_free(C); if (A) expr_free(A); expr_free(nn); return NULL;
    }
    Expr* val = Times_(C, mk_fn2("BesselK", cp(nn), cp(A)));
    expr_free(A); expr_free(nn);
    return simp2(val, as);
}

/* ---- (3) Bessel-K (exp): x^(nu-1) E^(-A x - B/x) ------------------------ */
static Expr* rec_besselk_exp(Expr* f, Expr* x, Expr* as) {
    Expr* fv[16]; size_t n = collect_factors(f, fv, 16);
    Expr* C = mk_int(1);
    Expr* g = NULL;
    Expr* rho = mk_int(0);        /* exponent of x (nu-1); default 0 */
    bool have_pow = false;
    for (size_t i = 0; i < n; i++) {
        Expr* fac = fv[i];
        if (free_of_x(fac, x)) { C = Times_(C, cp(fac)); continue; }
        Expr* e = exp_arg(fac);
        if (e && !g) { g = e; continue; }
        if (is_symbol(fac, x) && !have_pow) { expr_free(rho); rho = mk_int(1); have_pow = true; continue; }
        if (head_name_is(fac, "Power") && fac->data.function.arg_count == 2 &&
            is_symbol(fac->data.function.args[0], x) &&
            free_of_x(fac->data.function.args[1], x) && !have_pow) {
            expr_free(rho); rho = cp(fac->data.function.args[1]); have_pow = true; continue;
        }
        expr_free(C); expr_free(rho); return NULL;
    }
    if (!g) { expr_free(C); expr_free(rho); return NULL; }
    /* g = -A x - B/x.  Multiply by x: g x = -A x^2 - B  (no x^1, no x^3, ...). */
    Expr* gx = ev1("Expand", Times_(cp(g), cp(x)));
    Expr* c2 = coeff(gx, x, 2), *c1 = coeff(gx, x, 1), *c0 = coeff(gx, x, 0);
    Expr* A = simp(Times_(mk_int(-1), cp(c2)));
    Expr* B = simp(Times_(mk_int(-1), cp(c0)));
    bool shape = is_zero_simp(cp(c1)) &&
                 is_zero_simp(Plus_(cp(gx), Plus_(Times_(cp(A), Pow_(cp(x), mk_int(2))), cp(B))));
    expr_free(gx); expr_free(c2); expr_free(c1); expr_free(c0);
    if (!shape || !prove_pos(A, as) || !prove_pos(B, as)) {
        expr_free(C); expr_free(rho); expr_free(A); expr_free(B); return NULL;
    }
    Expr* nu = simp(Plus_(cp(rho), mk_int(1)));
    expr_free(rho);
    /* value = C 2 (B/A)^(nu/2) BesselK[nu, 2 Sqrt[A B]]. */
    Expr* val = Times_(C, Times_(mk_int(2), Times_(
        Pow_(Times_(cp(B), Pow_(cp(A), mk_int(-1))), Times_(mk_rat(1,2), cp(nu))),
        mk_fn2("BesselK", cp(nu), Times_(mk_int(2), mk_fn1("Sqrt", Times_(cp(A), cp(B))))))));
    expr_free(A); expr_free(B); expr_free(nu);
    /* Simplify, then TrigToExp as the LAST step (no trailing Simplify): for a
     * half-integer nu, Simplify renders BesselK as a Cosh/Sinh form and — a
     * Simplify complexity-metric quirk — would re-fold e^(-z) back into Cosh-Sinh
     * if run again, so TrigToExp has the final say and lands the elementary
     * Sqrt[Pi/A] e^(-2 Sqrt[A B]).  For a generic nu, TrigToExp is a no-op and the
     * simplified BesselK form survives. */
    return ev1("TrigToExp", simp2(val, as));
}

/* ---- (4) Airy: Cos[p x^3 + q x] ----------------------------------------- */
static Expr* rec_airy(Expr* f, Expr* x, Expr* as) {
    Expr* fv[16]; size_t n = collect_factors(f, fv, 16);
    Expr* C = mk_int(1);
    Expr* arg = NULL;
    for (size_t i = 0; i < n; i++) {
        Expr* fac = fv[i];
        if (free_of_x(fac, x)) { C = Times_(C, cp(fac)); continue; }
        if (head_name_is(fac, "Cos") && fac->data.function.arg_count == 1 && !arg) {
            arg = fac->data.function.args[0]; continue;
        }
        expr_free(C); return NULL;
    }
    if (!arg) { expr_free(C); return NULL; }
    /* arg = p x^3 + q x (no x^2, no constant). */
    Expr* p = coeff(arg, x, 3), *c2 = coeff(arg, x, 2), *q = coeff(arg, x, 1), *c0 = coeff(arg, x, 0);
    bool shape = is_zero_simp(cp(c2)) && is_zero_simp(cp(c0)) &&
                 is_zero_simp(Plus_(cp(arg), Times_(mk_int(-1),
                     Plus_(Times_(cp(p), Pow_(cp(x), mk_int(3))), Times_(cp(q), cp(x))))));
    expr_free(c2); expr_free(c0);
    if (!shape || !prove_pos(p, as) || !prove_real(q, as)) {
        expr_free(C); expr_free(p); expr_free(q); return NULL;
    }
    /* value = C Pi (3p)^(-1/3) AiryAi[q (3p)^(-1/3)]. */
    Expr* s = Pow_(Times_(mk_int(3), cp(p)), mk_rat(-1,3));
    Expr* val = Times_(C, Times_(mk_sym("Pi"), Times_(cp(s),
                    mk_fn1("AiryAi", Times_(cp(q), cp(s))))));
    expr_free(p); expr_free(q); expr_free(s);
    return simp2(val, as);
}

/* ---- entry -------------------------------------------------------------- */

Expr* integrate_intrep_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!is_zero_expr(a) || !is_pos_inf(b)) return NULL;   /* half line [0,Inf) only */
    if (!contains_symbol(f, x)) return NULL;

    Expr* v;
    if ((v = rec_laplace_bessel(f, x, assumptions))) return v;
    if ((v = rec_besselk_cosh(f, x, assumptions)))   return v;
    if ((v = rec_besselk_exp(f, x, assumptions)))    return v;
    if ((v = rec_airy(f, x, assumptions)))           return v;
    return NULL;
}

/* ---- builtin ------------------------------------------------------------ */

Expr* builtin_integrate_intrep(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];
    if (!head_name_is(spec, "List") || spec->data.function.arg_count != 3) return NULL;
    Expr* x = spec->data.function.args[0];
    Expr* a = spec->data.function.args[1];
    Expr* b = spec->data.function.args[2];
    if (x->type != EXPR_SYMBOL) return NULL;
    Expr* as = NULL;
    for (size_t t = 2; t < res->data.function.arg_count; t++) {
        Expr* opt = res->data.function.args[t];
        if (opt->type == EXPR_FUNCTION && opt->data.function.arg_count == 2 &&
            opt->data.function.head->type == EXPR_SYMBOL &&
            (opt->data.function.head->data.symbol.name == SYM_Rule ||
             opt->data.function.head->data.symbol.name == SYM_RuleDelayed) &&
            opt->data.function.args[0]->type == EXPR_SYMBOL &&
            strcmp(opt->data.function.args[0]->data.symbol.name, "Assumptions") == 0) {
            as = opt->data.function.args[1]; continue;
        }
        return NULL;
    }
    return integrate_intrep_try(f, x, a, b, as);
}

void integrate_intrep_init(void) {
    symtab_add_builtin("Integrate`IntegralRepresentation", builtin_integrate_intrep);
    symtab_get_def("Integrate`IntegralRepresentation")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Integrate`IntegralRepresentation",
        "Integrate`IntegralRepresentation[f, {x, 0, Infinity}] evaluates a "
        "half-line integral that is a classical integral representation of a "
        "special function: E^(-p x) BesselJ[nu, q x] -> Laplace-Bessel; "
        "E^(-A Cosh[x]) Cosh[n x] -> BesselK[n, A]; x^(nu-1) E^(-A x - B/x) -> "
        "2 (B/A)^(nu/2) BesselK[nu, 2 Sqrt[A B]]; Cos[p x^3 + q x] -> "
        "Pi (3p)^(-1/3) AiryAi[q (3p)^(-1/3)].  Each is gated on its convergence "
        "condition (p, A, B > 0) proved from Assumptions.  Returns unevaluated "
        "when the integrand or interval is not of a recognised form.");
}

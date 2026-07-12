/* risch_canonical.c — splitting factorization (Bronstein §3.5).
 *
 * See risch_canonical.h.  All polynomial arithmetic is over the field k = C(x,
 * lower monomials): gcd and exact division are the field operations from
 * risch_field.c, so pure-lower-field factors behave as units.
 */

#include "risch_canonical.h"

#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_intern.h"
#include "risch_field.h"

#include <stdio.h>
#include <stdlib.h>

/* const-accepting deep copy (expr_copy takes a mutable pointer). */
static Expr* rc_cp(const Expr* e) { return expr_copy((Expr*)e); }

/* ------------------------------------------------------------------ */
/* Local evaluation helpers (mirror risch_field.c).                    */
/* ------------------------------------------------------------------ */

static Expr* rc_eval_adopt(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}
static Expr* rc_fn(const char* head, Expr** args, size_t n) {
    return rc_eval_adopt(expr_new_function(expr_new_symbol(head), args, n));
}
static Expr* rc_call1(const char* head, Expr* a) {
    return rc_fn(head, (Expr*[]){ a }, 1);
}
static Expr* rc_call2(const char* head, Expr* a, Expr* b) {
    return rc_fn(head, (Expr*[]){ a, b }, 2);
}
static Expr* rc_times(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol("Times"), (Expr*[]){ a, b }, 2);
}
static Expr* rc_pow(Expr* a, long n) {
    return expr_new_function(expr_new_symbol("Power"),
        (Expr*[]){ a, expr_new_integer(n) }, 2);
}
static Expr* rc_neg(Expr* a) {
    return expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_new_integer(-1), a }, 2);
}
/* d p / d t (ordinary partial derivative w.r.t. the monomial variable). */
static Expr* rc_ddt(const Expr* p, const Expr* t) {
    return rc_call2("D", rc_cp(p), rc_cp(t));
}
/* Expand[a - b], adopting neither (copies). */
static Expr* rc_sub_expand(const Expr* a, const Expr* b) {
    Expr* diff = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ rc_cp(a), rc_neg(rc_cp(b)) }, 2);
    return rc_call1("Expand", rc_eval_adopt(diff));
}

/* ------------------------------------------------------------------ */
/* SplitFactor (Bronstein §3.5, p.100):                                */
/*   S <- gcd(p, Dp) / gcd(p, dp/dt)                                    */
/*   if deg(S) = 0: return (p, 1)                                       */
/*   (q_n, q_s) <- SplitFactor(p/S); return (q_n, S q_s)               */
/* ------------------------------------------------------------------ */

void risch_split_factor(const Expr* p, const Expr* t, const RischDeriv* d,
                        Expr** p_n, Expr** p_s) {
    Expr* Dp = risch_field_deriv(p, d);
    Expr* dpdt = rc_ddt(p, t);
    Expr* g1 = risch_field_gcd_t(p, Dp, t);
    Expr* g2 = risch_field_gcd_t(p, dpdt, t);
    expr_free(Dp); expr_free(dpdt);
    Expr* S = risch_field_divexact_t(g1, g2, t);
    expr_free(g1); expr_free(g2);

    long dS = S ? risch_field_degree_t(S, t) : -1;
    if (dS <= 0) {                 /* p has no special factor: p_n = p, p_s = 1 */
        if (S) expr_free(S);
        *p_n = rc_cp(p);
        *p_s = expr_new_integer(1);
        return;
    }

    Expr* pOverS = risch_field_divexact_t(p, S, t);
    if (!pOverS) {                 /* defensive: should not happen (S | p) */
        expr_free(S);
        *p_n = rc_cp(p);
        *p_s = expr_new_integer(1);
        return;
    }

    Expr* qn = NULL; Expr* qs = NULL;
    risch_split_factor(pOverS, t, d, &qn, &qs);
    expr_free(pOverS);

    *p_n = qn;
    *p_s = rc_call1("Expand", rc_times(S, qs));  /* S * q_s (adopts S, qs) */
}

/* ------------------------------------------------------------------ */
/* CanonicalRepresentation (Bronstein §3.5, p.103):                    */
/*   (a, d) <- (numerator(f), denominator(f)), d monic                  */
/*   (q, r) <- PolyDivide(a, d)                                         */
/*   (d_n, d_s) <- SplitFactor(d)                                       */
/*   (b, c) <- Diophantine(d_n, d_s, r), deg(b) < deg(d_s)             */
/*   f = q + b/d_s + c/d_n                                              */
/* ------------------------------------------------------------------ */

void risch_canonical_representation(const Expr* f, const Expr* t, const RischDeriv* d,
                                    Expr** f_p, Expr** f_s, Expr** f_n) {
    Expr* a = NULL; Expr* den = NULL;
    risch_field_num_den_t(f, t, &a, &den);

    if (risch_field_degree_t(den, t) <= 0) {   /* denominator = 1: f is polynomial */
        expr_free(den);
        *f_p = a; *f_s = expr_new_integer(0); *f_n = expr_new_integer(0);
        return;
    }

    Expr* q = NULL; Expr* r = NULL;
    risch_field_divmod_t(a, den, t, &q, &r);   /* a = q den + r */
    expr_free(a);

    Expr* dn = NULL; Expr* ds = NULL;
    risch_split_factor(den, t, d, &dn, &ds);   /* den = dn ds */
    expr_free(den);

    if (risch_field_degree_t(r, t) < 0) {      /* r == 0: no fractional part */
        expr_free(r); expr_free(dn); expr_free(ds);
        *f_p = q; *f_s = expr_new_integer(0); *f_n = expr_new_integer(0);
        return;
    }

    Expr* b = NULL; Expr* c = NULL;
    risch_field_diophantine_t(dn, ds, r, t, &b, &c);  /* b dn + c ds = r */
    expr_free(r);

    *f_p = q;
    *f_s = rc_call1("Cancel", rc_times(b, rc_pow(ds, -1)));  /* b / d_s */
    *f_n = rc_call1("Cancel", rc_times(c, rc_pow(dn, -1)));  /* c / d_n */
}

/* ------------------------------------------------------------------ */
/* Squarefree factorization over k[t] (Yun's algorithm, using d/dt).   */
/* Returns a malloc'd array *out of `*count` owned factors p_1..p_m     */
/* with p = prod_i p_i^i, each p_i squarefree and pairwise coprime.     */
/* Returns NULL / count 0 for a constant p.                            */
/* ------------------------------------------------------------------ */

static Expr** risch_squarefree_t(const Expr* p, const Expr* t, size_t* count) {
    *count = 0;
    long dp = risch_field_degree_t(p, t);
    if (dp <= 0) return NULL;               /* constant in t: no factors */

    size_t cap = (size_t)dp + 1, n = 0;
    Expr** out = malloc(cap * sizeof(Expr*));

    Expr* pp = rc_ddt(p, t);                        /* p' = dp/dt */
    Expr* g  = risch_field_gcd_t(p, pp, t);         /* g  = gcd(p, p') */
    Expr* c  = risch_field_divexact_t(p, g, t);     /* c1 = p / g */
    Expr* ppg = risch_field_divexact_t(pp, g, t);   /* p'/g */
    expr_free(pp); expr_free(g);
    /* dd = p'/g - c1'  */
    Expr* c1d = rc_ddt(c, t);
    Expr* dd = rc_sub_expand(ppg, c1d);
    expr_free(ppg); expr_free(c1d);

    while (n < cap && risch_field_degree_t(c, t) > 0) {
        Expr* pi = risch_field_gcd_t(c, dd, t);     /* p_i = gcd(c, dd) */
        out[n++] = rc_cp(pi);
        Expr* cnew = risch_field_divexact_t(c, pi, t);
        Expr* ddq  = risch_field_divexact_t(dd, pi, t);
        expr_free(pi);
        if (!cnew || !ddq) {                        /* defensive: division failed */
            if (cnew) expr_free(cnew);
            if (ddq) expr_free(ddq);
            break;
        }
        Expr* cnewd = rc_ddt(cnew, t);
        Expr* ddnew = rc_sub_expand(ddq, cnewd);
        expr_free(ddq); expr_free(cnewd);
        expr_free(c); expr_free(dd);
        c = cnew; dd = ddnew;
    }
    expr_free(c); expr_free(dd);
    *count = n;
    if (n == 0) { free(out); return NULL; }
    return out;
}

/* ------------------------------------------------------------------ */
/* Builtins.                                                           */
/* ------------------------------------------------------------------ */

/* Shared argument decode for the (p, t, deriv) builtins. */
static bool rc_decode3(Expr* res, const Expr** p, const Expr** t, RischDeriv* d) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return false;
    *p = res->data.function.args[0];
    *t = res->data.function.args[1];
    if ((*t)->type != EXPR_SYMBOL) return false;
    return risch_deriv_from_rules(res->data.function.args[2], d);
}

/* Risch`Derivation[p, deriv] -> D[p]. */
static Expr* builtin_risch_derivation(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    RischDeriv d;
    if (!risch_deriv_from_rules(res->data.function.args[1], &d)) return NULL;
    Expr* r = risch_field_deriv(res->data.function.args[0], &d);
    risch_deriv_free(&d);
    return r;
}

/* Risch`SplitFactor[p, t, deriv] -> {p_n, p_s}. */
static Expr* builtin_risch_splitfactor(Expr* res) {
    const Expr *p, *t; RischDeriv d;
    if (!rc_decode3(res, &p, &t, &d)) return NULL;
    Expr* pn = NULL; Expr* ps = NULL;
    risch_split_factor(p, t, &d, &pn, &ps);
    risch_deriv_free(&d);
    return expr_new_function(expr_new_symbol("List"), (Expr*[]){ pn, ps }, 2);
}

/* Risch`SplitSquarefreeFactor[p, t, deriv] -> {{N_1,...}, {S_1,...}}.
 * From the squarefree factors p_i of p, S_i = gcd(p_i, D p_i) and
 * N_i = p_i / S_i (Bronstein §3.5, p.102). */
static Expr* builtin_risch_splitsquarefree(Expr* res) {
    const Expr *p, *t; RischDeriv d;
    if (!rc_decode3(res, &p, &t, &d)) return NULL;

    size_t m = 0;
    Expr** facs = risch_squarefree_t(p, t, &m);
    Expr** Ns = (m ? malloc(m * sizeof(Expr*)) : NULL);
    Expr** Ss = (m ? malloc(m * sizeof(Expr*)) : NULL);
    for (size_t i = 0; i < m; i++) {
        Expr* Dpi = risch_field_deriv(facs[i], &d);
        Ss[i] = risch_field_gcd_t(facs[i], Dpi, t);
        expr_free(Dpi);
        Expr* Ni = risch_field_divexact_t(facs[i], Ss[i], t);
        Ns[i] = Ni ? Ni : rc_cp(facs[i]);
        expr_free(facs[i]);
    }
    free(facs);
    risch_deriv_free(&d);

    Expr* Nlist = expr_new_function(expr_new_symbol("List"), Ns, m);
    Expr* Slist = expr_new_function(expr_new_symbol("List"), Ss, m);
    free(Ns); free(Ss);
    return expr_new_function(expr_new_symbol("List"),
        (Expr*[]){ Nlist, Slist }, 2);
}

/* Risch`CanonicalRepresentation[f, t, deriv] -> {f_p, f_s, f_n}. */
static Expr* builtin_risch_canonicalrep(Expr* res) {
    const Expr *f, *t; RischDeriv d;
    if (!rc_decode3(res, &f, &t, &d)) return NULL;
    Expr* fp = NULL; Expr* fs = NULL; Expr* fn = NULL;
    risch_canonical_representation(f, t, &d, &fp, &fs, &fn);
    risch_deriv_free(&d);
    return expr_new_function(expr_new_symbol("List"), (Expr*[]){ fp, fs, fn }, 3);
}

/* Risch`PolyDivide[a, b, t] -> {q, r} with a = q b + r, deg_t(r) < deg_t(b). */
static Expr* builtin_risch_polydivide(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    const Expr* a = res->data.function.args[0];
    const Expr* b = res->data.function.args[1];
    const Expr* t = res->data.function.args[2];
    if (t->type != EXPR_SYMBOL) return NULL;
    Expr* q = NULL; Expr* r = NULL;
    if (!risch_field_divmod_t(a, b, t, &q, &r)) return NULL;   /* b == 0 */
    return expr_new_function(expr_new_symbol("List"), (Expr*[]){ q, r }, 2);
}

/* Risch`NormalQ / Risch`SpecialQ. */
static Expr* builtin_risch_normalq(Expr* res) {
    const Expr *p, *t; RischDeriv d;
    if (!rc_decode3(res, &p, &t, &d)) return NULL;
    bool r = risch_field_is_normal(p, t, &d);
    risch_deriv_free(&d);
    return expr_new_symbol(r ? "True" : "False");
}
static Expr* builtin_risch_specialq(Expr* res) {
    const Expr *p, *t; RischDeriv d;
    if (!rc_decode3(res, &p, &t, &d)) return NULL;
    bool r = risch_field_is_special(p, t, &d);
    risch_deriv_free(&d);
    return expr_new_symbol(r ? "True" : "False");
}

/* ------------------------------------------------------------------ */
/* Registration.                                                       */
/* ------------------------------------------------------------------ */

static void rc_install(const char* name, Expr* (*fn)(Expr*), const char* doc) {
    symtab_add_builtin(name, fn);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    if (doc) symtab_set_docstring(name, doc);
}

void integrate_risch_canonical_init(void) {
    rc_install("Risch`Derivation", builtin_risch_derivation,
        "Risch`Derivation[p, deriv] applies the monomial derivation D given by\n"
        "deriv = {var -> Dvar, ...} to p: D[p] = Sum Dvar D[p, var].");
    rc_install("Risch`SplitFactor", builtin_risch_splitfactor,
        "Risch`SplitFactor[p, t, deriv] returns {pn, ps}, the splitting\n"
        "factorization of the polynomial p in the monomial t: p = pn ps with ps\n"
        "special and every squarefree factor of pn normal w.r.t. the derivation\n"
        "deriv = {var -> Dvar, ...} (Bronstein, Symbolic Integration I, 3.5).");
    rc_install("Risch`SplitSquarefreeFactor", builtin_risch_splitsquarefree,
        "Risch`SplitSquarefreeFactor[p, t, deriv] returns {{N1,...}, {S1,...}},\n"
        "squarefree factorizations of the normal and special parts of p, so that\n"
        "p = (N1 N2^2 ...)(S1 S2^2 ...) is a splitting factorization.");
    rc_install("Risch`CanonicalRepresentation", builtin_risch_canonicalrep,
        "Risch`CanonicalRepresentation[f, t, deriv] returns {fp, fs, fn} with\n"
        "f = fp + fs + fn: fp the polynomial part, fs = b/ds the special part\n"
        "(special denominator, proper), fn = c/dn the normal part (normal\n"
        "denominator, proper) w.r.t. the derivation deriv = {var -> Dvar, ...}.");
    rc_install("Risch`PolyDivide", builtin_risch_polydivide,
        "Risch`PolyDivide[a, b, t] returns {q, r} with a = q b + r and\n"
        "deg_t(r) < deg_t(b), dividing in k[t] over the field of the other\n"
        "variables (the coefficient ring C(x, ...) is treated as a field).");
    rc_install("Risch`NormalQ", builtin_risch_normalq,
        "Risch`NormalQ[p, t, deriv] tests whether the polynomial p is normal\n"
        "w.r.t. the derivation deriv, i.e. gcd(p, D[p]) has degree 0 in t.");
    rc_install("Risch`SpecialQ", builtin_risch_specialq,
        "Risch`SpecialQ[p, t, deriv] tests whether the polynomial p is special\n"
        "w.r.t. the derivation deriv, i.e. p divides D[p] in k[t].");
}

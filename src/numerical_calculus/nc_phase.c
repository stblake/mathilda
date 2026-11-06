/*
 * nc_phase.c -- linear-phase splitting shared by NSum and NIntegrate.
 *
 * A summand or integrand f(x) = sum_j F_j(x) e^(lam_j x) with phase-free
 * amplitudes F_j is summed (Abel-Plana) or integrated (contour rotation) group
 * by group, each group's oscillation turned into exponential decay by moving
 * off the real axis.  This file holds the symbolic side: finding a phase
 * (nc_has_phase), stripping the explicit factors b^(c x + d) from a product
 * term (nc_strip_phase), and making the phase of a Bessel function explicit
 * through the scaled Hankel functions (nc_bessel_to_hankel / nc_hankel_expand).
 *
 * Memory: every function returns a fresh tree (or NULL) and never edits its
 * input in place (nodes are shared by reference count).
 */

#include "nc_phase.h"

#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include "eval.h"
#include "numeric.h"
#include "sym_intern.h"
#include "sym_names.h"

/* small builders (adopt their arguments) */
static Expr* nc_mk1(const char* h, Expr* a) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a }, 1);
}
static Expr* nc_mk2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}
static bool nc_has_head(const Expr* e, const char* h) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == h;
}
/* Numericalise a variable-free constant (consumed) to a machine complex. */
static bool nc_const_machine(Expr* e, double _Complex* out) {
    Expr* raw = eval_and_free(e);
    Expr* num = raw ? numericalize(raw, numeric_machine_spec()) : NULL;
    expr_free(raw);
    bool ok = false;
    if (num) {
        if (num->type == EXPR_REAL) { *out = num->data.real; ok = true; }
        else if (num->type == EXPR_INTEGER) { *out = (double)num->data.integer; ok = true; }
        else {
            Expr *re, *im;
            if (is_complex(num, &re, &im)) {
                double r, i;
                Expr* nr = numericalize(re, numeric_machine_spec());
                Expr* ni = numericalize(im, numeric_machine_spec());
                ok = nr && ni && nr->type == EXPR_REAL && ni->type == EXPR_REAL;
                if (ok) { r = nr->data.real; i = ni->data.real; *out = r + I * i; }
                expr_free(nr); expr_free(ni);
            }
        }
    }
    expr_free(num);
    return ok;
}

bool nc_free_of(const Expr* e, const char* var) {
    if (!e) return true;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name != var;
    if (e->type != EXPR_FUNCTION) return true;
    if (!nc_free_of(e->data.function.head, var)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!nc_free_of(e->data.function.args[i], var)) return false;
    return true;
}

/* A var-dependent phase anywhere in `e`: b^ex with b free of var and ex not,
 * or a circular function of var. */
bool nc_has_phase(const Expr* e, const char* var) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (nc_has_head(e, SYM_Power) && e->data.function.arg_count == 2
        && nc_free_of(e->data.function.args[0], var)
        && !nc_free_of(e->data.function.args[1], var))
        return true;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if ((h == SYM_Sin || h == SYM_Cos || h == SYM_Tan || h == SYM_Cot
             || h == SYM_Sec || h == SYM_Csc) && !nc_free_of(e, var))
            return true;
        /* J_nu(z), Y_nu(z) of a real argument oscillate like e^(+-i z): a hidden
         * linear phase when z is linear in the index (nc_bessel_to_hankel). */
        if ((h == SYM_BesselJ || h == SYM_BesselY) && e->data.function.arg_count == 2
            && !nc_free_of(e->data.function.args[1], var))
            return true;
    }
    if (nc_has_phase(e->data.function.head, var)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (nc_has_phase(e->data.function.args[i], var)) return true;
    return false;
}

/* Split one evaluated product term into F(var)·e^(lam·var), F phase-free.
 * Returns F (new) and sets *lam (new, exact), or NULL when the term's phase
 * is not linear in var or is not a top-level factor. */
Expr* nc_strip_phase(const Expr* term, const char* var, Expr** lam) {
    Expr* const* fs;
    Expr* one[1];
    size_t n;
    if (nc_has_head(term, SYM_Times)) {
        fs = term->data.function.args; n = term->data.function.arg_count;
    } else {
        one[0] = (Expr*)term; fs = one; n = 1;
    }
    Expr** keep = malloc(sizeof(Expr*) * n);
    if (!keep) return NULL;
    size_t nk = 0;
    Expr* L = expr_new_integer(0);
    bool ok = true;
    for (size_t i = 0; i < n && ok; i++) {
        Expr* f = fs[i];
        if (nc_has_head(f, SYM_Power) && f->data.function.arg_count == 2
            && nc_free_of(f->data.function.args[0], var)
            && !nc_free_of(f->data.function.args[1], var)) {
            Expr* b = f->data.function.args[0];
            Expr* ex = f->data.function.args[1];
            Expr* c = eval_and_free(nc_mk2("D", expr_copy(ex), expr_new_symbol(var)));
            if (!c || !nc_free_of(c, var)) { expr_free(c); ok = false; break; }
            L = eval_and_free(nc_mk2("Plus", L,
                                     nc_mk2("Times", c, nc_mk1("Log", expr_copy(b)))));
            /* b^(c var + d) = b^d · e^(c Log[b] var) for the principal power. */
            Expr* ex0 = eval_and_free(nc_mk2("ReplaceAll", expr_copy(ex),
                                             nc_mk2("Rule", expr_new_symbol(var),
                                                    expr_new_integer(0))));
            keep[nk++] = nc_mk2("Power", expr_copy(b), ex0);
        } else {
            keep[nk++] = expr_copy(f);
        }
    }
    if (!ok) {
        for (size_t i = 0; i < nk; i++) expr_free(keep[i]);
        free(keep); expr_free(L);
        return NULL;
    }
    Expr* s = expr_new_function(expr_new_symbol(SYM_Times), keep, nk);
    free(keep);
    s = eval_and_free(s);
    if (!s || nc_has_phase(s, var)) { expr_free(s); expr_free(L); return NULL; }
    *lam = L;
    return s;
}

/* Bessel functions as phase-carrying factors.
 *
 * J_nu(z) and Y_nu(z) of a real argument oscillate like cos, sin(z - phase):
 * Sum J_0(n) is Sin[n]/Sqrt[n]-like, but the phase is hidden inside the
 * function, so neither the symbolic split above nor any partial-sum model sees
 * it (the extrapolators returned 0.3348 for Sum_{n>=1} J_0(n) = 1/2).  The phase
 * is made explicit through the Hankel functions,
 *     J = (H1 + H2)/2,   Y = (H1 - H2)/(2 i),   H1,2(z) = e^(+-i z) S1,2(z),
 * where the scaled Hankel functions S1(z) = H1(z) e^(-i z) and S2(z) = H2(z) e^(i z)
 * are analytic and non-oscillatory for Re z > 0 (S ~ Sqrt[2/(Pi z)] e^(-+i(nu Pi/2
 * + Pi/4))).  With z = c n + d, c > 0, each half is a linear-phase group
 * e^(+-i c n) F(n), which the Abel-Plana rotation sums with exponentially
 * decaying kernels.  On the rotated contour z = c(n0 +- i y) + d the products
 * H1 e^(-i z) are evaluated without cancellation through the modified Bessel
 * function, H1_nu(w) = (2/Pi) i^(-nu-1) K_nu(-i w) and
 * H2_nu(w) = (2/Pi) i^(nu+1) K_nu(i w), valid for -Pi/2 < arg w <= Pi and
 * -Pi < arg w <= Pi/2 respectively -- the whole right half-plane.  (J + i Y
 * there would cancel e^(2 Im w) to leave e^(-Im w).)
 *
 * During the split the scaled functions are inert placeholder heads, so the
 * phase stripping sees only the explicit e^(+-i z) factors; they are expanded to
 * the K_nu form once the groups are formed. */
static const char* nc_hs_head(int kind) {
    static const char* h1 = NULL;
    static const char* h2 = NULL;
    if (!h1) { h1 = intern_symbol("NSum`HankelScaled1"); h2 = intern_symbol("NSum`HankelScaled2"); }
    return kind == 1 ? h1 : h2;
}

/* Copy of e with every BesselJ/BesselY whose argument depends on var replaced
 * by its Hankel decomposition.  *bad is set when such a factor is not of the
 * supported form (order depending on var, argument not real-linear in var with
 * a positive slope, or not in the right half-plane from n_check on). */
Expr* nc_bessel_to_hankel(const Expr* e, const char* var, double n_check, bool* bad) {
    if (*bad) return NULL;
    if (!e || e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    const Expr* hd = e->data.function.head;
    size_t n = e->data.function.arg_count;
    if (hd->type == EXPR_SYMBOL
        && (hd->data.symbol.name == SYM_BesselJ || hd->data.symbol.name == SYM_BesselY)
        && n == 2 && !nc_free_of(e->data.function.args[1], var)) {
        Expr* nu = e->data.function.args[0];
        Expr* z = e->data.function.args[1];
        double _Complex cd, zc;
        if (!nc_free_of(nu, var)) { *bad = true; return NULL; }
        Expr* c = eval_and_free(nc_mk2("D", expr_copy(z), expr_new_symbol(var)));
        bool lin = c && nc_free_of(c, var) && nc_const_machine(expr_copy(c), &cd)
                   && cimag(cd) == 0.0 && creal(cd) > 0.0;
        expr_free(c);
        if (lin) {
            Expr* zn = nc_mk2("ReplaceAll", expr_copy(z),
                              nc_mk2("Rule", expr_new_symbol(var), expr_new_real(n_check)));
            lin = nc_const_machine(zn, &zc) && cimag(zc) == 0.0 && creal(zc) > 0.0;
        }
        if (!lin) { *bad = true; return NULL; }
        bool is_y = (hd->data.symbol.name == SYM_BesselY);
        Expr* iz = nc_mk2("Times", expr_new_symbol(SYM_I), expr_copy(z));
        Expr* miz = nc_mk2("Times", expr_new_integer(-1), expr_copy(iz));
        Expr* t1 = nc_mk2("Times", nc_mk2("Power", expr_new_symbol(SYM_E), iz),
                          nc_mk2(nc_hs_head(1), expr_copy(nu), expr_copy(z)));
        Expr* t2 = nc_mk2("Times", nc_mk2("Power", expr_new_symbol(SYM_E), miz),
                          nc_mk2(nc_hs_head(2), expr_copy(nu), expr_copy(z)));
        if (is_y) t2 = nc_mk2("Times", expr_new_integer(-1), t2);
        /* J: 1/2 (t1 + t2);  Y: 1/(2 i) (t1 - t2) = -i/2 (t1 - t2) */
        Expr* pref = is_y ? nc_mk2("Times", nc_mk2("Rational", expr_new_integer(-1), expr_new_integer(2)),
                                   expr_new_symbol(SYM_I))
                          : nc_mk2("Rational", expr_new_integer(1), expr_new_integer(2));
        return nc_mk2("Times", pref, nc_mk2("Plus", t1, t2));
    }
    Expr* nh = nc_bessel_to_hankel(hd, var, n_check, bad);
    Expr** v = malloc(sizeof(Expr*) * (n ? n : 1));
    if (!v) { expr_free(nh); *bad = true; return NULL; }
    for (size_t i = 0; i < n; i++) {
        v[i] = nc_bessel_to_hankel(e->data.function.args[i], var, n_check, bad);
        if (*bad) {
            for (size_t j = 0; j < i; j++) expr_free(v[j]);
            free(v); expr_free(nh);
            return NULL;
        }
    }
    Expr* out = expr_new_function(nh, v, n);
    free(v);
    return out;
}

/* Copy of e with the placeholder heads replaced by their modified-Bessel form:
 *   S1_nu(z) = (2/Pi) I^(-nu-1) K_nu(-I z) E^(-I z),
 *   S2_nu(z) = (2/Pi) I^(nu+1)  K_nu(I z)  E^(I z).
 * (A fresh tree: nodes are shared by reference count, never edited in place.) */
bool nc_has_hankel(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* hd = e->data.function.head;
    if (hd->type == EXPR_SYMBOL
        && (hd->data.symbol.name == nc_hs_head(1) || hd->data.symbol.name == nc_hs_head(2)))
        return true;
    if (nc_has_hankel(hd)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (nc_has_hankel(e->data.function.args[i])) return true;
    return false;
}

Expr* nc_hankel_expand(const Expr* e) {
    if (!nc_has_hankel(e)) return expr_copy((Expr*)e);
    const Expr* hd = e->data.function.head;
    size_t n = e->data.function.arg_count;
    int kind = 0;
    if (hd->type == EXPR_SYMBOL && n == 2)
        kind = hd->data.symbol.name == nc_hs_head(1) ? 1
             : hd->data.symbol.name == nc_hs_head(2) ? 2 : 0;
    if (kind) {
        const Expr* nu = e->data.function.args[0];
        const Expr* z = e->data.function.args[1];
        int sg = (kind == 1) ? -1 : 1;
        Expr* isz = nc_mk2("Times", nc_mk2("Times", expr_new_integer(sg), expr_new_symbol(SYM_I)),
                           expr_copy((Expr*)z));
        Expr* ipow = nc_mk2("Power", expr_new_symbol(SYM_I),
                            nc_mk2("Times", expr_new_integer(sg),
                                   nc_mk2("Plus", expr_copy((Expr*)nu), expr_new_integer(1))));
        return nc_mk2("Times",
                      nc_mk2("Times", nc_mk2("Times", expr_new_integer(2),
                                             nc_mk2("Power", expr_new_symbol(SYM_Pi),
                                                    expr_new_integer(-1))), ipow),
                      nc_mk2("Times", nc_mk2(SYM_BesselK, expr_copy((Expr*)nu), expr_copy(isz)),
                             nc_mk2("Power", expr_new_symbol(SYM_E), isz)));
    }
    Expr** v = malloc(sizeof(Expr*) * (n ? n : 1));
    if (!v) return expr_copy((Expr*)e);
    for (size_t i = 0; i < n; i++) v[i] = nc_hankel_expand(e->data.function.args[i]);
    Expr* out = expr_new_function(nc_hankel_expand(hd), v, n);
    free(v);
    return out;
}

